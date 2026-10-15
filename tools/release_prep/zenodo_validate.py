"""Systematically validate a cropped Zenodo product against its full-history source.

Usage: python -I zenodo_validate.py <role> <source.nc> <cropped.nc> <report.json>

Checks (no filename assumptions):
  * cropped time axis == exactly the native Apr-Oct 2021 window, uniform step;
  * every dimension except time unchanged;
  * every variable present with identical dtype, dimensions, attributes and
    fill value;
  * every variable's raw values bit-identical to the source at the window
    indices (auto mask/scale disabled; NaNs compared as equal);
  * global attributes identical except exactly {publication_time_start,
    publication_time_end, history} where updates were documented.
"""
import json
import sys
from datetime import datetime, timedelta
import numpy as np
import netCDF4

START = datetime(2021, 4, 1)
END_EXCLUSIVE = datetime(2021, 11, 1)
CONTRACT = {
    "raw": ("time", 1, 308160), "qc": ("time", 1, 308160),
    "rain": ("time_5min", 5, 61632), "radar": ("time", 5, 61632),
    "gauge": ("time", 10, 30816), "temperature": ("time", 10, 30816),
    "humidity": ("time", 10, 30816),
}
ALLOWED_GLOBAL_DIFFS = {"publication_time_start", "publication_time_end", "history"}


def attrs_equal(x, y):
    """NaN-aware attribute equality (a NaN _FillValue must equal itself)."""
    ax, ay = np.asarray(x), np.asarray(y)
    if ax.shape != ay.shape:
        return False
    if ax.dtype.kind == "f" and ay.dtype.kind == "f":
        return bool(np.array_equal(ax, ay, equal_nan=True))
    try:
        return bool(np.array_equal(ax, ay))
    except Exception:
        return x == y


def arrays_identical(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    if a.dtype.kind == "f":
        return bool(np.array_equal(a, b, equal_nan=True))
    return bool(np.array_equal(a, b))


def main():
    role, src_path, crop_path, report_out = sys.argv[1:5]
    tdim, step_minutes, expected = CONTRACT[role]
    problems = []
    r = {"role": role, "source": src_path, "cropped": crop_path}
    with netCDF4.Dataset(src_path) as src, netCDF4.Dataset(crop_path) as cr:
        # window indices in source
        tv = src.variables[tdim]
        units, calendar = tv.units, getattr(tv, "calendar", "standard")
        x = np.asarray(tv[:], dtype="f8")
        low = netCDF4.date2num(START, units, calendar=calendar)
        high = netCDF4.date2num(END_EXCLUSIVE, units, calendar=calendar)
        i0, i1 = int(np.searchsorted(x, low)), int(np.searchsorted(x, high))

        # 1) time axis of cropped file
        tc = cr.variables[tdim]
        y = np.asarray(tc[:], dtype="f8")
        step = netCDF4.date2num(START + timedelta(minutes=step_minutes), units, calendar=calendar) - low
        r["n_time_cropped"] = int(y.size)
        if tc.units != units:
            problems.append("time units changed")
        if y.size != expected:
            problems.append(f"time axis has {y.size} slots, expected {expected}")
        else:
            if not abs(y[0] - low) < 1e-6:
                problems.append("cropped axis does not start 2021-04-01T00:00")
            if not abs(y[-1] - (high - step)) < 1e-6:
                problems.append("cropped axis does not end at last native slot before 2021-11-01")
            if not np.all(np.abs(np.diff(y) - step) < 1e-6):
                problems.append("cropped axis not uniform at native step")
        r["cropped_start"] = str(netCDF4.num2date(float(y[0]), units, calendar=calendar))
        r["cropped_end"] = str(netCDF4.num2date(float(y[-1]), units, calendar=calendar))

        # 2) dimensions
        for name, dim in src.dimensions.items():
            want = expected if name == tdim else len(dim)
            have = len(cr.dimensions[name]) if name in cr.dimensions else None
            if have != want:
                problems.append(f"dimension {name}: {have} != {want}")
        for name in cr.dimensions:
            if name not in src.dimensions:
                problems.append(f"new dimension {name}")

        # 3) global attributes
        g_src = {k: src.getncattr(k) for k in src.ncattrs()}
        g_cr = {k: cr.getncattr(k) for k in cr.ncattrs()}
        diffs = {k for k in set(g_src) | set(g_cr)
                 if not attrs_equal(g_src.get(k), g_cr.get(k))}
        r["global_attr_diffs"] = sorted(diffs)
        if not diffs <= ALLOWED_GLOBAL_DIFFS:
            problems.append(f"unexpected global attribute changes: {sorted(diffs - ALLOWED_GLOBAL_DIFFS)}")

        # 4) variables
        if set(src.variables) != set(cr.variables):
            problems.append(f"variable set changed: {sorted(set(src.variables) ^ set(cr.variables))}")
        n_checked = 0
        for name in src.variables:
            if name not in cr.variables:
                continue
            v, w = src.variables[name], cr.variables[name]
            if v.dimensions != w.dimensions:
                problems.append(f"{name}: dimensions {w.dimensions} != {v.dimensions}")
                continue
            if str(v.dtype) != str(w.dtype):
                problems.append(f"{name}: dtype {w.dtype} != {v.dtype}")
            va = {k: v.getncattr(k) for k in v.ncattrs()}
            wa = {k: w.getncattr(k) for k in w.ncattrs()}
            adiff = {k for k in set(va) | set(wa)
                     if not attrs_equal(va.get(k), wa.get(k))}
            if adiff:
                problems.append(f"{name}: attribute diffs {sorted(adiff)}")
            try:
                v.set_auto_maskandscale(False)
                w.set_auto_maskandscale(False)
            except Exception:
                pass
            if tdim in v.dimensions:
                taxis = v.dimensions.index(tdim)
                other = int(np.prod([len(src.dimensions[d]) for d in v.dimensions if d != tdim], dtype=np.int64))
                itemsize = v.dtype.itemsize if hasattr(v.dtype, "itemsize") else 8
                chunk = max(1, 200_000_000 // max(1, other * itemsize))
                ok = True
                for lo in range(0, expected, chunk):
                    hi = min(lo + chunk, expected)
                    sl_s = [slice(None)] * v.ndim
                    sl_s[taxis] = slice(i0 + lo, i0 + hi)
                    sl_c = [slice(None)] * v.ndim
                    sl_c[taxis] = slice(lo, hi)
                    if not arrays_identical(v[tuple(sl_s)], w[tuple(sl_c)]):
                        ok = False
                        break
                if not ok:
                    problems.append(f"{name}: values differ inside the window")
            else:
                if v.ndim == 0:
                    same = arrays_identical(v[...], w[...])
                else:
                    same = arrays_identical(v[:], w[:])
                if not same:
                    problems.append(f"{name}: static values differ")
            n_checked += 1
        r["n_variables_checked"] = n_checked

    r["problems"] = problems
    r["ok"] = not problems
    with open(report_out, "w", encoding="utf-8") as f:
        json.dump(r, f, indent=2)
    print(role, "OK" if not problems else f"PROBLEMS: {problems}")
    sys.exit(0 if not problems else 1)


if __name__ == "__main__":
    main()
