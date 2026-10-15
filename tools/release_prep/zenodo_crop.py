"""Physically crop one NetCDF product to 2021-04-01 00:00 .. 2021-10-31 23:59 (inclusive)
at its native time resolution, preserving dimensions, dtypes, fill values, compression,
chunking, variable attributes and all other metadata.

Usage: python -I zenodo_crop.py <role> <source.nc> <dest.nc>
Roles and their time dimension / expected cropped slots mirror src/essd_cz/catalog.py.

The ONLY metadata changes, both documented in the dry-run report:
  * global attributes `publication_time_start`/`publication_time_end` (when present)
    are updated to the cropped window so they do not mislabel the content;
  * a `history` line is appended describing the crop.
Values inside the window are copied bit-identically; nothing else is altered.
"""
import sys
from datetime import datetime, timedelta, timezone
import numpy as np
import netCDF4

START = datetime(2021, 4, 1)
END_EXCLUSIVE = datetime(2021, 11, 1)
CONTRACT = {
    "raw": ("time", 1, 308160),
    "qc": ("time", 1, 308160),
    "rain": ("time_5min", 5, 61632),
    "radar": ("time", 5, 61632),
    "gauge": ("time", 10, 30816),
    "temperature": ("time", 10, 30816),
    "humidity": ("time", 10, 30816),
}


def main():
    role, src_path, dst_path = sys.argv[1], sys.argv[2], sys.argv[3]
    tdim, step_minutes, expected = CONTRACT[role]
    with netCDF4.Dataset(src_path) as src:
        tv = src.variables[tdim]
        units, calendar = tv.units, getattr(tv, "calendar", "standard")
        x = np.asarray(tv[:], dtype="f8")
        low = netCDF4.date2num(START, units, calendar=calendar)
        high = netCDF4.date2num(END_EXCLUSIVE, units, calendar=calendar)
        i0, i1 = int(np.searchsorted(x, low)), int(np.searchsorted(x, high))
        n = i1 - i0
        if n != expected:
            raise SystemExit(f"{role}: selected {n} slots, expected {expected}")
        step = netCDF4.date2num(START + timedelta(minutes=step_minutes), units, calendar=calendar) - low
        sel = x[i0:i1]
        if not (abs(sel[0] - low) < 1e-6 and abs(sel[-1] - (high - step)) < 1e-6
                and np.all(np.abs(np.diff(sel) - step) < 1e-6)):
            raise SystemExit(f"{role}: April-October axis not exact/continuous")

        fmt = src.data_model
        with netCDF4.Dataset(dst_path, "w", format=fmt) as dst:
            # dimensions
            for name, dim in src.dimensions.items():
                size = n if name == tdim else (None if dim.isunlimited() else len(dim))
                dst.createDimension(name, size)
            # global attributes (ordered), with the two documented updates
            for key in src.ncattrs():
                dst.setncattr(key, src.getncattr(key))
            if "publication_time_start" in src.ncattrs():
                dst.setncattr("publication_time_start", START.strftime("%Y-%m-%dT%H:%M:%S"))
            if "publication_time_end" in src.ncattrs():
                dst.setncattr("publication_time_end",
                              (END_EXCLUSIVE - timedelta(minutes=step_minutes)).strftime("%Y-%m-%dT%H:%M:%S"))
            stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            note = (f"{stamp} physically cropped along {tdim} to 2021-04-01T00:00..2021-10-31T23:5"
                    f"{9 if step_minutes == 1 else 60 - step_minutes} UTC inclusive at native "
                    f"{step_minutes}-min resolution; values bit-identical to source "
                    f"{src_path.split('/')[-1]}; no other changes")
            old_hist = src.getncattr("history") if "history" in src.ncattrs() else None
            dst.setncattr("history", f"{note}\n{old_hist}" if old_hist else note)

            # variables
            for name, v in src.variables.items():
                fill = getattr(v, "_FillValue", None)
                kwargs = {}
                try:
                    filters = v.filters() or {}
                except Exception:
                    filters = {}
                if filters.get("zlib"):
                    kwargs["zlib"] = True
                    kwargs["complevel"] = filters.get("complevel", 4)
                    kwargs["shuffle"] = filters.get("shuffle", False)
                ch = None
                try:
                    ch = v.chunking()
                except Exception:
                    pass
                if isinstance(ch, list):
                    sizes = [n if d == tdim else len(src.dimensions[d]) for d in v.dimensions]
                    kwargs["chunksizes"] = [min(c, s) for c, s in zip(ch, sizes)]
                out = dst.createVariable(name, v.datatype, v.dimensions, fill_value=fill, **kwargs)
                for key in v.ncattrs():
                    if key == "_FillValue":
                        continue
                    out.setncattr(key, v.getncattr(key))
                out.set_auto_maskandscale(False)
                v.set_auto_maskandscale(False)
                if tdim in v.dimensions:
                    taxis = v.dimensions.index(tdim)
                    chunk = max(1, 200_000_000 // max(1, int(np.prod(
                        [len(src.dimensions[d]) for d in v.dimensions if d != tdim], dtype=np.int64)) * max(v.dtype.itemsize if hasattr(v.dtype, 'itemsize') else 8, 1)))
                    for lo in range(0, n, chunk):
                        hi = min(lo + chunk, n)
                        sl_in = [slice(None)] * v.ndim
                        sl_in[taxis] = slice(i0 + lo, i0 + hi)
                        sl_out = [slice(None)] * v.ndim
                        sl_out[taxis] = slice(lo, hi)
                        out[tuple(sl_out)] = v[tuple(sl_in)]
                else:
                    if v.ndim == 0:
                        out[...] = v[...]
                    else:
                        out[:] = v[:]
    print(f"{role}: wrote {n} native slots -> {dst_path}")


if __name__ == "__main__":
    main()
