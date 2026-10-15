"""Audit gate: compare the HPC-recomputed Graf/ITU2005/Pastorek rainfall against the
published final rainfall product, variable by variable, over their common grid.

Usage: python -I audit_rainfall_identity.py <published.nc> <recomputed.nc> <report.json>

Streams in time chunks; never modifies inputs.
"""
import json
import sys
import numpy as np
import netCDF4


def ids(ds, name="cml_id"):
    v = ds.variables[name][:]
    out = []
    for item in np.asarray(v):
        if isinstance(item, bytes):
            s = item.decode("utf-8", "replace")
        elif isinstance(item, np.ndarray):  # char array
            s = b"".join(item.tolist()).decode("utf-8", "replace") if item.dtype.kind == "S" else "".join(item.tolist())
        else:
            s = str(item)
        out.append(s.rstrip("\x00").strip())
    return out


def main():
    pub_path, rec_path, report_out = sys.argv[1], sys.argv[2], sys.argv[3]
    report = {"published": pub_path, "recomputed": rec_path, "variables": {}}
    with netCDF4.Dataset(pub_path) as pub, netCDF4.Dataset(rec_path) as rec:
        report["published_dims"] = {k: len(v) for k, v in pub.dimensions.items()}
        report["recomputed_dims"] = {k: len(v) for k, v in rec.dimensions.items()}
        ids_p, ids_r = ids(pub), ids(rec)
        report["n_cml_published"] = len(ids_p)
        report["n_cml_recomputed"] = len(ids_r)
        common = [c for c in ids_p if c in set(ids_r)]
        report["n_cml_common"] = len(common)
        idx_p = {c: i for i, c in enumerate(ids_p)}
        idx_r = {c: i for i, c in enumerate(ids_r)}
        sel_p = np.array([idx_p[c] for c in common])
        sel_r = np.array([idx_r[c] for c in common])
        same_order = bool(np.all(sel_p == sel_r))
        report["cml_order_identical"] = same_order

        tp = pub.variables["time_5min"]
        tr = rec.variables["time_5min"]
        t_p, t_r = tp[:], tr[:]
        report["time_units_identical"] = bool(getattr(tp, "units", None) == getattr(tr, "units", None))
        report["n_time_published"] = int(t_p.size)
        report["n_time_recomputed"] = int(t_r.size)
        report["time_axes_identical"] = bool(t_p.size == t_r.size and np.array_equal(np.asarray(t_p), np.asarray(t_r)))
        if not report["time_axes_identical"]:
            raise SystemExit("time axes differ; chunked comparison not implemented for that case")

        data_vars = [v for v in pub.variables if v in rec.variables
                     and pub.variables[v].dimensions == rec.variables[v].dimensions
                     and "time_5min" in pub.variables[v].dimensions
                     and v != "time_5min"]
        report["compared_variables"] = data_vars
        nt = t_p.size
        step = 20000
        for name in data_vars:
            vp, vr = pub.variables[name], rec.variables[name]
            dims = vp.dimensions
            taxis = dims.index("time_5min")
            n_shared = 0
            n_p_only = 0
            n_r_only = 0
            n_diff = 0
            max_abs = 0.0
            sum_sq = 0.0
            for lo in range(0, nt, step):
                hi = min(lo + step, nt)
                sl = [slice(None)] * len(dims)
                sl[taxis] = slice(lo, hi)
                a = np.ma.filled(np.ma.asarray(vp[tuple(sl)]).astype("f8"), np.nan)
                b = np.ma.filled(np.ma.asarray(vr[tuple(sl)]).astype("f8"), np.nan)
                if not same_order:
                    caxis = dims.index("cml_id")
                    a = np.take(a, sel_p, axis=caxis)
                    b = np.take(b, sel_r, axis=caxis)
                fa, fb = np.isfinite(a), np.isfinite(b)
                both = fa & fb
                n_shared += int(both.sum())
                n_p_only += int((fa & ~fb).sum())
                n_r_only += int((fb & ~fa).sum())
                d = np.abs(a[both] - b[both])
                if d.size:
                    n_diff += int((d > 0).sum())
                    max_abs = max(max_abs, float(d.max()))
                    sum_sq += float((d * d).sum())
            report["variables"][name] = {
                "shared_finite_points": n_shared,
                "finite_only_in_published": n_p_only,
                "finite_only_in_recomputed": n_r_only,
                "points_with_any_difference": n_diff,
                "max_abs_difference": max_abs,
                "rmse": (sum_sq / n_shared) ** 0.5 if n_shared else None,
            }
            print(name, report["variables"][name], flush=True)

    identical = all(v["points_with_any_difference"] == 0 and v["finite_only_in_published"] == 0
                    and v["finite_only_in_recomputed"] == 0
                    for v in report["variables"].values())
    report["values_and_masks_identical_on_common_grid"] = identical
    with open(report_out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print("RESULT:", "IDENTICAL" if identical else "DIFFERENCES PRESENT")


if __name__ == "__main__":
    main()
