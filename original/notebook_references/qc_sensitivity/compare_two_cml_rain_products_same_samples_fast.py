#!/usr/bin/env python3
"""
Fast fair comparison of TWO CML rainfall products against the SAME radar samples.

Optimized version:
- never loads the full common link x time matrices into RAM
- streams over link/time blocks
- computes exact metrics incrementally
- uses the exact same mask for new and old products
- keeps a capped random sample only for plotting

Both products are evaluated on identical:
  * cml_id
  * timestamps
  * finite positions
  * rain-rate range
  * radar-wet positions (optional)
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

try:
    from tqdm import tqdm
except Exception:
    def tqdm(x, **kwargs):
        return x


# ============================================================
# Helpers
# ============================================================
def find_dim(obj, candidates: List[str], required: bool = True) -> Optional[str]:
    dims = list(obj.dims)
    for cand in candidates:
        if cand in dims:
            return cand
    for d in dims:
        dl = d.lower()
        if any(c.lower() in dl for c in candidates):
            return d
    if required:
        raise KeyError(f"Cannot find dimension matching {candidates}. Available dims: {dims}")
    return None


def coord_as_str(ds: xr.Dataset, dim: str) -> np.ndarray:
    if dim in ds.coords:
        return np.asarray(ds[dim].values).astype(str)
    return np.arange(ds.sizes[dim]).astype(str)


def resolve_coord_value(ds: xr.Dataset, dim: str, value: str) -> int:
    vals = [str(v) for v in ds[dim].values]
    if value not in vals:
        raise ValueError(
            f"Requested {value!r} not found on dimension {dim!r}. Available: {vals}"
        )
    return vals.index(value)


def select_combo(ds: xr.Dataset, var_name: str, args) -> xr.DataArray:
    if var_name not in ds:
        raise KeyError(f"{var_name!r} not found. Available data variables: {list(ds.data_vars)}")
    da = ds[var_name]
    sel = {}
    for candidates, value in [
        (["baseline_method", "baseline"], args.baseline_method),
        (["kr_model", "kr", "itu"], args.kr_model),
        (["waa_method", "waa"], args.waa_method),
    ]:
        dim = find_dim(da, candidates, required=False)
        if dim is not None:
            sel[dim] = resolve_coord_value(ds, dim, value)
    return da.isel(sel) if sel else da


def init_stats():
    return {
        "N": 0,
        "sum_c": 0.0,
        "sum_r": 0.0,
        "sum_c2": 0.0,
        "sum_r2": 0.0,
        "sum_cr": 0.0,
        "sum_diff": 0.0,
        "sum_abs_diff": 0.0,
        "sum_diff2": 0.0,
    }


def update_stats(stats, c, r):
    if c.size == 0:
        return
    c = c.astype("float64", copy=False)
    r = r.astype("float64", copy=False)
    d = c - r

    stats["N"] += int(c.size)
    stats["sum_c"] += float(np.sum(c))
    stats["sum_r"] += float(np.sum(r))
    stats["sum_c2"] += float(np.sum(c * c))
    stats["sum_r2"] += float(np.sum(r * r))
    stats["sum_cr"] += float(np.sum(c * r))
    stats["sum_diff"] += float(np.sum(d))
    stats["sum_abs_diff"] += float(np.sum(np.abs(d)))
    stats["sum_diff2"] += float(np.sum(d * d))


def finalize_stats(stats):
    n = int(stats["N"])
    if n == 0:
        return {
            "N": 0, "MBE": np.nan, "MAE": np.nan, "RMSE": np.nan,
            "PCC": np.nan, "mean_cml": np.nan, "mean_radar": np.nan
        }

    mean_c = stats["sum_c"] / n
    mean_r = stats["sum_r"] / n
    mbe = stats["sum_diff"] / n
    mae = stats["sum_abs_diff"] / n
    rmse = math.sqrt(max(stats["sum_diff2"] / n, 0.0))

    cov = stats["sum_cr"] / n - mean_c * mean_r
    var_c = stats["sum_c2"] / n - mean_c * mean_c
    var_r = stats["sum_r2"] / n - mean_r * mean_r
    pcc = cov / math.sqrt(var_c * var_r) if n >= 2 and var_c > 0 and var_r > 0 else np.nan

    return {
        "N": n,
        "MBE": mbe,
        "MAE": mae,
        "RMSE": rmse,
        "PCC": pcc,
        "mean_cml": mean_c,
        "mean_radar": mean_r,
    }


def plot_pair(r, c, metrics, title, out_png, out_pdf, args):
    fig, ax = plt.subplots(figsize=(5.0, 4.5))

    if args.plot_type == "scatter":
        ax.scatter(r, c, s=args.scatter_size, alpha=args.scatter_alpha, linewidths=0)
    else:
        hb = ax.hexbin(
            r, c,
            gridsize=args.gridsize,
            mincnt=1,
            bins="log",
            extent=(args.min_rain, args.max_rain, args.min_rain, args.max_rain),
        )
        fig.colorbar(hb, ax=ax, label="Count")

    ax.plot(
        [args.min_rain, args.max_rain],
        [args.min_rain, args.max_rain],
        ls="--", lw=1.0, color="0.25",
    )
    ax.set_xlim(args.min_rain, args.max_rain)
    ax.set_ylim(args.min_rain, args.max_rain)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"Radar rain rate [mm h$^{-1}$]")
    ax.set_ylabel(r"CML rain rate [mm h$^{-1}$]")
    ax.set_title(title)

    txt = (
        f"N={metrics['N']:,}\n"
        f"RMSE={metrics['RMSE']:.2f}\n"
        f"MBE={metrics['MBE']:.2f}\n"
        f"PCC={metrics['PCC']:.2f}"
    )
    ax.text(
        0.97, 0.03, txt,
        transform=ax.transAxes,
        ha="right", va="bottom", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.75", alpha=0.9),
    )

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(
        description="Fast fair same-sample comparison of two CML rain products against radar."
    )

    p.add_argument("--new-rain-nc", required=True)
    p.add_argument("--old-rain-nc", required=True)
    p.add_argument("--out-dir", required=True)

    p.add_argument("--cml-var", default="R_mean_5min")
    p.add_argument("--radar-var", default="R_radar")
    p.add_argument("--radar-source", choices=["new", "old"], default="new")

    p.add_argument("--baseline-method", default="pycml_constant")
    p.add_argument("--kr-model", default="ITU_2003")
    p.add_argument("--waa-method", default="pastorek_2021")

    p.add_argument("--new-label", default="Strict QC V2")
    p.add_argument("--old-label", default="Original product")

    p.add_argument("--wet-only", dest="wet_only", action="store_true", default=True)
    p.add_argument("--all-pairs", dest="wet_only", action="store_false")
    p.add_argument("--radar-wet-thr", type=float, default=0.1)

    p.add_argument("--min-rain", type=float, default=0.0)
    p.add_argument("--max-rain", type=float, default=100.0)

    p.add_argument("--tstart", default=None)
    p.add_argument("--tend", default=None)

    p.add_argument("--link-chunk", type=int, default=32)
    p.add_argument("--time-chunk", type=int, default=2016,
                   help="Number of 5-min timesteps per block. 2016 = 7 days.")

    p.add_argument("--max-plot-points", type=int, default=500000)
    p.add_argument("--random-seed", type=int, default=42)
    p.add_argument("--plot-type", choices=["scatter", "hexbin"], default="hexbin")
    p.add_argument("--gridsize", type=int, default=50)
    p.add_argument("--scatter-size", type=float, default=2.0)
    p.add_argument("--scatter-alpha", type=float, default=0.15)

    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Opening datasets lazily...")
    ds_new = xr.open_dataset(args.new_rain_nc)
    ds_old = xr.open_dataset(args.old_rain_nc)

    try:
        new_da = select_combo(ds_new, args.cml_var, args)
        old_da = select_combo(ds_old, args.cml_var, args)

        new_cdim = find_dim(new_da, ["cml_id", "cml"])
        old_cdim = find_dim(old_da, ["cml_id", "cml"])
        new_tdim = find_dim(new_da, ["time_5min", "time"])
        old_tdim = find_dim(old_da, ["time_5min", "time"])

        # This fast script is intentionally link-level.
        if find_dim(new_da, ["sublink_id", "sublink"], required=False) is not None:
            raise ValueError("Use --cml-var R_mean_5min for this fast comparison.")
        if find_dim(old_da, ["sublink_id", "sublink"], required=False) is not None:
            raise ValueError("Use --cml-var R_mean_5min for this fast comparison.")

        new_ids = coord_as_str(ds_new, new_cdim)
        old_ids = coord_as_str(ds_old, old_cdim)

        old_id_to_pos = {x: i for i, x in enumerate(old_ids)}
        common_ids = [x for x in new_ids if x in old_id_to_pos]
        if not common_ids:
            raise RuntimeError("No common cml_id values.")

        new_id_to_pos = {x: i for i, x in enumerate(new_ids)}
        new_cpos = np.asarray([new_id_to_pos[x] for x in common_ids], dtype=np.int64)
        old_cpos = np.asarray([old_id_to_pos[x] for x in common_ids], dtype=np.int64)

        new_times = pd.DatetimeIndex(ds_new[new_tdim].values)
        old_times = pd.DatetimeIndex(ds_old[old_tdim].values)

        old_time_to_pos = {np.datetime64(t): i for i, t in enumerate(old_times.values)}
        common_times = [np.datetime64(t) for t in new_times.values if np.datetime64(t) in old_time_to_pos]

        if args.tstart:
            ts = np.datetime64(pd.Timestamp(args.tstart))
            common_times = [t for t in common_times if t >= ts]
        if args.tend:
            te = np.datetime64(pd.Timestamp(args.tend))
            common_times = [t for t in common_times if t <= te]

        if not common_times:
            raise RuntimeError("No common timestamps.")

        new_time_to_pos = {np.datetime64(t): i for i, t in enumerate(new_times.values)}
        new_tpos = np.asarray([new_time_to_pos[t] for t in common_times], dtype=np.int64)
        old_tpos = np.asarray([old_time_to_pos[t] for t in common_times], dtype=np.int64)

        radar_ds = ds_new if args.radar_source == "new" else ds_old
        radar_var = radar_ds[args.radar_var]
        radar_cdim = find_dim(radar_var, ["cml_id", "cml"])
        radar_tdim = find_dim(radar_var, ["time_5min", "time"])

        radar_ids = coord_as_str(radar_ds, radar_cdim)
        radar_id_to_pos = {x: i for i, x in enumerate(radar_ids)}
        missing_radar_ids = [x for x in common_ids if x not in radar_id_to_pos]
        if missing_radar_ids:
            raise RuntimeError(f"{len(missing_radar_ids)} common CML ids missing from radar source.")

        radar_cpos = np.asarray([radar_id_to_pos[x] for x in common_ids], dtype=np.int64)

        radar_times = pd.DatetimeIndex(radar_ds[radar_tdim].values)
        radar_time_to_pos = {np.datetime64(t): i for i, t in enumerate(radar_times.values)}
        missing_radar_times = [t for t in common_times if t not in radar_time_to_pos]
        if missing_radar_times:
            raise RuntimeError(f"{len(missing_radar_times)} common timestamps missing from radar source.")
        radar_tpos = np.asarray([radar_time_to_pos[t] for t in common_times], dtype=np.int64)

        n_links = len(common_ids)
        n_times = len(common_times)

        print(f"Common links      : {n_links:,}")
        print(f"Common timestamps : {n_times:,}")
        print(f"Time range        : {common_times[0]} -> {common_times[-1]}")
        print(f"Link chunk        : {args.link_chunk}")
        print(f"Time chunk        : {args.time_chunk}")
        print(f"Wet-only          : {args.wet_only}")
        if args.wet_only:
            print(f"Radar wet mask    : {args.radar_var} > {args.radar_wet_thr:g} mm h-1")

        stats_new = init_stats()
        stats_old = init_stats()

        rng = np.random.default_rng(args.random_seed)
        sample_r = []
        sample_new = []
        sample_old = []
        sample_seen = 0
        sample_cap = int(args.max_plot_points)

        # Stream over blocks.
        total_blocks = math.ceil(n_links / args.link_chunk) * math.ceil(n_times / args.time_chunk)
        pbar = tqdm(total=total_blocks, desc="Fair comparison blocks")

        for c0 in range(0, n_links, args.link_chunk):
            c1 = min(c0 + args.link_chunk, n_links)

            nc = new_cpos[c0:c1]
            oc = old_cpos[c0:c1]
            rc = radar_cpos[c0:c1]

            for t0 in range(0, n_times, args.time_chunk):
                t1 = min(t0 + args.time_chunk, n_times)

                nt = new_tpos[t0:t1]
                ot = old_tpos[t0:t1]
                rt = radar_tpos[t0:t1]

                # Use orthogonal indexing via sequential isel to avoid building global arrays.
                new_block = (
                    new_da.isel({new_cdim: nc})
                    .isel({new_tdim: nt})
                    .transpose(new_cdim, new_tdim)
                    .values
                    .astype("float32", copy=False)
                )
                old_block = (
                    old_da.isel({old_cdim: oc})
                    .isel({old_tdim: ot})
                    .transpose(old_cdim, old_tdim)
                    .values
                    .astype("float32", copy=False)
                )
                radar_block = (
                    radar_var.isel({radar_cdim: rc})
                    .isel({radar_tdim: rt})
                    .transpose(radar_cdim, radar_tdim)
                    .values
                    .astype("float32", copy=False)
                )

                mask = (
                    np.isfinite(new_block)
                    & np.isfinite(old_block)
                    & np.isfinite(radar_block)
                    & (new_block >= args.min_rain) & (new_block <= args.max_rain)
                    & (old_block >= args.min_rain) & (old_block <= args.max_rain)
                    & (radar_block >= args.min_rain) & (radar_block <= args.max_rain)
                )

                if args.wet_only:
                    mask &= radar_block > args.radar_wet_thr

                if mask.any():
                    nv = new_block[mask]
                    ov = old_block[mask]
                    rv = radar_block[mask]

                    update_stats(stats_new, nv, rv)
                    update_stats(stats_old, ov, rv)

                    # Reservoir-style probabilistic plot sample.
                    nblock = rv.size
                    sample_seen += nblock

                    if sample_cap > 0:
                        keep_prob = min(1.0, sample_cap / max(sample_seen, 1))
                        keep = rng.random(nblock) < keep_prob
                        if keep.any():
                            sample_r.append(rv[keep].copy())
                            sample_new.append(nv[keep].copy())
                            sample_old.append(ov[keep].copy())

                            current = sum(x.size for x in sample_r)
                            if current > 2 * sample_cap:
                                rr = np.concatenate(sample_r)
                                nn = np.concatenate(sample_new)
                                oo = np.concatenate(sample_old)
                                take = rng.choice(rr.size, size=sample_cap, replace=False)
                                sample_r = [rr[take]]
                                sample_new = [nn[take]]
                                sample_old = [oo[take]]

                pbar.update(1)

        pbar.close()

        m_new = finalize_stats(stats_new)
        m_old = finalize_stats(stats_old)

        if m_new["N"] == 0:
            raise RuntimeError("No fair common samples remained.")

        # Final plotting sample cap.
        if sample_r:
            rplot = np.concatenate(sample_r)
            nplot = np.concatenate(sample_new)
            oplot = np.concatenate(sample_old)
        else:
            rplot = np.array([], dtype="float32")
            nplot = np.array([], dtype="float32")
            oplot = np.array([], dtype="float32")

        if sample_cap > 0 and rplot.size > sample_cap:
            take = rng.choice(rplot.size, size=sample_cap, replace=False)
            rplot = rplot[take]
            nplot = nplot[take]
            oplot = oplot[take]

        print("\nFAIR COMMON SAMPLE RESULTS")
        print(f"Exact common N: {m_new['N']:,}")
        print("\nNEW PRODUCT")
        print(pd.Series(m_new).to_string())
        print("\nOLD PRODUCT")
        print(pd.Series(m_old).to_string())

        metrics_df = pd.DataFrame([
            {
                "product": args.new_label,
                "file": args.new_rain_nc,
                "baseline_method": args.baseline_method,
                "kr_model": args.kr_model,
                "waa_method": args.waa_method,
                "wet_only": args.wet_only,
                "radar_wet_thr": args.radar_wet_thr if args.wet_only else np.nan,
                **m_new,
            },
            {
                "product": args.old_label,
                "file": args.old_rain_nc,
                "baseline_method": args.baseline_method,
                "kr_model": args.kr_model,
                "waa_method": args.waa_method,
                "wet_only": args.wet_only,
                "radar_wet_thr": args.radar_wet_thr if args.wet_only else np.nan,
                **m_old,
            },
        ])
        metrics_df.to_csv(out_dir / "fair_common_metrics.csv", index=False)

        summary_df = pd.DataFrame([{
            "common_cml_links": n_links,
            "common_timestamps": n_times,
            "common_time_start": str(common_times[0]),
            "common_time_end": str(common_times[-1]),
            "final_identical_sample_positions": int(m_new["N"]),
            "wet_only": args.wet_only,
            "radar_wet_thr": args.radar_wet_thr if args.wet_only else np.nan,
            "min_rain": args.min_rain,
            "max_rain": args.max_rain,
            "link_chunk": args.link_chunk,
            "time_chunk": args.time_chunk,
            "plot_sample_points": int(rplot.size),
        }])
        summary_df.to_csv(out_dir / "fair_common_sample_summary.csv", index=False)

        plot_pair(
            rplot, nplot, m_new, args.new_label,
            out_dir / "scatter_new_vs_radar.png",
            out_dir / "scatter_new_vs_radar.pdf",
            args,
        )
        plot_pair(
            rplot, oplot, m_old, args.old_label,
            out_dir / "scatter_old_vs_radar.png",
            out_dir / "scatter_old_vs_radar.pdf",
            args,
        )

        print("\nSaved:")
        print(out_dir / "scatter_new_vs_radar.png")
        print(out_dir / "scatter_old_vs_radar.png")
        print(out_dir / "fair_common_metrics.csv")
        print(out_dir / "fair_common_sample_summary.csv")

    finally:
        ds_new.close()
        ds_old.close()


if __name__ == "__main__":
    main()
