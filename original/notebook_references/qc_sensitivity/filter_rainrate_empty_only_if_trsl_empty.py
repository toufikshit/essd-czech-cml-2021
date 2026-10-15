#!/usr/bin/env python3
"""
Filter estimated rain-rate links ONLY when raw CML data is also empty
=====================================================================

Why this script exists
----------------------
A rain-rate estimate that is all 0 does NOT automatically mean the CML link is
bad. It may simply mean the year/period was dry for that link/reference mask.

Therefore, this script checks BOTH:

1. Estimated rain rate:
      selected R estimate is all NaN / all 0 / only NaN+0

2. Raw CML signal:
      TRSL for the same cml_id is all NaN / all 0 / only NaN+0

Default drop rule:
      drop only if rain estimate is empty AND raw TRSL is empty

This avoids deleting valid dry links where the rain estimate is all zero but
the CML signal exists.

Recommended count-only check
----------------------------
python filter_rainrate_empty_only_if_trsl_empty.py \
  --rain-nc cml_rainrate_pycml_constant_kr_waa_methods.nc \
  --cml-nc cml_preprocessed_v3.nc \
  --rain-var R_mean_5min \
  --baseline-method pycml_constant \
  --kr-model ITU_2003 \
  --waa-method pastorek_2021 \
  --count-only

Create filtered rain-rate NetCDF
--------------------------------
python filter_rainrate_empty_only_if_trsl_empty.py \
  --rain-nc cml_rainrate_pycml_constant_kr_waa_methods.nc \
  --cml-nc cml_preprocessed_v3.nc \
  --out-nc cml_rainrate_pycml_itu03_pastorek_trsl_checked.nc \
  --rain-var R_mean_5min \
  --baseline-method pycml_constant \
  --kr-model ITU_2003 \
  --waa-method pastorek_2021 \
  --overwrite

Filter and plot against gauge
-----------------------------
python filter_rainrate_empty_only_if_trsl_empty.py \
  --rain-nc cml_rainrate_pycml_constant_kr_waa_methods.nc \
  --cml-nc cml_preprocessed_v3.nc \
  --out-nc cml_rainrate_pycml_itu03_pastorek_trsl_checked.nc \
  --rain-var R_mean_5min \
  --baseline-method pycml_constant \
  --kr-model ITU_2003 \
  --waa-method pastorek_2021 \
  --run-plot \
  --reference gauge \
  --plot-script "cml_compare_single_combo_pycml_itu03_pastorek_hexbin_oneplot_reference (1).py" \
  --gauge-nc CHMI_10MIN_v1.nc \
  --cml-coords-nc cml_daily_nc_v2_merged.nc \
  --max-gauge-dist-km 1 \
  --overwrite
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# ============================================================
# Basic helpers
# ============================================================
def as_array(x):
    if np.ma.isMaskedArray(x):
        return x.filled(np.nan)
    return np.asarray(x)


def values_to_str(arr) -> List[str]:
    arr = as_array(arr)
    if arr.dtype.kind == "S":
        return [x.decode("utf-8", errors="ignore") for x in arr]
    return [str(x) for x in arr]


def find_dim_from_candidates(dims: Sequence[str], candidates: Sequence[str], required: bool = True) -> Optional[str]:
    dims = list(dims)
    for cand in candidates:
        for d in dims:
            if d == cand:
                return d
    for d in dims:
        dl = d.lower()
        if any(c.lower() in dl for c in candidates):
            return d
    if required:
        raise KeyError(f"Cannot find dimension matching {candidates}. Available dims: {dims}")
    return None


def coord_index(nc, dim: str, value: str) -> int:
    if dim not in nc.variables:
        raise KeyError(f"Coordinate variable {dim!r} not found.")
    vals = values_to_str(nc.variables[dim][:])
    if value not in vals:
        raise ValueError(f"Value {value!r} not found on {dim!r}. Available: {vals}")
    return vals.index(value)


def build_selection_and_dims(
    dims: Sequence[str],
    cml_dim: str,
    cml_sel,
    time_dim: str,
    time_sel,
    selected_indices: Dict[str, int],
):
    sel = []
    out_dims = []
    for d in dims:
        if d == cml_dim:
            sel.append(cml_sel)
            out_dims.append(d)
        elif d == time_dim:
            sel.append(time_sel)
            out_dims.append(d)
        elif d in selected_indices:
            sel.append(int(selected_indices[d]))
            # Integer indexing removes this dimension.
        else:
            sel.append(slice(None))
            out_dims.append(d)
    return tuple(sel), tuple(out_dims)


def has_finite_nonzero_per_cml(arr, out_dims: Sequence[str], cml_dim: str, zero_eps: float) -> np.ndarray:
    arr = as_array(arr)
    cml_axis = list(out_dims).index(cml_dim)
    good = np.isfinite(arr) & (np.abs(arr) > float(zero_eps))
    axes = tuple(i for i in range(good.ndim) if i != cml_axis)
    if axes:
        return np.any(good, axis=axes)
    return good.astype(bool)


# ============================================================
# Rain-rate empty scan
# ============================================================
def scan_rainrate_has_signal(
    rain_nc: Path,
    rain_var: str,
    cml_dim: str,
    time_dim: str,
    baseline_method: str,
    kr_model: str,
    waa_method: str,
    filter_mode: str,
    chunk_cml: int,
    chunk_time: int,
    zero_eps: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return rain cml_ids and has_rain_signal boolean per rain cml_id."""
    import netCDF4

    with netCDF4.Dataset(rain_nc, "r") as nc:
        if cml_dim not in nc.dimensions:
            raise ValueError(f"{cml_dim!r} dimension not found in rain file.")
        if time_dim not in nc.dimensions:
            raise ValueError(f"{time_dim!r} dimension not found in rain file.")
        if rain_var not in nc.variables:
            raise ValueError(f"{rain_var!r} variable not found in rain file.")

        var = nc.variables[rain_var]
        dims = var.dimensions
        if cml_dim not in dims:
            raise ValueError(f"{rain_var!r} does not contain {cml_dim!r}.")
        if time_dim not in dims:
            raise ValueError(f"{rain_var!r} does not contain {time_dim!r}.")

        n_cml = len(nc.dimensions[cml_dim])
        n_time = len(nc.dimensions[time_dim])
        cml_ids = np.asarray(values_to_str(nc.variables[cml_dim][:]), dtype=str)

        selected_indices = {}
        selected_report = {}

        if filter_mode == "selected_combo":
            bdim = find_dim_from_candidates(dims, ["baseline_method", "baseline"], required=False)
            kdim = find_dim_from_candidates(dims, ["kr_model", "kr", "itu"], required=False)
            wdim = find_dim_from_candidates(dims, ["waa_method", "waa"], required=False)

            if bdim is not None:
                selected_indices[bdim] = coord_index(nc, bdim, baseline_method)
                selected_report[bdim] = baseline_method
            if kdim is not None:
                selected_indices[kdim] = coord_index(nc, kdim, kr_model)
                selected_report[kdim] = kr_model
            if wdim is not None:
                selected_indices[wdim] = coord_index(nc, wdim, waa_method)
                selected_report[wdim] = waa_method

        print("\nScanning estimated rain rate")
        print(f"  file         : {rain_nc}")
        print(f"  variable     : {rain_var}")
        print(f"  dims         : {dims}")
        print(f"  cml links    : {n_cml:,}")
        print(f"  time steps   : {n_time:,}")
        print(f"  filter mode  : {filter_mode}")
        if selected_report:
            for d, v in selected_report.items():
                print(f"  selected {d}: {v}")

        has_rain = np.zeros(n_cml, dtype=bool)

        for c0 in range(0, n_cml, chunk_cml):
            c1 = min(c0 + chunk_cml, n_cml)
            local_has = np.zeros(c1 - c0, dtype=bool)

            for t0 in range(0, n_time, chunk_time):
                t1 = min(t0 + chunk_time, n_time)

                sel, out_dims = build_selection_and_dims(
                    dims=dims,
                    cml_dim=cml_dim,
                    cml_sel=slice(c0, c1),
                    time_dim=time_dim,
                    time_sel=slice(t0, t1),
                    selected_indices=selected_indices,
                )

                arr = var[sel]
                local_has |= has_finite_nonzero_per_cml(arr, out_dims, cml_dim, zero_eps)

                if np.all(local_has):
                    break

            has_rain[c0:c1] = local_has
            print(
                f"  rain cml_id {c0:6d}:{c1:6d} | "
                f"has rain signal={int(local_has.sum()):5d} | "
                f"rain empty={int((~local_has).sum()):5d}"
            )

    return cml_ids, has_rain


# ============================================================
# Raw TRSL empty scan aligned to rain cml_ids
# ============================================================
def scan_trsl_has_signal_for_rain_ids(
    cml_nc: Path,
    rain_cml_ids: np.ndarray,
    cml_dim: str,
    raw_time_dim: str,
    trsl_var: str,
    chunk_cml: int,
    chunk_time: int,
    zero_eps: float,
) -> np.ndarray:
    """Return has_trsl_signal boolean aligned to rain_cml_ids."""
    import netCDF4

    rain_id_to_pos = {str(cid): i for i, cid in enumerate(rain_cml_ids)}
    has_trsl = np.zeros(len(rain_cml_ids), dtype=bool)
    found_in_raw = np.zeros(len(rain_cml_ids), dtype=bool)

    with netCDF4.Dataset(cml_nc, "r") as nc:
        if cml_dim not in nc.dimensions:
            raise ValueError(f"{cml_dim!r} dimension not found in raw CML file.")
        if raw_time_dim not in nc.dimensions:
            raise ValueError(f"{raw_time_dim!r} dimension not found in raw CML file.")
        if trsl_var not in nc.variables:
            raise ValueError(f"{trsl_var!r} variable not found in raw CML file.")

        raw_ids = np.asarray(values_to_str(nc.variables[cml_dim][:]), dtype=str)
        var = nc.variables[trsl_var]
        dims = var.dimensions

        if cml_dim not in dims:
            raise ValueError(f"{trsl_var!r} does not contain {cml_dim!r}.")
        if raw_time_dim not in dims:
            raise ValueError(f"{trsl_var!r} does not contain {raw_time_dim!r}.")

        n_raw_cml = len(nc.dimensions[cml_dim])
        n_time = len(nc.dimensions[raw_time_dim])

        print("\nScanning raw TRSL")
        print(f"  file       : {cml_nc}")
        print(f"  variable   : {trsl_var}")
        print(f"  dims       : {dims}")
        print(f"  raw links  : {n_raw_cml:,}")
        print(f"  time steps : {n_time:,}")

        for c0 in range(0, n_raw_cml, chunk_cml):
            c1 = min(c0 + chunk_cml, n_raw_cml)
            raw_chunk_ids = raw_ids[c0:c1]

            # Positions in rain file for raw links in this chunk.
            rain_pos = np.array(
                [rain_id_to_pos.get(str(cid), -1) for cid in raw_chunk_ids],
                dtype=np.int64,
            )
            match_mask = rain_pos >= 0
            if not np.any(match_mask):
                continue

            found_in_raw[rain_pos[match_mask]] = True
            local_has = np.zeros(c1 - c0, dtype=bool)

            for t0 in range(0, n_time, chunk_time):
                t1 = min(t0 + chunk_time, n_time)

                sel, out_dims = build_selection_and_dims(
                    dims=dims,
                    cml_dim=cml_dim,
                    cml_sel=slice(c0, c1),
                    time_dim=raw_time_dim,
                    time_sel=slice(t0, t1),
                    selected_indices={},
                )

                arr = var[sel]
                local_has |= has_finite_nonzero_per_cml(arr, out_dims, cml_dim, zero_eps)

                # Stop once all matched links in this chunk have raw signal.
                if np.all(local_has[match_mask]):
                    break

            matched_rain_pos = rain_pos[match_mask]
            has_trsl[matched_rain_pos] = local_has[match_mask]

            print(
                f"  raw cml_id {c0:6d}:{c1:6d} | "
                f"matched rain links={int(match_mask.sum()):5d} | "
                f"TRSL has signal={int(local_has[match_mask].sum()):5d} | "
                f"TRSL empty={int((~local_has[match_mask]).sum()):5d}"
            )

    missing_in_raw = ~found_in_raw
    if np.any(missing_in_raw):
        print("\nWARNING")
        print(f"  Rain file contains {int(missing_in_raw.sum()):,} cml_id values not found in raw CML file.")
        print("  These are treated as TRSL empty because raw CML evidence is missing.")
        has_trsl[missing_in_raw] = False

    return has_trsl


# ============================================================
# NetCDF subset writing
# ============================================================
def create_variable_like(src, dst, name: str, complevel: int, chunk_time: int, chunk_cml: int):
    src_var = src.variables[name]
    dims = src_var.dimensions
    fill_value = getattr(src_var, "_FillValue", None)

    is_numeric = np.issubdtype(src_var.dtype, np.number)
    kwargs = {}

    if is_numeric:
        kwargs["zlib"] = True
        kwargs["complevel"] = int(complevel)

        chunksizes = []
        for d in dims:
            n = len(dst.dimensions[d])
            dl = d.lower()
            if "time" in dl:
                chunksizes.append(min(n, chunk_time))
            elif "cml" in dl:
                chunksizes.append(min(n, chunk_cml))
            elif "sublink" in dl:
                chunksizes.append(min(n, 2))
            elif "baseline" in dl:
                chunksizes.append(min(n, 1))
            elif "kr" in dl:
                chunksizes.append(min(n, 2))
            elif "waa" in dl:
                chunksizes.append(min(n, 4))
            else:
                chunksizes.append(min(n, 1000))

        kwargs["chunksizes"] = tuple(chunksizes)

    if fill_value is not None:
        kwargs["fill_value"] = fill_value

    try:
        dst_var = dst.createVariable(name, src_var.datatype, dims, **kwargs)
    except TypeError:
        kwargs = {}
        if fill_value is not None:
            kwargs["fill_value"] = fill_value
        dst_var = dst.createVariable(name, src_var.datatype, dims, **kwargs)

    for attr in src_var.ncattrs():
        if attr == "_FillValue":
            continue
        try:
            dst_var.setncattr(attr, src_var.getncattr(attr))
        except Exception:
            dst_var.setncattr(attr, str(src_var.getncattr(attr)))

    return dst_var


def copy_variable_subset(src, dst, name: str, keep_idx: np.ndarray,
                         cml_dim: str, time_dim: str,
                         chunk_cml: int, chunk_time: int):
    src_var = src.variables[name]
    dst_var = dst.variables[name]
    dims = src_var.dimensions

    if cml_dim not in dims:
        dst_var[:] = src_var[:]
        return

    has_time = time_dim in dims
    n_keep = len(keep_idx)
    n_time = len(src.dimensions[time_dim]) if has_time else 1

    for ko0 in range(0, n_keep, chunk_cml):
        ko1 = min(ko0 + chunk_cml, n_keep)
        src_cml_idx = keep_idx[ko0:ko1]

        if has_time:
            for t0 in range(0, n_time, chunk_time):
                t1 = min(t0 + chunk_time, n_time)

                src_sel = []
                dst_sel = []
                for d in dims:
                    if d == cml_dim:
                        src_sel.append(src_cml_idx)
                        dst_sel.append(slice(ko0, ko1))
                    elif d == time_dim:
                        src_sel.append(slice(t0, t1))
                        dst_sel.append(slice(t0, t1))
                    else:
                        src_sel.append(slice(None))
                        dst_sel.append(slice(None))

                dst_var[tuple(dst_sel)] = src_var[tuple(src_sel)]
        else:
            src_sel = []
            dst_sel = []
            for d in dims:
                if d == cml_dim:
                    src_sel.append(src_cml_idx)
                    dst_sel.append(slice(ko0, ko1))
                else:
                    src_sel.append(slice(None))
                    dst_sel.append(slice(None))
            dst_var[tuple(dst_sel)] = src_var[tuple(src_sel)]


def write_filtered_rain_netcdf(
    in_nc: Path,
    out_nc: Path,
    keep_idx: np.ndarray,
    cml_dim: str,
    time_dim: str,
    chunk_cml: int,
    chunk_time: int,
    complevel: int,
    attrs: Dict[str, object],
):
    import netCDF4

    t0 = time.perf_counter()

    with netCDF4.Dataset(in_nc, "r") as src:
        with netCDF4.Dataset(out_nc, "w", format=getattr(src, "file_format", "NETCDF4")) as dst:
            for dname, dim in src.dimensions.items():
                if dname == cml_dim:
                    dst.createDimension(dname, len(keep_idx))
                else:
                    dst.createDimension(dname, None if dim.isunlimited() else len(dim))

            for attr in src.ncattrs():
                try:
                    dst.setncattr(attr, src.getncattr(attr))
                except Exception:
                    dst.setncattr(attr, str(src.getncattr(attr)))

            for k, v in attrs.items():
                dst.setncattr(k, v)

            for name in src.variables:
                create_variable_like(src, dst, name, complevel, chunk_time, chunk_cml)

            for name in src.variables:
                print(f"  copying variable: {name}")
                copy_variable_subset(src, dst, name, keep_idx, cml_dim, time_dim, chunk_cml, chunk_time)

    print(f"Copy elapsed: {time.perf_counter() - t0:.2f} s")


# ============================================================
# Plot command
# ============================================================
def run_plot_script(args, filtered_nc: Path):
    plot_script = Path(args.plot_script)
    if not plot_script.exists():
        raise FileNotFoundError(f"Plot script not found: {plot_script}")

    cmd = [
        sys.executable,
        str(plot_script),
        "--rain-nc", str(filtered_nc),
        "--out-dir", str(args.plot_out_dir),
        "--reference", args.reference,
        "--cml-var", args.rain_var,
        "--baseline-method", args.baseline_method,
        "--kr-model", args.kr_model,
        "--waa-method", args.waa_method,
        "--max-rain", str(args.max_rain),
        "--chunk-cml", str(args.plot_chunk_cml),
        "--gridsize", str(args.gridsize),
    ]

    if args.reference == "gauge":
        cmd.extend([
            "--gauge-nc", args.gauge_nc,
            "--cml-coords-nc", args.cml_coords_nc,
            "--max-gauge-dist-km", str(args.max_gauge_dist_km),
        ])

    if args.wet_only:
        cmd.append("--wet-only")

    print("\nRunning plot script:")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)


# ============================================================
# Main
# ============================================================
def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--rain-nc", default="cml_rainrate_pycml_constant_kr_waa_methods.nc")
    ap.add_argument("--cml-nc", default="cml_preprocessed_v3.nc",
                    help="Raw/preprocessed CML file containing TRSL.")
    ap.add_argument("--out-nc", default="cml_rainrate_trsl_checked_no_empty_links.nc")

    ap.add_argument("--rain-var", default="R_mean_5min",
                    help="R_mean_5min or R_sublink_5min.")
    ap.add_argument("--trsl-var", default="trsl")

    ap.add_argument("--rain-cml-dim", default="cml_id")
    ap.add_argument("--rain-time-dim", default="time_5min")
    ap.add_argument("--raw-cml-dim", default="cml_id")
    ap.add_argument("--raw-time-dim", default="time")

    ap.add_argument("--baseline-method", default="pycml_constant")
    ap.add_argument("--kr-model", default="ITU_2003")
    ap.add_argument("--waa-method", default="pastorek_2021")

    ap.add_argument("--filter-mode", choices=["selected_combo", "all_combos"], default="selected_combo")
    ap.add_argument("--drop-condition", choices=["both_empty", "trsl_empty", "rain_empty"], default="both_empty",
                    help="Default both_empty: drop only if rain estimate is empty AND raw TRSL is empty.")
    ap.add_argument("--zero-eps", type=float, default=0.0)

    ap.add_argument("--chunk-cml", type=int, default=128)
    ap.add_argument("--chunk-time-rain", type=int, default=10080)
    ap.add_argument("--chunk-time-trsl", type=int, default=10080)
    ap.add_argument("--complevel", type=int, default=4)

    ap.add_argument("--count-only", action="store_true")
    ap.add_argument("--overwrite", action="store_true")

    ap.add_argument("--run-plot", action="store_true")
    ap.add_argument("--plot-script", default="cml_compare_single_combo_pycml_itu03_pastorek_hexbin_oneplot_reference (1).py")
    ap.add_argument("--plot-out-dir", default="single_combo_trsl_checked_plot")
    ap.add_argument("--reference", choices=["radar", "gauge"], default="radar")
    ap.add_argument("--max-rain", type=float, default=100.0)
    ap.add_argument("--wet-only", action="store_true")
    ap.add_argument("--plot-chunk-cml", type=int, default=64)
    ap.add_argument("--gridsize", type=int, default=50)

    # Gauge passthrough options for plotting.
    ap.add_argument("--gauge-nc", default="CHMI_10MIN_v1.nc")
    ap.add_argument("--cml-coords-nc", default="cml_daily_nc_v2_merged.nc")
    ap.add_argument("--max-gauge-dist-km", type=float, default=1.0)

    args = ap.parse_args()

    rain_nc = Path(args.rain_nc)
    cml_nc = Path(args.cml_nc)
    out_nc = Path(args.out_nc)

    if not rain_nc.exists():
        raise FileNotFoundError(f"Rain NetCDF not found: {rain_nc}")
    if not cml_nc.exists():
        raise FileNotFoundError(f"Raw CML NetCDF not found: {cml_nc}")

    if out_nc.exists() and not args.overwrite and not args.count_only:
        raise FileExistsError(f"Output exists: {out_nc}. Use --overwrite.")

    if out_nc.exists() and args.overwrite and not args.count_only:
        out_nc.unlink()

    out_dir = out_nc.parent if str(out_nc.parent) != "." else Path(".")
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()

    rain_cml_ids, has_rain = scan_rainrate_has_signal(
        rain_nc=rain_nc,
        rain_var=args.rain_var,
        cml_dim=args.rain_cml_dim,
        time_dim=args.rain_time_dim,
        baseline_method=args.baseline_method,
        kr_model=args.kr_model,
        waa_method=args.waa_method,
        filter_mode=args.filter_mode,
        chunk_cml=args.chunk_cml,
        chunk_time=args.chunk_time_rain,
        zero_eps=args.zero_eps,
    )

    has_trsl = scan_trsl_has_signal_for_rain_ids(
        cml_nc=cml_nc,
        rain_cml_ids=rain_cml_ids,
        cml_dim=args.raw_cml_dim,
        raw_time_dim=args.raw_time_dim,
        trsl_var=args.trsl_var,
        chunk_cml=args.chunk_cml,
        chunk_time=args.chunk_time_trsl,
        zero_eps=args.zero_eps,
    )

    rain_empty = ~has_rain
    trsl_empty = ~has_trsl

    if args.drop_condition == "both_empty":
        drop = rain_empty & trsl_empty
    elif args.drop_condition == "trsl_empty":
        drop = trsl_empty
    else:
        drop = rain_empty

    keep = ~drop
    keep_idx = np.flatnonzero(keep)
    drop_idx = np.flatnonzero(drop)

    # Diagnostic categories.
    both_empty = rain_empty & trsl_empty
    rain_empty_trsl_ok = rain_empty & has_trsl
    trsl_empty_rain_ok = trsl_empty & has_rain
    both_ok = has_rain & has_trsl

    print("\nDiagnostic table")
    print(f"  Total rain-file links                     : {len(rain_cml_ids):,}")
    print(f"  Rain estimate empty                       : {int(rain_empty.sum()):,}")
    print(f"  Raw TRSL empty                            : {int(trsl_empty.sum()):,}")
    print(f"  BOTH rain empty AND TRSL empty             : {int(both_empty.sum()):,}")
    print(f"  Rain empty BUT TRSL exists -> NOT dropped by default : {int(rain_empty_trsl_ok.sum()):,}")
    print(f"  TRSL empty BUT rain exists -> WARNING category       : {int(trsl_empty_rain_ok.sum()):,}")
    print(f"  Both rain and TRSL have signal             : {int(both_ok.sum()):,}")

    print("\nFinal drop rule")
    print(f"  drop-condition : {args.drop_condition}")
    print(f"  Kept links     : {len(keep_idx):,}")
    print(f"  Dropped links  : {len(drop_idx):,}")

    stem = out_nc.stem
    report_csv = out_dir / f"{stem}_rain_trsl_empty_report.csv"
    kept_csv = out_dir / f"{stem}_kept_links.csv"
    dropped_csv = out_dir / f"{stem}_dropped_links.csv"

    report = pd.DataFrame({
        "cml_index_in_rain_file": np.arange(len(rain_cml_ids), dtype=int),
        "cml_id": rain_cml_ids,
        "has_rain_estimate_signal": has_rain,
        "rain_estimate_empty": rain_empty,
        "has_raw_trsl_signal": has_trsl,
        "raw_trsl_empty": trsl_empty,
        "both_empty": both_empty,
        "rain_empty_but_trsl_exists": rain_empty_trsl_ok,
        "trsl_empty_but_rain_exists": trsl_empty_rain_ok,
        "drop": drop,
        "drop_condition": args.drop_condition,
        "baseline_method": args.baseline_method,
        "kr_model": args.kr_model,
        "waa_method": args.waa_method,
    })
    report.to_csv(report_csv, index=False)
    report.loc[keep, ["cml_index_in_rain_file", "cml_id"]].to_csv(kept_csv, index=False)
    report.loc[drop, ["cml_index_in_rain_file", "cml_id", "drop_condition", "both_empty", "rain_estimate_empty", "raw_trsl_empty"]].to_csv(dropped_csv, index=False)

    print("\nSaved CSV")
    print(f"  {report_csv}")
    print(f"  {kept_csv}")
    print(f"  {dropped_csv}")

    if args.count_only:
        print("\nCount-only mode: no filtered NetCDF written.")
        print(f"Total elapsed: {time.perf_counter() - t0:.2f} s")
        return

    print("\nWriting filtered rain-rate NetCDF")
    print(f"  output: {out_nc}")

    attrs = {
        "rain_trsl_empty_filter": "Links filtered using estimated rain-rate emptiness and raw TRSL emptiness.",
        "rain_trsl_empty_filter_drop_condition": args.drop_condition,
        "rain_trsl_empty_filter_default_meaning": "both_empty means drop only if selected rain estimate is empty AND raw TRSL is empty.",
        "rain_trsl_empty_filter_rain_file": str(rain_nc),
        "rain_trsl_empty_filter_raw_cml_file": str(cml_nc),
        "rain_trsl_empty_filter_rain_variable": args.rain_var,
        "rain_trsl_empty_filter_trsl_variable": args.trsl_var,
        "rain_trsl_empty_filter_baseline_method": args.baseline_method,
        "rain_trsl_empty_filter_kr_model": args.kr_model,
        "rain_trsl_empty_filter_waa_method": args.waa_method,
        "rain_trsl_empty_filter_kept_links": int(len(keep_idx)),
        "rain_trsl_empty_filter_dropped_links": int(len(drop_idx)),
        "rain_trsl_empty_filter_report_csv": str(report_csv),
    }

    write_filtered_rain_netcdf(
        in_nc=rain_nc,
        out_nc=out_nc,
        keep_idx=keep_idx,
        cml_dim=args.rain_cml_dim,
        time_dim=args.rain_time_dim,
        chunk_cml=args.chunk_cml,
        chunk_time=args.chunk_time_rain,
        complevel=args.complevel,
        attrs=attrs,
    )

    if args.run_plot:
        run_plot_script(args, out_nc)

    print("\nDone")
    print(f"  filtered file: {out_nc}")
    print(f"Total elapsed: {time.perf_counter() - t0:.2f} s")


if __name__ == "__main__":
    main()
