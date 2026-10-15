#!/usr/bin/env python3
"""
Single-combination CML-vs-reference rain-rate comparison
========================================================

Makes ONE high-density hexbin plot for one selected rainfall-estimation
combination, using only fair matched wet-event pairs.

The final evaluation/plot mask is always built from the same paired samples:
    * CML estimate is finite,
    * selected comparison reference (radar or gauge) is finite,
    * both plotted values are inside --min-rain/--max-rain when
      --mask-mode range is used,
    * radar at the same CML/time position is strictly greater than
      --radar-wet-thr (default: radar > 0.1 mm h-1).

For a gauge comparison, the gauge remains the plotted/evaluated reference;
radar is used only to identify rainy instances. This ensures that metrics and
the hexbin use exactly the same fair matched sample pairs.

Default selected combo:
    baseline_method = pycml_constant
    kr_model        = ITU_2003
    waa_method      = pastorek_2021

Default plot density:
    gridsize = 50

Supported references
--------------------
--reference radar
    Compare CML against the file variable R_radar (or the variable provided by
    --radar-var).

--reference gauge
    Build a CML-aligned gauge reference field R_gauge(cml_id, time_5min):
      * compute each link midpoint from site_a/site_b coordinates,
      * assign the nearest gauge within --max-gauge-dist-km,
      * convert gauge accumulation (default 10-min amount) to mm h-1,
      * forward-fill each gauge value only within its own accumulation window,
      * compare the selected CML estimate against this aligned gauge reference.

Outputs
-------
For --reference radar:
    <out-dir>/hexbin_<combo>_vs_radar.png
    <out-dir>/hexbin_<combo>_vs_radar.pdf
    <out-dir>/global_metrics_<combo>_vs_radar.csv
    <out-dir>/per_sensor_metrics_<combo>_vs_radar.csv

For --reference gauge:
    <out-dir>/hexbin_<combo>_vs_gauge.png
    <out-dir>/hexbin_<combo>_vs_gauge.pdf
    <out-dir>/global_metrics_<combo>_vs_gauge.csv
    <out-dir>/per_sensor_metrics_<combo>_vs_gauge.csv
    <out-dir>/gauge_link_mapping.csv
"""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, BoundaryNorm, ListedColormap
from matplotlib.ticker import AutoMinorLocator

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    def tqdm(x, **kwargs):
        return x


# ============================================================
# Style
# ============================================================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

REFERENCE_HEXBIN_COLORS = [
    "#000080",
    "#0015d9",
    "#0057ff",
    "#00a5ff",
    "#00d5d5",
    "#36ef9a",
    "#8cff66",
    "#d9ff32",
    "#ffff00",
    "#ffd000",
    "#ff9a00",
    "#ff5a00",
    "#ff1a00",
    "#d00000",
    "#7f0000",
]
REFERENCE_CMAP = ListedColormap(REFERENCE_HEXBIN_COLORS, name="reference_blue_yellow_red")


# ============================================================
# General helpers
# ============================================================
def get_count_cmap(name: str):
    if str(name).lower() in {"reference", "same", "shown", "figure"}:
        return REFERENCE_CMAP
    return plt.get_cmap(name)


def find_dim_by_candidates(obj, candidates: List[str], required: bool = True) -> Optional[str]:
    dims = list(obj.dims)
    for cand in candidates:
        for d in dims:
            if d == cand:
                return d
    for d in dims:
        dl = d.lower()
        if any(c.lower() in dl for c in candidates):
            return d
    if required:
        raise KeyError(f"Cannot find dimension matching {candidates}. Dims: {dims}")
    return None


def coord_values_as_str(ds: xr.Dataset, dim: str) -> List[str]:
    if dim not in ds.coords:
        return [str(i) for i in range(ds.sizes[dim])]
    return [str(v) for v in ds[dim].values]


def resolve_coord_index(ds: xr.Dataset, dim: str, value: str) -> int:
    vals = coord_values_as_str(ds, dim)
    if value not in vals:
        raise ValueError(
            f"Requested value {value!r} not found on dimension {dim!r}. "
            f"Available values: {vals}"
        )
    return vals.index(value)


def nice_step(vmin: float, vmax: float, target_n: int = 5) -> float:
    span = vmax - vmin
    if span <= 0:
        return 1.0
    raw = span / target_n
    mag = 10 ** np.floor(np.log10(raw))
    norm = raw / mag
    if norm <= 1:
        return 1 * mag
    if norm <= 2:
        return 2 * mag
    if norm <= 5:
        return 5 * mag
    return 10 * mag


def make_nice_ticks(vmin: float, vmax: float, target_n: int = 5) -> np.ndarray:
    step = nice_step(vmin, vmax, target_n)
    tick_min = np.floor(vmin / step) * step
    tick_max = np.ceil(vmax / step) * step
    ticks = np.arange(tick_min, tick_max + 0.5 * step, step)
    return np.unique(np.round(ticks, 10))


def format_tick_labels(ticks: np.ndarray) -> List[str]:
    labels = []
    for t in ticks:
        if np.isclose(t, round(t)):
            labels.append(f"{int(round(t))}")
        else:
            labels.append(f"{t:g}")
    return labels


def style_axis(ax, lim_max: float, xlabel: str, ylabel: str):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_position(("outward", 6))
    ax.spines["left"].set_position(("outward", 6))

    ax.set_xlim(0.0, lim_max)
    ax.set_ylim(0.0, lim_max)

    ticks = make_nice_ticks(0.0, lim_max, target_n=5)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(format_tick_labels(ticks), fontsize=10)
    ax.set_yticklabels(format_tick_labels(ticks), fontsize=10)

    ax.spines["bottom"].set_bounds(float(ticks[0]), float(ticks[-1]))
    ax.spines["left"].set_bounds(float(ticks[0]), float(ticks[-1]))

    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.tick_params(axis="both", which="both", direction="out", top=False, right=False, bottom=True, left=True)
    ax.tick_params(axis="x", which="major", length=6, width=0.8, labelsize=10, pad=6)
    ax.tick_params(axis="y", which="major", length=6, width=0.8, labelsize=10, pad=6)
    ax.tick_params(axis="both", which="minor", length=3, width=0.5)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.grid(False)


def valid_range(a: np.ndarray, min_v: float, max_v: float) -> np.ndarray:
    return np.isfinite(a) & (a >= min_v) & (a <= max_v)


def parse_time_window(ds: xr.Dataset, time_dim: str, args) -> Tuple[slice, np.ndarray]:
    ntime = ds.sizes[time_dim]
    if time_dim in ds.coords:
        times = pd.DatetimeIndex(ds[time_dim].values)
    else:
        return slice(0, ntime), np.arange(ntime)

    start = 0
    end = ntime
    if args.tstart:
        start = int(np.searchsorted(times.values, np.datetime64(pd.Timestamp(args.tstart)), "left"))
    if args.tend:
        end = int(np.searchsorted(times.values, np.datetime64(pd.Timestamp(args.tend)), "right"))

    start = max(0, min(start, ntime))
    end = max(start, min(end, ntime))
    if end <= start:
        raise ValueError(
            f"Empty selected time window: --tstart={args.tstart!r}, --tend={args.tend!r}. "
            f"Available time range is {times[0]} -> {times[-1]}"
        )

    print(f"  Time window indices : {start}:{end} / {ntime}")
    print(f"  Time window used    : {times[start]} -> {times[end - 1]}")
    return slice(start, end), np.arange(start, end)


def make_log_count_boundaries(vmax: float, n_levels: int) -> np.ndarray:
    vmax = max(1.0, float(vmax))
    n_levels = max(4, int(n_levels))
    raw = np.geomspace(1.0, vmax, n_levels + 1)
    boundaries = np.unique(np.rint(raw).astype(int)).astype(float)
    if boundaries[0] > 1:
        boundaries = np.r_[1.0, boundaries]
    if boundaries[-1] < vmax:
        boundaries = np.r_[boundaries, vmax]
    if boundaries.size < 2:
        boundaries = np.array([1.0, max(2.0, vmax)], dtype=float)
    return boundaries


def log_count_ticks(vmax: float) -> List[float]:
    if not np.isfinite(vmax) or vmax <= 1:
        return [1.0]
    p1 = int(np.ceil(np.log10(max(1.0, float(vmax)))))
    return [10.0 ** p for p in range(0, p1 + 1)]


def format_count_tick(x: float, mode: str = "power") -> str:
    if mode == "power":
        if x <= 0:
            return ""
        p = int(round(np.log10(x)))
        return rf"$10^{{{p}}}$"
    if x >= 1000:
        return f"{int(round(x)):,}"
    return f"{int(round(x))}"


# ============================================================
# Statistics helpers
# ============================================================
def empty_stats() -> Dict[str, float]:
    return {
        "N": 0.0,
        "sum_c": 0.0,
        "sum_r": 0.0,
        "sum_c2": 0.0,
        "sum_r2": 0.0,
        "sum_cr": 0.0,
        "sum_diff": 0.0,
        "sum_abs_diff": 0.0,
        "sum_diff2": 0.0,
    }


def update_stats(stats: Dict[str, float], c: np.ndarray, r: np.ndarray, mask: np.ndarray):
    n = int(np.sum(mask))
    if n == 0:
        return

    cv = c[mask].astype("float64", copy=False)
    rv = r[mask].astype("float64", copy=False)
    d = cv - rv

    stats["N"] += n
    stats["sum_c"] += float(np.sum(cv))
    stats["sum_r"] += float(np.sum(rv))
    stats["sum_c2"] += float(np.sum(cv * cv))
    stats["sum_r2"] += float(np.sum(rv * rv))
    stats["sum_cr"] += float(np.sum(cv * rv))
    stats["sum_diff"] += float(np.sum(d))
    stats["sum_abs_diff"] += float(np.sum(np.abs(d)))
    stats["sum_diff2"] += float(np.sum(d * d))


def finalize_stats(stats: Dict[str, float]) -> Dict[str, float]:
    n = int(stats["N"])
    if n == 0:
        return {
            "N": 0,
            "MBE": np.nan,
            "MAE": np.nan,
            "RMSE": np.nan,
            "PCC": np.nan,
            "mean_ref": np.nan,
            "mean_cml": np.nan,
        }

    mean_c = stats["sum_c"] / n
    mean_r = stats["sum_r"] / n
    mbe = stats["sum_diff"] / n
    mae = stats["sum_abs_diff"] / n
    rmse = math.sqrt(max(stats["sum_diff2"] / n, 0.0))

    cov = stats["sum_cr"] / n - mean_c * mean_r
    var_c = stats["sum_c2"] / n - mean_c * mean_c
    var_r = stats["sum_r2"] / n - mean_r * mean_r
    if n >= 2 and var_c > 0 and var_r > 0:
        pcc = cov / math.sqrt(var_c * var_r)
    else:
        pcc = np.nan

    return {
        "N": n,
        "MBE": mbe,
        "MAE": mae,
        "RMSE": rmse,
        "PCC": pcc,
        "mean_ref": mean_r,
        "mean_cml": mean_c,
    }


def row_metrics(c2: np.ndarray, r2: np.ndarray, m2: np.ndarray) -> Tuple[np.ndarray, ...]:
    c = c2.astype("float64", copy=False)
    r = r2.astype("float64", copy=False)
    m = m2.astype(bool, copy=False)

    n = m.sum(axis=1).astype("float64")
    n_safe = np.where(n > 0, n, 1.0)

    c0 = np.where(m, c, 0.0)
    r0 = np.where(m, r, 0.0)
    d0 = c0 - r0

    s_c = c0.sum(axis=1)
    s_r = r0.sum(axis=1)
    s_c2 = (c0 * c0).sum(axis=1)
    s_r2 = (r0 * r0).sum(axis=1)
    s_cr = (c0 * r0).sum(axis=1)
    s_diff = d0.sum(axis=1)
    s_abs = np.abs(np.where(m, c - r, 0.0)).sum(axis=1)
    s_diff2 = (d0 * d0).sum(axis=1)

    mean_c = s_c / n_safe
    mean_r = s_r / n_safe
    mbe = s_diff / n_safe
    mae = s_abs / n_safe
    rmse = np.sqrt(np.maximum(s_diff2 / n_safe, 0.0))

    cov = s_cr / n_safe - mean_c * mean_r
    var_c = s_c2 / n_safe - mean_c * mean_c
    var_r = s_r2 / n_safe - mean_r * mean_r
    denom = np.sqrt(np.maximum(var_c, 0.0) * np.maximum(var_r, 0.0))
    ok_pcc = (n >= 2) & (var_c > 0) & (var_r > 0) & (denom > 0)
    pcc = np.full_like(denom, np.nan, dtype="float64")
    np.divide(cov, denom, out=pcc, where=ok_pcc)

    ok = n > 0
    return (
        n.astype("int64"),
        np.where(ok, mbe, np.nan),
        np.where(ok, mae, np.nan),
        np.where(ok, rmse, np.nan),
        pcc,
        np.where(ok, mean_r, np.nan),
        np.where(ok, mean_c, np.nan),
    )


# ============================================================
# CML selection
# ============================================================
def select_cml_combo(ds: xr.Dataset, cml_var: str, args) -> xr.DataArray:
    da = ds[cml_var]
    sel = {}
    for dim_candidates, value, label in [
        (["baseline_method", "baseline"], args.baseline_method, "baseline method"),
        (["kr_model", "kr", "itu"], args.kr_model, "k-R model"),
        (["waa_method", "waa"], args.waa_method, "WAA method"),
    ]:
        dim = find_dim_by_candidates(da, dim_candidates, required=False)
        if dim is not None:
            idx = resolve_coord_index(ds, dim, value)
            sel[dim] = idx
            print(f"  Selected {label:<15}: {value} on dim {dim!r}")
    if not sel:
        print("  WARNING: no baseline/kr/waa dimensions found. Using variable as-is.")
    return da.isel(sel) if sel else da


# ============================================================
# Gauge reference helpers
# ============================================================
def _series_from_time_da(da: xr.DataArray, time_name: str = "time") -> pd.Series:
    return pd.Series(np.asarray(da.values), index=pd.DatetimeIndex(da[time_name].values))


def haversine_km_vec(lon1, lat1, lon2, lat2):
    R = 6371.0
    lon1, lat1 = np.deg2rad(lon1), np.deg2rad(lat1)
    lon2 = np.deg2rad(np.asarray(lon2, dtype=float))
    lat2 = np.deg2rad(np.asarray(lat2, dtype=float))
    dlon, dlat = lon2 - lon1, lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def mm_amount_to_mmph(values, minutes):
    return values * (60.0 / float(minutes))


def _find_non_time_dim(da: xr.DataArray) -> str:
    for d in da.dims:
        if "time" not in d.lower():
            return d
    raise KeyError(f"Could not find a station dimension in gauge var dims {da.dims}")


def _link_midpoints_in_order(coords_ds: xr.Dataset, cml_ids: np.ndarray, args) -> Tuple[np.ndarray, np.ndarray]:
    a_lat = coords_ds[args.site_a_lat_var]
    cdim = find_dim_by_candidates(a_lat, ["cml_id", "cml"])
    src_ids = np.asarray(coords_ds[cdim].values, dtype=str)

    la = np.asarray(coords_ds[args.site_a_lat_var].values, dtype=float)
    lb = np.asarray(coords_ds[args.site_b_lat_var].values, dtype=float)
    loa = np.asarray(coords_ds[args.site_a_lon_var].values, dtype=float)
    lob = np.asarray(coords_ds[args.site_b_lon_var].values, dtype=float)

    mid_lon = 0.5 * (loa + lob)
    mid_lat = 0.5 * (la + lb)

    s_lon = pd.Series(mid_lon, index=src_ids).reindex(cml_ids)
    s_lat = pd.Series(mid_lat, index=src_ids).reindex(cml_ids)
    return s_lon.to_numpy(dtype=float), s_lat.to_numpy(dtype=float)


def build_gauge_reference(ds: xr.Dataset, args, dims: Dict[str, Optional[str]]) -> Tuple[xr.DataArray, pd.DataFrame]:
    cml_dim = dims["cml_dim"]
    time_dim = dims["time_dim"]
    cml_ids = np.asarray(ds[cml_dim].values, dtype=str)

    if args.site_a_lat_var in ds.variables:
        coords_ds = ds
        close_coords = False
        print("  Gauge: using CML site coordinates from the rain file.")
    else:
        if not args.cml_coords_nc:
            raise ValueError(
                f"CML site coordinate variable {args.site_a_lat_var!r} is not in the rain file. "
                "Pass --cml-coords-nc pointing to the CML metadata NetCDF."
            )
        coords_ds = xr.open_dataset(args.cml_coords_nc)
        if cml_dim in coords_ds.coords:
            coords_ds = coords_ds.assign_coords({cml_dim: coords_ds[cml_dim].astype(str)})
        close_coords = True
        print(f"  Gauge: using CML site coordinates from {args.cml_coords_nc}.")

    try:
        mid_lon, mid_lat = _link_midpoints_in_order(coords_ds, cml_ids, args)
    finally:
        if close_coords:
            coords_ds.close()

    ds_g = xr.open_dataset(args.gauge_nc)
    try:
        g_rr = ds_g[args.gauge_rr_var]
        g_dim = _find_non_time_dim(g_rr)
        g_time_name = find_dim_by_candidates(g_rr, ["time"], required=True)

        g_lon_all = np.asarray(ds_g[args.gauge_lon_var].values, dtype=float)
        g_lat_all = np.asarray(ds_g[args.gauge_lat_var].values, dtype=float)
        g_ok = np.isfinite(g_lon_all) & np.isfinite(g_lat_all)
        g_pos = np.where(g_ok)[0]
        g_lon, g_lat = g_lon_all[g_ok], g_lat_all[g_ok]
        if g_pos.size == 0:
            raise ValueError("No gauges with finite coordinates found in the gauge file.")

        if g_dim in ds_g.coords:
            g_station_ids = np.asarray(ds_g[g_dim].values).astype(str)
        else:
            g_station_ids = np.arange(g_lon_all.size).astype(str)

        cml_time = pd.DatetimeIndex(ds[time_dim].values)
        if cml_time.size >= 2:
            dt_cml_min = float(np.median(np.diff(cml_time.values).astype("timedelta64[s]").astype(float)) / 60.0)
        else:
            dt_cml_min = float(args.gauge_minutes)
        dt_cml_min = dt_cml_min if dt_cml_min > 0 else float(args.gauge_minutes)
        ffill_limit = max(0, int(round(float(args.gauge_minutes) / dt_cml_min)) - 1)

        n_cml = cml_ids.size
        n_time = cml_time.size
        R_gauge = np.full((n_cml, n_time), np.nan, dtype="float32")

        assigned_gauge_pos = np.full(n_cml, -1, dtype="int64")
        assigned_dist = np.full(n_cml, np.nan, dtype="float64")
        finite_mid = np.isfinite(mid_lon) & np.isfinite(mid_lat)
        for i in range(n_cml):
            if not finite_mid[i]:
                continue
            d = haversine_km_vec(mid_lon[i], mid_lat[i], g_lon, g_lat)
            if not np.isfinite(d).any():
                continue
            j = int(np.nanargmin(d))
            if d[j] <= args.max_gauge_dist_km:
                assigned_gauge_pos[i] = int(g_pos[j])
                assigned_dist[i] = float(d[j])

        used = np.unique(assigned_gauge_pos[assigned_gauge_pos >= 0])
        print(f"  Gauge: {len(used)} unique gauges assigned to "
              f"{int(np.sum(assigned_gauge_pos >= 0))}/{n_cml} links "
              f"(<= {args.max_gauge_dist_km:g} km).")

        gauge_on_cml_cache: Dict[int, np.ndarray] = {}
        for gpos in used:
            g_amt = _series_from_time_da(g_rr.isel({g_dim: int(gpos)}), time_name=g_time_name)
            g_mmph = mm_amount_to_mmph(g_amt, args.gauge_minutes)
            g_aligned = g_mmph.reindex(cml_time, method="ffill", limit=ffill_limit)
            gauge_on_cml_cache[int(gpos)] = g_aligned.to_numpy(dtype="float32")

        for i in range(n_cml):
            gpos = assigned_gauge_pos[i]
            if gpos >= 0:
                R_gauge[i, :] = gauge_on_cml_cache[int(gpos)]

        mapping = pd.DataFrame({
            "cml_id": cml_ids,
            "gauge_pos": assigned_gauge_pos,
            "gauge_id": [g_station_ids[p] if p >= 0 else "" for p in assigned_gauge_pos],
            "gauge_dist_km": assigned_dist,
            "mid_lon": mid_lon,
            "mid_lat": mid_lat,
        })
    finally:
        ds_g.close()

    R_gauge_da = xr.DataArray(
        R_gauge,
        dims=[cml_dim, time_dim],
        coords={cml_dim: ds[cml_dim].values, time_dim: ds[time_dim].values},
        name="R_gauge_ref",
    )
    return R_gauge_da, mapping


# ============================================================
# Data access and masks
# ============================================================
def get_block(
    cml_da_sel: xr.DataArray,
    ref_da: xr.DataArray,
    radar_da: xr.DataArray,
    cml_slice: slice,
    time_slice: slice,
    dims: Dict[str, Optional[str]],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read matched CML, comparison-reference, and radar wet-mask blocks."""
    cml_isel = {dims["cml_dim"]: cml_slice, dims["time_dim"]: time_slice}
    cml_sel = cml_da_sel.isel(cml_isel)

    if dims.get("sublink_dim") is not None and dims["sublink_dim"] in cml_sel.dims:
        cml_arr = cml_sel.transpose(
            dims["cml_dim"], dims["sublink_dim"], dims["time_dim"]
        ).values.astype("float32", copy=False)
    else:
        cml_arr = cml_sel.transpose(
            dims["cml_dim"], dims["time_dim"]
        ).values.astype("float32", copy=False)

    ref_arr = ref_da.isel({
        dims["ref_cml_dim"]: cml_slice,
        dims["ref_time_dim"]: time_slice,
    }).transpose(
        dims["ref_cml_dim"], dims["ref_time_dim"]
    ).values.astype("float32", copy=False)

    radar_arr = radar_da.isel({
        dims["radar_cml_dim"]: cml_slice,
        dims["radar_time_dim"]: time_slice,
    }).transpose(
        dims["radar_cml_dim"], dims["radar_time_dim"]
    ).values.astype("float32", copy=False)

    return cml_arr, ref_arr, radar_arr


def make_masks(
    cml_block: np.ndarray,
    ref_block: np.ndarray,
    radar_block: np.ndarray,
    min_rain: float,
    max_rain: float,
    wet_only: bool,
    radar_wet_thr: float,
    mask_mode: str = "range",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return broadcast references and the exact fair-pair evaluation mask.

    Wet means radar > radar_wet_thr (strict inequality). The plotted reference
    can independently be radar or gauge.
    """
    if cml_block.ndim == 2:
        ref_like = ref_block
        radar_like = radar_block
    elif cml_block.ndim == 3:
        ref_like = np.broadcast_to(ref_block[:, None, :], cml_block.shape)
        radar_like = np.broadcast_to(radar_block[:, None, :], cml_block.shape)
    else:
        raise ValueError(f"Unexpected cml_block ndim: {cml_block.ndim}")

    # Fair comparison: the selected reference and CML estimate must both exist.
    finite_pair = np.isfinite(ref_like) & np.isfinite(cml_block)

    if mask_mode == "finite":
        mask = finite_pair
    elif mask_mode == "range":
        mask = (
            finite_pair
            & valid_range(ref_like, min_rain, max_rain)
            & valid_range(cml_block, min_rain, max_rain)
        )
    else:
        raise ValueError("mask_mode must be 'range' or 'finite'")

    if wet_only:
        # Radar must also exist because it defines the rainy sample population.
        mask &= np.isfinite(radar_like) & (radar_like > radar_wet_thr)

    return ref_like, radar_like, mask


# ============================================================
# Diagnostics
# ============================================================
def _finite_quantiles_from_samples(samples: List[np.ndarray]) -> Dict[str, float]:
    if not samples:
        return {"min": np.nan, "p01": np.nan, "p50": np.nan, "p99": np.nan, "max": np.nan}
    a = np.concatenate(samples).astype("float64", copy=False)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {"min": np.nan, "p01": np.nan, "p50": np.nan, "p99": np.nan, "max": np.nan}
    q = np.nanpercentile(a, [1, 50, 99])
    return {"min": float(np.nanmin(a)), "p01": float(q[0]), "p50": float(q[1]), "p99": float(q[2]), "max": float(np.nanmax(a))}


def diagnose_pairs(ds, cml_da_sel, ref_da, radar_da, args, dims, time_slice) -> Dict[str, float]:
    rng = np.random.default_rng(args.random_seed)
    n_cml = ds.sizes[dims["cml_dim"]]
    total_positions = 0
    cml_finite = 0
    ref_finite = 0
    finite_overlap = 0
    cml_in_range = 0
    ref_in_range = 0
    range_overlap = 0
    wet_overlap = 0
    radar_finite_overlap = 0
    radar_wet_positions = 0

    cml_samples: List[np.ndarray] = []
    ref_samples: List[np.ndarray] = []
    max_samples = max(1, int(args.diagnostic_samples))

    for s in tqdm(range(0, n_cml, args.chunk_cml), desc="Diagnostic chunks"):
        e = min(s + args.chunk_cml, n_cml)
        cml_block, ref_block, radar_block = get_block(
            cml_da_sel, ref_da, radar_da, slice(s, e), time_slice, dims
        )
        if cml_block.ndim == 2:
            ref_like = ref_block
            radar_like = radar_block
        else:
            ref_like = np.broadcast_to(ref_block[:, None, :], cml_block.shape)
            radar_like = np.broadcast_to(radar_block[:, None, :], cml_block.shape)

        cf = np.isfinite(cml_block)
        rf = np.isfinite(ref_like)
        radar_finite = np.isfinite(radar_like)
        fo = cf & rf
        cr = valid_range(cml_block, args.min_rain, args.max_rain)
        rr = valid_range(ref_like, args.min_rain, args.max_rain)
        ro = fo & cr & rr
        wo = ro & radar_finite & (radar_like > args.radar_wet_thr)

        total_positions += int(cml_block.size)
        cml_finite += int(cf.sum())
        ref_finite += int(rf.sum())
        finite_overlap += int(fo.sum())
        cml_in_range += int((cf & cr).sum())
        ref_in_range += int((rf & rr).sum())
        range_overlap += int(ro.sum())
        radar_finite_overlap += int((ro & radar_finite).sum())
        radar_wet_positions += int((radar_finite & (radar_like > args.radar_wet_thr)).sum())
        wet_overlap += int(wo.sum())

        nfo = int(fo.sum())
        if nfo:
            keep_prob = min(1.0, max_samples / max(float(finite_overlap), 1.0))
            idx_keep = fo
            if keep_prob < 1.0:
                idx_keep = fo & (rng.random(fo.shape) < keep_prob)
            if idx_keep.any():
                cml_samples.append(cml_block[idx_keep].astype("float32", copy=False))
                ref_samples.append(ref_like[idx_keep].astype("float32", copy=False))
                if sum(x.size for x in cml_samples) > 2 * max_samples:
                    ca = np.concatenate(cml_samples)
                    ra = np.concatenate(ref_samples)
                    take = rng.choice(ca.size, size=max_samples, replace=False)
                    cml_samples = [ca[take]]
                    ref_samples = [ra[take]]

    cq = _finite_quantiles_from_samples(cml_samples)
    rq = _finite_quantiles_from_samples(ref_samples)
    ref_name = args.reference.capitalize()

    print("\nDiagnostics before final mask")
    print(f"  Total CML/reference positions     : {total_positions:,}")
    print(f"  CML finite values                 : {cml_finite:,}")
    print(f"  {ref_name} finite values              : {ref_finite:,}")
    print(f"  Finite CML & {args.reference} overlap  : {finite_overlap:,}")
    print(f"  CML values in [{args.min_rain:g}, {args.max_rain:g}]       : {cml_in_range:,}")
    print(f"  {ref_name} values in [{args.min_rain:g}, {args.max_rain:g}]   : {ref_in_range:,}")
    print(f"  Overlap after range filter        : {range_overlap:,}")
    print(f"  Range overlap with finite radar   : {radar_finite_overlap:,}")
    print(f"  Radar-wet positions (> {args.radar_wet_thr:g})    : {radar_wet_positions:,}")
    print(f"  FINAL fair wet pairs              : {wet_overlap:,}")
    print("  CML finite-overlap quantiles      : "
          f"min={cq['min']:.4g}, p01={cq['p01']:.4g}, p50={cq['p50']:.4g}, p99={cq['p99']:.4g}, max={cq['max']:.4g}")
    print(f"  {ref_name} finite-overlap quantiles    : "
          f"min={rq['min']:.4g}, p01={rq['p01']:.4g}, p50={rq['p50']:.4g}, p99={rq['p99']:.4g}, max={rq['max']:.4g}")

    return {
        "total_positions": total_positions,
        "cml_finite": cml_finite,
        "ref_finite": ref_finite,
        "finite_overlap": finite_overlap,
        "range_overlap": range_overlap,
        "radar_finite_overlap": radar_finite_overlap,
        "radar_wet_positions": radar_wet_positions,
        "wet_overlap": wet_overlap,
    }


# ============================================================
# Metrics and sampling
# ============================================================
def compute_metrics_and_rows(ds, cml_da_sel, ref_da, radar_da, args, dims, time_slice):
    n_cml = ds.sizes[dims["cml_dim"]]
    cml_ids = np.asarray(ds[dims["cml_dim"]].values, dtype=str)

    sublink_values = None
    if dims.get("sublink_dim") is not None and dims["sublink_dim"] in cml_da_sel.dims:
        sublink_values = np.asarray(ds[dims["sublink_dim"]].values)

    stats = empty_stats()
    rows = []

    for s in tqdm(range(0, n_cml, args.chunk_cml), desc="Metrics chunks"):
        e = min(s + args.chunk_cml, n_cml)
        cml_block, ref_block, radar_block = get_block(
            cml_da_sel, ref_da, radar_da, slice(s, e), time_slice, dims
        )
        ref_like, _radar_like, mask = make_masks(
            cml_block, ref_block, radar_block,
            min_rain=args.min_rain,
            max_rain=args.max_rain,
            wet_only=args.wet_only,
            radar_wet_thr=args.radar_wet_thr,
            mask_mode=args.mask_mode,
        )

        update_stats(stats, cml_block, ref_like, mask)

        if cml_block.ndim == 2:
            rows_c = cml_block
            rows_r = ref_like
            rows_m = mask
            row_ids = cml_ids[s:e]
            sl_ids = None
        else:
            C, S, T = cml_block.shape
            rows_c = cml_block.reshape(C * S, T)
            rows_r = ref_like.reshape(C * S, T)
            rows_m = mask.reshape(C * S, T)
            row_ids = np.repeat(cml_ids[s:e], S)
            sl_ids = np.tile(sublink_values, e - s)

        N, MBE, MAE, RMSE, PCC, mean_ref, mean_cml = row_metrics(rows_c, rows_r, rows_m)
        df = pd.DataFrame({
            "cml_id": row_ids,
            "N": N,
            "MBE": MBE,
            "MAE": MAE,
            "RMSE": RMSE,
            "PCC": PCC,
            f"mean_{args.reference}": mean_ref,
            "mean_cml": mean_cml,
        })
        if sl_ids is not None:
            df.insert(1, "sublink_id", sl_ids)
        rows.append(df)

    global_metrics = finalize_stats(stats)
    row_df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return global_metrics, row_df


def collect_plot_samples(ds, cml_da_sel, ref_da, radar_da, args, dims, time_slice, total_N: int):
    rng = np.random.default_rng(args.random_seed)
    n_cml = ds.sizes[dims["cml_dim"]]

    cap = int(args.max_plot_points)
    prob = 1.0 if cap == 0 else min(1.0, cap / max(float(total_N), 1.0))

    samples_r = []
    samples_c = []

    for s in tqdm(range(0, n_cml, args.chunk_cml), desc="Plot-sample chunks"):
        e = min(s + args.chunk_cml, n_cml)
        cml_block, ref_block, radar_block = get_block(
            cml_da_sel, ref_da, radar_da, slice(s, e), time_slice, dims
        )
        ref_like, _radar_like, mask = make_masks(
            cml_block, ref_block, radar_block,
            min_rain=args.min_rain,
            max_rain=args.max_rain,
            wet_only=args.wet_only,
            radar_wet_thr=args.radar_wet_thr,
            mask_mode=args.mask_mode,
        )

        n = int(np.sum(mask))
        if n == 0:
            continue

        rv = ref_like[mask].astype("float32", copy=False)
        cv = cml_block[mask].astype("float32", copy=False)

        if prob < 1.0:
            keep = rng.random(n) < prob
            rv = rv[keep]
            cv = cv[keep]

        if rv.size:
            samples_r.append(rv)
            samples_c.append(cv)

    if not samples_r:
        return np.array([], dtype="float32"), np.array([], dtype="float32")

    r = np.concatenate(samples_r)
    c = np.concatenate(samples_c)

    if cap > 0 and r.size > cap:
        idx = rng.choice(r.size, size=cap, replace=False)
        r = r[idx]
        c = c[idx]

    return r, c


# ============================================================
# Plotting
# ============================================================
def metrics_text(metrics: Dict[str, float]) -> str:
    return (
        f"N={int(metrics['N']):,}\n"
        f"RMSE={metrics['RMSE']:.2f}\n"
        f"MBE={metrics['MBE']:.2f}\n"
        f"PCC={metrics['PCC']:.2f}"
    )


def make_single_hexbin_plot(r, c, metrics, args, out_png: Path, out_pdf: Path):
    if r.size == 0:
        raise RuntimeError(
            "No plot samples were collected. No hexbin figure will be saved. "
            "Check the diagnostics, time window, reference assignment, or the mask settings."
        )

    fig, ax = plt.subplots(1, 1, figsize=(4.6, 4.2))

    hb = ax.hexbin(
        r,
        c,
        gridsize=args.gridsize,
        mincnt=args.mincnt,
        cmap=get_count_cmap(args.cmap),
        linewidths=0.0,
        norm=LogNorm() if args.log_counts else None,
        extent=(args.min_rain, args.max_rain, args.min_rain, args.max_rain),
    )

    ax.plot([0, args.max_rain], [0, args.max_rain], color="0.25", lw=0.9, ls="--", zorder=4)

    if args.metrics_box:
        ax.text(
            args.metrics_x, args.metrics_y, metrics_text(metrics),
            transform=ax.transAxes,
            ha="left", va="bottom",
            fontsize=args.metrics_fontsize,
            fontfamily="DejaVu Sans Mono",
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.75", alpha=0.90),
            zorder=10,
        )

    if args.title_mode == "combo":
        ax.set_title(
            f"{args.baseline_method} + {args.kr_model} + {args.waa_method}",
            fontsize=args.title_fontsize,
            fontweight="bold",
            pad=6,
        )
    elif args.title_mode == "metrics":
        ax.set_title(
            f"RMSE={metrics['RMSE']:.2f}, MBE={metrics['MBE']:.2f}, PCC={metrics['PCC']:.2f}",
            fontsize=args.title_fontsize,
            fontweight="normal",
            pad=5,
        )

    ref_label = "gauge" if args.reference == "gauge" else "radar"
    style_axis(
        ax,
        lim_max=args.max_rain,
        xlabel=rf"RR$_{{{ref_label}}}$ [mm h$^{{-1}}$]",
        ylabel=r"RR$_{CML}$ [mm h$^{-1}$]",
    )

    max_count = float(np.nanmax(hb.get_array())) if hb.get_array().size else 1.0
    base_cmap = get_count_cmap(args.cmap)

    if args.log_counts:
        vmax = float(args.count_vmax) if args.count_vmax > 0 else 10.0 ** np.ceil(np.log10(max(1.0, max_count)))
        boundaries = make_log_count_boundaries(vmax=vmax, n_levels=args.count_levels)
        n_extend = {"neither": 0, "min": 1, "max": 1, "both": 2}[args.colorbar_extend]
        n_regions = max(1, len(boundaries) - 1 + n_extend)
        cmap = mpl.colors.ListedColormap(base_cmap(np.linspace(0, 1, n_regions)), name=f"{base_cmap.name}_{n_regions}_bins")
        common_norm = BoundaryNorm(boundaries, ncolors=cmap.N, clip=(args.colorbar_extend == "neither"), extend=args.colorbar_extend)
        hb.set_cmap(cmap)
        hb.set_norm(common_norm)
        hb.set_clim(1, vmax)

        fig.subplots_adjust(right=0.84)
        cax = fig.add_axes([0.865, 0.20, 0.030, 0.62])
        cbar = fig.colorbar(hb, cax=cax, boundaries=boundaries, spacing="uniform", extend=args.colorbar_extend, drawedges=False)
        ticks = log_count_ticks(vmax)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([format_count_tick(t, args.count_tick_mode) for t in ticks])
        cbar.ax.minorticks_off()
    else:
        vmax = float(args.count_vmax) if args.count_vmax > 0 else max_count
        hb.set_clim(0, vmax)
        fig.subplots_adjust(right=0.84)
        cax = fig.add_axes([0.865, 0.20, 0.030, 0.62])
        cbar = fig.colorbar(hb, cax=cax)

    cbar.set_label("Count", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


# ============================================================
# CLI
# ============================================================
def parse_args():
    p = argparse.ArgumentParser(description="Single-combination CML-vs-reference hexbin/scatter plot.")

    p.add_argument("--rain-nc", default="cml_rainrate_baseline_kr_waa_methods.nc")
    p.add_argument("--out-dir", default="single_combo_reference_compare")
    p.add_argument("--cml-var", default="R_mean_5min",
                   help="Use R_mean_5min for link mean or R_sublink_5min for independent sublinks.")

    p.add_argument("--reference", choices=["radar", "gauge"], default="radar")
    p.add_argument("--radar-var", default="R_radar")

    p.add_argument("--baseline-method", default="pycml_constant")
    p.add_argument("--kr-model", default="ITU_2003")
    p.add_argument("--waa-method", default="pastorek_2021")

    p.add_argument("--tstart", default=None, help="Optional start time, e.g. '2021-06-01 00:00:00'.")
    p.add_argument("--tend", default=None, help="Optional inclusive end time, e.g. '2021-06-30 23:55:00'.")

    p.add_argument("--min-rain", type=float, default=0.0)
    p.add_argument("--max-rain", type=float, default=100.0)
    p.add_argument(
        "--wet-only", dest="wet_only", action="store_true", default=True,
        help="Keep only radar-defined wet pairs (default and recommended).",
    )
    p.add_argument(
        "--all-pairs", dest="wet_only", action="store_false",
        help="Disable the radar wet filter and use all fair finite/range pairs.",
    )
    p.add_argument(
        "--radar-wet-thr", type=float, default=0.1,
        help="Strict radar wet threshold: retained samples satisfy radar > threshold.",
    )
    p.add_argument("--mask-mode", choices=["range", "finite"], default="range")
    p.add_argument("--diagnostics-only", action="store_true")
    p.add_argument("--diagnostic-samples", type=int, default=300000)

    p.add_argument("--chunk-cml", type=int, default=64)
    p.add_argument("--max-plot-points", type=int, default=3_000_000,
                   help="Random points used only for plotting. Metrics are exact. Use 0 to plot all valid points.")
    p.add_argument("--random-seed", type=int, default=42)

    p.add_argument("--gridsize", type=int, default=50)
    p.add_argument("--mincnt", type=int, default=1)
    p.add_argument("--cmap", default="reference")
    p.add_argument("--log-counts", action="store_true", default=True)
    p.add_argument("--no-log-counts", dest="log_counts", action="store_false")
    p.add_argument("--count-vmax", type=float, default=0.0)
    p.add_argument("--count-levels", type=int, default=16)
    p.add_argument("--count-tick-mode", choices=["power", "plain"], default="power")
    p.add_argument("--colorbar-extend", choices=["neither", "min", "max", "both"], default="neither")

    p.add_argument("--title-mode", choices=["none", "combo", "metrics"], default="none")
    p.add_argument("--title-fontsize", type=float, default=10)
    p.add_argument("--metrics-box", action="store_true", default=True)
    p.add_argument("--no-metrics-box", dest="metrics_box", action="store_false")
    p.add_argument("--metrics-fontsize", type=float, default=8.6)
    p.add_argument("--metrics-x", type=float, default=0.64)
    p.add_argument("--metrics-y", type=float, default=0.03)

    # Gauge options
    p.add_argument("--gauge-nc", default="CHMI_10MIN_v1.nc")
    p.add_argument("--gauge-rr-var", default="rr")
    p.add_argument("--gauge-lon-var", default="lon")
    p.add_argument("--gauge-lat-var", default="lat")
    p.add_argument("--gauge-minutes", type=float, default=10.0,
                   help="Gauge accumulation window in minutes. Converted to mm h-1.")
    p.add_argument("--max-gauge-dist-km", type=float, default=5.0,
                   help="Maximum distance from CML midpoint to assigned gauge.")
    p.add_argument("--cml-coords-nc", default=None,
                   help="Optional CML metadata file with site_a/site_b coordinates if they are not in the rain file.")
    p.add_argument("--site-a-lat-var", default="site_a_latitude")
    p.add_argument("--site-a-lon-var", default="site_a_longitude")
    p.add_argument("--site-b-lat-var", default="site_b_latitude")
    p.add_argument("--site-b-lon-var", default="site_b_longitude")

    return p.parse_args()


# ============================================================
# Main
# ============================================================
def main():
    args = parse_args()
    t0 = time.perf_counter()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    wet_tag = f"radar_wet_gt_{args.radar_wet_thr:g}".replace(".", "p") if args.wet_only else "all_pairs"
    suffix = f"vs_{args.reference}_{wet_tag}"
    combo_key = f"{args.baseline_method}_{args.kr_model}_{args.waa_method}"
    out_png = out_dir / f"hexbin_{combo_key}_{suffix}.png"
    out_pdf = out_dir / f"hexbin_{combo_key}_{suffix}.pdf"
    global_csv = out_dir / f"global_metrics_{combo_key}_{suffix}.csv"
    row_csv = out_dir / f"per_sensor_metrics_{combo_key}_{suffix}.csv"
    mapping_csv = out_dir / "gauge_link_mapping.csv"

    print("\nLoading dataset")
    print(f"  file: {args.rain_nc}")
    ds = xr.open_dataset(args.rain_nc)

    try:
        if args.cml_var not in ds.data_vars:
            raise KeyError(f"CML variable {args.cml_var!r} not found. Available data variables: {list(ds.data_vars)}")

        print("\nSelecting CML estimate")
        cml_da = select_cml_combo(ds, args.cml_var, args)

        dims = {
            "cml_dim": find_dim_by_candidates(cml_da, ["cml_id", "cml"]),
            "time_dim": find_dim_by_candidates(cml_da, ["time_5min", "time"]),
            "sublink_dim": find_dim_by_candidates(cml_da, ["sublink_id", "sublink"], required=False),
            "ref_cml_dim": None,
            "ref_time_dim": None,
            "radar_cml_dim": None,
            "radar_time_dim": None,
        }

        # Radar is always loaded because it defines wet instances, even when
        # the plotted comparison reference is a rain gauge.
        if args.radar_var not in ds.data_vars:
            raise KeyError(
                f"Radar wet-mask variable {args.radar_var!r} not found. "
                f"Available: {list(ds.data_vars)}"
            )
        radar_da = ds[args.radar_var]
        dims["radar_cml_dim"] = find_dim_by_candidates(radar_da, ["cml_id", "cml"])
        dims["radar_time_dim"] = find_dim_by_candidates(radar_da, ["time_5min", "time"])

        # Select the plotted/evaluated reference.
        gauge_mapping = None
        if args.reference == "radar":
            ref_da = radar_da
            dims["ref_cml_dim"] = dims["radar_cml_dim"]
            dims["ref_time_dim"] = dims["radar_time_dim"]
        else:
            print(f"  Gauge file        : {args.gauge_nc}")
            print(f"  Gauge var         : {args.gauge_rr_var}  ({args.gauge_minutes:g}-min amounts -> mm h-1)")
            print(f"  Max gauge dist    : {args.max_gauge_dist_km:g} km")
            ref_da, gauge_mapping = build_gauge_reference(ds, args, dims)
            dims["ref_cml_dim"] = find_dim_by_candidates(ref_da, ["cml_id", "cml"])
            dims["ref_time_dim"] = find_dim_by_candidates(ref_da, ["time_5min", "time"])
            gauge_mapping.to_csv(mapping_csv, index=False)

        if ds.sizes[dims["cml_dim"]] != ref_da.sizes[dims["ref_cml_dim"]]:
            raise ValueError("CML and reference cml dimensions have different sizes.")
        if ds.sizes[dims["time_dim"]] != ref_da.sizes[dims["ref_time_dim"]]:
            raise ValueError("CML and reference time dimensions have different sizes. Use matched variables.")
        if ds.sizes[dims["cml_dim"]] != radar_da.sizes[dims["radar_cml_dim"]]:
            raise ValueError("CML and radar wet-mask cml dimensions have different sizes.")
        if ds.sizes[dims["time_dim"]] != radar_da.sizes[dims["radar_time_dim"]]:
            raise ValueError("CML and radar wet-mask time dimensions have different sizes.")

        print("\nInput")
        print(f"  CML variable      : {args.cml_var}")
        print(f"  Reference         : {args.reference}")
        if args.reference == "radar":
            print(f"  Radar variable    : {args.radar_var}")
        else:
            print(f"  Gauge variable    : {args.gauge_rr_var}")
        print(f"  Selected combo    : {args.baseline_method} + {args.kr_model} + {args.waa_method}")
        print(f"  Dimensions        : {dict(cml_da.sizes)}")
        print(f"  Cutoff            : {args.min_rain:g} to {args.max_rain:g} mm h-1")
        print(f"  Mask mode         : {args.mask_mode}")
        print(f"  Fair pairs only   : finite CML + finite {args.reference}")
        print(f"  Wet only          : {args.wet_only}")
        if args.wet_only:
            print(f"  Wet selector      : {args.radar_var} > {args.radar_wet_thr:g} mm h-1 (strict)")
        print(f"  Hexbin gridsize   : {args.gridsize}")
        print(f"  Output directory  : {out_dir}")

        time_slice, _ = parse_time_window(ds, dims["time_dim"], args)

        diag = diagnose_pairs(ds, cml_da, ref_da, radar_da, args, dims, time_slice)
        if args.diagnostics_only:
            print("\nDiagnostics-only mode: stopping before metrics/plotting.")
            return

        print("\nComputing exact metrics")
        metrics, row_df = compute_metrics_and_rows(ds, cml_da, ref_da, radar_da, args, dims, time_slice)
        metrics_df = pd.DataFrame([{
            "baseline_method": args.baseline_method,
            "kr_model": args.kr_model,
            "waa_method": args.waa_method,
            "reference": args.reference,
            **metrics,
            "min_rain": args.min_rain,
            "max_rain": args.max_rain,
            "wet_only": args.wet_only,
            "wet_mask_source": args.radar_var if args.wet_only else "none",
            "radar_wet_operator": ">" if args.wet_only else "none",
            "radar_wet_thr": args.radar_wet_thr if args.wet_only else np.nan,
            "mask_mode": args.mask_mode,
            "diagnostic_finite_overlap": diag.get("finite_overlap", np.nan),
            "diagnostic_range_overlap": diag.get("range_overlap", np.nan),
            "diagnostic_radar_finite_overlap": diag.get("radar_finite_overlap", np.nan),
            "diagnostic_final_wet_pairs": diag.get("wet_overlap", np.nan),
            "gridsize": args.gridsize,
        }])

        metrics_df.to_csv(global_csv, index=False)
        row_df.to_csv(row_csv, index=False)

        print("\nGlobal metrics")
        print(metrics_df.to_string(index=False))

        print("\nCollecting plot samples")
        r, c = collect_plot_samples(
            ds, cml_da, ref_da, radar_da, args, dims, time_slice, total_N=int(metrics["N"])
        )
        print(f"  Plot samples: {r.size:,}")

        print("\nPlotting")
        make_single_hexbin_plot(r, c, metrics, args, out_png, out_pdf)

        print("\nSaved")
        print(f"  {out_png}")
        print(f"  {out_pdf}")
        print(f"  {global_csv}")
        print(f"  {row_csv}")
        if gauge_mapping is not None:
            print(f"  {mapping_csv}")
        print(f"\nTotal elapsed: {time.perf_counter() - t0:.2f} s")

    finally:
        ds.close()


if __name__ == "__main__":
    main()
