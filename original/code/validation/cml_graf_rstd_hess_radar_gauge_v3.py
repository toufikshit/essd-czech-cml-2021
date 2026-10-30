#!/usr/bin/env python3
"""
cml_graf_rstd_hess_radar_gauge.py
=================================

Apply the established Graf/RSTD wet-dry detector to preprocessed 1-min CML
data and evaluate it independently against:

1. path-averaged radar at native 5-min resolution
2. nearest rain gauge within 1 km of the CML link midpoint at native 10-min
   resolution

No directional-sublink averaging is performed. Every sublink is an
independent CML observation in the contingency counts.

Graf/RSTD detector
------------------
For every directional sublink independently:

    RSTD = centred rolling sample standard deviation over 60 min
    Q80  = 80th percentile of RSTD computed where the complete 60-sample
           threshold window exists
    threshold = 1.12 * Q80
    wet = RSTD_detection > threshold

Defaults:
    rstd_window             = 60
    graf_quantile           = 0.80
    graf_factor             = 1.12
    threshold_min_periods   = 60
    detection_min_periods   = 45

Reference-resolution comparison
-------------------------------
The Graf/RSTD detector itself is computed at the native 1-min CML resolution.
Only AFTER detection, the 1-min CML wet/dry flags are transformed to the
native resolution of each independent reference:

    1-min CML Graf flags -> exact native 5-min radar intervals
    1-min CML Graf flags -> exact native 10-min gauge intervals

No CML signal interpolation, reference filling, or forward filling is used.
This preserves the 60-min Graf detector definition while ensuring that every
contingency comparison is performed at the resolution of the corresponding
reference.

Radar:
    native resolution       = normally 5 min
    Czech rr amount -> rate = raw * 12
    wet                     = rate > 0.1 mm h-1
    timestamp label         = start

Gauge:
    native resolution       = normally 10 min
    CHMI rr amount -> rate  = raw * 6
    wet                     = rate > 0.1 mm h-1
    closest gauge           = nearest to CML LINK MIDPOINT
    maximum distance        = 1.0 km
    timestamp label         = end

The 1-min CML detector is aggregated independently to each native reference
interval. By default, ANY wet 1-min detector value makes that reference
interval CML-wet, but ALL expected CML minutes in the interval must have a
valid detector result (--min-cml-valid-fraction 1.0).

There is no reference filling or interpolation. Missing CML detector values,
missing reference values, incomplete CML intervals, and unmapped gauges are
excluded, never interpreted as dry.

Outputs
-------
<out-dir>/
    radar/
        graf_rstd_vs_radar_contingency_hess.png
        graf_rstd_vs_radar_contingency_hess.pdf
        graf_rstd_vs_radar_contingency.csv
        graf_rstd_vs_radar_metrics.csv
        graf_rstd_vs_radar_per_sublink.csv
        graf_rstd_vs_radar_contingency_table.tex

    gauge/
        graf_rstd_vs_gauge_contingency_hess.png
        graf_rstd_vs_gauge_contingency_hess.pdf
        graf_rstd_vs_gauge_contingency.csv
        graf_rstd_vs_gauge_metrics.csv
        graf_rstd_vs_gauge_per_sublink.csv
        graf_rstd_vs_gauge_contingency_table.tex
        closest_gauge_assignment.csv

    graf_rstd_radar_gauge_contingency_hess.png
    graf_rstd_radar_gauge_contingency_hess.pdf
    graf_rstd_thresholds.csv
    run_summary.json

Figures:
    vector PDF + 600-dpi PNG
    clean HESS-oriented contingency heatmap layout
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# Utilities
# =============================================================================
def ordered_intersection(first: Sequence[str], second: Sequence[str]) -> List[str]:
    s = set(str(x) for x in second)
    return [str(x) for x in first if str(x) in s]


def detect_trsl_source(ds: xr.Dataset, requested: Optional[str]) -> Optional[str]:
    if requested:
        if requested not in ds.variables:
            raise KeyError(
                f"Requested TRSL variable {requested!r} not found. "
                f"Available: {list(ds.variables)}"
            )
        return requested

    for name in ("trsl", "trsl_qc"):
        if name in ds.variables:
            return name

    if "tx" in ds.variables and "rx" in ds.variables:
        return None

    raise KeyError("No trsl/trsl_qc and cannot reconstruct from tx-rx.")


def detect_var(
    ds: xr.Dataset,
    requested: Optional[str],
    candidates: Sequence[str],
    role: str,
) -> str:
    if requested:
        if requested in ds.variables:
            return requested
        raise KeyError(
            f"Requested {role} variable {requested!r} not found. "
            f"Available: {list(ds.variables)}"
        )

    for name in candidates:
        if name in ds.variables:
            return name

    raise KeyError(
        f"Could not detect {role}. Tried {list(candidates)}. "
        f"Available: {list(ds.variables)}"
    )


def detect_time_dim(da: xr.DataArray) -> str:
    for d in da.dims:
        if "time" in d.lower():
            return d
    raise KeyError(f"No time dimension in {da.dims}")


def detect_station_dim(da: xr.DataArray, time_dim: str) -> str:
    others = [d for d in da.dims if d != time_dim]
    if len(others) != 1:
        raise ValueError(
            f"Expected one station dimension besides {time_dim!r}; got {da.dims}"
        )
    return others[0]


def median_native_step_ns(times: pd.DatetimeIndex, fallback_minutes: float) -> int:
    if len(times) < 2:
        return int(round(fallback_minutes * 60.0 * 1e9))

    diffs = np.diff(times.view("int64"))
    diffs = diffs[diffs > 0]
    if diffs.size == 0:
        return int(round(fallback_minutes * 60.0 * 1e9))

    values, counts = np.unique(diffs, return_counts=True)
    return int(values[counts == counts.max()].min())


def haversine_matrix_km(
    lon1: np.ndarray,
    lat1: np.ndarray,
    lon2: np.ndarray,
    lat2: np.ndarray,
) -> np.ndarray:
    """
    Distance matrix shape (len(lon1), len(lon2)).
    """
    radius = 6371.0

    lon1r = np.radians(np.asarray(lon1, dtype=np.float64))[:, None]
    lat1r = np.radians(np.asarray(lat1, dtype=np.float64))[:, None]
    lon2r = np.radians(np.asarray(lon2, dtype=np.float64))[None, :]
    lat2r = np.radians(np.asarray(lat2, dtype=np.float64))[None, :]

    dlon = lon2r - lon1r
    dlat = lat2r - lat1r

    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2.0) ** 2
    )
    return 2.0 * radius * np.arcsin(np.sqrt(a))


# =============================================================================
# Graf/RSTD detector
# =============================================================================
def centred_window_bounds(n: int, window: int) -> Tuple[np.ndarray, np.ndarray]:
    i = np.arange(n, dtype=np.int64)
    left = window // 2
    right = window - left
    return np.maximum(0, i - left), np.minimum(n, i + right)


def graf_rstd_matrix(
    values: np.ndarray,
    window: int,
    quantile: float,
    factor: float,
    threshold_min_periods: int,
    detection_min_periods: int,
):
    """
    Vectorized detector for a time x independent_sublink matrix.
    """
    x = np.asarray(values, dtype=np.float64)
    n_time, n_obs = x.shape

    finite = np.isfinite(x)
    clean = np.where(finite, x, 0.0)

    csum = np.empty((n_time + 1, n_obs), dtype=np.float64)
    csum2 = np.empty((n_time + 1, n_obs), dtype=np.float64)
    ccount = np.empty((n_time + 1, n_obs), dtype=np.int32)

    csum[0] = 0.0
    csum2[0] = 0.0
    ccount[0] = 0

    np.cumsum(clean, axis=0, out=csum[1:])
    np.cumsum(clean * clean, axis=0, out=csum2[1:])
    np.cumsum(finite.astype(np.int32), axis=0, out=ccount[1:])

    start, end = centred_window_bounds(n_time, window)

    sums = csum[end] - csum[start]
    sums2 = csum2[end] - csum2[start]
    counts = ccount[end] - ccount[start]

    denom = np.maximum(counts - 1, 1)
    variance = (
        sums2 - (sums * sums) / np.maximum(counts, 1)
    ) / denom
    variance = np.maximum(variance, 0.0)

    std_detection = np.sqrt(variance)
    std_detection[counts < detection_min_periods] = np.nan

    std_threshold = std_detection.copy()
    std_threshold[counts < threshold_min_periods] = np.nan

    with np.errstate(all="ignore"):
        q_values = np.nanquantile(std_threshold, quantile, axis=0)

    thresholds = factor * q_values

    valid = (
        np.isfinite(std_detection)
        & np.isfinite(thresholds[None, :])
    )
    wet = valid & (std_detection > thresholds[None, :])

    finite_input = finite.sum(axis=0, dtype=np.int64)

    return wet, valid, thresholds, q_values, finite_input


# =============================================================================
# Native interval handling
# =============================================================================
def build_reference_intervals(
    cml_times: pd.DatetimeIndex,
    reference_times: pd.DatetimeIndex,
    timestamp_label: str,
    fallback_minutes: float,
) -> Dict[str, np.ndarray | int]:
    """
    Build exact native reference intervals on the CML time grid.

    No fill. Missing reference timestamps do not widen neighbouring intervals.
    """
    step_ns = median_native_step_ns(reference_times, fallback_minutes)
    cml_step_ns = median_native_step_ns(cml_times, 1.0)
    expected = max(1, int(round(step_ns / cml_step_ns)))

    ref_ns = reference_times.view("int64")

    if timestamp_label == "start":
        interval_start = ref_ns
        interval_end = ref_ns + step_ns
    elif timestamp_label == "end":
        interval_start = ref_ns - step_ns
        interval_end = ref_ns
    elif timestamp_label == "center":
        half = step_ns // 2
        interval_start = ref_ns - half
        interval_end = interval_start + step_ns
    else:
        raise ValueError(timestamp_label)

    cml_ns = cml_times.view("int64")

    start_idx = np.searchsorted(
        cml_ns, interval_start, side="left"
    ).astype(np.int64)
    end_idx = np.searchsorted(
        cml_ns, interval_end, side="left"
    ).astype(np.int64)

    observed = end_idx - start_idx

    cml_end_exclusive = int(cml_ns[-1] + cml_step_ns)

    inside = (
        (interval_start >= int(cml_ns[0]))
        & (interval_end <= cml_end_exclusive)
    )

    complete = inside & (observed == expected)

    return {
        "start_idx": start_idx,
        "end_idx": end_idx,
        "complete": complete.astype(bool),
        "expected": int(expected),
        "step_ns": int(step_ns),
    }


def aggregate_to_reference(
    wet_1min: np.ndarray,
    valid_1min: np.ndarray,
    mapping: Dict[str, np.ndarray | int],
    min_valid_fraction: float,
    temporal_aggregation: str,
):
    n_time, n_obs = wet_1min.shape

    cwet = np.empty((n_time + 1, n_obs), dtype=np.int32)
    cvalid = np.empty((n_time + 1, n_obs), dtype=np.int32)
    cwet[0] = 0
    cvalid[0] = 0

    np.cumsum(
        wet_1min & valid_1min,
        axis=0,
        dtype=np.int32,
        out=cwet[1:],
    )
    np.cumsum(
        valid_1min,
        axis=0,
        dtype=np.int32,
        out=cvalid[1:],
    )

    i0 = np.asarray(mapping["start_idx"], dtype=np.int64)
    i1 = np.asarray(mapping["end_idx"], dtype=np.int64)

    wet_count = cwet[i1] - cwet[i0]
    valid_count = cvalid[i1] - cvalid[i0]

    required = max(
        1,
        int(
            math.ceil(
                int(mapping["expected"]) * float(min_valid_fraction)
            )
        ),
    )

    valid_ref = (
        np.asarray(mapping["complete"], dtype=bool)[:, None]
        & (valid_count >= required)
    )

    if temporal_aggregation == "any":
        pred_wet = wet_count > 0
    elif temporal_aggregation == "majority":
        pred_wet = wet_count > (valid_count / 2.0)
    elif temporal_aggregation == "all":
        pred_wet = (
            (valid_count > 0)
            & (wet_count == valid_count)
        )
    else:
        raise ValueError(temporal_aggregation)

    pred_wet &= valid_ref

    return pred_wet, valid_ref


# =============================================================================
# Gauge mapping
# =============================================================================
def prepare_gauge_mapping(
    cml: xr.Dataset,
    cml_ids: np.ndarray,
    gauge: xr.Dataset,
    gauge_var: str,
    args,
):
    """
    Nearest gauge to CML LINK MIDPOINT; accept only <= max distance.
    """
    needed = [
        args.site_a_lon,
        args.site_a_lat,
        args.site_b_lon,
        args.site_b_lat,
    ]
    missing = [name for name in needed if name not in cml.variables]
    if missing:
        raise KeyError(
            f"CML endpoint coordinates missing: {missing}. "
            "Use the publication preprocessed CML file containing site_a/site_b coordinates."
        )

    source_ids = np.asarray(cml["cml_id"].values).astype(str)
    source_pos = {v: i for i, v in enumerate(source_ids)}
    order = np.asarray(
        [source_pos[str(x)] for x in cml_ids],
        dtype=np.int64,
    )

    a_lon = np.asarray(
        cml[args.site_a_lon].isel(cml_id=order).values,
        dtype=np.float64,
    ).reshape(-1)
    a_lat = np.asarray(
        cml[args.site_a_lat].isel(cml_id=order).values,
        dtype=np.float64,
    ).reshape(-1)
    b_lon = np.asarray(
        cml[args.site_b_lon].isel(cml_id=order).values,
        dtype=np.float64,
    ).reshape(-1)
    b_lat = np.asarray(
        cml[args.site_b_lat].isel(cml_id=order).values,
        dtype=np.float64,
    ).reshape(-1)

    mid_lon = 0.5 * (a_lon + b_lon)
    mid_lat = 0.5 * (a_lat + b_lat)

    rr = gauge[gauge_var]
    gauge_time_dim = detect_time_dim(rr)
    station_dim = detect_station_dim(rr, gauge_time_dim)

    lon_name = detect_var(
        gauge,
        args.gauge_lon,
        ("lon", "longitude", "station_lon", "station_longitude", "LON"),
        "gauge longitude",
    )
    lat_name = detect_var(
        gauge,
        args.gauge_lat,
        ("lat", "latitude", "station_lat", "station_latitude", "LAT"),
        "gauge latitude",
    )

    gauge_lon = np.asarray(
        gauge[lon_name].values,
        dtype=np.float64,
    ).reshape(-1)
    gauge_lat = np.asarray(
        gauge[lat_name].values,
        dtype=np.float64,
    ).reshape(-1)

    if gauge_lon.size != gauge.sizes[station_dim]:
        raise ValueError(
            f"Gauge longitude size {gauge_lon.size} != station dimension "
            f"{station_dim}={gauge.sizes[station_dim]}"
        )
    if gauge_lat.size != gauge.sizes[station_dim]:
        raise ValueError(
            f"Gauge latitude size {gauge_lat.size} != station dimension "
            f"{station_dim}={gauge.sizes[station_dim]}"
        )

    finite_station = np.isfinite(gauge_lon) & np.isfinite(gauge_lat)
    valid_station_pos = np.flatnonzero(finite_station)

    if valid_station_pos.size == 0:
        raise RuntimeError("No finite gauge coordinates.")

    distances = haversine_matrix_km(
        mid_lon,
        mid_lat,
        gauge_lon[valid_station_pos],
        gauge_lat[valid_station_pos],
    )

    nearest_local = np.argmin(distances, axis=1)
    nearest_dist = distances[
        np.arange(len(cml_ids)),
        nearest_local,
    ]
    nearest_station_pos = valid_station_pos[nearest_local]

    accepted = (
        np.isfinite(mid_lon)
        & np.isfinite(mid_lat)
        & np.isfinite(nearest_dist)
        & (nearest_dist <= float(args.max_gauge_dist_km))
    )

    assigned = np.where(
        accepted,
        nearest_station_pos,
        -1,
    ).astype(np.int64)

    # Read station IDs if possible.
    id_candidates = [
        args.gauge_id,
        "WSI",
        "station_id",
        "station",
        "gauge",
        "GH_ID",
        station_dim,
    ]
    id_name = next(
        (
            name for name in id_candidates
            if name and name in gauge.variables
            and gauge[name].size == gauge.sizes[station_dim]
        ),
        None,
    )

    if id_name is not None:
        station_ids = np.asarray(
            gauge[id_name].values
        ).astype(str).reshape(-1)
    else:
        station_ids = np.asarray(
            [f"gauge_{i:04d}" for i in range(gauge.sizes[station_dim])],
            dtype=str,
        )

    assigned_ids = np.asarray(
        [
            station_ids[p] if p >= 0 else ""
            for p in assigned
        ],
        dtype=str,
    )

    mapping_df = pd.DataFrame(
        {
            "cml_id": cml_ids,
            "link_midpoint_longitude": mid_lon,
            "link_midpoint_latitude": mid_lat,
            "gauge_position": assigned,
            "gauge_id": assigned_ids,
            "distance_km": nearest_dist,
            "within_1km": accepted.astype(np.uint8),
        }
    )

    return {
        "rr": rr,
        "time_dim": gauge_time_dim,
        "station_dim": station_dim,
        "assigned_positions": assigned,
        "mapping": mapping_df,
    }


# =============================================================================
# Contingency metrics
# =============================================================================
def counts_for_column(
    pred: np.ndarray,
    cml_valid: np.ndarray,
    ref_wet: np.ndarray,
    ref_valid: np.ndarray,
):
    common = cml_valid & ref_valid

    tp = int(np.sum(common & pred & ref_wet))
    tn = int(np.sum(common & ~pred & ~ref_wet))
    fp = int(np.sum(common & pred & ~ref_wet))
    fn = int(np.sum(common & ~pred & ref_wet))

    return tp, tn, fp, fn


def metrics_from_counts(
    tp: int,
    tn: int,
    fp: int,
    fn: int,
):
    total = tp + tn + fp + fn

    pod = tp / (tp + fn) if (tp + fn) else np.nan
    far = fp / (tp + fp) if (tp + fp) else np.nan
    specificity = tn / (tn + fp) if (tn + fp) else np.nan
    accuracy = (tp + tn) / total if total else np.nan
    csi = tp / (tp + fp + fn) if (tp + fp + fn) else np.nan
    frequency_bias = (tp + fp) / (tp + fn) if (tp + fn) else np.nan

    denominator = math.sqrt(
        (tp + fp)
        * (tp + fn)
        * (tn + fp)
        * (tn + fn)
    )
    mcc = (
        ((tp * tn) - (fp * fn)) / denominator
        if denominator
        else np.nan
    )

    return {
        "N": int(total),
        "POD": float(pod),
        "FAR": float(far),
        "specificity": float(specificity),
        "accuracy": float(accuracy),
        "CSI": float(csi),
        "frequency_bias": float(frequency_bias),
        "MCC": float(mcc),
    }


# =============================================================================
# HESS publication figures/tables
# =============================================================================
def contingency_matrix(counts: Dict[str, int]) -> np.ndarray:
    """
    Row order   : CML wet, CML dry
    Column order: Reference wet, Reference dry
    """
    return np.array(
        [
            [counts["TP"], counts["FP"]],
            [counts["FN"], counts["TN"]],
        ],
        dtype=np.float64,
    )


def draw_contingency_plot_on_axis(
    ax,
    counts: Dict[str, int],
    metrics: Dict[str, float],
    reference_label: str,
    panel_label: Optional[str] = None,
):
    """
    HESS-ready 2x2 contingency heatmap.
    """
    mat = contingency_matrix(counts)
    plot_values = np.log10(mat + 1.0)

    vmax = float(np.nanmax(plot_values))
    if not np.isfinite(vmax) or vmax <= 0.0:
        vmax = 1.0

    im = ax.imshow(
        plot_values,
        origin="upper",
        cmap="Blues",
        vmin=0.0,
        vmax=vmax,
        interpolation="nearest",
        aspect="equal",
    )

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(
        [f"{reference_label} wet", f"{reference_label} dry"],
        fontsize=9,
    )
    ax.set_yticklabels(
        ["CML wet", "CML dry"],
        fontsize=9,
    )
    ax.tick_params(axis="x", bottom=False, top=False, pad=6)
    ax.tick_params(axis="y", left=False, right=False, pad=6)

    ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 2, 1), minor=True)
    ax.grid(which="minor", color="black", linestyle="-", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    for spine in ax.spines.values():
        spine.set_visible(False)

    tp = counts["TP"]
    fp = counts["FP"]
    fn = counts["FN"]
    tn = counts["TN"]

    ref_wet = tp + fn
    ref_dry = fp + tn

    cell_names = np.array(
        [
            ["Hit (TP)", "False alarm (FP)"],
            ["Miss (FN)", "Correct negative (TN)"],
        ],
        dtype=object,
    )

    def pct(value: float, denom: float) -> str:
        return f"{100.0 * value / denom:.1f}%" if denom else "–"

    ref_denoms = np.array(
        [
            [ref_wet, ref_dry],
            [ref_wet, ref_dry],
        ],
        dtype=np.float64,
    )

    threshold = vmax * 0.55

    for i in range(2):
        for j in range(2):
            value = int(mat[i, j])
            color = "white" if plot_values[i, j] >= threshold else "black"
            ax.text(
                j,
                i,
                f"{cell_names[i, j]}\n{value:,}\n{pct(value, ref_denoms[i, j])}",
                ha="center",
                va="center",
                fontsize=9.2,
                color=color,
            )

    metric_text = (
        f"N = {int(metrics['N']):,}    "
        f"POD = {metrics['POD']:.3f}    "
        f"FAR = {metrics['FAR']:.3f}    "
        f"MCC = {metrics['MCC']:.3f}    "
        f"CSI = {metrics['CSI']:.3f}"
    )

    ax.text(
        0.5,
        -0.23,
        metric_text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=9.0,
    )

    ax.text(
        0.5,
        -0.33,
        "Cell shading uses log10(count + 1).",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=8.0,
    )

    if panel_label is not None:
        ax.text(
            -0.18,
            1.06,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=10.0,
            fontweight="bold",
        )

    return im


def save_single_figure(
    counts,
    metrics,
    reference_label,
    out_png,
    out_pdf,
):
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    im = draw_contingency_plot_on_axis(
        ax,
        counts,
        metrics,
        reference_label,
    )

    cbar = fig.colorbar(
        im,
        ax=ax,
        fraction=0.046,
        pad=0.04,
    )
    cbar.set_label("log10(count + 1)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fig.subplots_adjust(
        left=0.18,
        right=0.88,
        bottom=0.28,
        top=0.96,
    )

    fig.savefig(
        out_pdf,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    fig.savefig(
        out_png,
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    plt.close(fig)


def save_combined_figure(
    radar_counts,
    radar_metrics,
    gauge_counts,
    gauge_metrics,
    out_png,
    out_pdf,
):
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11.4, 5.0),
    )

    draw_contingency_plot_on_axis(
        axes[0],
        radar_counts,
        radar_metrics,
        "Radar",
        panel_label="(a)",
    )
    draw_contingency_plot_on_axis(
        axes[1],
        gauge_counts,
        gauge_metrics,
        "Rain gauge",
        panel_label="(b)",
    )

    fig.subplots_adjust(
        left=0.06,
        right=0.98,
        bottom=0.28,
        top=0.95,
        wspace=0.32,
    )

    fig.savefig(
        out_pdf,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    fig.savefig(
        out_png,
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    plt.close(fig)


def write_latex_table(
    counts,
    metrics,
    reference_label,
    path,
):
    text = rf"""\begin{{table}}[t]
\centering
\caption{{Contingency table for Graf/RSTD rain detection against the {reference_label.lower()} reference.}}
\label{{tab:graf_rstd_{reference_label.lower().replace(' ', '_')}}}
\begin{{tabular}}{{lrr}}
\hline
 & Reference wet & Reference dry \\
\hline
CML wet & {counts['TP']:,} (TP) & {counts['FP']:,} (FP) \\
CML dry & {counts['FN']:,} (FN) & {counts['TN']:,} (TN) \\
\hline
\end{{tabular}}

\vspace{{0.5em}}
POD = {metrics['POD']:.3f},
FAR = {metrics['FAR']:.3f},
MCC = {metrics['MCC']:.3f},
CSI = {metrics['CSI']:.3f},
$N$ = {int(metrics['N']):,}.
\end{{table}}
"""
    path.write_text(text, encoding="utf-8")


# =============================================================================
# Main experiment
# =============================================================================
def run(args):
    started = time.perf_counter()

    out_root = Path(args.out_dir)
    radar_dir = out_root / "radar"
    gauge_dir = out_root / "gauge"

    radar_dir.mkdir(parents=True, exist_ok=True)
    gauge_dir.mkdir(parents=True, exist_ok=True)

    cml = xr.open_dataset(args.cml_nc)
    radar = xr.open_dataset(args.radar_nc)
    gauge = xr.open_dataset(args.gauge_nc)

    try:
        if "cml_id" not in cml.coords:
            raise KeyError("CML file has no cml_id coordinate.")
        if "sublink_id" not in cml.dims:
            raise KeyError("CML file has no sublink_id dimension.")
        if "time" not in cml.coords:
            raise KeyError("CML file has no time coordinate.")

        cml = cml.assign_coords(
            cml_id=cml["cml_id"].astype(str)
        )
        radar = radar.assign_coords(
            cml_id=radar["cml_id"].astype(str)
        )

        trsl_var = detect_trsl_source(
            cml,
            args.trsl_var,
        )
        radar_var = detect_var(
            radar,
            args.radar_var,
            ("rr", "radar_mmph", "R_radar"),
            "radar rain",
        )
        gauge_var = detect_var(
            gauge,
            args.gauge_var,
            ("rr", "rainfall", "precipitation"),
            "gauge rain",
        )

        cml_ids = np.asarray(
            cml["cml_id"].values
        ).astype(str)
        radar_ids = np.asarray(
            radar["cml_id"].values
        ).astype(str)

        # Radar uses all CMLs that have path-associated radar.
        radar_cml_ids = ordered_intersection(
            cml_ids,
            radar_ids,
        )

        if not radar_cml_ids:
            raise RuntimeError(
                "No common cml_id values between CML and radar."
            )

        # Gauge assignment is made for ALL CMLs, independently of radar
        # availability. Thus radar-unavailable links are not lost from gauge
        # evaluation.
        gauge_info = prepare_gauge_mapping(
            cml=cml,
            cml_ids=cml_ids,
            gauge=gauge,
            gauge_var=gauge_var,
            args=args,
        )

        gauge_mapping = gauge_info["mapping"]
        gauge_mapping.to_csv(
            gauge_dir / "closest_gauge_assignment.csv",
            index=False,
        )

        n_gauge_links = int(
            (gauge_info["assigned_positions"] >= 0).sum()
        )

        cml_times = pd.DatetimeIndex(
            cml["time"].values
        )

        # ---------------------------------------------------------------------
        # Radar time axis
        # ---------------------------------------------------------------------
        radar_da = radar[radar_var]
        radar_time_dim = detect_time_dim(radar_da)
        all_radar_times = pd.DatetimeIndex(
            radar[radar_time_dim].values
        )

        radar_time_keep = np.ones(
            len(all_radar_times),
            dtype=bool,
        )

        if args.tstart:
            radar_time_keep &= (
                all_radar_times >= pd.Timestamp(args.tstart)
            )
        if args.tend:
            radar_time_keep &= (
                all_radar_times <= pd.Timestamp(args.tend)
            )

        radar_times = all_radar_times[
            radar_time_keep
        ]
        radar_time_pos = np.flatnonzero(
            radar_time_keep
        )

        if len(radar_times) == 0:
            raise RuntimeError(
                "No radar timestamps remain after time filtering."
            )

        radar_interval_map = build_reference_intervals(
            cml_times,
            radar_times,
            args.radar_timestamp_label,
            args.radar_fallback_minutes,
        )

        # ---------------------------------------------------------------------
        # Gauge time axis
        # ---------------------------------------------------------------------
        gauge_rr = gauge_info["rr"]
        gauge_time_dim = gauge_info["time_dim"]
        gauge_station_dim = gauge_info["station_dim"]

        all_gauge_times = pd.DatetimeIndex(
            gauge[gauge_time_dim].values
        )

        gauge_time_keep = np.ones(
            len(all_gauge_times),
            dtype=bool,
        )

        if args.tstart:
            gauge_time_keep &= (
                all_gauge_times >= pd.Timestamp(args.tstart)
            )
        if args.tend:
            gauge_time_keep &= (
                all_gauge_times <= pd.Timestamp(args.tend)
            )

        gauge_times = all_gauge_times[
            gauge_time_keep
        ]
        gauge_time_pos = np.flatnonzero(
            gauge_time_keep
        )

        if len(gauge_times) == 0:
            raise RuntimeError(
                "No gauge timestamps remain after time filtering."
            )

        gauge_interval_map = build_reference_intervals(
            cml_times,
            gauge_times,
            args.gauge_timestamp_label,
            args.gauge_fallback_minutes,
        )

        # Cache only gauges actually assigned within 1 km.
        assigned = gauge_info["assigned_positions"]
        unique_gauges = np.unique(
            assigned[assigned >= 0]
        ).astype(np.int64)

        gauge_cache = {}
        for station_pos in unique_gauges:
            arr = gauge_rr.isel(
                {
                    gauge_station_dim: int(station_pos),
                    gauge_time_dim: gauge_time_pos,
                }
            ).values
            gauge_cache[int(station_pos)] = (
                np.asarray(arr, dtype=np.float64)
                * float(args.gauge_scale)
            )

        # ---------------------------------------------------------------------
        # CML index maps
        # ---------------------------------------------------------------------
        cml_pos = {
            str(v): i
            for i, v in enumerate(cml_ids)
        }
        radar_pos = {
            str(v): i
            for i, v in enumerate(radar_ids)
        }

        n_sublinks = int(
            cml.sizes["sublink_id"]
        )
        sublink_values = np.asarray(
            cml["sublink_id"].values
        )

        # Global results.
        radar_global = {
            "TP": 0,
            "TN": 0,
            "FP": 0,
            "FN": 0,
        }
        gauge_global = {
            "TP": 0,
            "TN": 0,
            "FP": 0,
            "FN": 0,
        }

        radar_rows = []
        gauge_rows = []
        threshold_rows = []

        print("=" * 100)
        print("GRAF/RSTD — RADAR + RAIN-GAUGE BENCHMARK")
        print("=" * 100)
        print(
            f"CML                         : {args.cml_nc}"
        )
        print(
            f"TRSL                        : {trsl_var or 'tx-rx'}"
        )
        print(
            f"Radar                       : {args.radar_nc}"
        )
        print(
            f"Gauge                       : {args.gauge_nc}"
        )
        print(
            f"All preprocessed CML links  : {len(cml_ids):,}"
        )
        print(
            f"Links with path radar       : {len(radar_cml_ids):,}"
        )
        print(
            f"Links with gauge <=1 km     : {n_gauge_links:,}"
        )
        print(
            f"Independent sublinks/link   : {n_sublinks}"
        )
        print(
            f"Radar wet threshold         : > {args.radar_wet_threshold} mm h-1"
        )
        print(
            f"Gauge wet threshold         : > {args.gauge_wet_threshold} mm h-1"
        )
        print(
            f"Gauge assignment            : nearest to link midpoint, <= {args.max_gauge_dist_km:g} km"
        )
        print(
            f"Radar native interval       : {float(radar_interval_map['step_ns'])/60e9:.1f} min"
        )
        print(
            f"Gauge native interval       : {float(gauge_interval_map['step_ns'])/60e9:.1f} min"
        )
        print("=" * 100)

        # Process ALL CMLs. Radar and gauge mappings are handled independently.
        for b0 in range(
            0,
            len(cml_ids),
            args.link_batch,
        ):
            b1 = min(
                b0 + args.link_batch,
                len(cml_ids),
            )

            ids_batch = cml_ids[b0:b1]

            print(
                f"links {b0+1:,}-{b1:,}/{len(cml_ids):,}",
                flush=True,
            )

            cpos = np.arange(
                b0,
                b1,
                dtype=np.int64,
            )

            if trsl_var is not None:
                da = cml[trsl_var].isel(
                    cml_id=cpos
                )
            else:
                da = (
                    cml["tx"].isel(cml_id=cpos)
                    - cml["rx"].isel(cml_id=cpos)
                )

            # IMPORTANT:
            # Use the underlying xarray Variable rather than DataArray.transpose().
            #
            # The publication NetCDF contains auxiliary byte/string coordinates
            # (for example polarization).  DataArray.transpose() also attempts to
            # transpose attached coordinates, and xarray's StackedBytesArray does
            # not implement transpose.  The TRSL data itself is numeric and can be
            # transposed safely at the Variable level without touching auxiliary
            # coordinates.
            trsl3 = da.variable.transpose(
                "time",
                "cml_id",
                "sublink_id",
            ).values

            values = np.asarray(
                trsl3,
                dtype=np.float64,
            ).reshape(
                len(cml_times),
                (b1 - b0) * n_sublinks,
            )

            (
                wet_1min,
                valid_1min,
                thresholds,
                q_values,
                finite_input,
            ) = graf_rstd_matrix(
                values=values,
                window=args.rstd_window,
                quantile=args.graf_quantile,
                factor=args.graf_factor,
                threshold_min_periods=args.threshold_min_periods,
                detection_min_periods=args.detection_min_periods,
            )

            pred_radar, valid_radar_pred = (
                aggregate_to_reference(
                    wet_1min,
                    valid_1min,
                    radar_interval_map,
                    args.min_cml_valid_fraction,
                    args.temporal_aggregation,
                )
            )

            pred_gauge, valid_gauge_pred = (
                aggregate_to_reference(
                    wet_1min,
                    valid_1min,
                    gauge_interval_map,
                    args.min_cml_valid_fraction,
                    args.temporal_aggregation,
                )
            )

            # -------------------------------------------------------------
            # References per physical link
            # -------------------------------------------------------------
            for local_link, cid in enumerate(
                ids_batch
            ):
                global_cml_pos = b0 + local_link

                # -------------------------
                # RADAR reference
                # -------------------------
                ridx = radar_pos.get(
                    str(cid),
                    None,
                )

                if ridx is not None:
                    radar_raw = radar_da.isel(
                        {
                            "cml_id": int(ridx),
                            radar_time_dim: radar_time_pos,
                        }
                    ).values

                    radar_rate = (
                        np.asarray(
                            radar_raw,
                            dtype=np.float64,
                        )
                        * float(args.radar_scale)
                    )

                    radar_ref_valid = (
                        np.isfinite(radar_rate)
                        & (radar_rate >= 0)
                    )
                    radar_ref_wet = (
                        radar_ref_valid
                        & (
                            radar_rate
                            > float(
                                args.radar_wet_threshold
                            )
                        )
                    )
                else:
                    radar_ref_valid = None
                    radar_ref_wet = None

                # -------------------------
                # GAUGE reference
                # -------------------------
                gauge_position = int(
                    assigned[
                        global_cml_pos
                    ]
                )

                if gauge_position >= 0:
                    gauge_rate = gauge_cache[
                        gauge_position
                    ]

                    gauge_ref_valid = (
                        np.isfinite(gauge_rate)
                        & (gauge_rate >= 0)
                    )
                    gauge_ref_wet = (
                        gauge_ref_valid
                        & (
                            gauge_rate
                            > float(
                                args.gauge_wet_threshold
                            )
                        )
                    )
                else:
                    gauge_ref_valid = None
                    gauge_ref_wet = None

                # -------------------------
                # Each sublink independent
                # -------------------------
                for s in range(
                    n_sublinks
                ):
                    obs = (
                        local_link
                        * n_sublinks
                        + s
                    )

                    threshold_rows.append(
                        {
                            "cml_id": str(cid),
                            "sublink_id": sublink_values[s],
                            "graf_q80_db": float(
                                q_values[obs]
                            ),
                            "graf_threshold_db": float(
                                thresholds[obs]
                            ),
                            "finite_input_samples": int(
                                finite_input[obs]
                            ),
                            "valid_detector_samples_1min": int(
                                valid_1min[:, obs].sum()
                            ),
                            "wet_detector_samples_1min": int(
                                wet_1min[:, obs].sum()
                            ),
                        }
                    )

                    # Radar
                    if radar_ref_valid is not None:
                        tp, tn, fp, fn = (
                            counts_for_column(
                                pred_radar[:, obs],
                                valid_radar_pred[:, obs],
                                radar_ref_wet,
                                radar_ref_valid,
                            )
                        )
                    else:
                        tp = tn = fp = fn = 0

                    m = metrics_from_counts(
                        tp,
                        tn,
                        fp,
                        fn,
                    )

                    radar_rows.append(
                        {
                            "cml_id": str(cid),
                            "sublink_id": sublink_values[s],
                            "has_radar_reference": int(
                                ridx is not None
                            ),
                            "TP": tp,
                            "TN": tn,
                            "FP": fp,
                            "FN": fn,
                            **m,
                        }
                    )

                    radar_global["TP"] += tp
                    radar_global["TN"] += tn
                    radar_global["FP"] += fp
                    radar_global["FN"] += fn

                    # Gauge
                    if gauge_ref_valid is not None:
                        tp, tn, fp, fn = (
                            counts_for_column(
                                pred_gauge[:, obs],
                                valid_gauge_pred[:, obs],
                                gauge_ref_wet,
                                gauge_ref_valid,
                            )
                        )
                    else:
                        tp = tn = fp = fn = 0

                    m = metrics_from_counts(
                        tp,
                        tn,
                        fp,
                        fn,
                    )

                    gauge_rows.append(
                        {
                            "cml_id": str(cid),
                            "sublink_id": sublink_values[s],
                            "gauge_position": gauge_position,
                            "gauge_within_1km": int(
                                gauge_position >= 0
                            ),
                            "TP": tp,
                            "TN": tn,
                            "FP": fp,
                            "FN": fn,
                            **m,
                        }
                    )

                    gauge_global["TP"] += tp
                    gauge_global["TN"] += tn
                    gauge_global["FP"] += fp
                    gauge_global["FN"] += fn

            del (
                trsl3,
                values,
                wet_1min,
                valid_1min,
                pred_radar,
                valid_radar_pred,
                pred_gauge,
                valid_gauge_pred,
            )

        radar_metrics = metrics_from_counts(
            radar_global["TP"],
            radar_global["TN"],
            radar_global["FP"],
            radar_global["FN"],
        )

        gauge_metrics = metrics_from_counts(
            gauge_global["TP"],
            gauge_global["TN"],
            gauge_global["FP"],
            gauge_global["FN"],
        )

        # ---------------------------------------------------------------------
        # Save tables
        # ---------------------------------------------------------------------
        pd.DataFrame(
            [radar_global]
        ).to_csv(
            radar_dir
            / "graf_rstd_vs_radar_contingency.csv",
            index=False,
        )

        pd.DataFrame(
            [{**radar_global, **radar_metrics}]
        ).to_csv(
            radar_dir
            / "graf_rstd_vs_radar_metrics.csv",
            index=False,
        )

        pd.DataFrame(
            radar_rows
        ).to_csv(
            radar_dir
            / "graf_rstd_vs_radar_per_sublink.csv",
            index=False,
        )

        pd.DataFrame(
            [gauge_global]
        ).to_csv(
            gauge_dir
            / "graf_rstd_vs_gauge_contingency.csv",
            index=False,
        )

        pd.DataFrame(
            [{**gauge_global, **gauge_metrics}]
        ).to_csv(
            gauge_dir
            / "graf_rstd_vs_gauge_metrics.csv",
            index=False,
        )

        pd.DataFrame(
            gauge_rows
        ).to_csv(
            gauge_dir
            / "graf_rstd_vs_gauge_per_sublink.csv",
            index=False,
        )

        pd.DataFrame(
            threshold_rows
        ).to_csv(
            out_root
            / "graf_rstd_thresholds.csv",
            index=False,
        )

        # ---------------------------------------------------------------------
        # HESS figures
        # ---------------------------------------------------------------------
        radar_png = (
            radar_dir
            / "graf_rstd_vs_radar_contingency_hess.png"
        )
        radar_pdf = (
            radar_dir
            / "graf_rstd_vs_radar_contingency_hess.pdf"
        )
        gauge_png = (
            gauge_dir
            / "graf_rstd_vs_gauge_contingency_hess.png"
        )
        gauge_pdf = (
            gauge_dir
            / "graf_rstd_vs_gauge_contingency_hess.pdf"
        )

        save_single_figure(
            radar_global,
            radar_metrics,
            "Radar",
            radar_png,
            radar_pdf,
        )

        save_single_figure(
            gauge_global,
            gauge_metrics,
            "Rain gauge",
            gauge_png,
            gauge_pdf,
        )

        combined_png = (
            out_root
            / "graf_rstd_radar_gauge_contingency_hess.png"
        )
        combined_pdf = (
            out_root
            / "graf_rstd_radar_gauge_contingency_hess.pdf"
        )

        save_combined_figure(
            radar_global,
            radar_metrics,
            gauge_global,
            gauge_metrics,
            combined_png,
            combined_pdf,
        )

        write_latex_table(
            radar_global,
            radar_metrics,
            "Radar",
            radar_dir
            / "graf_rstd_vs_radar_contingency_table.tex",
        )

        write_latex_table(
            gauge_global,
            gauge_metrics,
            "Rain gauge",
            gauge_dir
            / "graf_rstd_vs_gauge_contingency_table.tex",
        )

        # ---------------------------------------------------------------------
        # Run summary
        # ---------------------------------------------------------------------
        summary = {
            "cml_file": args.cml_nc,
            "radar_file": args.radar_nc,
            "gauge_file": args.gauge_nc,
            "trsl_source": trsl_var or "tx-rx",
            "all_preprocessed_cml_links": int(
                len(cml_ids)
            ),
            "links_with_path_radar": int(
                len(radar_cml_ids)
            ),
            "links_with_gauge_within_1km": int(
                n_gauge_links
            ),
            "independent_sublinks_per_link": int(
                n_sublinks
            ),
            "sublinks_averaged": False,
            "rstd_window_minutes": int(
                args.rstd_window
            ),
            "graf_quantile": float(
                args.graf_quantile
            ),
            "graf_factor": float(
                args.graf_factor
            ),
            "threshold_min_periods": int(
                args.threshold_min_periods
            ),
            "detection_min_periods": int(
                args.detection_min_periods
            ),
            "temporal_aggregation": str(
                args.temporal_aggregation
            ),
            "min_cml_valid_fraction": float(
                args.min_cml_valid_fraction
            ),
            "reference_filling": "none",
            "paired_valid_only": True,
            "invalid_cml_treated_as_dry": False,
            "invalid_reference_treated_as_dry": False,
            "radar": {
                "variable": radar_var,
                "scale_to_mmph": float(
                    args.radar_scale
                ),
                "wet_threshold_mmph": float(
                    args.radar_wet_threshold
                ),
                "timestamp_label": str(
                    args.radar_timestamp_label
                ),
                "native_minutes": float(
                    radar_interval_map["step_ns"]
                ) / 60e9,
                "counts": radar_global,
                "metrics": radar_metrics,
            },
            "gauge": {
                "variable": gauge_var,
                "scale_to_mmph": float(
                    args.gauge_scale
                ),
                "wet_threshold_mmph": float(
                    args.gauge_wet_threshold
                ),
                "timestamp_label": str(
                    args.gauge_timestamp_label
                ),
                "native_minutes": float(
                    gauge_interval_map["step_ns"]
                ) / 60e9,
                "assignment": (
                    "nearest gauge to CML link midpoint"
                ),
                "max_distance_km": float(
                    args.max_gauge_dist_km
                ),
                "counts": gauge_global,
                "metrics": gauge_metrics,
            },
            "elapsed_seconds": float(
                time.perf_counter() - started
            ),
        }

        (
            out_root
            / "run_summary.json"
        ).write_text(
            json.dumps(
                summary,
                indent=2,
                allow_nan=True,
            ),
            encoding="utf-8",
        )

        print("\n" + "=" * 100)
        print("RADAR")
        print("=" * 100)
        print(
            f"TP={radar_global['TP']:,}  "
            f"FP={radar_global['FP']:,}  "
            f"FN={radar_global['FN']:,}  "
            f"TN={radar_global['TN']:,}"
        )
        print(
            f"N={radar_metrics['N']:,}  "
            f"POD={radar_metrics['POD']:.4f}  "
            f"FAR={radar_metrics['FAR']:.4f}  "
            f"MCC={radar_metrics['MCC']:.4f}  "
            f"CSI={radar_metrics['CSI']:.4f}"
        )

        print("\n" + "=" * 100)
        print("RAIN GAUGE <= 1 KM FROM LINK MIDPOINT")
        print("=" * 100)
        print(
            f"TP={gauge_global['TP']:,}  "
            f"FP={gauge_global['FP']:,}  "
            f"FN={gauge_global['FN']:,}  "
            f"TN={gauge_global['TN']:,}"
        )
        print(
            f"N={gauge_metrics['N']:,}  "
            f"POD={gauge_metrics['POD']:.4f}  "
            f"FAR={gauge_metrics['FAR']:.4f}  "
            f"MCC={gauge_metrics['MCC']:.4f}  "
            f"CSI={gauge_metrics['CSI']:.4f}"
        )

        print("\nSaved:")
        print(f"  {combined_png}")
        print(f"  {combined_pdf}")
        print(f"  {radar_png}")
        print(f"  {gauge_png}")
        print(
            f"Elapsed: "
            f"{(time.perf_counter()-started)/60.0:.1f} min"
        )

    finally:
        cml.close()
        radar.close()
        gauge.close()


# =============================================================================
# CLI
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Graf/RSTD wet-dry detection compared independently "
            "against path radar and nearest <=1-km rain gauges."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument(
        "--cml-nc",
        required=True,
    )
    p.add_argument(
        "--radar-nc",
        required=True,
    )
    p.add_argument(
        "--gauge-nc",
        required=True,
    )
    p.add_argument(
        "--out-dir",
        required=True,
    )

    p.add_argument(
        "--trsl-var",
        default=None,
        help=(
            "Auto-detect order: trsl, trsl_qc, then tx-rx."
        ),
    )

    # Detector.
    p.add_argument(
        "--rstd-window",
        type=int,
        default=60,
    )
    p.add_argument(
        "--graf-quantile",
        type=float,
        default=0.80,
    )
    p.add_argument(
        "--graf-factor",
        type=float,
        default=1.12,
    )
    p.add_argument(
        "--threshold-min-periods",
        type=int,
        default=60,
    )
    p.add_argument(
        "--detection-min-periods",
        type=int,
        default=45,
    )

    p.add_argument(
        "--temporal-aggregation",
        choices=[
            "any",
            "majority",
            "all",
        ],
        default="any",
    )
    p.add_argument(
        "--min-cml-valid-fraction",
        type=float,
        default=1.0,
    )

    # Shared optional time window.
    p.add_argument(
        "--tstart",
        default=None,
    )
    p.add_argument(
        "--tend",
        default=None,
    )

    # Radar.
    p.add_argument(
        "--radar-var",
        default=None,
    )
    p.add_argument(
        "--radar-scale",
        type=float,
        default=12.0,
    )
    p.add_argument(
        "--radar-wet-threshold",
        type=float,
        default=0.1,
    )
    p.add_argument(
        "--radar-timestamp-label",
        choices=[
            "start",
            "end",
            "center",
        ],
        default="start",
    )
    p.add_argument(
        "--radar-fallback-minutes",
        type=float,
        default=5.0,
    )

    # Gauge.
    p.add_argument(
        "--gauge-var",
        default=None,
    )
    p.add_argument(
        "--gauge-scale",
        type=float,
        default=6.0,
    )
    p.add_argument(
        "--gauge-wet-threshold",
        type=float,
        default=0.1,
    )
    p.add_argument(
        "--gauge-timestamp-label",
        choices=[
            "start",
            "end",
            "center",
        ],
        default="end",
        help=(
            "CHMI publication dataset timestamps are treated as "
            "the end of the preceding 10-min accumulation interval."
        ),
    )
    p.add_argument(
        "--gauge-fallback-minutes",
        type=float,
        default=10.0,
    )
    p.add_argument(
        "--max-gauge-dist-km",
        type=float,
        default=1.0,
    )
    p.add_argument(
        "--gauge-lon",
        default=None,
    )
    p.add_argument(
        "--gauge-lat",
        default=None,
    )
    p.add_argument(
        "--gauge-id",
        default=None,
    )

    # CML endpoint coordinates.
    p.add_argument(
        "--site-a-lon",
        default="site_a_longitude",
    )
    p.add_argument(
        "--site-a-lat",
        default="site_a_latitude",
    )
    p.add_argument(
        "--site-b-lon",
        default="site_b_longitude",
    )
    p.add_argument(
        "--site-b-lat",
        default="site_b_latitude",
    )

    # Memory control.
    p.add_argument(
        "--link-batch",
        type=int,
        default=8,
    )

    args = p.parse_args()

    if not 0.0 < args.graf_quantile < 1.0:
        p.error(
            "--graf-quantile must be in (0,1)"
        )
    if not 0.0 < args.min_cml_valid_fraction <= 1.0:
        p.error(
            "--min-cml-valid-fraction must be in (0,1]"
        )
    if (
        args.threshold_min_periods
        > args.rstd_window
    ):
        p.error(
            "--threshold-min-periods cannot exceed --rstd-window"
        )
    if (
        args.detection_min_periods
        > args.rstd_window
    ):
        p.error(
            "--detection-min-periods cannot exceed --rstd-window"
        )
    if args.max_gauge_dist_km <= 0:
        p.error(
            "--max-gauge-dist-km must be >0"
        )
    if args.link_batch < 1:
        p.error(
            "--link-batch must be >=1"
        )

    return args


if __name__ == "__main__":
    run(parse_args())
