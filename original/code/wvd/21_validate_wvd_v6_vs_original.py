#!/usr/bin/env python3
"""
21_validate_wvd_v6_vs_original.py
=================================

Strict CHMI validation of the original fixed internal-anchor WVD retrieval
against the V6 causal progressive-baseline retrieval.

PRIMARY PRINCIPLE
-----------------
The comparison is PAIRED. Original and V6 are evaluated on EXACTLY the same
CML/CHMI 10-min observations for each directional sublink and wet/dry branch.
This prevents missing V6 values from making V6 look artificially better.

Validation population
---------------------
- E-band directional sublinks already present in the retrieval products:
      71-76 GHz OR 81-86 GHz
- physical path length 2-10 km
- nearest CHMI station with BOTH H and T available <= 5 km from CML path
- exact CHMI 10-min timestamps; no RH or T interpolation
- CML 5-min -> 10-min requires exactly two finite 5-min WVD values AND
  retrieval_status == 1 for BOTH samples
- the 10-min CML value is the arithmetic mean of those two 5-min values
- CML pair [t, t+5 min] is labelled t+10 min and matched exactly to CHMI

CHMI water-vapour-density reference
-----------------------------------
    es = 6.112 * exp(17.67*T / (T + 243.5))       [hPa]
    e  = RH/100 * es                              [hPa]
    rho_v = 216.7 * e / (T + 273.15)             [g m-3]

IMPORTANT
---------
CHMI H/T are used ONLY here for validation. They are not used to calibrate,
fit, tune, reset, or otherwise modify either retrieval.

Outputs
-------
<out-dir>/
  station_assignments.csv
  pooled_metrics_common.csv
  coverage_summary.csv
  per_sublink_metrics_common.csv
  stratified_summary.csv
  min_length_sensitivity.csv

  pooled/
    graf_original_vs_v6_hexbin.png/pdf
    gauge_original_vs_v6_hexbin.png/pdf

  per_sublink/
    graf_metrics_original_vs_v6.png/pdf
    gauge_metrics_original_vs_v6.png/pdf

  by_length/
    graf_metrics_by_length.png/pdf
    gauge_metrics_by_length.png/pdf

  by_frequency/
    graf_metrics_by_frequency.png/pdf
    gauge_metrics_by_frequency.png/pdf

  sensitivity/
    graf_min_length_sensitivity.png/pdf
    gauge_min_length_sensitivity.png/pdf

The length/frequency figures are diagnostics. Do NOT choose/drop individual
links using CHMI performance metrics. Any final population restriction should
be defined from independent CML properties / physical criteria.

Example
-------
python 21_validate_wvd_v6_vs_original.py \
  --original-nc CZ_CML_WVD_60_100GHz_5MIN_APR_OCT_2021.nc \
  --original-diagnostics-csv CZ_CML_WVD_60_100GHz_DIAGNOSTICS.csv \
  --v6-nc CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc \
  --v6-diagnostics-csv CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv \
  --raw-cml-nc CZ_CML_RAW_1MIN_20200101_20211231.nc \
  --humidity-nc CHMI_10MIN_HUMIDITY_2021.nc \
  --temperature-nc CHMI_10MIN_TEMPERATURE_2021.nc \
  --out-dir WVD_V6_VALIDATION
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import xarray as xr


METHODS = ("graf", "gauge")
MODEL_LABELS = ("Original", "V6")
LENGTH_BINS = [2.0, 4.0, 6.0, 8.0, 10.000001]
LENGTH_LABELS = ["2–4", "4–6", "6–8", "8–10"]
FREQ_BANDS = ("71–76 GHz", "81–86 GHz")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--original-nc", required=True)
    p.add_argument("--original-diagnostics-csv", required=True)
    p.add_argument("--v6-nc", required=True)
    p.add_argument("--v6-diagnostics-csv", required=True)
    p.add_argument("--raw-cml-nc", required=True)
    p.add_argument("--humidity-nc", required=True)
    p.add_argument("--temperature-nc", required=True)
    p.add_argument("--out-dir", type=Path, required=True)

    p.add_argument("--max-station-distance-km", type=float, default=5.0)
    p.add_argument("--min-pairs-per-sublink", type=int, default=100)
    p.add_argument("--scatter-max-gm3", type=float, default=60.0)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def to_strings(values):
    a = np.asarray(values)
    out = []
    for x in a.reshape(-1):
        if isinstance(x, (bytes, np.bytes_)):
            out.append(x.decode("utf-8", errors="replace").strip())
        else:
            out.append(str(x).strip())
    return np.asarray(out, dtype=object)


def point_to_segment_distance_km(lat_p, lon_p, lat1, lon1, lat2, lon2):
    """Local equirectangular point-to-segment distance; sufficient at <~10 km."""
    vals = np.asarray([lat_p, lon_p, lat1, lon1, lat2, lon2], dtype=float)
    if not np.all(np.isfinite(vals)):
        return np.nan

    lat0 = float(np.mean([lat_p, lat1, lat2]))
    R = 6371.0
    c = math.cos(math.radians(lat0))

    def xy(lat, lon):
        return np.array([
            math.radians(lon) * R * c,
            math.radians(lat) * R,
        ], dtype=float)

    p = xy(lat_p, lon_p)
    a = xy(lat1, lon1)
    b = xy(lat2, lon2)
    ab = b - a
    den = float(np.dot(ab, ab))
    if den <= 0:
        return float(np.linalg.norm(p - a))
    u = float(np.dot(p - a, ab) / den)
    u = min(1.0, max(0.0, u))
    q = a + u * ab
    return float(np.linalg.norm(p - q))


def metric_block(ref, est):
    ref = np.asarray(ref, dtype=float)
    est = np.asarray(est, dtype=float)
    m = np.isfinite(ref) & np.isfinite(est)
    ref = ref[m]
    est = est[m]
    n = int(ref.size)
    if n == 0:
        return dict(N=0, PCC=np.nan, RMSE=np.nan, Bias=np.nan, MAE=np.nan)
    err = est - ref
    pcc = (
        float(np.corrcoef(ref, est)[0, 1])
        if n >= 2 and np.std(ref) > 0 and np.std(est) > 0
        else np.nan
    )
    return dict(
        N=n,
        PCC=pcc,
        RMSE=float(np.sqrt(np.mean(err ** 2))),
        Bias=float(np.mean(err)),
        MAE=float(np.mean(np.abs(err))),
    )


def rh_t_to_wvd(rh_pct, t_c):
    rh = np.asarray(rh_pct, dtype=float)
    t = np.asarray(t_c, dtype=float)
    es = 6.112 * np.exp(17.67 * t / (t + 243.5))
    e = (rh / 100.0) * es
    rho = 216.7 * e / (t + 273.15)
    bad = (
        ~np.isfinite(rh) | ~np.isfinite(t)
        | (rh < 0) | (rh > 105)
        | (t <= -100) | (t >= 70)
    )
    rho[bad] = np.nan
    return rho


def decode_retrieval_time(ds):
    for name in ("time", "time_5min"):
        if name in ds.coords or name in ds.variables:
            return pd.DatetimeIndex(pd.to_datetime(ds[name].values))
    raise RuntimeError("No time/time_5min coordinate in retrieval NetCDF")


def strict_5_to_10(time5, wvd5, status5):
    """
    Pair exactly two consecutive 5-min samples:
       00:00 + 00:05 -> label 00:10
       00:10 + 00:15 -> label 00:20
    Both samples must be finite and status==1.
    """
    t = pd.DatetimeIndex(time5)
    x = np.asarray(wvd5, dtype=float)
    s = np.asarray(status5)

    if len(t) != len(x) or len(t) != len(s):
        raise RuntimeError("5-min value/status/time length mismatch")
    if len(t) < 2:
        return pd.DatetimeIndex([]), np.empty(0, float)

    # We require the product to start on a 10-min boundary and be regular.
    dt = np.diff(t.view("i8"))
    expected = pd.Timedelta(minutes=5).value
    if np.any(dt != expected):
        raise RuntimeError("Retrieval time grid is not an exact regular 5-min grid")
    if t[0].minute % 10 != 0 or t[0].second != 0:
        raise RuntimeError(
            f"5-min grid starts at {t[0]}, not an exact 10-min boundary"
        )

    n = (len(t) // 2) * 2
    x0, x1 = x[:n:2], x[1:n:2]
    s0, s1 = s[:n:2], s[1:n:2]
    ok = (
        np.isfinite(x0) & np.isfinite(x1)
        & (s0 == 1) & (s1 == 1)
    )

    y = np.full(n // 2, np.nan, dtype=float)
    y[ok] = 0.5 * (x0[ok] + x1[ok])
    t10 = t[:n:2] + pd.Timedelta(minutes=10)
    return pd.DatetimeIndex(t10), y


def length_bin(length_km):
    x = float(length_km)
    for lo, hi, lab in zip(LENGTH_BINS[:-1], LENGTH_BINS[1:], LENGTH_LABELS):
        if lo <= x < hi:
            return lab
    return np.nan


def frequency_band(f):
    f = float(f)
    if 71.0 <= f <= 76.0:
        return "71–76 GHz"
    if 81.0 <= f <= 86.0:
        return "81–86 GHz"
    return np.nan


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 4,
        "ytick.major.size": 4,
    })


# ---------------------------------------------------------------------------
# Metadata / station assignment
# ---------------------------------------------------------------------------

def load_diag(path):
    d = pd.read_csv(path)
    required = {"series_index", "cml_id", "sublink_id", "frequency_ghz", "length_km"}
    missing = required - set(d.columns)
    if missing:
        raise RuntimeError(f"{path} missing columns: {sorted(missing)}")
    d = d.copy()
    d["series_index"] = pd.to_numeric(d["series_index"], errors="raise").astype(int)
    d["cml_id"] = d["cml_id"].astype(str)
    d["sublink_id"] = d["sublink_id"].astype(str)
    if d.duplicated(["cml_id", "sublink_id"]).any():
        raise RuntimeError(f"Duplicate cml_id/sublink_id rows in {path}")
    return d


def raw_geometry(raw_path):
    ds = xr.open_dataset(raw_path, cache=False)
    try:
        cids = to_strings(ds["cml_id"].values)
        out = {}
        for i, cid in enumerate(cids):
            out[str(cid)] = {
                "site_a_latitude": float(ds["site_a_latitude"].values[i]),
                "site_a_longitude": float(ds["site_a_longitude"].values[i]),
                "site_b_latitude": float(ds["site_b_latitude"].values[i]),
                "site_b_longitude": float(ds["site_b_longitude"].values[i]),
            }
        return out
    finally:
        ds.close()


def common_met_stations(hum, temp):
    for ds, label in ((hum, "humidity"), (temp, "temperature")):
        for name in ("station", "lat", "lon", "time"):
            if name not in ds.coords and name not in ds.variables:
                raise RuntimeError(f"{label} file missing `{name}`")

    h_ids = to_strings(hum["station"].values)
    t_ids = to_strings(temp["station"].values)
    h_map = {x: i for i, x in enumerate(h_ids)}
    t_map = {x: i for i, x in enumerate(t_ids)}
    ids = [x for x in h_ids if x in t_map]

    rows = []
    for sid in ids:
        hi = h_map[sid]
        ti = t_map[sid]
        hlat = float(hum["lat"].values[hi])
        hlon = float(hum["lon"].values[hi])
        tlat = float(temp["lat"].values[ti])
        tlon = float(temp["lon"].values[ti])

        # Same station ID should be spatially identical. Use humidity coords
        # for path matching but retain mismatch as a QC diagnostic.
        coord_diff = point_to_segment_distance_km(
            tlat, tlon, hlat, hlon, hlat, hlon
        )
        rows.append({
            "station_id": sid,
            "hum_index": hi,
            "temp_index": ti,
            "lat": hlat,
            "lon": hlon,
            "temp_coord_offset_km": coord_diff,
        })
    return pd.DataFrame(rows)


def assign_stations(cml_ids, geom, stations, max_km):
    rows = []
    slat = stations["lat"].to_numpy(float)
    slon = stations["lon"].to_numpy(float)
    sids = stations["station_id"].astype(str).to_numpy()

    for cid in cml_ids:
        g = geom[str(cid)]
        dist = np.asarray([
            point_to_segment_distance_km(
                la, lo,
                g["site_a_latitude"], g["site_a_longitude"],
                g["site_b_latitude"], g["site_b_longitude"],
            )
            for la, lo in zip(slat, slon)
        ])
        if not np.any(np.isfinite(dist)):
            best_i = -1
            best_d = np.nan
            best_id = ""
        else:
            best_i = int(np.nanargmin(dist))
            best_d = float(dist[best_i])
            best_id = str(sids[best_i])

        rows.append({
            "cml_id": str(cid),
            "station_id": best_id,
            "station_distance_to_path_km": best_d,
            "within_validation_radius": bool(np.isfinite(best_d) and best_d <= max_km),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def savefig(fig, stem: Path, dpi):
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_pooled_hexbin(method, pooled, metrics_rows, out_dir, scatter_max, dpi):
    """
    Two panels on the EXACT SAME paired observations:
      Original vs CHMI | V6 vs CHMI
    """
    ref = pooled[method]["ref"]
    old = pooled[method]["old"]
    v6 = pooled[method]["v6"]

    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4.45), sharex=True, sharey=True)
    hbs = []
    for ax, est, model in zip(axes, (old, v6), MODEL_LABELS):
        hb = ax.hexbin(
            ref, est,
            gridsize=95,
            extent=(0, scatter_max, 0, scatter_max),
            mincnt=1,
            bins="log",
            cmap="magma",
            linewidths=0,
        )
        hbs.append(hb)
        ax.plot([0, scatter_max], [0, scatter_max], "--", lw=1.1, color="0.35")
        ax.set_xlim(0, scatter_max)
        ax.set_ylim(0, scatter_max)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.17, linewidth=0.6)

        mr = next(
            r for r in metrics_rows
            if r["method"] == method and r["model"] == model
        )
        txt = (
            f"N = {mr['N']:,}\n"
            f"PCC = {mr['PCC']:.3f}\n"
            f"RMSE = {mr['RMSE']:.2f}\n"
            f"Bias = {mr['Bias']:+.2f}"
        )
        ax.text(
            0.04, 0.96, txt, transform=ax.transAxes,
            ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.92),
        )
        ax.set_title(model)

    axes[0].set_ylabel(r"CML WVD [g m$^{-3}$]")
    for ax in axes:
        ax.set_xlabel(r"CHMI WVD [g m$^{-3}$]")

    fig.suptitle(
        f"{method.upper()} — paired strict 10-min validation",
        fontsize=12.5, fontweight="bold"
    )
    # Synchronize color norm for direct visual comparison.
    vmax = max(float(np.nanmax(h.get_array())) for h in hbs if h.get_array().size)
    norm = LogNorm(vmin=1, vmax=max(2.0, vmax))
    for h in hbs:
        h.set_norm(norm)
    cb = fig.colorbar(hbs[-1], ax=axes, fraction=0.035, pad=0.025)
    cb.set_label("Count")
    fig.tight_layout(rect=[0, 0, 0.95, 0.94])
    savefig(fig, out_dir / f"{method}_original_vs_v6_hexbin", dpi)


def jitter(n, width=0.075, seed=0):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, int(n))


def plot_overall_metric_boxes(method, df, out_dir, dpi):
    d = df[(df["method"] == method) & (df["eligible_metrics"] == 1)].copy()

    metric_specs = [
        ("PCC", "PCC"),
        ("RMSE", r"RMSE [g m$^{-3}$]"),
        ("Bias", r"Bias [g m$^{-3}$]"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(11.3, 4.1))
    for ai, (metric, ylabel) in enumerate(metric_specs):
        ax = axes[ai]
        vals_old = d[f"original_{metric.lower()}"].to_numpy(float)
        vals_v6 = d[f"v6_{metric.lower()}"].to_numpy(float)
        arrays = [vals_old[np.isfinite(vals_old)], vals_v6[np.isfinite(vals_v6)]]

        bp = ax.boxplot(
            arrays, positions=[1, 2], widths=0.48,
            showfliers=False, patch_artist=True,
            medianprops=dict(color="black", linewidth=1.3),
        )
        for patch, fc in zip(bp["boxes"], ("#d9d9d9", "#b9d7ea")):
            patch.set_facecolor(fc)
            patch.set_edgecolor("0.3")

        # Paired thin lines retain link identity.
        both = np.isfinite(vals_old) & np.isfinite(vals_v6)
        for a, b in zip(vals_old[both], vals_v6[both]):
            ax.plot([1, 2], [a, b], color="0.8", lw=0.45, alpha=0.35, zorder=1)

        ax.scatter(
            1 + jitter(np.isfinite(vals_old).sum(), seed=10 + ai),
            vals_old[np.isfinite(vals_old)],
            s=10, alpha=0.42, edgecolors="none", color="#555555", zorder=2
        )
        ax.scatter(
            2 + jitter(np.isfinite(vals_v6).sum(), seed=20 + ai),
            vals_v6[np.isfinite(vals_v6)],
            s=10, alpha=0.42, edgecolors="none", color="#2471a3", zorder=2
        )
        ax.set_xticks([1, 2], MODEL_LABELS)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.18, linewidth=0.6)
        if metric == "PCC":
            ax.axhline(0, color="0.7", lw=0.7)
        if metric == "Bias":
            ax.axhline(0, color="0.45", lw=0.8)

    fig.suptitle(
        f"{method.upper()} — per-sublink paired performance "
        f"(N≥{int(d['min_pairs_rule'].iloc[0]) if len(d) else ''})",
        fontsize=12.3, fontweight="bold"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    savefig(fig, out_dir / f"{method}_metrics_original_vs_v6", dpi)


def grouped_positions(n_groups):
    centers = np.arange(1, n_groups + 1, dtype=float)
    return centers - 0.18, centers + 0.18, centers


def plot_stratified(method, df, category_col, categories, out_stem, xlabel, dpi):
    d = df[(df["method"] == method) & (df["eligible_metrics"] == 1)].copy()
    specs = [
        ("PCC", "PCC"),
        ("RMSE", r"RMSE [g m$^{-3}$]"),
        ("Bias", r"Bias [g m$^{-3}$]"),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(8.8, 9.0), sharex=True)
    po, pv, centers = grouped_positions(len(categories))

    for ai, (metric, ylabel) in enumerate(specs):
        ax = axes[ai]
        for gi, cat in enumerate(categories):
            g = d[d[category_col] == cat]
            old = g[f"original_{metric.lower()}"].to_numpy(float)
            v6 = g[f"v6_{metric.lower()}"].to_numpy(float)
            old = old[np.isfinite(old)]
            v6 = v6[np.isfinite(v6)]

            if old.size:
                bp = ax.boxplot(
                    [old], positions=[po[gi]], widths=0.28,
                    showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", lw=1.1),
                )
                bp["boxes"][0].set_facecolor("#d9d9d9")
                bp["boxes"][0].set_edgecolor("0.35")
                ax.scatter(
                    po[gi] + jitter(old.size, 0.045, 1000 + ai*100 + gi),
                    old, s=8, alpha=0.38, color="#555555", edgecolors="none"
                )

            if v6.size:
                bp = ax.boxplot(
                    [v6], positions=[pv[gi]], widths=0.28,
                    showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", lw=1.1),
                )
                bp["boxes"][0].set_facecolor("#b9d7ea")
                bp["boxes"][0].set_edgecolor("0.35")
                ax.scatter(
                    pv[gi] + jitter(v6.size, 0.045, 2000 + ai*100 + gi),
                    v6, s=8, alpha=0.38, color="#2471a3", edgecolors="none"
                )

        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.18, linewidth=0.6)
        if metric in ("PCC", "Bias"):
            ax.axhline(0, color="0.6", lw=0.7)

    axes[-1].set_xticks(centers, categories)
    axes[-1].set_xlabel(xlabel)
    handles = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#d9d9d9",
               markeredgecolor="0.35", label="Original", markersize=9),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#b9d7ea",
               markeredgecolor="0.35", label="V6", markersize=9),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False)
    fig.suptitle(
        f"{method.upper()} — per-sublink metrics by {xlabel.lower()}",
        fontsize=12.3, fontweight="bold", y=0.985
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    savefig(fig, out_stem, dpi)


def plot_length_sensitivity(method, sens, out_dir, dpi):
    d = sens[sens["method"] == method].copy()
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.9))

    specs = [
        ("median_PCC", "Median per-sublink PCC"),
        ("median_RMSE", r"Median per-sublink RMSE [g m$^{-3}$]"),
        ("median_abs_Bias", r"Median |bias| [g m$^{-3}$]"),
    ]
    for ax, (col, ylabel) in zip(axes, specs):
        for model, marker in (("Original", "o"), ("V6", "s")):
            g = d[d["model"] == model]
            ax.plot(
                g["min_length_km"], g[col],
                marker=marker, lw=1.4, ms=4.5, label=model
            )
        ax.set_xlabel("Minimum retained path length [km]")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.18, linewidth=0.6)

    axes[0].legend(frameon=False)
    fig.suptitle(
        f"{method.upper()} — length-threshold sensitivity (diagnostic only)",
        fontsize=12.2, fontweight="bold"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    savefig(fig, out_dir / f"{method}_min_length_sensitivity", dpi)


# ---------------------------------------------------------------------------
# Main validation
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    style()

    out_dir = args.out_dir
    pooled_dir = out_dir / "pooled"
    per_dir = out_dir / "per_sublink"
    len_dir = out_dir / "by_length"
    freq_dir = out_dir / "by_frequency"
    sens_dir = out_dir / "sensitivity"
    for p in (out_dir, pooled_dir, per_dir, len_dir, freq_dir, sens_dir):
        p.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("STRICT PAIRED WVD VALIDATION — ORIGINAL vs V6")
    print("=" * 100)
    print("CHMI humidity/temperature are VALIDATION ONLY; no calibration/tuning.")
    print(f"Station radius           : <= {args.max_station_distance_km:g} km to CML path")
    print("CML aggregation          : exact two valid 5-min values -> one 10-min value")
    print("retrieval_status         : both 5-min samples must equal 1")
    print("CHMI matching            : exact 10-min timestamp; no interpolation")
    print(f"Per-sublink metric rule  : N >= {args.min_pairs_per_sublink}")
    print("Primary comparison       : paired common Original/V6 observations")
    print()

    old_diag = load_diag(args.original_diagnostics_csv)
    v6_diag = load_diag(args.v6_diagnostics_csv)

    old_map = {
        (r.cml_id, r.sublink_id): int(r.series_index)
        for r in old_diag.itertuples()
    }
    v6_map = {
        (r.cml_id, r.sublink_id): int(r.series_index)
        for r in v6_diag.itertuples()
    }
    common_keys = [
        (r.cml_id, r.sublink_id)
        for r in v6_diag.sort_values("series_index").itertuples()
        if (r.cml_id, r.sublink_id) in old_map
    ]
    if not common_keys:
        raise RuntimeError("No common cml_id/sublink_id between original and V6")

    meta = v6_diag.set_index(["cml_id", "sublink_id"])
    geom = raw_geometry(args.raw_cml_nc)

    hum = xr.open_dataset(args.humidity_nc, cache=False)
    temp = xr.open_dataset(args.temperature_nc, cache=False)
    old = xr.open_dataset(args.original_nc, cache=False, decode_times=True)
    v6 = xr.open_dataset(args.v6_nc, cache=False, decode_times=True)

    try:
        if "H" not in hum:
            raise RuntimeError(f"`H` missing from humidity file; vars={list(hum.data_vars)}")
        if "T" not in temp:
            raise RuntimeError(f"`T` missing from temperature file; vars={list(temp.data_vars)}")

        stations = common_met_stations(hum, temp)
        if stations.empty:
            raise RuntimeError("No station IDs common to humidity and temperature files")

        physical_cmls = list(dict.fromkeys(k[0] for k in common_keys))
        missing_geom = [x for x in physical_cmls if x not in geom]
        if missing_geom:
            raise RuntimeError(f"Raw geometry missing CML IDs, e.g. {missing_geom[:5]}")

        assignments = assign_stations(
            physical_cmls, geom, stations, args.max_station_distance_km
        )
        assignments.to_csv(out_dir / "station_assignments.csv", index=False)
        assign_map = assignments.set_index("cml_id")

        print(f"Common directional sublinks : {len(common_keys):,}")
        print(f"Physical CMLs               : {len(physical_cmls):,}")
        print(
            f"CMLs with common H/T <=5 km : "
            f"{int(assignments['within_validation_radius'].sum()):,}"
        )

        # Retrieval grids must be identical at the 5-min level for a strict paired comparison.
        old_t5 = decode_retrieval_time(old)
        v6_t5 = decode_retrieval_time(v6)
        if not old_t5.equals(v6_t5):
            raise RuntimeError("Original and V6 5-min time coordinates differ")
        t10, _dummy = strict_5_to_10(
            old_t5, np.zeros(len(old_t5)), np.ones(len(old_t5), dtype=np.uint8)
        )

        # CHMI exact-time arrays are cached by station because both directions of a
        # physical CML share the same validation station.
        truth_cache = {}

        def truth_for_station(sid):
            sid = str(sid)
            if sid in truth_cache:
                return truth_cache[sid]

            h = hum["H"].sel(station=sid)
            tt = temp["T"].sel(station=sid)

            hs = pd.Series(
                np.asarray(h.values, dtype=float),
                index=pd.DatetimeIndex(pd.to_datetime(h["time"].values))
            )
            ts = pd.Series(
                np.asarray(tt.values, dtype=float),
                index=pd.DatetimeIndex(pd.to_datetime(tt["time"].values))
            )

            # EXACT matching only. No interpolation / nearest-time filling.
            hv = hs.reindex(t10).to_numpy(float)
            tv = ts.reindex(t10).to_numpy(float)
            rho = rh_t_to_wvd(hv, tv)
            truth_cache[sid] = rho
            return rho

        # Verify dimensions/semantics.
        for ds, label in ((old, "original"), (v6, "V6")):
            if "wvd_internal_anchor" not in ds:
                raise RuntimeError(f"{label} missing wvd_internal_anchor")
            if "retrieval_status" not in ds:
                raise RuntimeError(f"{label} missing retrieval_status")
            if ds["wvd_internal_anchor"].dims != ("method", "series", "time"):
                raise RuntimeError(
                    f"{label} unexpected WVD dims: {ds['wvd_internal_anchor'].dims}"
                )
            if ds["retrieval_status"].dims != ("method", "series", "time"):
                raise RuntimeError(
                    f"{label} unexpected status dims: {ds['retrieval_status'].dims}"
                )

        pooled = {
            m: {"ref": [], "old": [], "v6": []}
            for m in METHODS
        }
        per_rows = []
        coverage_rows = []

        total = len(common_keys)
        for kk, key in enumerate(common_keys, start=1):
            cid, sid = key
            ar = assign_map.loc[cid]
            if not bool(ar["within_validation_radius"]):
                continue

            station_id = str(ar["station_id"])
            truth = truth_for_station(station_id)

            mr = meta.loc[(cid, sid)]
            fghz = float(mr["frequency_ghz"])
            L = float(mr["length_km"])
            lbin = length_bin(L)
            fband = frequency_band(fghz)

            oi = old_map[key]
            vi = v6_map[key]

            for mi, method in enumerate(METHODS):
                old5 = np.asarray(
                    old["wvd_internal_anchor"].isel(method=mi, series=oi).values,
                    dtype=float
                )
                olds = np.asarray(
                    old["retrieval_status"].isel(method=mi, series=oi).values
                )
                v65 = np.asarray(
                    v6["wvd_internal_anchor"].isel(method=mi, series=vi).values,
                    dtype=float
                )
                v6s = np.asarray(
                    v6["retrieval_status"].isel(method=mi, series=vi).values
                )

                old_t10, old10 = strict_5_to_10(old_t5, old5, olds)
                v6_t10, v610 = strict_5_to_10(v6_t5, v65, v6s)
                if not old_t10.equals(t10) or not v6_t10.equals(t10):
                    raise RuntimeError("Internal 10-min time mismatch")

                truth_ok = np.isfinite(truth)
                old_ok = np.isfinite(old10) & truth_ok
                v6_ok = np.isfinite(v610) & truth_ok
                common = old_ok & v6_ok

                coverage_rows.append({
                    "method": method,
                    "cml_id": cid,
                    "sublink_id": sid,
                    "frequency_ghz": fghz,
                    "length_km": L,
                    "station_id": station_id,
                    "station_distance_to_path_km": float(ar["station_distance_to_path_km"]),
                    "n_truth": int(truth_ok.sum()),
                    "n_original": int(old_ok.sum()),
                    "n_v6": int(v6_ok.sum()),
                    "n_common": int(common.sum()),
                    "n_original_only": int((old_ok & ~v6_ok).sum()),
                    "n_v6_only": int((v6_ok & ~old_ok).sum()),
                })

                ref_c = truth[common]
                old_c = old10[common]
                v6_c = v610[common]

                if ref_c.size:
                    pooled[method]["ref"].append(ref_c.astype(np.float32))
                    pooled[method]["old"].append(old_c.astype(np.float32))
                    pooled[method]["v6"].append(v6_c.astype(np.float32))

                mo = metric_block(ref_c, old_c)
                mv = metric_block(ref_c, v6_c)
                eligible = int(ref_c.size >= args.min_pairs_per_sublink)

                per_rows.append({
                    "method": method,
                    "cml_id": cid,
                    "sublink_id": sid,
                    "frequency_ghz": fghz,
                    "frequency_band": fband,
                    "length_km": L,
                    "length_bin": lbin,
                    "station_id": station_id,
                    "station_distance_to_path_km": float(ar["station_distance_to_path_km"]),
                    "n_common": int(ref_c.size),
                    "eligible_metrics": eligible,
                    "min_pairs_rule": int(args.min_pairs_per_sublink),
                    "original_pcc": mo["PCC"],
                    "original_rmse": mo["RMSE"],
                    "original_bias": mo["Bias"],
                    "original_mae": mo["MAE"],
                    "v6_pcc": mv["PCC"],
                    "v6_rmse": mv["RMSE"],
                    "v6_bias": mv["Bias"],
                    "v6_mae": mv["MAE"],
                    "delta_pcc_v6_minus_original": mv["PCC"] - mo["PCC"],
                    "delta_rmse_v6_minus_original": mv["RMSE"] - mo["RMSE"],
                    "delta_abs_bias_v6_minus_original": (
                        abs(mv["Bias"]) - abs(mo["Bias"])
                        if np.isfinite(mv["Bias"]) and np.isfinite(mo["Bias"])
                        else np.nan
                    ),
                })

            if kk % 100 == 0 or kk == total:
                print(f"[{kk:4d}/{total:4d}] validated directional sublinks")

        # Concatenate pooled pairs.
        for method in METHODS:
            for k in ("ref", "old", "v6"):
                pooled[method][k] = (
                    np.concatenate(pooled[method][k])
                    if pooled[method][k]
                    else np.empty(0, dtype=np.float32)
                )

        per = pd.DataFrame(per_rows)
        coverage = pd.DataFrame(coverage_rows)
        per.to_csv(out_dir / "per_sublink_metrics_common.csv", index=False)
        coverage.to_csv(out_dir / "coverage_summary.csv", index=False)

        # Pooled metrics on the SAME common pairs.
        pooled_rows = []
        for method in METHODS:
            ref = pooled[method]["ref"]
            oldx = pooled[method]["old"]
            v6x = pooled[method]["v6"]
            for model, est in (("Original", oldx), ("V6", v6x)):
                mm = metric_block(ref, est)
                pooled_rows.append({
                    "method": method,
                    "model": model,
                    **mm,
                })
        pooled_df = pd.DataFrame(pooled_rows)
        pooled_df.to_csv(out_dir / "pooled_metrics_common.csv", index=False)

        # Stratified median / IQR per-sublink metrics.
        strat_rows = []
        eligible = per[per["eligible_metrics"] == 1].copy()

        def add_strat(kind, col, categories):
            for method in METHODS:
                dm = eligible[eligible["method"] == method]
                for cat in categories:
                    g = dm[dm[col] == cat]
                    for model, prefix in (("Original", "original"), ("V6", "v6")):
                        row = {
                            "stratification": kind,
                            "group": cat,
                            "method": method,
                            "model": model,
                            "n_sublinks": int(len(g)),
                            "n_pairs_sum": int(g["n_common"].sum()),
                        }
                        for metric in ("pcc", "rmse", "bias"):
                            x = g[f"{prefix}_{metric}"].to_numpy(float)
                            x = x[np.isfinite(x)]
                            row[f"{metric}_median"] = float(np.median(x)) if x.size else np.nan
                            row[f"{metric}_q25"] = float(np.quantile(x, .25)) if x.size else np.nan
                            row[f"{metric}_q75"] = float(np.quantile(x, .75)) if x.size else np.nan
                        strat_rows.append(row)

        add_strat("length_bin_km", "length_bin", LENGTH_LABELS)
        add_strat("frequency_band", "frequency_band", list(FREQ_BANDS))
        strat = pd.DataFrame(strat_rows)
        strat.to_csv(out_dir / "stratified_summary.csv", index=False)

        # Length-threshold sensitivity, using paired per-link metrics.
        sens_rows = []
        for method in METHODS:
            dm = eligible[eligible["method"] == method]
            for threshold in (2, 3, 4, 5, 6, 7, 8):
                g = dm[dm["length_km"] >= threshold]
                for model, prefix in (("Original", "original"), ("V6", "v6")):
                    pcc = g[f"{prefix}_pcc"].to_numpy(float)
                    rmse = g[f"{prefix}_rmse"].to_numpy(float)
                    bias = g[f"{prefix}_bias"].to_numpy(float)
                    sens_rows.append({
                        "method": method,
                        "model": model,
                        "min_length_km": threshold,
                        "n_sublinks": int(len(g)),
                        "median_PCC": float(np.nanmedian(pcc)) if len(g) else np.nan,
                        "median_RMSE": float(np.nanmedian(rmse)) if len(g) else np.nan,
                        "median_abs_Bias": float(np.nanmedian(np.abs(bias))) if len(g) else np.nan,
                    })
        sens = pd.DataFrame(sens_rows)
        sens.to_csv(out_dir / "min_length_sensitivity.csv", index=False)

        # Figures.
        for method in METHODS:
            plot_pooled_hexbin(
                method, pooled, pooled_rows, pooled_dir,
                args.scatter_max_gm3, args.dpi
            )
            plot_overall_metric_boxes(method, per, per_dir, args.dpi)
            plot_stratified(
                method, per, "length_bin", LENGTH_LABELS,
                len_dir / f"{method}_metrics_by_length",
                "Path length [km]", args.dpi
            )
            plot_stratified(
                method, per, "frequency_band", list(FREQ_BANDS),
                freq_dir / f"{method}_metrics_by_frequency",
                "Frequency band", args.dpi
            )
            plot_length_sensitivity(method, sens, sens_dir, args.dpi)

        print("\n" + "=" * 100)
        print("POOLED PAIRED METRICS")
        print("=" * 100)
        print(pooled_df.to_string(index=False))

        print("\nPer-sublink eligible counts:")
        print(
            eligible.groupby("method").size().rename("n_sublinks").to_string()
        )

        print("\nOutput:", out_dir.resolve())
        print("\nIMPORTANT: length/frequency stratification is diagnostic only.")
        print("Do not select individual CMLs using CHMI PCC/RMSE/bias.")

    finally:
        old.close()
        v6.close()
        hum.close()
        temp.close()


if __name__ == "__main__":
    main()
