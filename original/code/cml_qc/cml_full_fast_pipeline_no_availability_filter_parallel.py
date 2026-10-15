#!/usr/bin/env python3
"""
Full CML rainfall pipeline: diagnostics, Graf wet/dry classification,
Blettner filtering, rainfall retrieval experiments, gridding, publication plots,
GIF animation, and hexbin comparison.

Designed for the Czech CML/radar/gauge dataset.

Main outputs
------------
output_dir/
  diagnostics/
    link_diagnostics.csv
    sublink_diagnostics.csv
    stage_summary.csv
    CML_filtered_has_data_and_radar_reference.nc
    CML_filtered_blettner_full.nc
  wet_dry/
    graf_wet_detection_5min.nc              (optional with --save-wet-flags)
    graf_vs_radar_contingency.csv
    graf_vs_radar_metrics.csv
    graf_vs_radar_contingency.png
  rainfall/
    cml_rainrate_all_experiments.nc         (optional all combinations)
    cml_rainrate_selected_pastorek_itu03.nc
    cml_rainfall_estimation_log.csv
  gridded/
    rainfall_maps_selected_event.nc
  plots/
    snapshot_cml_vs_native_radolan.png/.pdf
    animation_cml_vs_native_radolan.gif
    event_frame_metrics.csv
    hexbin_cml_vs_radar.png/.pdf
    hexbin_metrics.csv

Notes
-----
* Radar path-average file is used for wet/dry reference and rainfall retrieval.
* Native RADOLAN/RADKLIM daily files are used only for the realistic field plot.
* Graf wet/dry uses exactly the rolling-standard-deviation logic requested:
      wet_threshold = 1.12 * rolling_std_60min.quantile(0.8)
      wet = rolling_std_60min_min45 > wet_threshold
* Rainfall retrieval uses radar wet/dry, not Graf wet/dry, unless --use-graf-wet-for-estimation.
* Blettner/Graf filter is used in two ways:
      1) diagnostics and link filtering;
      2) optional cleaned TRSL variant for the filtered experiment.
* In this optimized version, availability is kept as a diagnostic by default,
  but it is NOT used to reject links unless --use-availability-filter is set.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import os
# Keep BLAS/numexpr/OpenMP libraries from using many threads inside each worker.
# The pipeline parallelizes by CML link, so one numeric thread per process gives
# better CPU utilization and avoids oversubscription.
for _env_name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_env_name, "1")
import re
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LogNorm, ListedColormap
from matplotlib.ticker import AutoMinorLocator

try:
    from scipy.spatial import cKDTree
except Exception as exc:  # pragma: no cover
    raise ImportError("scipy is required. Install with: conda install -c conda-forge scipy") from exc

try:
    import imageio.v2 as imageio
except Exception:  # pragma: no cover
    imageio = None

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    def tqdm(x, **kwargs):
        return x

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.minor.size": 2.5,
    "ytick.minor.size": 2.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

RAIN_BOUNDS = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 50, 80]
ACCUM_BOUNDS = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 50, 80]
RADAR_VAR_CANDIDATES = ["YW", "RW", "RR", "rainfall_rate", "rain_rate", "R", "rainrate", "precip", "__xarray_dataarray_variable__"]
HEXBIN_COLORS = [
    "#000080", "#0015d9", "#0057ff", "#00a5ff", "#00d5d5",
    "#36ef9a", "#8cff66", "#d9ff32", "#ffff00", "#ffd000",
    "#ff9a00", "#ff5a00", "#ff1a00", "#d00000", "#7f0000",
]
HEXBIN_CMAP = ListedColormap(HEXBIN_COLORS, name="reference_blue_yellow_red")

# ---------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------
def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def find_dim(ds: xr.Dataset | xr.DataArray, candidates: Sequence[str], required: bool = True) -> Optional[str]:
    dims = list(ds.dims)
    for cand in candidates:
        if cand in dims:
            return cand
    for d in dims:
        dl = d.lower()
        for cand in candidates:
            if cand.lower() in dl:
                return d
    if required:
        raise KeyError(f"Cannot find dimension among {candidates}. Available dims: {dims}")
    return None


def pick_var(ds: xr.Dataset, candidates: Sequence[str], required: bool = True) -> Optional[str]:
    for c in candidates:
        if c in ds.data_vars or c in ds.variables:
            return c
    if required:
        raise KeyError(f"None of {list(candidates)} found. Available: {list(ds.variables)}")
    return None


def as_datetime_index(coord) -> pd.DatetimeIndex:
    vals = coord.values if hasattr(coord, "values") else coord
    return pd.DatetimeIndex(vals)


def time_slice_indices(times: pd.DatetimeIndex, tstart: pd.Timestamp, tend: pd.Timestamp) -> Tuple[int, int]:
    vals = times.values
    i0 = int(np.searchsorted(vals, np.datetime64(tstart), "left"))
    i1 = int(np.searchsorted(vals, np.datetime64(tend), "right"))
    return max(0, i0), min(len(times), i1)


def pol_to_HV(x) -> str:
    if isinstance(x, (bytes, np.bytes_)):
        x = x.decode("utf-8", errors="ignore")
    s = str(x).strip().lower()
    if "v" in s and "h" not in s:
        return "V"
    return "H"


def scalar_float(x) -> float:
    try:
        arr = np.asarray(x, dtype="float64")
        if arr.size == 0:
            return float("nan")
        return float(arr.reshape(-1)[0])
    except Exception:
        return float("nan")


def safe_values(da: xr.DataArray) -> np.ndarray:
    return np.asarray(da.values)


def to_hz(freq_raw: float) -> float:
    f = float(freq_raw)
    if not np.isfinite(f):
        return np.nan
    # Most CML metadata store MHz. Some derived files store GHz.
    if f < 200:      # GHz
        return f * 1e9
    if f < 200000:   # MHz
        return f * 1e6
    return f         # already Hz


def to_ghz(freq_raw: float) -> float:
    return to_hz(freq_raw) / 1e9


def finite_count(a: np.ndarray) -> int:
    return int(np.isfinite(a).sum())


def radar_values_to_mm_h(raw: np.ndarray, scale: float) -> np.ndarray:
    return np.asarray(raw, dtype="float64") * float(scale)


def encode_times_seconds(times: pd.DatetimeIndex) -> np.ndarray:
    return (pd.DatetimeIndex(times).view("int64") // 10**9).astype("int64")


def write_dataset(ds: xr.Dataset, path: Path, complevel: int = 4):
    path.parent.mkdir(parents=True, exist_ok=True)
    encoding = {}
    for name, da in ds.data_vars.items():
        if np.issubdtype(da.dtype, np.floating):
            encoding[name] = {"zlib": True, "complevel": int(complevel), "_FillValue": np.float32(np.nan)}
        elif np.issubdtype(da.dtype, np.integer):
            encoding[name] = {"zlib": True, "complevel": int(complevel)}
    try:
        ds.to_netcdf(path, encoding=encoding)
    except Exception:
        # scipy backend fallback does not support zlib/complevel
        ds.to_netcdf(path)

# ---------------------------------------------------------------------
# Blettner/Graf processing functions
# ---------------------------------------------------------------------
def detect_plateaus_np(rsl: np.ndarray, time_span: int = 3, max_thld: float = -85, std_thld: float = 0.5, pad: int = 5) -> np.ndarray:
    s = pd.Series(np.asarray(rsl, dtype="float64"))
    low = s.rolling(time_span, center=True, min_periods=1).max().to_numpy() < max_thld
    steady = s.rolling(time_span, center=True, min_periods=1).std().to_numpy() < std_thld
    sane = ~(low & steady)
    # pad rejected samples using rolling min of sane mask
    sane2 = pd.Series(sane.astype(float)).rolling(pad, center=True, min_periods=1).min().to_numpy() > 0
    return sane2.astype(bool)


def get_blackout_start_end(rsl: np.ndarray, rsl_threshold: float = -65) -> Tuple[np.ndarray, np.ndarray]:
    rsl = np.asarray(rsl, dtype="float64")
    below = rsl < rsl_threshold
    isn = np.isnan(rsl)
    gap_start = np.roll(isn, -1) & below
    gap_end = np.roll(isn, 1) & below
    if gap_start.size:
        gap_start[-1] = False
        gap_end[0] = False
    return gap_start, gap_end


def blackout_mask_from_markers(rsl: np.ndarray, gap_start: np.ndarray, gap_end: np.ndarray, max_gap_length: int = 60) -> np.ndarray:
    rsl = np.asarray(rsl, dtype="float64")
    mask = np.zeros(rsl.shape, dtype=bool)
    in_gap = False
    cnt = 0
    for i in range(len(rsl)):
        if bool(gap_start[i]):
            in_gap = True
            cnt = 0
        if (not bool(gap_start[i])) and np.isfinite(rsl[i]):
            in_gap = False
        if bool(gap_end[i]):
            in_gap = False
        if in_gap and np.isnan(rsl[i]) and cnt < max_gap_length:
            mask[i] = True
            cnt += 1
            if cnt > max_gap_length:
                mask[max(0, i - max_gap_length):i + 1] = False
                in_gap = False
    return mask


def blackout_filling_np(tsl: np.ndarray, rsl: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    tsl_o = np.asarray(tsl, dtype="float64").copy()
    rsl_o = np.asarray(rsl, dtype="float64").copy()
    gs, ge = get_blackout_start_end(rsl_o, -65)
    mf = blackout_mask_from_markers(rsl_o, gs, ge)
    mr = blackout_mask_from_markers(rsl_o[::-1], ge[::-1], gs[::-1])[::-1]
    mask = mf | mr
    tsl_max = np.nanmax(tsl_o) if np.isfinite(tsl_o).any() else np.nan
    rsl_min = np.nanmin(rsl_o) if np.isfinite(rsl_o).any() else np.nan
    tsl_o[mask] = tsl_max
    rsl_o[mask] = rsl_min
    return tsl_o, rsl_o, mask


def interpolate_short_gaps(values: np.ndarray, times: pd.DatetimeIndex, max_gap: str = "5min") -> np.ndarray:
    s = pd.Series(np.asarray(values, dtype="float64"), index=times)
    # time interpolation with limit equal to number of 1-minute gaps; robust for 1-min CML
    step_min = 1.0
    if len(times) > 1:
        diffs = np.diff(times.view("int64")) / 1e9 / 60.0
        diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
        if diffs.size:
            step_min = float(np.nanmedian(diffs))
    limit = max(1, int(round(pd.Timedelta(max_gap).total_seconds() / 60.0 / step_min)))
    return s.interpolate(method="time", limit=limit, limit_area="inside").to_numpy(dtype="float64")


def temporal_sanity_check_np(trsl: np.ndarray, window: int, thld: float = 2.0, perc: float = 0.1) -> bool:
    s = pd.Series(np.asarray(trsl, dtype="float64"))
    rs = s.rolling(window, center=True, min_periods=max(2, min(window, 10))).std().to_numpy()
    high = np.isfinite(rs) & (rs > thld)
    low = np.isfinite(rs) & (rs <= thld)
    denom = int(high.sum() + low.sum())
    if denom == 0:
        return False
    return (float(high.sum()) / float(denom)) <= perc


def availability_check_np(trsl: np.ndarray, perc_thld: float = 0.5) -> bool:
    trsl = np.asarray(trsl, dtype="float64")
    if trsl.size == 0:
        return False
    return (float(np.isnan(trsl).sum()) / float(trsl.size)) < perc_thld


def get_itu_ab(f_ghz: float, pol: str, model: str) -> Tuple[float, float]:
    try:
        import pycomlink as pycml
        try:
            a, b = pycml.processing.k_R_relation.a_b(f_GHz=f_ghz, pol=pol, approx_type=model)
        except TypeError:
            a, b = pycml.processing.k_R_relation.a_b(f_GHz=f_ghz, pol=pol, a_b_approximation=model)
        return scalar_float(a), scalar_float(b)
    except Exception:
        # Fallback approximation for demo/syntax testing. Real runs should use pycomlink.
        # Values are intentionally approximate, not a replacement for ITU tables.
        f = float(f_ghz)
        if f < 10:
            a, b = 0.01, 1.25
        elif f < 40:
            a, b = 0.08, 1.08
        else:
            a, b = 0.20, 0.95
        if str(model).upper().endswith("2005"):
            a *= 1.02
        return float(a), float(b)


def detection_limit_from_kR(length_km: float, freq_hz: float, pol: str, a_quantization_db: float = 0.33) -> float:
    a, _ = get_itu_ab(freq_hz / 1e9, pol, "ITU_2003")
    sens = a * float(length_km)
    if not np.isfinite(sens) or sens <= 0:
        return np.nan
    return float(a_quantization_db) / sens


def apply_blettner_signal_filter(tsl: np.ndarray, rsl: np.ndarray, times: pd.DatetimeIndex, proc_line: str) -> Tuple[np.ndarray, Dict[str, float]]:
    tsl2 = np.asarray(tsl, dtype="float64").copy()
    rsl2 = np.asarray(rsl, dtype="float64").copy()
    info = {"plateau_removed_samples": 0, "blackout_filled_samples": 0}
    if proc_line == "full":
        sane = detect_plateaus_np(rsl2)
        info["plateau_removed_samples"] = int((~sane).sum())
        tsl2[~sane] = np.nan
        rsl2[~sane] = np.nan
        tsl2, rsl2, mask_b = blackout_filling_np(tsl2, rsl2)
        info["blackout_filled_samples"] = int(mask_b.sum())
    trsl = tsl2 - rsl2
    if proc_line in ("graf_2020", "full"):
        trsl = interpolate_short_gaps(trsl, times, "5min")
    return trsl, info


def blettner_screen_sublink(tsl: np.ndarray, rsl: np.ndarray, times: pd.DatetimeIndex, length_km: float, freq_hz: float, pol: str, proc_line: str, min_valid_samples: int, use_availability_filter: bool = False, availability_perc_thld: float = 0.5) -> Dict[str, float]:
    info: Dict[str, float] = {}
    tsl = np.asarray(tsl, dtype="float64")
    rsl = np.asarray(rsl, dtype="float64")
    info["valid_tsl"] = finite_count(tsl)
    info["valid_rsl"] = finite_count(rsl)
    info["has_cml_data_sublink"] = int(info["valid_tsl"] >= min_valid_samples and info["valid_rsl"] >= min_valid_samples)
    info["plateau_removed_samples"] = 0
    info["blackout_filled_samples"] = 0
    info["detection_limit"] = np.nan
    info["pass_detection_limit"] = 0
    info["pass_long_std"] = 0
    info["pass_short_std"] = 0
    info["pass_availability"] = 0
    info["availability_nan_fraction"] = np.nan
    info["availability_used_as_filter"] = int(bool(use_availability_filter))
    info["pass_blettner"] = 0
    if not info["has_cml_data_sublink"]:
        return info
    trsl, clean_info = apply_blettner_signal_filter(tsl, rsl, times, proc_line)
    info.update(clean_info)
    det = detection_limit_from_kR(length_km, freq_hz, pol)
    info["detection_limit"] = float(det) if np.isfinite(det) else np.nan
    info["pass_detection_limit"] = int(np.isfinite(det) and det < 2.0)
    info["pass_long_std"] = int(temporal_sanity_check_np(trsl, window=60 * 24 * 7, thld=2.0, perc=0.1))
    info["pass_short_std"] = int(temporal_sanity_check_np(trsl, window=300, thld=2.0, perc=0.1))
    info["availability_nan_fraction"] = float(np.isnan(trsl).sum() / trsl.size) if trsl.size else np.nan
    info["pass_availability"] = int(availability_check_np(trsl, perc_thld=float(availability_perc_thld)))
    info["availability_used_as_filter"] = int(bool(use_availability_filter))
    if proc_line == "no_filter":
        info["pass_blettner"] = 1
    elif proc_line == "graf_2020":
        info["pass_blettner"] = int(info["pass_long_std"] and info["pass_short_std"])
    elif proc_line == "full":
        base_pass = int(info["pass_detection_limit"] and info["pass_long_std"] and info["pass_short_std"])
        if bool(use_availability_filter):
            info["pass_blettner"] = int(base_pass and info["pass_availability"])
        else:
            # Availability is still reported, but does NOT reject the link.
            info["pass_blettner"] = int(base_pass)
    else:
        raise ValueError("proc_line must be no_filter, graf_2020, or full")
    return info

# ---------------------------------------------------------------------
# Graf wet/dry detection
# ---------------------------------------------------------------------
def graf_wet_from_trsl(trsl_raw: np.ndarray) -> Tuple[np.ndarray, float]:
    trsl_raw = np.asarray(trsl_raw, dtype="float64")
    s = pd.Series(trsl_raw)
    std60_full = s.rolling(60, center=True, min_periods=60).std().to_numpy()
    q80 = np.nanquantile(std60_full, 0.8) if np.isfinite(std60_full).any() else np.nan
    wet_threshold = 1.12 * q80 if np.isfinite(q80) else np.nan
    std60 = s.rolling(60, center=True, min_periods=45).std().to_numpy()
    wet = std60 > wet_threshold if np.isfinite(wet_threshold) else np.zeros(trsl_raw.shape, dtype=bool)
    wet[~np.isfinite(std60)] = False
    return wet.astype(bool), float(wet_threshold) if np.isfinite(wet_threshold) else np.nan


def forward_fill_radar_wet_to_cml_time(radar_time: pd.DatetimeIndex, cml_time: pd.DatetimeIndex, wet_5min: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(radar_time.values, cml_time.values, "right") - 1
    pre = idx < 0
    idx = np.clip(idx, 0, max(0, len(wet_5min) - 1))
    out = np.asarray(wet_5min, dtype=bool)[idx]
    out[pre] = False
    return out


def confusion_counts(pred_wet: np.ndarray, ref_wet: np.ndarray, valid: Optional[np.ndarray] = None) -> Dict[str, int]:
    p = np.asarray(pred_wet, dtype=bool)
    r = np.asarray(ref_wet, dtype=bool)
    if valid is None:
        v = np.ones(p.shape, dtype=bool)
    else:
        v = np.asarray(valid, dtype=bool)
    return {
        "hit": int((v & p & r).sum()),
        "false_alarm": int((v & p & ~r).sum()),
        "miss": int((v & ~p & r).sum()),
        "correct_negative": int((v & ~p & ~r).sum()),
    }


def detection_metrics_from_counts(c: Dict[str, int]) -> Dict[str, float]:
    H = c.get("hit", 0)
    F = c.get("false_alarm", 0)
    M = c.get("miss", 0)
    C = c.get("correct_negative", 0)
    total = H + F + M + C
    pod = H / (H + M) if (H + M) else np.nan
    far = F / (H + F) if (H + F) else np.nan
    acc = (H + C) / total if total else np.nan
    denom = math.sqrt((H + F) * (H + M) * (C + F) * (C + M))
    mcc = ((H * C) - (F * M)) / denom if denom else np.nan
    return {"N": total, "POD": pod, "FAR": far, "ACC": acc, "MCC": mcc}


def plot_contingency(counts: Dict[str, int], out_png: Path):
    H, F, M, C = counts["hit"], counts["false_alarm"], counts["miss"], counts["correct_negative"]
    radar_wet = H + M
    radar_dry = F + C
    vals = np.array([[H, F], [M, C]], dtype=float)
    labels = [["Hit", "False alarm"], ["Miss", "Correct negative"]]
    subs = [
        [f"{100*H/radar_wet:.1f}% of radar wet" if radar_wet else "", f"{100*F/radar_dry:.1f}% of radar dry" if radar_dry else ""],
        [f"{100*M/radar_wet:.1f}% of radar wet" if radar_wet else "", f"{100*C/radar_dry:.1f}% of radar dry" if radar_dry else ""],
    ]
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    im = ax.imshow(vals, cmap="Blues")
    vmax = np.nanmax(vals)
    for i in range(2):
        for j in range(2):
            color = "white" if vals[i, j] > 0.55 * vmax else "#1a1a1a"
            ax.text(j, i - 0.25, labels[i][j], ha="center", va="center", fontsize=20, fontweight="bold", color=color)
            ax.text(j, i + 0.02, f"{int(vals[i,j]):,}", ha="center", va="center", fontsize=22, color=color)
            ax.text(j, i + 0.28, subs[i][j], ha="center", va="center", fontsize=14, color=color)
    ax.set_xticks([0, 1], ["Wet", "Dry"], fontsize=20)
    ax.set_yticks([0, 1], ["Wet", "Dry"], fontsize=20, rotation=90, va="center")
    ax.set_xlabel("Radar reference", fontsize=24, labelpad=16)
    ax.set_ylabel("CML rolling-std detection", fontsize=24, labelpad=16)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 2, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=4)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------------
# Rainfall estimation
# ---------------------------------------------------------------------
def baseline_last_dry(trsl: np.ndarray, wet: np.ndarray) -> np.ndarray:
    trsl = np.asarray(trsl, dtype="float64")
    wet = np.asarray(wet, dtype=bool)
    n = trsl.size
    dry_valid = (~wet) & np.isfinite(trsl)
    dry_idx = np.where(dry_valid, np.arange(n, dtype=np.int64), -1)
    last_idx = np.maximum.accumulate(dry_idx)
    out = np.full(n, np.nan, dtype="float64")
    ok = last_idx >= 0
    out[ok] = trsl[last_idx[ok]]
    return out


def baseline_constant_fallback(trsl: np.ndarray, wet: np.ndarray, n_average_last_dry: int = 1) -> np.ndarray:
    # Closest previous dry average over the last n dry samples. This fallback is used only when pycomlink is unavailable.
    trsl = np.asarray(trsl, dtype="float64")
    wet = np.asarray(wet, dtype=bool)
    out = np.full(trsl.shape, np.nan, dtype="float64")
    buf: List[float] = []
    n = max(1, int(n_average_last_dry))
    last = np.nan
    for i, val in enumerate(trsl):
        if (not wet[i]) and np.isfinite(val):
            buf.append(float(val))
            if len(buf) > n:
                buf = buf[-n:]
            last = float(np.nanmean(buf))
        out[i] = last
    return out


def compute_baseline(method: str, trsl: np.ndarray, wet: np.ndarray, times: pd.DatetimeIndex, n_average_last_dry: int) -> np.ndarray:
    if method == "last_dry":
        return baseline_last_dry(trsl, wet)
    if method == "pycml_constant":
        try:
            import pycomlink as pycml
            xr_trsl = xr.DataArray(trsl.astype("float64"), coords={"time": times.values}, dims=["time"])
            xr_wet = xr.DataArray(wet.astype(bool), coords={"time": times.values}, dims=["time"])
            b = pycml.processing.baseline.baseline_constant(xr_trsl, xr_wet, n_average_last_dry=n_average_last_dry)
            return np.asarray(b.values, dtype="float64")
        except Exception:
            return baseline_constant_fallback(trsl, wet, n_average_last_dry)
    raise ValueError(f"Unknown baseline method {method}")


def compute_waa(method: str, A_obs: np.ndarray, trsl: np.ndarray, baseline: np.ndarray, wet: np.ndarray, times: pd.DatetimeIndex, f_ghz: float, pol: str, length_km: float, args) -> np.ndarray:
    A_obs = np.asarray(A_obs, dtype="float64")
    if method == "none":
        return np.zeros_like(A_obs)
    try:
        import pycomlink as pycml
        xr_A = xr.DataArray(A_obs, coords={"time": times.values}, dims=["time"])
        xr_trsl = xr.DataArray(trsl.astype("float64"), coords={"time": times.values}, dims=["time"])
        xr_base = xr.DataArray(baseline.astype("float64"), coords={"time": times.values}, dims=["time"])
        xr_wet = xr.DataArray(wet.astype(bool), coords={"time": times.values}, dims=["time"])
        if method == "leijnse_2008":
            waa = pycml.processing.wet_antenna.waa_leijnse_2008_from_A_obs(
                A_obs=xr_A, f_Hz=f_ghz * 1e9, pol=pol, L_km=length_km,
                gamma=args.waa_gamma, l_antenna=args.waa_l_antenna, delta=args.waa_delta,
            )
        elif method == "pastorek_2021":
            waa = pycml.processing.wet_antenna.waa_pastorek_2021_from_A_obs(
                A_obs=xr_A, f_Hz=f_ghz * 1e9, pol=pol, L_km=length_km,
                A_max=args.pastorek_a_max, zeta=args.pastorek_zeta, d=args.pastorek_d,
            )
        elif method == "schleiss_2013":
            waa = pycml.processing.wet_antenna.waa_schleiss_2013(
                rsl=xr_trsl, baseline=xr_base, wet=xr_wet,
                waa_max=args.schleiss_waa_max, delta_t=args.schleiss_delta_t, tau=args.schleiss_tau,
            )
        else:
            raise ValueError(method)
        out = np.asarray(waa.values, dtype="float64")
    except Exception:
        # Simple fallback for demo. Real run uses pycomlink.
        if method == "pastorek_2021":
            out = np.minimum(args.pastorek_a_max, args.pastorek_zeta * np.power(np.maximum(A_obs, 0), 1.0) + args.pastorek_d)
        elif method == "leijnse_2008":
            out = np.minimum(2.5, 0.15 * np.maximum(A_obs, 0))
        elif method == "schleiss_2013":
            out = np.minimum(args.schleiss_waa_max, 0.10 * np.maximum(A_obs, 0))
        else:
            out = np.zeros_like(A_obs)
    out = np.where(np.isfinite(out), out, 0.0)
    out = np.clip(out, 0.0, np.where(np.isfinite(A_obs), A_obs, 0.0))
    return out


def calc_R_from_A(A_rain: np.ndarray, kr_model: str, f_ghz: float, pol: str, length_km: float) -> np.ndarray:
    A_rain = np.asarray(A_rain, dtype="float64")
    try:
        import pycomlink as pycml
        xr_A = xr.DataArray(A_rain, dims=["time"])
        try:
            R = pycml.processing.k_R_relation.calc_R_from_A(A=xr_A, L_km=length_km, f_GHz=f_ghz, pol=pol, a_b_approximation=kr_model)
        except TypeError:
            R = pycml.processing.k_R_relation.calc_R_from_A(A=xr_A, L_km=length_km, f_GHz=f_ghz, pol=pol, approx_type=kr_model)
        return np.asarray(R.values, dtype="float32")
    except Exception:
        a, b = get_itu_ab(f_ghz, pol, kr_model)
        gamma = np.maximum(A_rain, 0) / max(float(length_km), 1e-6)
        R = np.power(np.maximum(gamma, 0) / max(a, 1e-12), 1.0 / max(b, 1e-12))
        R[~np.isfinite(R)] = np.nan
        return R.astype("float32")


def resample_1min_rate_to_5min(values_1min: np.ndarray, time_1min: pd.DatetimeIndex, target_5min: pd.DatetimeIndex) -> np.ndarray:
    # Rate -> 1-minute depth -> 5-minute depth -> rate
    s = pd.Series(np.asarray(values_1min, dtype="float64"), index=time_1min)
    out = (s / 60.0).resample("5min", label="left", closed="left").sum(min_count=1) * 12.0
    return out.reindex(target_5min).to_numpy(dtype="float32")

# ---------------------------------------------------------------------
# Worker initialization and processing
# ---------------------------------------------------------------------
_WORKER: Dict[str, object] = {}


def _init_worker(cml_nc: str, radar_nc: str, cfg_json: str):
    global _WORKER
    cfg = json.loads(cfg_json)
    _WORKER = {
        "cfg": cfg,
        "ds_c": xr.open_dataset(cml_nc),
        "ds_r": xr.open_dataset(radar_nc),
    }
    ds_c = _WORKER["ds_c"]
    ds_r = _WORKER["ds_r"]
    if "cml_id" in ds_c.coords:
        _WORKER["cml_ids"] = ds_c["cml_id"].astype(str).values.astype(str)
    else:
        _WORKER["cml_ids"] = np.arange(ds_c.sizes["cml_id"]).astype(str)
    if "cml_id" in ds_r.coords:
        radar_ids = ds_r["cml_id"].astype(str).values.astype(str)
    else:
        radar_ids = np.arange(ds_r.sizes["cml_id"]).astype(str)
    _WORKER["radar_pos"] = {str(cid): i for i, cid in enumerate(radar_ids)}


def _read_link_arrays(ds_c: xr.Dataset, cidx: int, proc0: int, proc1: int, sublink_dim: str, tx_name: str, rx_name: str) -> Tuple[np.ndarray, np.ndarray]:
    tx_da = ds_c[tx_name].isel({"cml_id": cidx, "time": slice(proc0, proc1)})
    rx_da = ds_c[rx_name].isel({"cml_id": cidx, "time": slice(proc0, proc1)})
    if sublink_dim in tx_da.dims:
        tx = tx_da.transpose("time", sublink_dim).values
        rx = rx_da.transpose("time", sublink_dim).values
    else:
        tx = tx_da.transpose("time").values[:, None]
        rx = rx_da.transpose("time").values[:, None]
    return np.asarray(tx, dtype="float64"), np.asarray(rx, dtype="float64")


def _get_scalar_or_indexed(ds: xr.Dataset, name: str, cidx: int, sublink_dim: Optional[str] = None, sl: Optional[int] = None, default=np.nan):
    if name not in ds.variables:
        return default
    da = ds[name]
    indexer = {}
    if "cml_id" in da.dims:
        indexer["cml_id"] = cidx
    if sublink_dim and sl is not None and sublink_dim in da.dims:
        indexer[sublink_dim] = sl
    try:
        return da.isel(indexer).values
    except Exception:
        return default


def process_diagnostics_link(cidx: int) -> Tuple[Dict, List[Dict]]:
    ds_c: xr.Dataset = _WORKER["ds_c"]
    ds_r: xr.Dataset = _WORKER["ds_r"]
    cfg: Dict = _WORKER["cfg"]
    cml_ids = _WORKER["cml_ids"]
    radar_pos = _WORKER["radar_pos"]
    cid = str(cml_ids[cidx])
    sublink_dim = cfg["sublink_dim"]
    n_sublink = int(cfg["n_sublink"])
    proc_time = pd.DatetimeIndex(np.array(cfg["proc_time"], dtype="datetime64[ns]"))
    proc0, proc1 = int(cfg["proc0"]), int(cfg["proc1"])
    rproc0, rproc1 = int(cfg["rproc0"]), int(cfg["rproc1"])

    tx, rx = _read_link_arrays(ds_c, cidx, proc0, proc1, sublink_dim, cfg["tx_name"], cfg["rx_name"])
    length = scalar_float(_get_scalar_or_indexed(ds_c, cfg["length_var"], cidx, default=np.nan))

    rpos = radar_pos.get(cid)
    if rpos is None:
        radar_mm = np.full(max(1, rproc1-rproc0), np.nan)
    else:
        radar_raw = ds_r[cfg["radar_var"]].isel({"cml_id": rpos, "time": slice(rproc0, rproc1)}).values
        radar_mm = radar_values_to_mm_h(radar_raw, cfg["radar_scale_to_mmph"])
    has_ref_finite = bool(np.isfinite(radar_mm).any())
    has_ref_positive = bool(np.isfinite(radar_mm).any() and np.nanmax(radar_mm) > cfg["radar_wet_thr_5min"])

    sub_rows: List[Dict] = []
    link_has_cml_data = False
    link_pass_det = False
    link_pass_long = False
    link_pass_short = False
    link_pass_avail = False
    link_pass_blettner = False
    plateau_sum = 0
    blackout_sum = 0

    for sl in range(n_sublink):
        freq_raw = scalar_float(_get_scalar_or_indexed(ds_c, cfg["frequency_var"], cidx, sublink_dim, sl, default=np.nan))
        pol = pol_to_HV(_get_scalar_or_indexed(ds_c, cfg["polarization_var"], cidx, sublink_dim, sl, default="H"))
        info = blettner_screen_sublink(
            tx[:, sl], rx[:, sl], proc_time, length, to_hz(freq_raw), pol,
            cfg["blettner_proc_line"], int(cfg["min_valid_samples"]),
            bool(cfg.get("use_availability_filter", False)),
            float(cfg.get("availability_perc_thld", 0.5)),
        )
        row = {"cml_id": cid, "sublink_id": sl, "length_km": length, "frequency_raw": freq_raw, "polarization": pol}
        row.update(info)
        sub_rows.append(row)
        link_has_cml_data = link_has_cml_data or bool(info["has_cml_data_sublink"])
        link_pass_det = link_pass_det or bool(info["pass_detection_limit"])
        link_pass_long = link_pass_long or bool(info["pass_long_std"])
        link_pass_short = link_pass_short or bool(info["pass_short_std"])
        link_pass_avail = link_pass_avail or bool(info["pass_availability"])
        link_pass_blettner = link_pass_blettner or bool(info["pass_blettner"])
        plateau_sum += int(info.get("plateau_removed_samples", 0))
        blackout_sum += int(info.get("blackout_filled_samples", 0))

    link_row = {
        "cml_id": cid,
        "cml_index": cidx,
        "length_km": length,
        "has_cml_data": int(link_has_cml_data),
        "has_radar_reference": int(has_ref_finite),
        "has_positive_radar_reference": int(has_ref_positive),
        "has_cml_data_and_radar_reference": int(link_has_cml_data and has_ref_finite),
        "pass_detection_limit": int(link_pass_det),
        "pass_long_std": int(link_pass_long),
        "pass_short_std": int(link_pass_short),
        "pass_availability": int(link_pass_avail),
        "pass_blettner": int(link_pass_blettner),
        "final_blettner_link": int(link_has_cml_data and has_ref_finite and link_pass_blettner),
        "plateau_removed_samples_total": int(plateau_sum),
        "blackout_filled_samples_total": int(blackout_sum),
    }
    return link_row, sub_rows


def process_estimation_link(task: Tuple[int, str, bool]) -> Dict:
    cidx, cid, is_blettner_link = task
    ds_c: xr.Dataset = _WORKER["ds_c"]
    ds_r: xr.Dataset = _WORKER["ds_r"]
    cfg: Dict = _WORKER["cfg"]
    radar_pos = _WORKER["radar_pos"]
    sublink_dim = cfg["sublink_dim"]
    n_sublink = int(cfg["n_sublink"])
    proc0, proc1 = int(cfg["proc0"]), int(cfg["proc1"])
    out_rel0, out_rel1 = int(cfg["out_rel0"]), int(cfg["out_rel1"])
    rproc0, rproc1 = int(cfg["rproc0"]), int(cfg["rproc1"])
    rout0, rout1 = int(cfg["rout0"]), int(cfg["rout1"])
    proc_time = pd.DatetimeIndex(np.array(cfg["proc_time"], dtype="datetime64[ns]"))
    out_time = pd.DatetimeIndex(np.array(cfg["out_time"], dtype="datetime64[ns]"))
    radar_proc_time = pd.DatetimeIndex(np.array(cfg["radar_proc_time"], dtype="datetime64[ns]"))
    radar_out_time = pd.DatetimeIndex(np.array(cfg["radar_out_time"], dtype="datetime64[ns]"))

    tx, rx = _read_link_arrays(ds_c, cidx, proc0, proc1, sublink_dim, cfg["tx_name"], cfg["rx_name"])
    length = scalar_float(_get_scalar_or_indexed(ds_c, cfg["length_var"], cidx, default=np.nan))
    rpos = radar_pos.get(str(cid))
    radar_proc_raw = ds_r[cfg["radar_var"]].isel({"cml_id": rpos, "time": slice(rproc0, rproc1)}).values
    radar_proc = radar_values_to_mm_h(radar_proc_raw, cfg["radar_scale_to_mmph"])
    radar_out = radar_values_to_mm_h(ds_r[cfg["radar_var"]].isel({"cml_id": rpos, "time": slice(rout0, rout1)}).values, cfg["radar_scale_to_mmph"]).astype("float32")
    radar_wet_5 = np.isfinite(radar_proc) & (radar_proc >= float(cfg["radar_wet_thr_5min"]))
    radar_wet_1min = forward_fill_radar_wet_to_cml_time(radar_proc_time, proc_time, radar_wet_5)

    filter_variants = cfg["filter_variants"]
    baseline_methods = cfg["baseline_methods"]
    kr_models = cfg["kr_models"]
    waa_methods = cfg["waa_methods"]
    T5 = len(radar_out_time)
    R_mean = np.full((len(filter_variants), len(baseline_methods), len(kr_models), len(waa_methods), T5), np.nan, dtype="float32")
    R_sel = np.full(T5, np.nan, dtype="float32")
    graf_counts = {"hit": 0, "false_alarm": 0, "miss": 0, "correct_negative": 0}
    graf_wet_5_acc = np.zeros(T5, dtype=bool)
    thresholds: List[float] = []
    log_rows: List[Dict] = []

    # Process each variant independently. Raw always allowed. Blettner variant only for final links.
    for fv_i, fv in enumerate(filter_variants):
        if fv == "blettner_full" and not is_blettner_link:
            continue
        sub_R_store = np.full((n_sublink, len(baseline_methods), len(kr_models), len(waa_methods), T5), np.nan, dtype="float32")
        for sl in range(n_sublink):
            trsl_raw = tx[:, sl] - rx[:, sl]
            graf_wet, thr = graf_wet_from_trsl(trsl_raw)
            thresholds.append(thr)
            valid_det = np.isfinite(trsl_raw) & np.isfinite(radar_wet_1min)
            cc = confusion_counts(graf_wet, radar_wet_1min, valid_det)
            for k, v in cc.items():
                graf_counts[k] += int(v)
            # Save 5-min wet flag as any wet minute in the bin
            try:
                graf_5 = pd.Series(graf_wet.astype(int), index=proc_time).resample("5min", label="left", closed="left").max().reindex(radar_out_time).fillna(0).to_numpy(dtype=bool)
                graf_wet_5_acc |= graf_5
            except Exception:
                pass

            if fv == "raw":
                trsl = trsl_raw
            else:
                trsl, _ = apply_blettner_signal_filter(tx[:, sl], rx[:, sl], proc_time, cfg["blettner_proc_line"])

            freq_raw = scalar_float(_get_scalar_or_indexed(ds_c, cfg["frequency_var"], cidx, sublink_dim, sl, default=np.nan))
            f_ghz = to_ghz(freq_raw)
            pol = pol_to_HV(_get_scalar_or_indexed(ds_c, cfg["polarization_var"], cidx, sublink_dim, sl, default="H"))
            if not np.isfinite(length) or length <= 0 or not np.isfinite(f_ghz):
                continue
            wet_for_est = graf_wet if cfg.get("use_graf_wet_for_estimation", False) else radar_wet_1min
            for b_i, bmethod in enumerate(baseline_methods):
                baseline = compute_baseline(bmethod, trsl, wet_for_est, proc_time, int(cfg["n_average_last_dry"]))
                A_obs = np.maximum(trsl - baseline, 0.0)
                A_obs[~wet_for_est] = 0.0
                A_obs[~np.isfinite(A_obs)] = np.nan
                for w_i, wmethod in enumerate(waa_methods):
                    waa = compute_waa(wmethod, A_obs, trsl, baseline, wet_for_est, proc_time, f_ghz, pol, length, argparse.Namespace(**cfg))
                    A_rain = np.maximum(A_obs - waa, 0.0)
                    A_rain[~wet_for_est] = 0.0
                    for k_i, kr in enumerate(kr_models):
                        R_1 = calc_R_from_A(A_rain, kr, f_ghz, pol, length)
                        R_1 = R_1[out_rel0:out_rel1]
                        R_5 = resample_1min_rate_to_5min(R_1, out_time, radar_out_time)
                        sub_R_store[sl, b_i, k_i, w_i, :] = R_5
                        log_rows.append({"cml_id": cid, "sublink_id": sl, "filter_variant": fv, "baseline_method": bmethod, "kr_model": kr, "waa_method": wmethod, "status": "processed"})
        # link mean after sublink estimation
        with np.errstate(invalid="ignore"):
            R_mean[fv_i, :, :, :, :] = np.nanmean(sub_R_store, axis=0)

    # selected product: blettner_full if available, pycml_constant/ITU_2003/pastorek_2021
    try:
        fsel = filter_variants.index(cfg["selected_filter_variant"])
        bsel = baseline_methods.index(cfg["selected_baseline_method"])
        ksel = kr_models.index(cfg["selected_kr_model"])
        wsel = waa_methods.index(cfg["selected_waa_method"])
        R_sel = R_mean[fsel, bsel, ksel, wsel, :].astype("float32")
    except Exception:
        pass

    return {
        "cml_id": str(cid),
        "R_mean": R_mean,
        "R_selected": R_sel,
        "R_radar": radar_out,
        "graf_counts": graf_counts,
        "graf_wet_5": graf_wet_5_acc.astype("int8"),
        "graf_threshold_median": float(np.nanmedian(thresholds)) if thresholds else np.nan,
        "logs": log_rows,
    }

# ---------------------------------------------------------------------
# Diagnostics and stage summary
# ---------------------------------------------------------------------
def build_stage_summary(link_df: pd.DataFrame, n_selected: int, use_availability_filter: bool = False) -> pd.DataFrame:
    def cnt(col):
        return int(pd.to_numeric(link_df[col], errors="coerce").fillna(0).astype(int).sum()) if col in link_df.columns else 0
    with_cml = cnt("has_cml_data")
    with_ref = cnt("has_radar_reference")
    with_both = cnt("has_cml_data_and_radar_reference")
    det = int(((link_df.get("has_cml_data_and_radar_reference", 0).astype(int) == 1) & (link_df.get("pass_detection_limit", 0).astype(int) == 1)).sum()) if not link_df.empty else 0
    long_std = int(((link_df.get("has_cml_data_and_radar_reference", 0).astype(int) == 1) & (link_df.get("pass_detection_limit", 0).astype(int) == 1) & (link_df.get("pass_long_std", 0).astype(int) == 1)).sum()) if not link_df.empty else 0
    short_std = int(((link_df.get("has_cml_data_and_radar_reference", 0).astype(int) == 1) & (link_df.get("pass_detection_limit", 0).astype(int) == 1) & (link_df.get("pass_long_std", 0).astype(int) == 1) & (link_df.get("pass_short_std", 0).astype(int) == 1)).sum()) if not link_df.empty else 0
    avail = int(((link_df.get("has_cml_data_and_radar_reference", 0).astype(int) == 1) & (link_df.get("pass_detection_limit", 0).astype(int) == 1) & (link_df.get("pass_long_std", 0).astype(int) == 1) & (link_df.get("pass_short_std", 0).astype(int) == 1) & (link_df.get("pass_availability", 0).astype(int) == 1)).sum()) if not link_df.empty else 0
    final = cnt("final_blettner_link")
    rows = [
        ("all_cml_links_before_filtering", n_selected, 0),
        ("with_cml_data", with_cml, n_selected - with_cml),
        ("with_radar_reference", with_ref, n_selected - with_ref),
        ("with_cml_data_and_radar_reference", with_both, max(with_cml, with_ref) - with_both),
        ("after_detection_limit", det, with_both - det),
        ("after_long_std", long_std, det - long_std),
        ("after_short_std", short_std, long_std - short_std),
    ]
    if use_availability_filter:
        rows.append(("after_availability", avail, short_std - avail))
        rows.append(("final_blettner_links", final, with_both - final))
    else:
        rows.append(("availability_diagnostic_pass_count_not_filter", avail, short_std - avail))
        rows.append(("final_blettner_links_no_availability_filter", final, with_both - final))
    return pd.DataFrame(rows, columns=["stage", "n_links_remaining_or_count", "n_links_dropped_or_missing"])

# ---------------------------------------------------------------------
# Gridding
# ---------------------------------------------------------------------
def link_midpoints_from_dataset(ds: xr.Dataset, cml_ids: np.ndarray, lon_a_name: str, lat_a_name: str, lon_b_name: str, lat_b_name: str) -> Tuple[np.ndarray, np.ndarray, Tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]]:
    ds2 = ds.assign_coords(cml_id=ds["cml_id"].astype(str)) if "cml_id" in ds.coords else ds
    sub = ds2.sel(cml_id=cml_ids.astype(str))
    lon_a = np.asarray(sub[lon_a_name].values, dtype=float)
    lat_a = np.asarray(sub[lat_a_name].values, dtype=float)
    lon_b = np.asarray(sub[lon_b_name].values, dtype=float)
    lat_b = np.asarray(sub[lat_b_name].values, dtype=float)
    return 0.5*(lon_a+lon_b), 0.5*(lat_a+lat_b), (lon_a, lat_a, lon_b, lat_b)


def make_grid(lon: np.ndarray, lat: np.ndarray, resolution: float, pad: float = 0.05, extent: Optional[Sequence[float]] = None) -> Tuple[np.ndarray, np.ndarray]:
    if extent is None:
        ok = np.isfinite(lon) & np.isfinite(lat)
        extent = [float(lon[ok].min()-pad), float(lon[ok].max()+pad), float(lat[ok].min()-pad), float(lat[ok].max()+pad)]
    x = np.arange(extent[0], extent[1] + 0.5*resolution, resolution, dtype="float32")
    y = np.arange(extent[2], extent[3] + 0.5*resolution, resolution, dtype="float32")
    return x, y


def idw_interpolate_points(xobs: np.ndarray, yobs: np.ndarray, vals: np.ndarray, xg: np.ndarray, yg: np.ndarray, nnear: int, p: int, max_dist: float) -> np.ndarray:
    out = np.full((len(yg), len(xg)), np.nan, dtype="float32")
    ok = np.isfinite(xobs) & np.isfinite(yobs) & np.isfinite(vals)
    if ok.sum() < 3:
        return out
    pts = np.column_stack([xobs[ok], yobs[ok]])
    vv = np.asarray(vals[ok], dtype="float64")
    tree = cKDTree(pts)
    X, Y = np.meshgrid(xg, yg)
    targets = np.column_stack([X.ravel(), Y.ravel()])
    k = min(int(nnear), len(vv))
    d, ind = tree.query(targets, k=k, distance_upper_bound=float(max_dist), workers=-1)
    if k == 1:
        d = d[:, None]
        ind = ind[:, None]
    valid = np.isfinite(d) & (ind < len(vv))
    weights = np.zeros_like(d, dtype="float64")
    weights[valid] = 1.0 / np.maximum(d[valid], 1e-12) ** int(p)
    vals_near = np.zeros_like(d, dtype="float64")
    clipped = np.minimum(ind, max(0, len(vv)-1))
    vals_near[valid] = vv[clipped[valid]]
    sw = weights.sum(axis=1)
    pred = np.full(targets.shape[0], np.nan, dtype="float64")
    good = sw > 0
    pred[good] = (weights[good] * vals_near[good]).sum(axis=1) / sw[good]
    return pred.reshape(len(yg), len(xg)).astype("float32")


def grid_selected_product(selected_nc: Path, coords_nc: Path, gauge_nc: Optional[Path], out_nc: Path, args) -> Path:
    print("Gridding selected product ...")
    ds = xr.open_dataset(selected_nc)
    dsc = xr.open_dataset(coords_nc)
    try:
        cml_ids = ds["cml_id"].values.astype(str)
        times = as_datetime_index(ds["time_5min"] if "time_5min" in ds.coords else ds["time"])
        t0 = pd.Timestamp(args.grid_tstart) if args.grid_tstart else times[0]
        t1 = pd.Timestamp(args.grid_tend) if args.grid_tend else times[-1]
        i0, i1 = time_slice_indices(times, t0, t1)
        times_sel = times[i0:i1]
        R = ds["R_selected"].isel(time_5min=slice(i0, i1)).values.astype("float32")
        Rrad = ds["R_radar"].isel(time_5min=slice(i0, i1)).values.astype("float32")
        xobs, yobs, geom = link_midpoints_from_dataset(dsc, cml_ids, args.site_a_lon_var, args.site_a_lat_var, args.site_b_lon_var, args.site_b_lat_var)
        if args.plot_extent:
            extent = [float(v) for v in args.plot_extent]
        else:
            extent = None
        xg, yg = make_grid(xobs, yobs, args.resolution, extent=extent)
        R_grid = np.full((len(times_sel), len(yg), len(xg)), np.nan, dtype="float32")
        Rrad_grid = np.full_like(R_grid, np.nan)
        iterator = range(len(times_sel))
        for ti in tqdm(iterator, desc="IDW CML/radar", unit="time"):
            R_grid[ti] = idw_interpolate_points(xobs, yobs, R[:, ti], xg, yg, args.nnear, args.p, args.max_dist)
            Rrad_grid[ti] = idw_interpolate_points(xobs, yobs, Rrad[:, ti], xg, yg, args.nnear, args.p, args.max_dist)

        data_vars = {
            "R": (("time", "y", "x"), R_grid, {"units": "mm h-1", "description": "IDW gridded selected CML rainfall product"}),
            "R_radar_path": (("time", "y", "x"), Rrad_grid, {"units": "mm h-1", "description": "IDW gridded path-averaged radar values at CML midpoints"}),
        }

        if gauge_nc is not None and Path(gauge_nc).exists() and not args.no_gauge:
            try:
                dsg = xr.open_dataset(gauge_nc)
                gtime_name = find_dim(dsg, [args.gauge_time_dim, "time"], required=True)
                gst_dim = find_dim(dsg, [args.gauge_station_dim, "station", "WSI", "id"], required=True)
                gvar = args.gauge_var if args.gauge_var in dsg.variables else pick_var(dsg, ["rr", "RR", "rain", "precipitation"])
                glon_name = args.gauge_lon if args.gauge_lon in dsg.variables else pick_var(dsg, ["lon", "longitude"])
                glat_name = args.gauge_lat if args.gauge_lat in dsg.variables else pick_var(dsg, ["lat", "latitude"])
                gtimes = as_datetime_index(dsg[gtime_name])
                gi0, gi1 = time_slice_indices(gtimes, t0, t1)
                gtimes_sel = gtimes[gi0:gi1]
                gr = dsg[gvar].isel({gtime_name: slice(gi0, gi1)}).transpose(gtime_name, gst_dim).values.astype("float32")
                # default gauge is 10-min depth in mm; convert to mm/h
                if args.gauge_is_rate:
                    grate = gr
                else:
                    grate = gr * 6.0
                glon = np.asarray(dsg[glon_name].values, dtype=float)
                glat = np.asarray(dsg[glat_name].values, dtype=float)
                G_grid = np.full((len(gtimes_sel), len(yg), len(xg)), np.nan, dtype="float32")
                for ti in tqdm(range(len(gtimes_sel)), desc="IDW gauges", unit="time"):
                    G_grid[ti] = idw_interpolate_points(glon, glat, grate[ti, :], xg, yg, args.nnear, args.p, args.max_dist)
                data_vars["R_gauge"] = (("time_gauge", "y", "x"), G_grid, {"units": "mm h-1"})
                coords = {"time": times_sel.values, "time_gauge": gtimes_sel.values, "x": xg, "y": yg}
                dsg.close()
            except Exception as e:
                print(f"Gauge gridding skipped: {e}")
                coords = {"time": times_sel.values, "x": xg, "y": yg}
        else:
            coords = {"time": times_sel.values, "x": xg, "y": yg}

        out = xr.Dataset(data_vars=data_vars, coords=coords, attrs={"source_selected_nc": str(selected_nc), "idw_nnear": int(args.nnear), "idw_p": int(args.p), "idw_max_dist_deg": float(args.max_dist)})
        write_dataset(out, out_nc, args.complevel)
        print(f"Saved gridded NetCDF: {out_nc}")
        return out_nc
    finally:
        ds.close(); dsc.close()

# ---------------------------------------------------------------------
# Native RADOLAN/RADKLIM plotting
# ---------------------------------------------------------------------
def discover_native_radar_files(root: str, pattern: str, exclude_substring: str = "_CML_RR") -> Dict[object, str]:
    files = sorted(Path(root).glob("**/" + pattern))
    out = {}
    for p in files:
        if exclude_substring and exclude_substring in str(p):
            continue
        m = re.search(r"(\d{8})", p.name)
        if not m:
            continue
        date = pd.Timestamp(m.group(1)).date()
        out[date] = str(p)
    return out


def radolan_crs():
    try:
        import cartopy.crs as ccrs
        return ccrs.Stereographic(
            central_latitude=90.0,
            central_longitude=10.0,
            true_scale_latitude=60.0,
            globe=ccrs.Globe(ellipse=None, semimajor_axis=6370040, semiminor_axis=6370040),
        )
    except Exception:
        return None


def read_native_radar_field(path: str, target_time: pd.Timestamp, scale: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, pd.Timestamp]:
    ds = xr.open_dataset(path)
    try:
        var = pick_var(ds, RADAR_VAR_CANDIDATES, required=True)
        da = ds[var]
        tdim = next((d for d in da.dims if "time" in d.lower()), None)
        if tdim:
            sel = da.sel({tdim: target_time}, method="nearest")
            actual = pd.Timestamp(sel[tdim].values)
        else:
            sel = da
            actual = target_time
        field = np.asarray(sel.values, dtype="float32") * np.float32(scale)
        lon_name = next((c for c in ["lon", "longitude"] if c in ds.coords or c in ds.variables), None)
        lat_name = next((c for c in ["lat", "latitude"] if c in ds.coords or c in ds.variables), None)
        if lon_name and lat_name:
            lon = np.asarray(ds[lon_name].values, dtype="float32")
            lat = np.asarray(ds[lat_name].values, dtype="float32")
            if lon.ndim == 1 and lat.ndim == 1:
                lon, lat = np.meshgrid(lon, lat)
            return lon, lat, field, actual
        if "x" in ds.variables and "y" in ds.variables:
            X, Y = np.meshgrid(np.asarray(ds["x"].values), np.asarray(ds["y"].values))
            try:
                import cartopy.crs as ccrs
                pts = ccrs.PlateCarree().transform_points(radolan_crs(), X, Y)
                lon = pts[..., 0].astype("float32")
                lat = pts[..., 1].astype("float32")
                return lon, lat, field, actual
            except Exception as exc:
                raise RuntimeError("Native radar has x/y but cartopy is unavailable for projection conversion.") from exc
        raise RuntimeError(f"No lat/lon or projected x/y coordinates in native radar file: {path}")
    finally:
        ds.close()


def setup_map_axis(ax, extent: Sequence[float], show_xlabels: bool = True):
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
        proj = ccrs.PlateCarree()
        ax.set_extent(extent, crs=proj)
        ax.add_feature(cfeature.BORDERS, linewidth=0.55, edgecolor="0.25", zorder=4)
        ax.gridlines(crs=proj, draw_labels=False, linewidth=0.25, color="0.75", alpha=0.8, linestyle="--", zorder=0)
        step = 0.5 if max(extent[1]-extent[0], extent[3]-extent[2]) <= 3 else 1.0
        xticks = np.arange(np.ceil(extent[0]/step)*step, extent[1]+1e-9, step)
        yticks = np.arange(np.ceil(extent[2]/step)*step, extent[3]+1e-9, step)
        ax.set_xticks(xticks, crs=proj)
        ax.set_yticks(yticks, crs=proj)
        ax.xaxis.set_major_formatter(LongitudeFormatter(number_format=".1f", degree_symbol="°"))
        ax.yaxis.set_major_formatter(LatitudeFormatter(number_format=".1f", degree_symbol="°"))
    except Exception:
        ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis="both", which="major", labelsize=8.5, length=4, width=0.8)
    ax.tick_params(axis="both", which="minor", length=2.2, width=0.5)
    if not show_xlabels:
        ax.tick_params(axis="x", labelbottom=False)


def plot_links(ax, geom, extent):
    if geom is None:
        return
    try:
        import cartopy.crs as ccrs
        transform = ccrs.PlateCarree()
    except Exception:
        transform = None
    lon_a, lat_a, lon_b, lat_b = geom
    keep = np.isfinite(lon_a)&np.isfinite(lat_a)&np.isfinite(lon_b)&np.isfinite(lat_b)&(
        ((lon_a>=extent[0])&(lon_a<=extent[1])&(lat_a>=extent[2])&(lat_a<=extent[3]))|
        ((lon_b>=extent[0])&(lon_b<=extent[1])&(lat_b>=extent[2])&(lat_b<=extent[3]))
    )
    for i in np.where(keep)[0]:
        kwargs = {"color":"0.25", "linewidth":0.28, "alpha":0.35, "zorder":3}
        if transform is not None:
            kwargs["transform"] = transform
        ax.plot([lon_a[i], lon_b[i]], [lat_a[i], lat_b[i]], **kwargs)


def make_comparison_figure(cml_field, xg, yg, radar_field, rlon, rlat, extent, geom, out_png, out_pdf, args, panel_time: Optional[pd.Timestamp] = None):
    try:
        import cartopy.crs as ccrs
        proj_kw = {"projection": ccrs.PlateCarree()}
        trans = ccrs.PlateCarree()
    except Exception:
        proj_kw = {}
        trans = None
    bounds = np.asarray(args.bounds, dtype=float)
    norm = BoundaryNorm(bounds, ncolors=256, extend="both")
    cmap = mpl.colormaps[args.cmap].copy(); cmap.set_bad("white")
    fig, axes = plt.subplots(2, 1, figsize=(6.3, 8.0), subplot_kw=proj_kw, squeeze=False)
    ax1, ax2 = axes[:, 0]
    fig.subplots_adjust(left=0.12, right=0.86, bottom=0.08, top=0.94, hspace=0.16)
    cml = np.asarray(cml_field, dtype="float32").copy(); cml[cml < args.dry_threshold] = np.nan
    rad = np.asarray(radar_field, dtype="float32").copy(); rad[rad < args.dry_threshold] = np.nan
    setup_map_axis(ax1, extent, show_xlabels=False)
    plot_links(ax1, geom, extent)
    pc_kwargs = {"shading":"nearest", "cmap":cmap, "norm":norm, "zorder":2, "rasterized":True}
    if trans is not None: pc_kwargs["transform"] = trans
    pc = ax1.pcolormesh(xg, yg, cml, **pc_kwargs)
    ax1.text(0.018, 0.975, "(a) CML", transform=ax1.transAxes, ha="left", va="top", fontsize=10.5, fontweight="bold", bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="0.5", lw=0.55, alpha=0.92), zorder=6)
    setup_map_axis(ax2, extent, show_xlabels=True)
    plot_links(ax2, geom, extent)
    ax2.pcolormesh(rlon, rlat, rad, **pc_kwargs)
    ax2.text(0.018, 0.975, "(b) RADOLAN/RADKLIM", transform=ax2.transAxes, ha="left", va="top", fontsize=10.5, fontweight="bold", bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="0.5", lw=0.55, alpha=0.92), zorder=6)
    # No time subtitle by default. For debugging, user can enable --plot-time-title.
    if args.plot_time_title and panel_time is not None:
        fig.text(0.49, 0.965, panel_time.strftime("%d %B %Y, %H:%M UTC"), ha="center", va="top", fontsize=10)
    ptop, pbot = ax1.get_position(), ax2.get_position()
    cax = fig.add_axes([0.885, pbot.y0, 0.020, ptop.y1 - pbot.y0])
    cb = fig.colorbar(pc, cax=cax, extend="both", ticks=bounds)
    cb.set_label(r"Rain rate (mm h$^{-1}$)", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def plot_snapshot_and_gif(gridded_nc: Path, coords_nc: Path, out_dir: Path, args):
    print("Creating snapshot and GIF ...")
    out_dir = ensure_dir(out_dir)
    frames_dir = ensure_dir(out_dir / "frames_png")
    ds = xr.open_dataset(gridded_nc)
    dsc = xr.open_dataset(coords_nc)
    try:
        times = as_datetime_index(ds["time"])
        xg = ds["x"].values.astype("float32")
        yg = ds["y"].values.astype("float32")
        extent = args.plot_extent if args.plot_extent else [float(np.nanmin(xg)), float(np.nanmax(xg)), float(np.nanmin(yg)), float(np.nanmax(yg))]
        try:
            geom = (
                np.asarray(dsc[args.site_a_lon_var].values, dtype=float), np.asarray(dsc[args.site_a_lat_var].values, dtype=float),
                np.asarray(dsc[args.site_b_lon_var].values, dtype=float), np.asarray(dsc[args.site_b_lat_var].values, dtype=float),
            )
        except Exception:
            geom = None
        files_by_date = discover_native_radar_files(args.radar_root, args.radar_pattern, args.exclude_radar_substring)
        if not files_by_date:
            raise FileNotFoundError(f"No native RADOLAN files found below {args.radar_root} matching {args.radar_pattern}")
        target_times = pd.date_range(pd.Timestamp(args.grid_tstart), pd.Timestamp(args.grid_tend), freq=args.frame_freq)
        target_times = [t for t in target_times if times.min() <= t <= times.max()]
        if not target_times:
            target_times = list(times)
        metrics_rows = []
        pngs = []
        for t in tqdm(target_times, desc="render CML/native radar", unit="frame"):
            pos = times.get_indexer([t], method="nearest")[0]
            cml_time = times[pos]
            cml_field = np.asarray(ds["R"].isel(time=pos).values, dtype="float32")
            rfile = files_by_date.get(pd.Timestamp(t).date())
            if rfile is None:
                print(f"No native radar file for {t.date()}; skip frame")
                continue
            rlon, rlat, rfield, ractual = read_native_radar_field(rfile, pd.Timestamp(t), args.native_radar_scale_to_mmph)
            stem = pd.Timestamp(t).strftime("frame_%Y%m%d_%H%M")
            out_png = frames_dir / f"{stem}.png"
            out_pdf = out_dir / f"{stem}.pdf"
            make_comparison_figure(cml_field, xg, yg, rfield, rlon, rlat, extent, geom, out_png, out_pdf, args, panel_time=pd.Timestamp(t))
            pngs.append(out_png)
            metrics_rows.append({"time": t, "cml_time": cml_time, "native_radar_time": ractual, "cml_max": float(np.nanmax(cml_field)) if np.isfinite(cml_field).any() else np.nan, "radar_max": float(np.nanmax(rfield)) if np.isfinite(rfield).any() else np.nan, "png": out_png.name})
        pd.DataFrame(metrics_rows).to_csv(out_dir / "event_frame_metrics.csv", index=False)
        # Snapshot from requested snapshot time or first frame
        if pngs:
            snap_t = pd.Timestamp(args.snapshot_time) if args.snapshot_time else pd.Timestamp(target_times[0])
            pos = times.get_indexer([snap_t], method="nearest")[0]
            cml_field = np.asarray(ds["R"].isel(time=pos).values, dtype="float32")
            rfile = files_by_date.get(snap_t.date())
            if rfile:
                rlon, rlat, rfield, ractual = read_native_radar_field(rfile, snap_t, args.native_radar_scale_to_mmph)
                make_comparison_figure(cml_field, xg, yg, rfield, rlon, rlat, extent, geom, out_dir / "snapshot_cml_vs_native_radolan.png", out_dir / "snapshot_cml_vs_native_radolan.pdf", args, panel_time=snap_t)
        if args.make_gif and imageio is not None and pngs:
            imgs = [imageio.imread(p) for p in pngs]
            imageio.mimsave(out_dir / "animation_cml_vs_native_radolan.gif", imgs, duration=float(args.gif_duration))
    finally:
        ds.close(); dsc.close()

# ---------------------------------------------------------------------
# Hexbin CML vs radar
# ---------------------------------------------------------------------
def finalize_stats(c: np.ndarray, r: np.ndarray) -> Dict[str, float]:
    ok = np.isfinite(c) & np.isfinite(r)
    if not ok.any():
        return {"N": 0, "MBE": np.nan, "MAE": np.nan, "RMSE": np.nan, "PCC": np.nan}
    cc = c[ok].astype(float); rr = r[ok].astype(float); d = cc - rr
    pcc = np.corrcoef(cc, rr)[0,1] if len(cc) >= 3 and np.nanstd(cc) > 0 and np.nanstd(rr) > 0 else np.nan
    return {"N": int(len(cc)), "MBE": float(np.nanmean(d)), "MAE": float(np.nanmean(np.abs(d))), "RMSE": float(np.sqrt(np.nanmean(d*d))), "PCC": float(pcc) if np.isfinite(pcc) else np.nan, "mean_cml": float(np.nanmean(cc)), "mean_ref": float(np.nanmean(rr))}


def make_hexbin(selected_nc: Path, out_dir: Path, args):
    print("Creating CML-vs-radar hexbin ...")
    out_dir = ensure_dir(out_dir)
    ds = xr.open_dataset(selected_nc)
    try:
        c = ds["R_selected"].values.reshape(-1)
        r = ds["R_radar"].values.reshape(-1)
        ok = np.isfinite(c) & np.isfinite(r) & (c >= args.hexbin_min) & (r >= args.hexbin_min) & (c <= args.hexbin_max) & (r <= args.hexbin_max)
        c = c[ok]; r = r[ok]
        metrics = finalize_stats(c, r)
        pd.DataFrame([metrics]).to_csv(out_dir / "hexbin_metrics.csv", index=False)
        fig, ax = plt.subplots(figsize=(5.7, 5.3))
        if c.size:
            hb = ax.hexbin(r, c, gridsize=args.hexbin_gridsize, extent=(0, args.hexbin_max, 0, args.hexbin_max), mincnt=1, bins="log", cmap=HEXBIN_CMAP, linewidths=0)
            cb = fig.colorbar(hb, ax=ax, pad=0.02)
            cb.set_label("Number of samples", fontsize=10)
        ax.plot([0, args.hexbin_max], [0, args.hexbin_max], color="0.25", lw=1.0, ls="--")
        ax.set_xlim(0, args.hexbin_max); ax.set_ylim(0, args.hexbin_max)
        ax.set_xlabel(r"Radar reference (mm h$^{-1}$)", fontsize=12)
        ax.set_ylabel(r"CML rainfall estimate (mm h$^{-1}$)", fontsize=12)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.minorticks_on(); ax.xaxis.set_minor_locator(AutoMinorLocator(2)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        txt = f"N = {metrics['N']:,}\nRMSE = {metrics['RMSE']:.2f}\nMBE = {metrics['MBE']:.2f}\nPCC = {metrics['PCC']:.2f}"
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, ha="left", va="top", fontsize=9, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", alpha=0.9))
        fig.tight_layout()
        fig.savefig(out_dir / "hexbin_cml_vs_radar.png", dpi=600, bbox_inches="tight")
        fig.savefig(out_dir / "hexbin_cml_vs_radar.pdf", bbox_inches="tight")
        plt.close(fig)
    finally:
        ds.close()

# ---------------------------------------------------------------------
# Demo dataset
# ---------------------------------------------------------------------
def make_demo_dataset(out_dir: Path, n_cml: int = 80, n_gauge: int = 30):
    out_dir = ensure_dir(out_dir)
    rng = np.random.default_rng(42)
    time = pd.date_range("2021-06-20", "2021-06-23", freq="1min", inclusive="left")
    time5 = pd.date_range("2021-06-20", "2021-06-23", freq="5min", inclusive="left")
    cml_ids = np.array([f"CML_{i:04d}" for i in range(n_cml)], dtype=object)
    sub = np.array([0,1], dtype=int)
    lon_a = rng.uniform(12.0, 18.0, n_cml); lat_a = rng.uniform(49.0, 50.8, n_cml)
    lon_b = lon_a + rng.normal(0, 0.05, n_cml); lat_b = lat_a + rng.normal(0, 0.03, n_cml)
    length = np.sqrt((lon_b-lon_a)**2 + (lat_b-lat_a)**2)*90 + 0.5
    freq = rng.choice([18000, 23000, 38000, 74000], size=(n_cml,2)).astype("float32")
    pol = np.array([["H", "V"] for _ in range(n_cml)], dtype=object)
    base_tx = rng.normal(20, 1, (n_cml,2))
    base_rx = rng.normal(-45, 2, (n_cml,2))
    t_hours = (time - time[0]).total_seconds().to_numpy()/3600.0
    rain = 8*np.exp(-0.5*((t_hours-46)/3)**2) + 3*np.exp(-0.5*((t_hours-26)/2)**2)
    spatial = rng.uniform(0.5, 1.5, n_cml)
    tx = np.repeat(base_tx[None,:,:], len(time), axis=0).transpose(1,0,2)
    rx = np.repeat(base_rx[None,:,:], len(time), axis=0).transpose(1,0,2)
    for i in range(n_cml):
        att = rain * spatial[i] * length[i] * 0.025
        rx[i,:,:] -= att[:,None]
    tx += rng.normal(0,0.03,tx.shape); rx += rng.normal(0,0.07,rx.shape)
    ds_c = xr.Dataset(
        data_vars={"tx": (("cml_id","time","sublink_id"), tx.astype("float32")), "rx": (("cml_id","time","sublink_id"), rx.astype("float32")), "length": ("cml_id", length.astype("float32")), "frequency": (("cml_id","sublink_id"), freq), "polarization": (("cml_id","sublink_id"), pol), "site_a_longitude": ("cml_id", lon_a.astype("float32")), "site_a_latitude": ("cml_id", lat_a.astype("float32")), "site_b_longitude": ("cml_id", lon_b.astype("float32")), "site_b_latitude": ("cml_id", lat_b.astype("float32"))},
        coords={"cml_id": cml_ids, "time": time, "sublink_id": sub},
    )
    # path radar 5-min, raw convention requires /12 if scale=12
    rain5 = pd.Series(rain, index=time).resample("5min").mean().reindex(time5).to_numpy()
    rr = (rain5[None,:] * spatial[:,None] / 12.0).astype("float32")
    ds_r = xr.Dataset({"rr": (("cml_id","time"), rr)}, coords={"cml_id": cml_ids, "time": time5})
    st = np.array([f"G{i:03d}" for i in range(n_gauge)], dtype=object)
    glon = rng.uniform(12.0,18.0,n_gauge); glat = rng.uniform(49.0,50.8,n_gauge)
    timeg = pd.date_range("2021-06-20", "2021-06-23", freq="10min", inclusive="left")
    gra = pd.Series(rain, index=time).resample("10min").mean().reindex(timeg).to_numpy()/6.0
    gvals = np.repeat(gra[:,None], n_gauge, axis=1) * rng.uniform(0.6,1.4,n_gauge)[None,:]
    ds_g = xr.Dataset({"rr": (("time","station"), gvals.astype("float32")), "lon": ("station", glon.astype("float32")), "lat": ("station", glat.astype("float32"))}, coords={"time": timeg, "station": st})
    write_dataset(ds_c, out_dir / "demo_CML_S2021_1MIN_CZ.nc")
    write_dataset(ds_r, out_dir / "demo_RADAR_S2021_5MIN_CZ.nc")
    write_dataset(ds_g, out_dir / "demo_CHMI_10MIN_CZ.nc")
    print(f"Demo dataset saved in {out_dir}")

# ---------------------------------------------------------------------
# Main pipeline orchestration
# ---------------------------------------------------------------------
def run_pipeline(args):
    t0 = time.perf_counter()
    out = ensure_dir(args.out_dir)
    d_diag = ensure_dir(out / "diagnostics")
    d_wet = ensure_dir(out / "wet_dry")
    d_rain = ensure_dir(out / "rainfall")
    d_grid = ensure_dir(out / "gridded")
    d_plot = ensure_dir(out / "plots")

    print("Opening datasets for metadata ...")
    ds_c = xr.open_dataset(args.cml_nc)
    ds_r = xr.open_dataset(args.radar_nc)
    try:
        cml_dim = find_dim(ds_c, ["cml_id"])
        time_dim = find_dim(ds_c, ["time"])
        sublink_dim = find_dim(ds_c, ["sublink_id", "sublink"])
        r_time_dim = find_dim(ds_r, ["time"])
        r_cml_dim = find_dim(ds_r, ["cml_id"])
        tx_name = args.tx_var if args.tx_var in ds_c.variables else pick_var(ds_c, ["tx", "tsl", "TSL"])
        rx_name = args.rx_var if args.rx_var in ds_c.variables else pick_var(ds_c, ["rx", "rsl", "RSL"])
        radar_var = args.radar_var if args.radar_var in ds_r.variables else pick_var(ds_r, ["rr", "R_radar", "rain", "rainfall_rate", "R"])
        cml_ids_all = ds_c[cml_dim].astype(str).values.astype(str)
        radar_ids = ds_r[r_cml_dim].astype(str).values.astype(str)
        radar_id_set = set(radar_ids)
        full_time = as_datetime_index(ds_c[time_dim])
        radar_time = as_datetime_index(ds_r[r_time_dim])
        tstart = pd.Timestamp(args.tstart) if args.tstart else full_time[0]
        tend = pd.Timestamp(args.tend) if args.tend else full_time[-1]
        pstart = tstart - pd.Timedelta(hours=float(args.baseline_pad_hours))
        proc0, proc1 = time_slice_indices(full_time, pstart, tend)
        out0, out1 = time_slice_indices(full_time, tstart, tend)
        rproc0, rproc1 = time_slice_indices(radar_time, pstart, tend)
        rout0, rout1 = time_slice_indices(radar_time, tstart, tend)
        proc_time = full_time[proc0:proc1]
        out_time = full_time[out0:out1]
        radar_proc_time = radar_time[rproc0:rproc1]
        radar_out_time = radar_time[rout0:rout1]
        out_rel0 = out0 - proc0
        out_rel1 = out1 - proc0
        n_sublink = ds_c.sizes[sublink_dim]
        # Do NOT pre-filter links by radar-ID intersection or link length here.
        # Keep all CML links for diagnostics, then report/drop them explicitly in
        # the stage summary: CML data availability, radar reference availability,
        # detection limit, long/short std, availability, and final Blettner links.
        selected_indices = []
        for i, cid in enumerate(cml_ids_all):
            if args.only_cml_id and str(cid) != str(args.only_cml_id):
                continue
            selected_indices.append(i)
            if args.max_links and len(selected_indices) >= int(args.max_links):
                break
        print(f"Selected for diagnostics before filtering: {len(selected_indices)} / {len(cml_ids_all)}")
        cfg = {
            "proc0": proc0, "proc1": proc1, "out_rel0": out_rel0, "out_rel1": out_rel1,
            "rproc0": rproc0, "rproc1": rproc1, "rout0": rout0, "rout1": rout1,
            "proc_time": [str(t) for t in proc_time], "out_time": [str(t) for t in out_time],
            "radar_proc_time": [str(t) for t in radar_proc_time], "radar_out_time": [str(t) for t in radar_out_time],
            "sublink_dim": sublink_dim, "n_sublink": int(n_sublink), "tx_name": tx_name, "rx_name": rx_name,
            "radar_var": radar_var, "radar_scale_to_mmph": float(args.radar_scale_to_mmph), "radar_wet_thr_5min": float(args.radar_wet_thr_5min),
            "length_var": args.length_var, "frequency_var": args.frequency_var, "polarization_var": args.polarization_var,
            "min_valid_samples": int(args.min_valid_samples), "blettner_proc_line": args.blettner_proc_line,
            "use_availability_filter": bool(args.use_availability_filter), "availability_perc_thld": float(args.availability_threshold),
            "filter_variants": args.filter_variants, "baseline_methods": args.baseline_methods, "kr_models": args.kr_models, "waa_methods": args.waa_methods,
            "n_average_last_dry": int(args.n_average_last_dry),
            "waa_gamma": args.waa_gamma, "waa_l_antenna": args.waa_l_antenna, "waa_delta": args.waa_delta,
            "pastorek_a_max": args.pastorek_a_max, "pastorek_zeta": args.pastorek_zeta, "pastorek_d": args.pastorek_d,
            "schleiss_waa_max": args.schleiss_waa_max, "schleiss_delta_t": args.schleiss_delta_t, "schleiss_tau": args.schleiss_tau,
            "selected_filter_variant": args.selected_filter_variant, "selected_baseline_method": args.selected_baseline_method,
            "selected_kr_model": args.selected_kr_model, "selected_waa_method": args.selected_waa_method,
            "use_graf_wet_for_estimation": bool(args.use_graf_wet_for_estimation),
        }
    finally:
        ds_c.close(); ds_r.close()

    cfg_json = json.dumps(cfg)
    cpu_count = os.cpu_count() or 1
    workers = cpu_count if int(args.n_workers) <= 0 else int(args.n_workers)
    workers = max(1, workers)
    chunksize = max(1, int(args.mp_chunksize))
    print(f"Multiprocessing workers: {workers} requested / {cpu_count} logical CPUs available; chunksize={chunksize}")
    print("Running diagnostics with multiprocessing ...")
    diag_rows = []
    sub_rows = []
    with cf.ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(args.cml_nc, args.radar_nc, cfg_json)) as ex:
        iterator = ex.map(process_diagnostics_link, selected_indices, chunksize=chunksize)
        for lr, sr in tqdm(iterator, total=len(selected_indices), desc="diagnostics", unit="link"):
            diag_rows.append(lr); sub_rows.extend(sr)
    link_df = pd.DataFrame(diag_rows).sort_values("cml_index")
    sub_df = pd.DataFrame(sub_rows)
    stage_df = build_stage_summary(link_df, len(selected_indices), bool(args.use_availability_filter))
    link_df.to_csv(d_diag / "link_diagnostics.csv", index=False)
    sub_df.to_csv(d_diag / "sublink_diagnostics.csv", index=False)
    stage_df.to_csv(d_diag / "stage_summary.csv", index=False)
    print(stage_df.to_string(index=False))

    has_both_ids = link_df.loc[link_df["has_cml_data_and_radar_reference"].astype(int) == 1, "cml_id"].astype(str).to_numpy()
    final_ids = link_df.loc[link_df["final_blettner_link"].astype(int) == 1, "cml_id"].astype(str).to_numpy()
    print(f"Links with CML data + radar reference: {len(has_both_ids)}")
    print(f"Final Blettner links: {len(final_ids)}")

    # Save filtered CML datasets
    print("Saving filtered CML datasets ...")
    dsc = xr.open_dataset(args.cml_nc)
    try:
        dsc = dsc.assign_coords(cml_id=dsc["cml_id"].astype(str)) if "cml_id" in dsc.coords else dsc
        if len(has_both_ids):
            write_dataset(dsc.sel(cml_id=has_both_ids), d_diag / "CML_filtered_has_data_and_radar_reference.nc", args.complevel)
        if len(final_ids):
            write_dataset(dsc.sel(cml_id=final_ids), d_diag / "CML_filtered_blettner_full.nc", args.complevel)
    finally:
        dsc.close()

    if args.diagnostics_only:
        return

    # Estimation set: raw variant uses has_both, blettner variant uses final. Workers return NaNs for blettner if not final.
    id_to_idx = {str(row.cml_id): int(row.cml_index) for row in link_df.itertuples()}
    selected_for_est = has_both_ids
    tasks = [(id_to_idx[cid], cid, cid in set(final_ids)) for cid in selected_for_est]
    T5 = len(pd.DatetimeIndex(np.array(cfg["radar_out_time"], dtype="datetime64[ns]")))
    nC = len(selected_for_est)
    nF, nB, nK, nW = len(args.filter_variants), len(args.baseline_methods), len(args.kr_models), len(args.waa_methods)
    R_all_path = d_rain / "cml_rainrate_all_experiments.nc"
    R_sel_path = d_rain / "cml_rainrate_selected_pastorek_itu03.nc"

    R_all = np.full((nC, nF, nB, nK, nW, T5), np.nan, dtype="float32") if args.save_all_experiments else None
    R_sel = np.full((nC, T5), np.nan, dtype="float32")
    R_rad = np.full((nC, T5), np.nan, dtype="float32")
    graf_wet_5 = np.zeros((nC, T5), dtype="int8") if args.save_wet_flags else None
    all_logs: List[Dict] = []
    global_counts = {"hit": 0, "false_alarm": 0, "miss": 0, "correct_negative": 0}

    print("Running rainfall estimation + Graf wet/dry detection ...")
    with cf.ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(args.cml_nc, args.radar_nc, cfg_json)) as ex:
        iterator = ex.map(process_estimation_link, tasks, chunksize=chunksize)
        for i, outres in enumerate(tqdm(iterator, total=len(tasks), desc="estimation", unit="link")):
            if R_all is not None:
                R_all[i] = outres["R_mean"]
            R_sel[i] = outres["R_selected"]
            R_rad[i] = outres["R_radar"]
            if graf_wet_5 is not None:
                graf_wet_5[i] = outres["graf_wet_5"]
            for k in global_counts:
                global_counts[k] += int(outres["graf_counts"].get(k, 0))
            all_logs.extend(outres["logs"])
    pd.DataFrame(all_logs).to_csv(d_rain / "cml_rainfall_estimation_log.csv", index=False)

    # Wet/dry outputs
    pd.DataFrame([global_counts]).to_csv(d_wet / "graf_vs_radar_contingency.csv", index=False)
    metrics = detection_metrics_from_counts(global_counts)
    pd.DataFrame([metrics]).to_csv(d_wet / "graf_vs_radar_metrics.csv", index=False)
    plot_contingency(global_counts, d_wet / "graf_vs_radar_contingency.png")
    if graf_wet_5 is not None:
        ds_w = xr.Dataset({"graf_wet_5min": (("cml_id","time_5min"), graf_wet_5, {"description":"Graf rolling std wet flag resampled to 5-min; 1 wet, 0 dry"})}, coords={"cml_id": selected_for_est, "time_5min": pd.DatetimeIndex(np.array(cfg["radar_out_time"], dtype="datetime64[ns]")).values})
        write_dataset(ds_w, d_wet / "graf_wet_detection_5min.nc", args.complevel)

    coords_common = {"cml_id": selected_for_est, "time_5min": pd.DatetimeIndex(np.array(cfg["radar_out_time"], dtype="datetime64[ns]")).values}
    ds_sel = xr.Dataset({
        "R_selected": (("cml_id","time_5min"), R_sel, {"units":"mm h-1", "description":"Selected CML product"}),
        "R_radar": (("cml_id","time_5min"), R_rad, {"units":"mm h-1", "description":"Path averaged radar reference"}),
    }, coords=coords_common, attrs={"selected_filter_variant": args.selected_filter_variant, "selected_baseline_method": args.selected_baseline_method, "selected_kr_model": args.selected_kr_model, "selected_waa_method": args.selected_waa_method})
    write_dataset(ds_sel, R_sel_path, args.complevel)
    print(f"Saved selected rainfall product: {R_sel_path}")

    if R_all is not None:
        ds_all = xr.Dataset({
            "R_mean_5min": (("cml_id","filter_variant","baseline_method","kr_model","waa_method","time_5min"), R_all, {"units":"mm h-1"}),
            "R_radar": (("cml_id","time_5min"), R_rad, {"units":"mm h-1"}),
        }, coords={"cml_id": selected_for_est, "filter_variant": args.filter_variants, "baseline_method": args.baseline_methods, "kr_model": args.kr_models, "waa_method": args.waa_methods, "time_5min": coords_common["time_5min"]})
        write_dataset(ds_all, R_all_path, args.complevel)
        print(f"Saved all experiments rainfall product: {R_all_path}")

    gridded_path = d_grid / "rainfall_maps_selected_event.nc"
    if not args.skip_gridding:
        grid_selected_product(R_sel_path, d_diag / "CML_filtered_has_data_and_radar_reference.nc", Path(args.gauge_nc) if args.gauge_nc else None, gridded_path, args)
    if not args.skip_plots:
        plot_snapshot_and_gif(gridded_path, d_diag / "CML_filtered_has_data_and_radar_reference.nc", d_plot, args)
        make_hexbin(R_sel_path, d_plot, args)

    report = {
        "elapsed_seconds": time.perf_counter() - t0,
        "n_links_selected_basic": len(selected_indices),
        "n_links_with_cml_data_and_radar_reference": int(len(has_both_ids)),
        "n_final_blettner_links": int(len(final_ids)),
        "wet_dry_metrics": metrics,
        "paths": {
            "diagnostics": str(d_diag), "wet_dry": str(d_wet), "rainfall": str(d_rain), "gridded": str(d_grid), "plots": str(d_plot),
        },
    }
    (out / "pipeline_report.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"DONE. Total elapsed: {report['elapsed_seconds']:.1f} s")

# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Full fast CML rainfall pipeline with Blettner diagnostics, Graf wet/dry, rainfall experiments, gridding and plots.")
    p.add_argument("--cml-nc", default="CML_S2020E2021_1MIN_CZ.nc")
    p.add_argument("--radar-nc", default="RADAR_S2021E2021_5MIN_CZ.nc")
    p.add_argument("--gauge-nc", default="CHMI_10MIN_v1.nc")
    p.add_argument("--out-dir", default="full_cml_pipeline_output")
    p.add_argument("--tstart", default=None)
    p.add_argument("--tend", default=None)
    p.add_argument("--grid-tstart", default="2021-06-21 21:50:00")
    p.add_argument("--grid-tend", default="2021-06-22 04:50:00")
    p.add_argument("--baseline-pad-hours", type=float, default=240.0)
    p.add_argument("--min-link-km", type=float, default=0.5)
    p.add_argument("--min-valid-samples", type=int, default=100)
    p.add_argument("--availability-threshold", type=float, default=0.5,
                   help="Diagnostic availability threshold: pass if NaN fraction is below this value. Default 0.5.")
    p.add_argument("--use-availability-filter", action="store_true",
                   help="Old behavior: require the availability check for full Blettner filtering. Default is off; availability is diagnostic only.")
    p.add_argument("--radar-wet-thr-5min", type=float, default=0.00833)
    p.add_argument("--radar-scale-to-mmph", type=float, default=12.0)
    p.add_argument("--native-radar-scale-to-mmph", type=float, default=12.0)
    p.add_argument("--tx-var", default="tx")
    p.add_argument("--rx-var", default="rx")
    p.add_argument("--radar-var", default="rr")
    p.add_argument("--length-var", default="length")
    p.add_argument("--frequency-var", default="frequency")
    p.add_argument("--polarization-var", default="polarization")
    p.add_argument("--site-a-lon-var", default="site_a_longitude")
    p.add_argument("--site-a-lat-var", default="site_a_latitude")
    p.add_argument("--site-b-lon-var", default="site_b_longitude")
    p.add_argument("--site-b-lat-var", default="site_b_latitude")
    p.add_argument("--filter-variants", nargs="+", default=["raw", "blettner_full"], choices=["raw", "blettner_full"])
    p.add_argument("--blettner-proc-line", default="full", choices=["no_filter", "graf_2020", "full"])
    p.add_argument("--baseline-methods", nargs="+", default=["pycml_constant"], choices=["pycml_constant", "last_dry"])
    p.add_argument("--kr-models", nargs="+", default=["ITU_2003", "ITU_2005"], choices=["ITU_2003", "ITU_2005"])
    p.add_argument("--waa-methods", nargs="+", default=["none", "leijnse_2008", "pastorek_2021", "schleiss_2013"], choices=["none", "leijnse_2008", "pastorek_2021", "schleiss_2013"])
    p.add_argument("--selected-filter-variant", default="blettner_full")
    p.add_argument("--selected-baseline-method", default="pycml_constant")
    p.add_argument("--selected-kr-model", default="ITU_2003")
    p.add_argument("--selected-waa-method", default="pastorek_2021")
    p.add_argument("--n-average-last-dry", type=int, default=1)
    p.add_argument("--waa-gamma", type=float, default=(2.06e-05)/1.4)
    p.add_argument("--waa-l-antenna", type=float, default=0.0043)
    p.add_argument("--waa-delta", type=float, default=0.24*1.5)
    p.add_argument("--pastorek-a-max", type=float, default=14.0)
    p.add_argument("--pastorek-zeta", type=float, default=0.55)
    p.add_argument("--pastorek-d", type=float, default=0.1)
    p.add_argument("--schleiss-waa-max", type=float, default=2.0)
    p.add_argument("--schleiss-delta-t", type=float, default=1.0)
    p.add_argument("--schleiss-tau", type=float, default=15.0)
    p.add_argument("--use-graf-wet-for-estimation", action="store_true")
    p.add_argument("--save-all-experiments", action="store_true", help="Save large all-experiment rainfall NetCDF. Can be very large for full dataset.")
    p.add_argument("--save-wet-flags", action="store_true")
    p.add_argument("--diagnostics-only", action="store_true")
    p.add_argument("--skip-gridding", action="store_true")
    p.add_argument("--skip-plots", action="store_true")
    p.add_argument("--only-cml-id", default=None)
    p.add_argument("--max-links", type=int, default=None)
    p.add_argument("--n-workers", type=int, default=0,
                   help="Number of worker processes. Use 0 for all logical CPUs. Old script default was capped at 8.")
    p.add_argument("--mp-chunksize", type=int, default=2,
                   help="ProcessPool map chunksize. Increase to 4-8 to reduce scheduling overhead; use 1 for best load balancing.")
    p.add_argument("--complevel", type=int, default=4)
    # gridding/gauge
    p.add_argument("--resolution", type=float, default=0.01)
    p.add_argument("--nnear", type=int, default=15)
    p.add_argument("--p", type=int, default=2)
    p.add_argument("--max-dist", type=float, default=0.3)
    p.add_argument("--no-gauge", action="store_true")
    p.add_argument("--gauge-var", default="rr")
    p.add_argument("--gauge-time-dim", default="time")
    p.add_argument("--gauge-station-dim", default="station")
    p.add_argument("--gauge-lon", default="lon")
    p.add_argument("--gauge-lat", default="lat")
    p.add_argument("--gauge-is-rate", action="store_true")
    # native radar plotting
    p.add_argument("--radar-root", default="dwd_radolan/2021")
    p.add_argument("--radar-pattern", default="YW_2017.002_*.nc")
    p.add_argument("--exclude-radar-substring", default="_CML_RR")
    p.add_argument("--snapshot-time", default=None)
    p.add_argument("--frame-freq", default="10min")
    p.add_argument("--make-gif", action="store_true", default=True)
    p.add_argument("--no-gif", dest="make_gif", action="store_false")
    p.add_argument("--gif-duration", type=float, default=0.35)
    p.add_argument("--bounds", nargs="+", type=float, default=RAIN_BOUNDS)
    p.add_argument("--accum-bounds", nargs="+", type=float, default=ACCUM_BOUNDS)
    p.add_argument("--cmap", default="turbo")
    p.add_argument("--dry-threshold", type=float, default=0.1)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--plot-time-title", action="store_true", help="By default no time subtitle is drawn, matching requested publication style.")
    p.add_argument("--plot-extent", nargs=4, type=float, default=None)
    # hexbin
    p.add_argument("--hexbin-min", type=float, default=0.1)
    p.add_argument("--hexbin-max", type=float, default=50.0)
    p.add_argument("--hexbin-gridsize", type=int, default=50)
    # demo
    p.add_argument("--make-demo", action="store_true")
    p.add_argument("--demo-out-dir", default="demo_cml_pipeline")
    p.add_argument("--demo-n-cml", type=int, default=80)
    p.add_argument("--demo-n-gauge", type=int, default=30)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.make_demo:
        make_demo_dataset(Path(args.demo_out_dir), args.demo_n_cml, args.demo_n_gauge)
        return
    # Normalize p in case notebook passes 2.0 as string via shell
    args.p = int(float(args.p))
    run_pipeline(args)


if __name__ == "__main__":
    main()
