#!/usr/bin/env python3
"""
Optimized dual wet/dry CML rainfall retrieval for the Czech dataset.

Two compact NetCDF products are written in one streamed run:

1) GRAF product at native CML resolution
   * full preprocessing retained from Example_Processing
   * wet/dry classification read from graf_rstd_wet_1min and
     graf_rstd_valid_1min
   * Pastorek et al. (2021) wet-antenna correction
   * explicit ITU k-R inversion

2) Radar-conditioned product at the native radar resolution
   * the same processed CML TRSL is aggregated to exact radar intervals
   * baseline and rainfall retrieval are performed AFTER temporal aggregation
   * the path-associated radar rain rate supplies the wet/dry mask
   * no radar or CML temporal filling is performed

Only CML IDs that exist in the radar file and have at least one finite,
non-negative radar observation in the requested output period are written.
Each sublink is estimated independently; the link mean is calculated afterward.

Frequency and length conventions
--------------------------------
* Input frequency is converted once to Hz. Pastorek receives Hz.
* ITU coefficient lookup and inversion receive GHz = Hz / 1e9.
* Link length is converted once to km.
* calc_R_from_A returns mm h-1 directly; there is no multiplication by 60.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd
import xarray as xr

try:
    from numba import njit, prange
except Exception as exc:  # pragma: no cover
    raise ImportError("numba is required for optimized preprocessing") from exc

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    def tqdm(x: Iterable, **_: object):
        return x

SCRIPT_VERSION = "1.0"
ALLOWED_KR_MODELS = ("ITU_2003", "ITU_2005")


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------
def pol_to_HV(pol_raw: object) -> str:
    if isinstance(pol_raw, (bytes, np.bytes_)):
        pol_raw = bytes(pol_raw).decode("utf-8", errors="ignore")
    s = str(pol_raw).replace("\x00", "").strip().lower()
    if "v" in s and "h" not in s:
        return "V"
    if "h" in s and "v" not in s:
        return "H"
    return "H"


def _decode_one(value: object) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return bytes(value).decode("utf-8", errors="ignore").replace("\x00", "").strip()
    return str(value).replace("\x00", "").strip()


def _decode_string_coord(da: xr.DataArray, primary_dim: str) -> np.ndarray:
    """Decode one string per primary_dim, including fixed S1 character axes."""
    if primary_dim not in da.dims:
        raise ValueError(f"{da.name!r} does not contain dimension {primary_dim!r}")
    arr = np.asarray(da.values)
    dims = list(da.dims)
    paxis = dims.index(primary_dim)
    if paxis != 0:
        arr = np.moveaxis(arr, paxis, 0)
    n = arr.shape[0]
    if arr.ndim == 1:
        return np.asarray([_decode_one(v) for v in arr], dtype=object)
    flat = arr.reshape(n, -1)
    out = np.empty(n, dtype=object)
    for i, row in enumerate(flat):
        out[i] = "".join(_decode_one(v) for v in row).strip()
    return out


def _get_ids(ds: xr.Dataset, dim: str) -> np.ndarray:
    if dim in ds.variables:
        return _decode_string_coord(ds[dim], dim)
    if "cml_id" in ds.variables and dim in ds["cml_id"].dims:
        return _decode_string_coord(ds["cml_id"], dim)
    return np.asarray([str(i) for i in range(ds.sizes[dim])], dtype=object)


def _find_dim(dims: Sequence[str], candidates: Sequence[str], label: str) -> str:
    for name in candidates:
        if name in dims:
            return name
    raise KeyError(f"Could not find {label} dimension. Available dimensions: {tuple(dims)}")


def _find_var(ds: xr.Dataset, preferred: str, alternatives: Sequence[str], label: str) -> str:
    for name in (preferred, *alternatives):
        if name and name in ds.variables:
            return name
    raise KeyError(
        f"Could not find {label}. Tried {(preferred, *alternatives)}. "
        f"Available variables include {list(ds.variables)[:80]}"
    )


def _as_datetime_index(values: np.ndarray) -> pd.DatetimeIndex:
    idx = pd.DatetimeIndex(pd.to_datetime(np.asarray(values)))
    if idx.tz is not None:
        idx = idx.tz_convert("UTC").tz_localize(None)
    return idx


def _time_ns(values: Sequence[object]) -> np.ndarray:
    return _as_datetime_index(np.asarray(values)).to_numpy(dtype="datetime64[ns]").astype(np.int64)


def _median_step_ns(times: pd.DatetimeIndex, fallback_minutes: float) -> int:
    ns = _time_ns(times)
    d = np.diff(ns)
    d = d[d > 0]
    if d.size == 0:
        return int(round(float(fallback_minutes) * 60 * 1e9))
    return int(np.median(d))


def _frequency_to_hz(values: np.ndarray, unit: str) -> tuple[np.ndarray, str]:
    f = np.asarray(values, dtype=np.float64)
    if unit == "Hz":
        return f, "Hz"
    if unit == "MHz":
        return f * 1e6, "MHz"
    if unit == "GHz":
        return f * 1e9, "GHz"
    finite = np.abs(f[np.isfinite(f)])
    if finite.size == 0:
        raise ValueError("Frequency metadata contains no finite values")
    med = float(np.nanmedian(finite))
    if med > 1e8:
        return f, "Hz"
    if med > 1e3:
        return f * 1e6, "MHz"
    return f * 1e9, "GHz"


def _length_to_km(values: np.ndarray, unit: str) -> tuple[np.ndarray, str]:
    x = np.asarray(values, dtype=np.float64)
    if unit == "km":
        return x, "km"
    if unit == "m":
        return x / 1000.0, "m"
    finite = np.abs(x[np.isfinite(x)])
    if finite.size == 0:
        raise ValueError("Length metadata contains no finite values")
    inferred = "m" if float(np.nanmedian(finite)) > 100.0 else "km"
    return (x / 1000.0 if inferred == "m" else x), inferred


def _metadata_numeric_2d(
    da: xr.DataArray,
    sublink_dim: str,
    cml_dim: str,
    cml_positions: np.ndarray,
    n_sub: int,
) -> np.ndarray:
    work = da.isel({cml_dim: cml_positions}) if cml_dim in da.dims else da
    if sublink_dim not in work.dims:
        work = work.expand_dims({sublink_dim: np.arange(n_sub)})
    if cml_dim not in work.dims:
        work = work.expand_dims({cml_dim: np.arange(len(cml_positions))})
    extra = [d for d in work.dims if d not in {sublink_dim, cml_dim}]
    if extra:
        raise ValueError(f"Numeric metadata {da.name!r} has unsupported extra dimensions {extra}")
    arr = np.asarray(work.transpose(sublink_dim, cml_dim).values, dtype=np.float64)
    return np.broadcast_to(arr, (n_sub, len(cml_positions))).copy()


def _metadata_string_2d(
    da: xr.DataArray,
    sublink_dim: str,
    cml_dim: str,
    cml_positions: np.ndarray,
    n_sub: int,
) -> np.ndarray:
    work = da.isel({cml_dim: cml_positions}) if cml_dim in da.dims else da
    stored_dims = list(work.dims)
    arr = np.asarray(work.values)
    for dim, size in ((sublink_dim, n_sub), (cml_dim, len(cml_positions))):
        if dim not in stored_dims:
            arr = np.expand_dims(arr, axis=0)
            stored_dims.insert(0, dim)
    extras = [d for d in stored_dims if d not in {sublink_dim, cml_dim}]
    desired = [sublink_dim, cml_dim] + extras
    axes = [stored_dims.index(d) for d in desired]
    arr = np.transpose(arr, axes=axes)
    arr = np.broadcast_to(arr, (n_sub, len(cml_positions)) + arr.shape[2:])
    flat = arr.reshape(n_sub, len(cml_positions), -1)
    out = np.empty((n_sub, len(cml_positions)), dtype=object)
    for s in range(n_sub):
        for c in range(len(cml_positions)):
            out[s, c] = "".join(_decode_one(v) for v in flat[s, c]).strip()
    return out


def _read_signal_block(
    da: xr.DataArray,
    time_dim: str,
    sublink_dim: str,
    cml_dim: str,
    time_slice: slice,
    cml_positions: np.ndarray,
    dtype=np.float64,
) -> np.ndarray:
    work = da.isel({time_dim: time_slice, cml_dim: cml_positions}).transpose(
        sublink_dim, cml_dim, time_dim
    )
    arr = np.asarray(work.values, dtype=dtype)
    return arr.reshape(arr.shape[0] * arr.shape[1], arr.shape[2])


# -----------------------------------------------------------------------------
# Full preprocessing retained from v5
# -----------------------------------------------------------------------------
def detect_plateaus_vectorized(
    rsl: xr.DataArray,
    time_dim: str,
    time_span: int,
    max_thld: float,
    std_thld: float,
    pad: int,
) -> xr.DataArray:
    low = rsl.rolling({time_dim: int(time_span)}, center=True).max() < max_thld
    steady = rsl.rolling({time_dim: int(time_span)}, center=True).std() < std_thld
    sane = ~(low & steady)
    sane = ~(sane.rolling({time_dim: int(pad)}, center=True).min() == 0)
    return sane.fillna(True).astype(bool)


@njit(cache=True)
def _blackout_mask_one(rsl: np.ndarray, threshold: float, max_gap: int) -> np.ndarray:
    n = rsl.size
    below = np.zeros(n, dtype=np.bool_)
    isn = np.zeros(n, dtype=np.bool_)
    for i in range(n):
        isn[i] = np.isnan(rsl[i])
        below[i] = (not isn[i]) and (rsl[i] < threshold)
    gap_start = np.zeros(n, dtype=np.bool_)
    gap_end = np.zeros(n, dtype=np.bool_)
    for i in range(n):
        gap_start[i] = isn[(i + 1) % n] and below[i]
        gap_end[i] = isn[(i - 1 + n) % n] and below[i]
    mask = np.zeros(n, dtype=np.bool_)
    in_gap = False
    count = 0
    for i in range(n):
        if gap_start[i]:
            in_gap = True
            count = 0
        if (not gap_start[i]) and (not isn[i]):
            in_gap = False
        if gap_end[i]:
            in_gap = False
        if in_gap and isn[i] and count < max_gap:
            mask[i] = True
            count += 1
    return mask


@njit(cache=True, parallel=True)
def blackout_filling_parallel(
    tsl_in: np.ndarray,
    rsl_in: np.ndarray,
    tsl_max: np.ndarray,
    rsl_min: np.ndarray,
    threshold: float,
    max_gap: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_series, n_time = rsl_in.shape
    tsl_out = tsl_in.copy()
    rsl_out = rsl_in.copy()
    filled = np.zeros((n_series, n_time), dtype=np.bool_)
    for s in prange(n_series):
        forward = _blackout_mask_one(rsl_in[s], threshold, max_gap)
        reverse = _blackout_mask_one(rsl_in[s, ::-1].copy(), threshold, max_gap)
        for t in range(n_time):
            use = forward[t] or reverse[n_time - 1 - t]
            if use:
                filled[s, t] = True
                tsl_out[s, t] = tsl_max[s]
                rsl_out[s, t] = rsl_min[s]
    return tsl_out, rsl_out, filled


@njit(cache=True, parallel=True)
def interpolate_short_gaps_parallel(
    values: np.ndarray,
    time_ns: np.ndarray,
    max_gap_ns: int,
) -> np.ndarray:
    out = values.copy()
    n_series, n_time = out.shape
    for s in prange(n_series):
        i = 0
        while i < n_time:
            if not np.isnan(out[s, i]):
                i += 1
                continue
            start = i
            while i < n_time and np.isnan(out[s, i]):
                i += 1
            end = i
            left = start - 1
            right = end
            if left >= 0 and right < n_time:
                y0 = out[s, left]
                y1 = out[s, right]
                gap = time_ns[right] - time_ns[left]
                if not np.isnan(y0) and not np.isnan(y1) and 0 < gap <= max_gap_ns:
                    for j in range(start, end):
                        frac = float(time_ns[j] - time_ns[left]) / float(gap)
                        out[s, j] = y0 + frac * (y1 - y0)
    return out


def temporal_sanity_check_vectorized(
    trsl: xr.DataArray,
    time_dim: str,
    window: int,
    threshold: float,
    fraction: float,
) -> xr.DataArray:
    roll_std = trsl.rolling({time_dim: int(window)}, center=True).std()
    high = roll_std > threshold
    low = roll_std <= threshold
    denom = high.sum(time_dim) + low.sum(time_dim)
    ratio = high.sum(time_dim) / denom
    return (ratio <= fraction).fillna(False).astype(bool)


def availability_check_vectorized(
    trsl: xr.DataArray,
    time_dim: str,
    missing_fraction_threshold: float,
) -> xr.DataArray:
    return (
        (trsl.isnull().sum(time_dim) / trsl.sizes[time_dim]) < missing_fraction_threshold
    ).astype(bool)


# -----------------------------------------------------------------------------
# k-R and Pastorek retrieval
# -----------------------------------------------------------------------------
def _normalize_kr_model(model: str) -> str:
    model = "ITU_2005" if str(model).upper() == "ITU" else str(model)
    if model not in ALLOWED_KR_MODELS:
        raise ValueError(f"Unsupported k-R model {model!r}")
    return model


def _call_with_kr_model(fn, *, model: str, **kwargs):
    model = _normalize_kr_model(model)
    errors = []
    for keyword in ("a_b_approximation", "approx_type"):
        try:
            return fn(**kwargs, **{keyword: model})
        except (TypeError, ValueError) as exc:
            errors.append(f"{keyword}: {type(exc).__name__}: {exc}")
    raise RuntimeError(
        f"pycomlink rejected {model!r} with both keyword variants: {' | '.join(errors)}"
    )


def _ab_coefficients(pycml, f_ghz: float, pol: str, model: str) -> tuple[float, float]:
    a, b = _call_with_kr_model(
        pycml.processing.k_R_relation.a_b,
        model=model,
        f_GHz=float(f_ghz),
        pol=pol_to_HV(pol),
    )
    return float(np.asarray(a).squeeze()), float(np.asarray(b).squeeze())


def detection_limits(
    pycml,
    length_km: np.ndarray,
    frequency_hz: np.ndarray,
    polarization: np.ndarray,
    quantization_db: float,
    kr_model: str,
) -> np.ndarray:
    out = np.full(length_km.shape, np.nan, dtype=np.float64)
    flat_l = length_km.reshape(-1)
    flat_f = frequency_hz.reshape(-1)
    flat_p = polarization.reshape(-1)
    flat_o = out.reshape(-1)
    cache: dict[tuple[float, str], float] = {}
    for i in range(flat_o.size):
        L = float(flat_l[i])
        f_hz = float(flat_f[i])
        pol = pol_to_HV(flat_p[i])
        if not (np.isfinite(L) and L > 0 and np.isfinite(f_hz) and f_hz > 0):
            continue
        key = (round(f_hz / 1e9, 6), pol)
        if key not in cache:
            a, _ = _ab_coefficients(pycml, key[0], pol, kr_model)
            cache[key] = a
        sensitivity = cache[key] * L
        if np.isfinite(sensitivity) and sensitivity > 0:
            flat_o[i] = float(quantization_db) / sensitivity
    return out


def _calc_rain(pycml, A_rain: xr.DataArray, L_km: float, f_hz: float, pol: str, model: str):
    return _call_with_kr_model(
        pycml.processing.k_R_relation.calc_R_from_A,
        model=model,
        A=A_rain,
        L_km=float(L_km),
        f_GHz=float(f_hz) / 1e9,
        pol=pol_to_HV(pol),
        R_min=0,
    )


@dataclass
class Estimate:
    R: np.ndarray
    error: str = ""


def estimate_pastorek_one(
    pycml,
    trsl: np.ndarray,
    wet: np.ndarray,
    valid: np.ndarray,
    time_values: np.ndarray,
    length_km: float,
    frequency_hz: float,
    polarization: str,
    kr_model: str,
    n_average_last_dry: int,
    pastorek_a_max: float,
    pastorek_zeta: float,
    pastorek_d: float,
) -> Estimate:
    n = len(trsl)
    empty = np.full(n, np.nan, dtype=np.float32)
    try:
        valid = np.asarray(valid, dtype=bool) & np.isfinite(trsl)
        wet = np.asarray(wet, dtype=bool) & valid
        if valid.sum() < 2:
            return Estimate(empty, "fewer than two valid observations")

        trsl_use = np.where(valid, np.asarray(trsl, dtype=np.float64), np.nan)
        xr_trsl = xr.DataArray(trsl_use, dims=("time",), coords={"time": time_values})
        xr_wet = xr.DataArray(wet, dims=("time",), coords={"time": time_values})

        baseline = pycml.processing.baseline.baseline_constant(
            trsl=xr_trsl,
            wet=xr_wet,
            n_average_last_dry=int(n_average_last_dry),
        )
        A_obs = xr.where((xr_trsl - baseline) > 0, xr_trsl - baseline, 0.0)
        waa = pycml.processing.wet_antenna.waa_pastorek_2021_from_A_obs(
            A_obs=A_obs,
            f_Hz=float(frequency_hz),
            pol=pol_to_HV(polarization),
            L_km=float(length_km),
            A_max=float(pastorek_a_max),
            zeta=float(pastorek_zeta),
            d=float(pastorek_d),
        )
        # Same defensive clipping used in the user's tested multi-method script.
        waa = xr.where(waa < 0, 0.0, waa)
        waa = xr.where(waa > A_obs, A_obs, waa)
        A_rain = xr.where((A_obs - waa) > 0, A_obs - waa, 0.0)

        R = _calc_rain(
            pycml,
            A_rain=A_rain,
            L_km=float(length_km),
            f_hz=float(frequency_hz),
            pol=pol_to_HV(polarization),
            model=kr_model,
        )
        out = np.asarray(getattr(R, "values", R), dtype=np.float32)
        out[valid & ~wet] = 0.0
        out[~valid] = np.nan
        out[~np.isfinite(out)] = np.nan
        out[out < 0] = 0.0
        return Estimate(out)
    except Exception as exc:  # one bad series must not abort a network run
        return Estimate(empty, f"{type(exc).__name__}: {exc}")


# -----------------------------------------------------------------------------
# Exact temporal mapping to radar intervals
# -----------------------------------------------------------------------------
@dataclass
class RadarIntervalMap:
    radar_proc_indices: np.ndarray
    radar_proc_times: pd.DatetimeIndex
    radar_output_slice: slice
    cml_read_start: int
    cml_read_end: int
    starts_rel: np.ndarray
    ends_rel: np.ndarray
    actual_counts: np.ndarray
    expected_count: int


def _interval_bounds(times: pd.DatetimeIndex, step_ns: int, label: str) -> tuple[np.ndarray, np.ndarray]:
    t = _time_ns(times)
    if label == "start":
        return t, t + step_ns
    if label == "end":
        return t - step_ns, t
    half = step_ns // 2
    return t - half, t + (step_ns - half)


def build_radar_interval_map(
    cml_times: pd.DatetimeIndex,
    radar_times: pd.DatetimeIndex,
    output_start: pd.Timestamp,
    output_end: pd.Timestamp,
    baseline_pad_hours: float,
    timestamp_label: str,
    fallback_minutes: float,
) -> RadarIntervalMap:
    rns = _time_ns(radar_times)
    out0 = int(np.searchsorted(rns, output_start.to_datetime64().astype("datetime64[ns]").astype(np.int64), "left"))
    out1 = int(np.searchsorted(rns, output_end.to_datetime64().astype("datetime64[ns]").astype(np.int64), "right"))
    if out1 <= out0:
        raise ValueError("No radar timestamps fall inside the requested output period")
    proc_start = radar_times[out0] - pd.Timedelta(hours=float(baseline_pad_hours))
    proc0 = int(np.searchsorted(rns, proc_start.to_datetime64().astype("datetime64[ns]").astype(np.int64), "left"))
    proc_indices = np.arange(proc0, out1, dtype=np.int64)
    proc_times = radar_times[proc_indices]

    rstep = _median_step_ns(radar_times, fallback_minutes)
    cstep = _median_step_ns(cml_times, 1.0)
    expected = max(1, int(round(rstep / cstep)))
    left, right = _interval_bounds(proc_times, rstep, timestamp_label)
    cns = _time_ns(cml_times)
    starts_abs = np.searchsorted(cns, left, side="left").astype(np.int64)
    ends_abs = np.searchsorted(cns, right, side="left").astype(np.int64)
    read_start = int(starts_abs.min())
    read_end = int(ends_abs.max())
    if read_start < 0 or read_end > len(cml_times):
        raise ValueError("Computed CML read bounds are outside the CML time coordinate")
    return RadarIntervalMap(
        radar_proc_indices=proc_indices,
        radar_proc_times=proc_times,
        radar_output_slice=slice(out0 - proc0, out1 - proc0),
        cml_read_start=read_start,
        cml_read_end=read_end,
        starts_rel=starts_abs - read_start,
        ends_rel=ends_abs - read_start,
        actual_counts=ends_abs - starts_abs,
        expected_count=expected,
    )


def aggregate_complete_mean(
    values_series_time: np.ndarray,
    interval_map: RadarIntervalMap,
    min_valid_fraction: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return interval mean, validity, and finite minute count as (series,time)."""
    values = np.asarray(values_series_time, dtype=np.float64).T  # time, series
    finite = np.isfinite(values)
    safe = np.where(finite, values, 0.0)

    ps = np.empty((safe.shape[0] + 1, safe.shape[1]), dtype=np.float64)
    pc = np.empty((safe.shape[0] + 1, safe.shape[1]), dtype=np.int32)
    ps[0] = 0.0
    pc[0] = 0
    np.cumsum(safe, axis=0, out=ps[1:])
    np.cumsum(finite.astype(np.int32), axis=0, out=pc[1:])
    sums = ps[interval_map.ends_rel] - ps[interval_map.starts_rel]
    counts = pc[interval_map.ends_rel] - pc[interval_map.starts_rel]

    required = int(math.ceil(interval_map.expected_count * float(min_valid_fraction) - 1e-12))
    complete_grid = interval_map.actual_counts == interval_map.expected_count
    valid = complete_grid[:, None] & (counts >= required)
    mean = np.full(sums.shape, np.nan, dtype=np.float32)
    np.divide(sums, counts, out=mean, where=valid & (counts > 0))
    return mean.T, valid.T, np.minimum(counts.T, 255).astype(np.uint8)


# -----------------------------------------------------------------------------
# Radar link matching and availability
# -----------------------------------------------------------------------------
def _radar_block(
    radar_da: xr.DataArray,
    radar_time_dim: str,
    radar_cml_dim: str,
    time_indices: np.ndarray,
    cml_positions: np.ndarray,
    scale: float,
) -> np.ndarray:
    work = radar_da.isel({radar_time_dim: time_indices, radar_cml_dim: cml_positions}).transpose(
        radar_time_dim, radar_cml_dim
    )
    arr = np.asarray(work.values, dtype=np.float64) * float(scale)
    arr[~np.isfinite(arr) | (arr < 0)] = np.nan
    return arr.astype(np.float32)


def select_radar_available_links(
    ds_cml: xr.Dataset,
    ds_radar: xr.Dataset,
    cml_dim: str,
    radar_cml_dim: str,
    radar_time_dim: str,
    radar_var: str,
    radar_output_indices: np.ndarray,
    radar_scale: float,
    scan_batch: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cml_ids = _get_ids(ds_cml, cml_dim)
    radar_ids = _get_ids(ds_radar, radar_cml_dim)
    radar_pos = {str(cid): i for i, cid in enumerate(radar_ids)}

    src_candidates = []
    radar_candidates = []
    ids_candidates = []
    for i, cid in enumerate(cml_ids):
        rp = radar_pos.get(str(cid))
        if rp is not None:
            src_candidates.append(i)
            radar_candidates.append(rp)
            ids_candidates.append(str(cid))
    if not src_candidates:
        raise ValueError("No common CML IDs were found between CML and radar files")

    src_candidates = np.asarray(src_candidates, dtype=np.int64)
    radar_candidates = np.asarray(radar_candidates, dtype=np.int64)
    ids_candidates = np.asarray(ids_candidates, dtype=object)
    keep = np.zeros(len(src_candidates), dtype=bool)
    da = ds_radar[radar_var]
    for b0 in tqdm(range(0, len(src_candidates), int(scan_batch)), desc="radar availability", unit="batch"):
        b1 = min(b0 + int(scan_batch), len(src_candidates))
        arr = _radar_block(
            da,
            radar_time_dim,
            radar_cml_dim,
            radar_output_indices,
            radar_candidates[b0:b1],
            radar_scale,
        )
        keep[b0:b1] = np.any(np.isfinite(arr), axis=0)

    if not keep.any():
        raise ValueError("No common CML link has a finite radar observation in the requested period")
    return src_candidates[keep], radar_candidates[keep], ids_candidates[keep]


# -----------------------------------------------------------------------------
# NetCDF streamed output
# -----------------------------------------------------------------------------
def _compression(level: int, chunks: tuple[int, ...]) -> dict[str, object]:
    return {
        "zlib": int(level) > 0,
        "complevel": int(level),
        "shuffle": True,
        "chunksizes": tuple(max(1, int(x)) for x in chunks),
    }


def _write_time_var(nc, name: str, dim: str, times: pd.DatetimeIndex) -> object:
    from netCDF4 import date2num
    v = nc.createVariable(name, "i8", (dim,))
    v.units = "seconds since 1970-01-01 00:00:00 UTC"
    v.calendar = "standard"
    v[:] = np.asarray(date2num(times.to_pydatetime(), v.units, v.calendar), dtype=np.int64)
    return v


@dataclass
class OutputHandles:
    nc: object
    R: object
    R_mean: object
    sane: object
    sane_long: object
    sane_short: object
    sane_available: object
    sane_detection: object
    detection_limit: object
    plateau: Optional[object] = None
    blackout: Optional[object] = None
    wet: Optional[object] = None
    wet_valid: Optional[object] = None
    radar_rain: Optional[object] = None
    radar_wet: Optional[object] = None
    radar_valid: Optional[object] = None
    trsl_interval: Optional[object] = None
    interval_valid: Optional[object] = None
    valid_count: Optional[object] = None


def _create_base_metadata(
    nc,
    times: pd.DatetimeIndex,
    time_dim: str,
    cml_ids: np.ndarray,
    source_cml_positions: np.ndarray,
    source_radar_positions: np.ndarray,
    sublink_values: np.ndarray,
    length_km: np.ndarray,
    frequency_input: np.ndarray,
    frequency_hz: np.ndarray,
    polarization: np.ndarray,
    time_chunk: int,
    link_chunk: int,
    compression_level: int,
    title: str,
):
    ntime = len(times)
    nsub, nlink = length_km.shape
    nc.createDimension(time_dim, ntime)
    nc.createDimension("sublink_id", nsub)
    nc.createDimension("cml_id", nlink)
    _write_time_var(nc, time_dim, time_dim, times)

    vc = nc.createVariable("cml_id", str, ("cml_id",))
    vc[:] = np.asarray(cml_ids, dtype=object)
    vs = nc.createVariable("sublink_id", str, ("sublink_id",))
    vs[:] = np.asarray([str(x) for x in sublink_values], dtype=object)
    vsrc = nc.createVariable("source_cml_index", "i4", ("cml_id",))
    vsrc[:] = source_cml_positions.astype(np.int32)
    vradsrc = nc.createVariable("source_radar_index", "i4", ("cml_id",))
    vradsrc[:] = source_radar_positions.astype(np.int32)

    chunks2 = (1, min(nlink, int(link_chunk)))
    vlen = nc.createVariable("length_km", "f4", ("sublink_id", "cml_id"), **_compression(compression_level, chunks2))
    vlen.units = "km"
    vlen[:] = length_km.astype(np.float32)
    vf0 = nc.createVariable("frequency_input", "f8", ("sublink_id", "cml_id"), **_compression(compression_level, chunks2))
    vf0[:] = frequency_input.astype(np.float64)
    vfg = nc.createVariable("frequency_GHz", "f4", ("sublink_id", "cml_id"), **_compression(compression_level, chunks2))
    vfg.units = "GHz"
    vfg[:] = (frequency_hz / 1e9).astype(np.float32)
    vp = nc.createVariable("polarization", str, ("sublink_id", "cml_id"))
    vp[:] = np.asarray(polarization, dtype=object)

    nc.setncattr("title", title)
    nc.setncattr("script", Path(__file__).name)
    nc.setncattr("script_version", SCRIPT_VERSION)
    nc.setncattr("created_utc", pd.Timestamp.utcnow().isoformat())
    nc.setncattr("only_links_with_radar_available", "true")
    nc.setncattr("waa_model", "Pastorek_2021")
    nc.setncattr("rain_rate_units", "mm h-1")
    nc.setncattr("rain_rate_formula", "R=(A_rain/(a*L_km))**(1/b)")
    nc.setncattr("frequency_internal_units", "Hz for Pastorek; GHz for ITU k-R")


def create_graf_output(
    path: Path,
    args,
    times: pd.DatetimeIndex,
    cml_ids: np.ndarray,
    source_cml_positions: np.ndarray,
    source_radar_positions: np.ndarray,
    sublink_values: np.ndarray,
    length_km: np.ndarray,
    frequency_input: np.ndarray,
    frequency_hz: np.ndarray,
    polarization: np.ndarray,
) -> OutputHandles:
    from netCDF4 import Dataset
    nc = Dataset(path, "w", format="NETCDF4")
    _create_base_metadata(
        nc, times, "time", cml_ids, source_cml_positions, source_radar_positions,
        sublink_values, length_km, frequency_input, frequency_hz, polarization,
        args.time_chunk_1min, args.output_link_chunk, args.compression_level,
        "CML rainfall retrieval using GRAF wet/dry classification and Pastorek WAA",
    )
    nt, ns, nl = len(times), len(sublink_values), len(cml_ids)
    c3 = (min(nt, args.time_chunk_1min), 1, min(nl, args.output_link_chunk))
    c2 = (min(nt, args.time_chunk_1min), min(nl, args.output_link_chunk))
    cs = (1, min(nl, args.output_link_chunk))

    R = nc.createVariable("R_cml", "f4", ("time", "sublink_id", "cml_id"), fill_value=np.float32(np.nan), **_compression(args.compression_level, c3))
    R.units = "mm h-1"
    R.long_name = "CML rain rate estimated independently for each sublink"
    Rm = nc.createVariable("R_link_mean", "f4", ("time", "cml_id"), fill_value=np.float32(np.nan), **_compression(args.compression_level, c2))
    Rm.units = "mm h-1"
    wet = nc.createVariable("graf_wet", "u1", ("time", "sublink_id", "cml_id"), fill_value=np.uint8(255), **_compression(args.compression_level, c3))
    wet.flag_values = np.asarray([0, 1], dtype=np.uint8)
    wet.flag_meanings = "dry wet"
    wetv = nc.createVariable("graf_wet_valid", "u1", ("time", "sublink_id", "cml_id"), **_compression(args.compression_level, c3))
    wetv.flag_values = np.asarray([0, 1], dtype=np.uint8)
    wetv.flag_meanings = "invalid valid"
    plateau = nc.createVariable("isSane_plateaus_rsl", "u1", ("time", "sublink_id", "cml_id"), **_compression(args.compression_level, c3))
    blackout = nc.createVariable("flag_blackout_fill", "u1", ("time", "sublink_id", "cml_id"), **_compression(args.compression_level, c3))

    sane = nc.createVariable("isSane", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    slong = nc.createVariable("isSane_std_long", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    sshort = nc.createVariable("isSane_std_short", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    savail = nc.createVariable("isSane_available", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    sdetect = nc.createVariable("isSane_detection_limit", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    dlim = nc.createVariable("detection_limit", "f4", ("sublink_id", "cml_id"), fill_value=np.float32(np.nan), **_compression(args.compression_level, cs))
    dlim.units = "mm h-1"

    nc.setncattr("wet_dry_source", args.graf_wet_var)
    nc.setncattr("wet_dry_valid_source", args.graf_valid_var)
    nc.setncattr("temporal_resolution", "native CML resolution")
    nc.setncattr("temporal_resampling", "none")
    return OutputHandles(nc, R, Rm, sane, slong, sshort, savail, sdetect, dlim, plateau, blackout, wet, wetv)


def create_radar_output(
    path: Path,
    args,
    times: pd.DatetimeIndex,
    cml_ids: np.ndarray,
    source_cml_positions: np.ndarray,
    source_radar_positions: np.ndarray,
    sublink_values: np.ndarray,
    length_km: np.ndarray,
    frequency_input: np.ndarray,
    frequency_hz: np.ndarray,
    polarization: np.ndarray,
) -> OutputHandles:
    from netCDF4 import Dataset
    nc = Dataset(path, "w", format="NETCDF4")
    _create_base_metadata(
        nc, times, "time_5min", cml_ids, source_cml_positions, source_radar_positions,
        sublink_values, length_km, frequency_input, frequency_hz, polarization,
        args.time_chunk_radar, args.output_link_chunk, args.compression_level,
        "CML rainfall retrieval using radar wet/dry classification and Pastorek WAA",
    )
    nt, ns, nl = len(times), len(sublink_values), len(cml_ids)
    c3 = (min(nt, args.time_chunk_radar), 1, min(nl, args.output_link_chunk))
    c2 = (min(nt, args.time_chunk_radar), min(nl, args.output_link_chunk))
    cs = (1, min(nl, args.output_link_chunk))

    R = nc.createVariable("R_cml", "f4", ("time_5min", "sublink_id", "cml_id"), fill_value=np.float32(np.nan), **_compression(args.compression_level, c3))
    R.units = "mm h-1"
    Rm = nc.createVariable("R_link_mean", "f4", ("time_5min", "cml_id"), fill_value=np.float32(np.nan), **_compression(args.compression_level, c2))
    Rm.units = "mm h-1"
    rr = nc.createVariable("radar_rain", "f4", ("time_5min", "cml_id"), fill_value=np.float32(np.nan), **_compression(args.compression_level, c2))
    rr.units = "mm h-1"
    rw = nc.createVariable("radar_wet", "u1", ("time_5min", "cml_id"), fill_value=np.uint8(255), **_compression(args.compression_level, c2))
    rw.flag_values = np.asarray([0, 1], dtype=np.uint8)
    rw.flag_meanings = "dry wet"
    rv = nc.createVariable("radar_valid", "u1", ("time_5min", "cml_id"), **_compression(args.compression_level, c2))
    tr = nc.createVariable("trsl_5min", "f4", ("time_5min", "sublink_id", "cml_id"), fill_value=np.float32(np.nan), **_compression(args.compression_level, c3))
    tr.units = "dB"
    iv = nc.createVariable("cml_interval_valid", "u1", ("time_5min", "sublink_id", "cml_id"), **_compression(args.compression_level, c3))
    vc = nc.createVariable("cml_valid_minute_count", "u1", ("time_5min", "sublink_id", "cml_id"), **_compression(args.compression_level, c3))

    sane = nc.createVariable("isSane", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    slong = nc.createVariable("isSane_std_long", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    sshort = nc.createVariable("isSane_std_short", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    savail = nc.createVariable("isSane_available", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    sdetect = nc.createVariable("isSane_detection_limit", "u1", ("sublink_id", "cml_id"), **_compression(args.compression_level, cs))
    dlim = nc.createVariable("detection_limit", "f4", ("sublink_id", "cml_id"), fill_value=np.float32(np.nan), **_compression(args.compression_level, cs))
    dlim.units = "mm h-1"

    nc.setncattr("wet_dry_source", "path-associated radar rain rate")
    nc.setncattr("radar_scale", float(args.radar_scale))
    nc.setncattr("radar_wet_threshold_mmph", float(args.radar_wet_threshold))
    nc.setncattr("radar_timestamp_label", args.radar_timestamp_label)
    nc.setncattr("temporal_resolution", "native radar resolution")
    nc.setncattr("temporal_workflow", "processed CML TRSL aggregated first; baseline, WAA, and k-R inversion applied afterward")
    nc.setncattr("expected_cml_minutes_per_radar_interval", int(args._interval_expected_count))
    nc.setncattr("min_cml_valid_fraction", float(args.min_cml_valid_fraction))
    nc.setncattr("missing_policy", "no temporal filling; incomplete CML intervals and missing radar remain missing")
    return OutputHandles(
        nc, R, Rm, sane, slong, sshort, savail, sdetect, dlim,
        radar_rain=rr, radar_wet=rw, radar_valid=rv,
        trsl_interval=tr, interval_valid=iv, valid_count=vc,
    )


def _safe_nanmean_sublinks(values_time_sub_link: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values_time_sub_link)
    count = finite.sum(axis=1)
    total = np.nansum(values_time_sub_link, axis=1, dtype=np.float64)
    out = np.full(total.shape, np.nan, dtype=np.float32)
    np.divide(total, count, out=out, where=count > 0)
    return out


# -----------------------------------------------------------------------------
# Main run
# -----------------------------------------------------------------------------
def run(args: argparse.Namespace) -> None:
    try:
        import pycomlink as pycml
        from netCDF4 import Dataset  # noqa: F401
    except Exception as exc:
        raise RuntimeError("pycomlink and netCDF4 are required in the runtime environment") from exc

    started = time.perf_counter()
    input_path = Path(args.input_nc).resolve()
    radar_path = Path(args.radar_nc).resolve()
    graf_path = Path(args.output_graf_nc).resolve()
    radar_out_path = Path(args.output_radar_nc).resolve()
    if graf_path == radar_out_path:
        raise ValueError("The two output paths must differ")
    if input_path in {graf_path, radar_out_path}:
        raise ValueError("Outputs must not overwrite the input CML file")
    for p in (graf_path, radar_out_path):
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Add --overwrite to replace it")
        partial = Path(str(p) + ".partial")
        if partial.exists():
            partial.unlink()
        if p.exists() and args.overwrite:
            p.unlink()

    print(f"Running {Path(__file__).name} version {SCRIPT_VERSION}")
    print(f"Input CML             : {input_path}")
    print(f"Input radar           : {radar_path}")
    print("WAA model             : Pastorek 2021")
    print(f"k-R model             : {args.kr_model}")

    ds_cml = xr.open_dataset(input_path)
    ds_radar = xr.open_dataset(radar_path)
    graf_handles = None
    radar_handles = None
    graf_partial = Path(str(graf_path) + ".partial")
    radar_partial = Path(str(radar_out_path) + ".partial")
    log_rows: list[dict[str, object]] = []
    try:
        cml_dim = args.cml_dim if args.cml_dim in ds_cml.dims else _find_dim(ds_cml.dims, ("cml_id", "cml", "link"), "CML")
        sublink_dim = args.sublink_dim if args.sublink_dim in ds_cml.dims else _find_dim(ds_cml.dims, ("sublink_id", "sublink"), "sublink")
        time_dim = args.time_dim if args.time_dim in ds_cml.dims else _find_dim(ds_cml.dims, ("time",), "CML time")
        radar_time_dim = args.radar_time_dim if args.radar_time_dim in ds_radar.dims else _find_dim(ds_radar.dims, ("time", "time_5min"), "radar time")
        radar_cml_dim = args.radar_cml_dim if args.radar_cml_dim in ds_radar.dims else _find_dim(ds_radar.dims, ("cml_id", "cml", "link"), "radar CML")

        tsl_name = _find_var(ds_cml, args.tsl_var, ("tx", "tsl"), "TSL")
        rsl_name = _find_var(ds_cml, args.rsl_var, ("rx", "rsl"), "RSL")
        graf_wet_name = _find_var(ds_cml, args.graf_wet_var, ("graf_rstd_wet_1min",), "GRAF wet flag")
        graf_valid_name = _find_var(ds_cml, args.graf_valid_var, ("graf_rstd_valid_1min",), "GRAF validity flag")
        frequency_name = _find_var(ds_cml, args.frequency_var, ("frequency", "freq"), "frequency")
        length_name = _find_var(ds_cml, args.length_var, ("length",), "length")
        polarization_name = _find_var(ds_cml, args.polarization_var, ("polarization", "pol"), "polarization")
        radar_var_name = _find_var(ds_radar, args.radar_var, ("rr", "radar_mmph"), "radar rain")

        cml_times = _as_datetime_index(ds_cml[time_dim].values)
        radar_times = _as_datetime_index(ds_radar[radar_time_dim].values)
        overlap_start = max(cml_times[0], radar_times[0])
        overlap_end = min(cml_times[-1], radar_times[-1])
        output_start = pd.Timestamp(args.tstart) if args.tstart else overlap_start
        output_end = pd.Timestamp(args.tend) if args.tend else overlap_end
        if output_end < output_start:
            raise ValueError("Requested time range has no CML/radar overlap")

        interval_map = build_radar_interval_map(
            cml_times,
            radar_times,
            output_start,
            output_end,
            args.baseline_pad_hours,
            args.radar_timestamp_label,
            args.radar_fallback_minutes,
        )
        args._interval_expected_count = interval_map.expected_count
        radar_output_indices = interval_map.radar_proc_indices[interval_map.radar_output_slice]
        radar_output_times = interval_map.radar_proc_times[interval_map.radar_output_slice]

        # GRAF process window includes the same baseline pad and ends at output_end.
        cns = _time_ns(cml_times)
        graf_output_start_idx = int(np.searchsorted(cns, output_start.to_datetime64().astype("datetime64[ns]").astype(np.int64), "left"))
        graf_output_end_idx = int(np.searchsorted(cns, output_end.to_datetime64().astype("datetime64[ns]").astype(np.int64), "right"))
        graf_proc_start_time = output_start - pd.Timedelta(hours=float(args.baseline_pad_hours))
        graf_proc_start_idx = int(np.searchsorted(cns, graf_proc_start_time.to_datetime64().astype("datetime64[ns]").astype(np.int64), "left"))
        # One common CML read window must cover both branches.
        common_read_start = min(graf_proc_start_idx, interval_map.cml_read_start)
        common_read_end = max(graf_output_end_idx, interval_map.cml_read_end)
        common_times = cml_times[common_read_start:common_read_end]
        graf_output_slice = slice(graf_output_start_idx - common_read_start, graf_output_end_idx - common_read_start)
        radar_starts_common = interval_map.starts_rel + interval_map.cml_read_start - common_read_start
        radar_ends_common = interval_map.ends_rel + interval_map.cml_read_start - common_read_start
        interval_common = RadarIntervalMap(
            interval_map.radar_proc_indices,
            interval_map.radar_proc_times,
            interval_map.radar_output_slice,
            common_read_start,
            common_read_end,
            radar_starts_common,
            radar_ends_common,
            interval_map.actual_counts,
            interval_map.expected_count,
        )

        src_positions, radar_positions, selected_ids = select_radar_available_links(
            ds_cml,
            ds_radar,
            cml_dim,
            radar_cml_dim,
            radar_time_dim,
            radar_var_name,
            radar_output_indices,
            args.radar_scale,
            args.availability_scan_batch,
        )
        nlink = len(src_positions)
        nsub = ds_cml.sizes[sublink_dim]
        sublink_values = np.asarray(ds_cml[sublink_dim].values if sublink_dim in ds_cml.variables else np.arange(nsub))

        print(f"Requested period       : {output_start} -> {output_end}")
        print(f"CML source samples     : {len(common_times):,} including baseline pad")
        print(f"GRAF output samples    : {graf_output_end_idx-graf_output_start_idx:,}")
        print(f"Radar output samples   : {len(radar_output_times):,}")
        print(f"Expected CML/radar bin : {interval_map.expected_count} native CML samples")
        print(f"Radar-available links  : {nlink:,}")
        print(f"Sublinks               : {nsub}")

        length_input = _metadata_numeric_2d(ds_cml[length_name], sublink_dim, cml_dim, src_positions, nsub)
        length_km, inferred_length = _length_to_km(length_input, args.length_unit)
        frequency_input = _metadata_numeric_2d(ds_cml[frequency_name], sublink_dim, cml_dim, src_positions, nsub)
        frequency_hz, inferred_frequency = _frequency_to_hz(frequency_input, args.frequency_unit)
        polarization_raw = _metadata_string_2d(ds_cml[polarization_name], sublink_dim, cml_dim, src_positions, nsub)
        polarization = np.asarray([[pol_to_HV(v) for v in row] for row in polarization_raw], dtype=object)

        print(f"Frequency input unit   : {inferred_frequency}")
        print(f"Length input unit      : {inferred_length}")
        print(f"Polarization           : H={int(np.sum(polarization=='H')):,}, V={int(np.sum(polarization=='V')):,}")

        dlimit = detection_limits(
            pycml,
            length_km,
            frequency_hz,
            polarization,
            args.a_quantization_db,
            args.kr_model,
        )
        sane_detection_all = dlimit < args.detection_limit_threshold

        graf_handles = create_graf_output(
            graf_partial, args, cml_times[graf_output_start_idx:graf_output_end_idx],
            selected_ids, src_positions, radar_positions, sublink_values,
            length_km, frequency_input, frequency_hz, polarization,
        )
        radar_handles = create_radar_output(
            radar_partial, args, radar_output_times,
            selected_ids, src_positions, radar_positions, sublink_values,
            length_km, frequency_input, frequency_hz, polarization,
        )
        for h in (graf_handles, radar_handles):
            h.nc.setncattr("source_cml_file", str(input_path))
            h.nc.setncattr("source_radar_file", str(radar_path))
            h.nc.setncattr("frequency_input_unit", inferred_frequency)
            h.nc.setncattr("length_input_unit", inferred_length)
            h.nc.setncattr("kr_model", args.kr_model)
            h.nc.setncattr("baseline_method", "pycomlink baseline_constant")
            h.nc.setncattr("n_average_last_dry", int(args.n_average_last_dry))
            h.nc.setncattr("pastorek_A_max", float(args.pastorek_a_max))
            h.nc.setncattr("pastorek_zeta", float(args.pastorek_zeta))
            h.nc.setncattr("pastorek_d", float(args.pastorek_d))
            h.nc.setncattr("preprocessing_line", "full")

        common_time_ns = _time_ns(common_times)
        max_gap_ns = int(pd.Timedelta(args.interpolation_max_gap).value)
        batch_starts = range(0, nlink, int(args.link_batch))
        pbar = tqdm(batch_starts, total=math.ceil(nlink / int(args.link_batch)), desc="dual rainfall", unit="batch")
        for b0 in pbar:
            b1 = min(b0 + int(args.link_batch), nlink)
            src = src_positions[b0:b1]
            radpos = radar_positions[b0:b1]
            nb = b1 - b0
            batch_started = time.perf_counter()

            tsl = _read_signal_block(ds_cml[tsl_name], time_dim, sublink_dim, cml_dim, slice(common_read_start, common_read_end), src)
            rsl = _read_signal_block(ds_cml[rsl_name], time_dim, sublink_dim, cml_dim, slice(common_read_start, common_read_end), src)

            target = xr.DataArray(
                rsl.reshape(nsub, nb, len(common_times)),
                dims=(sublink_dim, cml_dim, time_dim),
                coords={
                    sublink_dim: np.arange(nsub),
                    cml_dim: np.arange(nb),
                    time_dim: common_times,
                },
            )
            plateau_da = detect_plateaus_vectorized(
                target,
                time_dim,
                args.plateau_time_span,
                args.plateau_max_threshold,
                args.plateau_std_threshold,
                args.plateau_pad,
            )
            plateau = np.asarray(plateau_da.values, dtype=bool).reshape(nsub * nb, len(common_times))
            tsl[~plateau] = np.nan
            rsl[~plateau] = np.nan

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                tsl_max = np.nanmax(tsl, axis=1)
                rsl_min = np.nanmin(rsl, axis=1)
            tsl, rsl, blackout = blackout_filling_parallel(
                tsl, rsl, tsl_max, rsl_min,
                args.blackout_rsl_threshold,
                args.blackout_max_gap,
            )
            trsl_raw = tsl - rsl
            trsl = interpolate_short_gaps_parallel(trsl_raw, common_time_ns, max_gap_ns)

            trsl_da = xr.DataArray(
                trsl.reshape(nsub, nb, len(common_times)),
                dims=(sublink_dim, cml_dim, time_dim),
                coords={
                    sublink_dim: np.arange(nsub),
                    cml_dim: np.arange(nb),
                    time_dim: common_times,
                },
            )
            sane_long = np.asarray(
                temporal_sanity_check_vectorized(
                    trsl_da, time_dim,
                    args.std_long_window,
                    args.std_long_threshold,
                    args.std_long_fraction,
                ).values,
                dtype=bool,
            ).reshape(nsub * nb)
            sane_short = np.asarray(
                temporal_sanity_check_vectorized(
                    trsl_da, time_dim,
                    args.std_short_window,
                    args.std_short_threshold,
                    args.std_short_fraction,
                ).values,
                dtype=bool,
            ).reshape(nsub * nb)
            sane_available = np.asarray(
                availability_check_vectorized(
                    trsl_da, time_dim, args.missing_fraction_threshold
                ).values,
                dtype=bool,
            ).reshape(nsub * nb)
            sane_detection = sane_detection_all[:, b0:b1].reshape(-1)
            sane = sane_long & sane_short & sane_available & sane_detection

            graf_wet_raw = _read_signal_block(
                ds_cml[graf_wet_name], time_dim, sublink_dim, cml_dim,
                slice(common_read_start, common_read_end), src, dtype=np.float32,
            )
            graf_valid_raw = _read_signal_block(
                ds_cml[graf_valid_name], time_dim, sublink_dim, cml_dim,
                slice(common_read_start, common_read_end), src, dtype=np.float32,
            )
            graf_valid = np.isfinite(graf_valid_raw) & (graf_valid_raw == 1) & np.isfinite(graf_wet_raw) & ((graf_wet_raw == 0) | (graf_wet_raw == 1))
            graf_wet = graf_wet_raw == 1

            trsl5, cml5_valid, valid_count5 = aggregate_complete_mean(
                trsl,
                interval_common,
                args.min_cml_valid_fraction,
            )
            radar_proc = _radar_block(
                ds_radar[radar_var_name],
                radar_time_dim,
                radar_cml_dim,
                interval_map.radar_proc_indices,
                radpos,
                args.radar_scale,
            )
            radar_valid = np.isfinite(radar_proc)
            radar_wet = radar_valid & (radar_proc >= args.radar_wet_threshold)

            Rg = np.full_like(trsl, np.nan, dtype=np.float32)
            Rr = np.full_like(trsl5, np.nan, dtype=np.float32)

            def task(si: int):
                sl = si // nb
                lc = si % nb
                g = Estimate(np.full(trsl.shape[1], np.nan, dtype=np.float32), "series failed sanity")
                r = Estimate(np.full(trsl5.shape[1], np.nan, dtype=np.float32), "series failed sanity")
                if sane[si]:
                    g = estimate_pastorek_one(
                        pycml, trsl[si], graf_wet[si], graf_valid[si], common_times.values,
                        length_km[sl, b0+lc], frequency_hz[sl, b0+lc], polarization[sl, b0+lc],
                        args.kr_model, args.n_average_last_dry,
                        args.pastorek_a_max, args.pastorek_zeta, args.pastorek_d,
                    )
                    rv = cml5_valid[si] & radar_valid[:, lc]
                    r = estimate_pastorek_one(
                        pycml, trsl5[si], radar_wet[:, lc], rv, interval_map.radar_proc_times.values,
                        length_km[sl, b0+lc], frequency_hz[sl, b0+lc], polarization[sl, b0+lc],
                        args.kr_model, args.n_average_last_dry,
                        args.pastorek_a_max, args.pastorek_zeta, args.pastorek_d,
                    )
                return si, sl, lc, g, r

            indices = list(range(nsub * nb))
            if args.n_workers <= 1:
                results = [task(i) for i in indices]
            else:
                results = []
                with ThreadPoolExecutor(max_workers=int(args.n_workers)) as ex:
                    futs = [ex.submit(task, i) for i in indices]
                    for fut in as_completed(futs):
                        results.append(fut.result())

            for si, sl, lc, ge, re in results:
                Rg[si] = ge.R
                Rr[si] = re.R
                log_rows.append({
                    "cml_id": str(selected_ids[b0+lc]),
                    "source_cml_index": int(src[lc]),
                    "source_radar_index": int(radpos[lc]),
                    "sublink": int(sl),
                    "sane": int(sane[si]),
                    "graf_error": ge.error,
                    "radar_error": re.error,
                })

            # Write GRAF output slice.
            gos = graf_output_slice
            Rg_out = Rg[:, gos].reshape(nsub, nb, -1).transpose(2, 0, 1)
            gw_out = graf_wet[:, gos].reshape(nsub, nb, -1).transpose(2, 0, 1)
            gv_out = graf_valid[:, gos].reshape(nsub, nb, -1).transpose(2, 0, 1)
            pl_out = plateau[:, gos].reshape(nsub, nb, -1).transpose(2, 0, 1)
            bo_out = blackout[:, gos].reshape(nsub, nb, -1).transpose(2, 0, 1)
            graf_handles.R[:, :, b0:b1] = Rg_out
            graf_handles.R_mean[:, b0:b1] = _safe_nanmean_sublinks(Rg_out)
            graf_handles.wet[:, :, b0:b1] = np.where(gv_out, gw_out.astype(np.uint8), np.uint8(255))
            graf_handles.wet_valid[:, :, b0:b1] = gv_out.astype(np.uint8)
            graf_handles.plateau[:, :, b0:b1] = pl_out.astype(np.uint8)
            graf_handles.blackout[:, :, b0:b1] = bo_out.astype(np.uint8)

            # Write radar output slice.
            ros = interval_map.radar_output_slice
            Rr_out = Rr[:, ros].reshape(nsub, nb, -1).transpose(2, 0, 1)
            tr5_out = trsl5[:, ros].reshape(nsub, nb, -1).transpose(2, 0, 1)
            iv_out = cml5_valid[:, ros].reshape(nsub, nb, -1).transpose(2, 0, 1)
            vc_out = valid_count5[:, ros].reshape(nsub, nb, -1).transpose(2, 0, 1)
            rr_out = radar_proc[ros]
            rv_out = radar_valid[ros]
            rw_out = radar_wet[ros]
            radar_handles.R[:, :, b0:b1] = Rr_out
            radar_handles.R_mean[:, b0:b1] = _safe_nanmean_sublinks(Rr_out)
            radar_handles.radar_rain[:, b0:b1] = rr_out
            radar_handles.radar_valid[:, b0:b1] = rv_out.astype(np.uint8)
            radar_handles.radar_wet[:, b0:b1] = np.where(rv_out, rw_out.astype(np.uint8), np.uint8(255))
            radar_handles.trsl_interval[:, :, b0:b1] = tr5_out
            radar_handles.interval_valid[:, :, b0:b1] = iv_out.astype(np.uint8)
            radar_handles.valid_count[:, :, b0:b1] = vc_out

            # Shared static flags.
            for h in (graf_handles, radar_handles):
                h.sane[:, b0:b1] = sane.reshape(nsub, nb).astype(np.uint8)
                h.sane_long[:, b0:b1] = sane_long.reshape(nsub, nb).astype(np.uint8)
                h.sane_short[:, b0:b1] = sane_short.reshape(nsub, nb).astype(np.uint8)
                h.sane_available[:, b0:b1] = sane_available.reshape(nsub, nb).astype(np.uint8)
                h.sane_detection[:, b0:b1] = sane_detection.reshape(nsub, nb).astype(np.uint8)
                h.detection_limit[:, b0:b1] = dlimit[:, b0:b1].astype(np.float32)

            if ((b0 // args.link_batch) + 1) % args.sync_every_batches == 0:
                graf_handles.nc.sync()
                radar_handles.nc.sync()
            pbar.set_postfix_str(
                f"links {b0}:{b1} | sane {int(sane.sum())}/{len(sane)} | {time.perf_counter()-batch_started:.1f}s"
            )

        graf_handles.nc.setncattr("elapsed_seconds", float(time.perf_counter() - started))
        radar_handles.nc.setncattr("elapsed_seconds", float(time.perf_counter() - started))
        graf_handles.nc.sync()
        radar_handles.nc.sync()
        graf_handles.nc.close()
        radar_handles.nc.close()
        graf_handles = None
        radar_handles = None

        Path(args.log_csv).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(log_rows).to_csv(args.log_csv, index=False)
        graf_partial.replace(graf_path)
        radar_partial.replace(radar_out_path)

        summary = {
            "script": Path(__file__).name,
            "version": SCRIPT_VERSION,
            "input_cml": str(input_path),
            "input_radar": str(radar_path),
            "output_graf": str(graf_path),
            "output_radar": str(radar_out_path),
            "n_links": int(nlink),
            "n_sublinks": int(nsub),
            "graf_time_samples": int(graf_output_end_idx - graf_output_start_idx),
            "radar_time_samples": int(len(radar_output_times)),
            "waa_model": "Pastorek_2021",
            "kr_model": args.kr_model,
            "elapsed_seconds": float(time.perf_counter() - started),
        }
        Path(args.summary_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.summary_json).write_text(json.dumps(summary, indent=2), encoding="utf-8")

        print(f"Saved GRAF product     : {graf_path}")
        print(f"Saved radar product    : {radar_out_path}")
        print(f"Saved processing log   : {args.log_csv}")
        print(f"Saved run summary      : {args.summary_json}")
        print(f"Total elapsed          : {time.perf_counter()-started:.1f} s")
    except Exception:
        if graf_handles is not None:
            try:
                graf_handles.nc.close()
            except Exception:
                pass
        if radar_handles is not None:
            try:
                radar_handles.nc.close()
            except Exception:
                pass
        for p in (graf_partial, radar_partial):
            if p.exists():
                p.unlink()
        raise
    finally:
        ds_cml.close()
        ds_radar.close()


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Streamed Pastorek rainfall retrieval with GRAF and radar wet/dry products",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--input-nc", required=True)
    ap.add_argument("--radar-nc", required=True)
    ap.add_argument("--output-graf-nc", required=True)
    ap.add_argument("--output-radar-nc", required=True)
    ap.add_argument("--log-csv", required=True)
    ap.add_argument("--summary-json", required=True)
    ap.add_argument("--overwrite", action="store_true")

    ap.add_argument("--time-dim", default="time")
    ap.add_argument("--cml-dim", default="cml_id")
    ap.add_argument("--sublink-dim", default="sublink_id")
    ap.add_argument("--radar-time-dim", default="time")
    ap.add_argument("--radar-cml-dim", default="cml_id")
    ap.add_argument("--tsl-var", default="tx")
    ap.add_argument("--rsl-var", default="rx")
    ap.add_argument("--graf-wet-var", default="graf_rstd_wet_1min")
    ap.add_argument("--graf-valid-var", default="graf_rstd_valid_1min")
    ap.add_argument("--frequency-var", default="frequency")
    ap.add_argument("--length-var", default="length")
    ap.add_argument("--polarization-var", default="polarization")
    ap.add_argument("--radar-var", default="rr")

    ap.add_argument("--tstart", default=None)
    ap.add_argument("--tend", default=None)
    ap.add_argument("--baseline-pad-hours", type=float, default=48.0)
    ap.add_argument("--radar-scale", type=float, default=12.0)
    ap.add_argument("--radar-wet-threshold", type=float, default=0.1)
    ap.add_argument("--radar-timestamp-label", choices=("start", "end", "center"), default="start")
    ap.add_argument("--radar-fallback-minutes", type=float, default=5.0)
    ap.add_argument("--min-cml-valid-fraction", type=float, default=1.0)

    ap.add_argument("--frequency-unit", choices=("Hz", "MHz", "GHz", "auto"), default="auto")
    ap.add_argument("--length-unit", choices=("m", "km", "auto"), default="auto")
    ap.add_argument("--kr-model", choices=ALLOWED_KR_MODELS, default="ITU_2005")
    ap.add_argument("--n-average-last-dry", type=int, default=1)
    ap.add_argument("--pastorek-a-max", type=float, default=14.0)
    ap.add_argument("--pastorek-zeta", type=float, default=0.55)
    ap.add_argument("--pastorek-d", type=float, default=0.1)

    # Full preprocessing parameters retained from v5.
    ap.add_argument("--a-quantization-db", type=float, default=0.33)
    ap.add_argument("--detection-limit-threshold", type=float, default=2.0)
    ap.add_argument("--plateau-time-span", type=int, default=3)
    ap.add_argument("--plateau-max-threshold", type=float, default=-85.0)
    ap.add_argument("--plateau-std-threshold", type=float, default=0.5)
    ap.add_argument("--plateau-pad", type=int, default=5)
    ap.add_argument("--blackout-rsl-threshold", type=float, default=-65.0)
    ap.add_argument("--blackout-max-gap", type=int, default=60)
    ap.add_argument("--interpolation-max-gap", default="5min")
    ap.add_argument("--std-long-window", type=int, default=300)
    ap.add_argument("--std-long-threshold", type=float, default=2.0)
    ap.add_argument("--std-long-fraction", type=float, default=0.1)
    ap.add_argument("--std-short-window", type=int, default=60)
    ap.add_argument("--std-short-threshold", type=float, default=0.8)
    ap.add_argument("--std-short-fraction", type=float, default=0.33)
    ap.add_argument("--missing-fraction-threshold", type=float, default=0.5)

    workers = max(1, min(8, os.cpu_count() or 1))
    ap.add_argument("--n-workers", type=int, default=workers)
    ap.add_argument("--link-batch", type=int, default=16)
    ap.add_argument("--availability-scan-batch", type=int, default=128)
    ap.add_argument("--sync-every-batches", type=int, default=4)
    ap.add_argument("--compression-level", type=int, default=1)
    ap.add_argument("--time-chunk-1min", type=int, default=1440)
    ap.add_argument("--time-chunk-radar", type=int, default=288)
    ap.add_argument("--output-link-chunk", type=int, default=8)

    args = ap.parse_args(argv)
    if args.n_workers < 1 or args.link_batch < 1 or args.availability_scan_batch < 1:
        ap.error("worker and batch sizes must be >= 1")
    if not 0 <= args.compression_level <= 9:
        ap.error("--compression-level must be between 0 and 9")
    if not 0 < args.min_cml_valid_fraction <= 1:
        ap.error("--min-cml-valid-fraction must be in (0,1]")
    return args


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
