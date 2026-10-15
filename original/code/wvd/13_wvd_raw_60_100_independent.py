#!/usr/bin/env python3
"""
13_wvd_raw_60_100_independent.py
================================

Build a 5-min CML water-vapour-density (WVD) product directly from the
publication RAW Czech CML file.

Scientific design
-----------------
* Input CML signal source:
    CZ_CML_RAW_1MIN_20200101_20211231.nc
  using TX and RX directly, with TL = TX - RX.

* Population:
    60 <= frequency <= 100 GHz
    2 <= path length <= 10 km
  Every directional sublink is processed independently.

* No external humidity is used anywhere in baseline estimation or retrieval.

* Two independent wet/dry screening branches are produced:
    1) "graf":
         centred 60-min rolling standard deviation of TL
         Q80 threshold
         threshold = 1.12 * Q80
    2) "gauge":
         nearest CHMI precipitation gauge to the CML path
         wet when the native 10-min accumulation is > 0 mm

* Following the Fencl-style dry-period protection discussed for this study,
  detected wet intervals are expanded by +/- 6 h before dry samples are used.

* Evaluation/retrieval period:
    2021-04-01 00:00 through 2021-10-31 23:59 UTC

* Calibration/reference candidate period:
    all available observations outside the evaluation period.
  For the GRAF branch this may include 2020.
  For the GAUGE branch only times covered by the gauge file are classifiable.

* Temperature:
    measured time-varying T(t) from the nearest CHMI temperature station.
    Native 10-min values are interpolated to the CML minute grid only across
    short station gaps, then averaged to the retrieval 5-min interval.
    Temperature is a physical P.676 input, not humidity calibration.

* Pressure:
    the only constant thermodynamic input (default 980 hPa), saved explicitly
    in metadata.

* Atmospheric model:
    ITU-R P.676 exact line-by-line specific gaseous attenuation through
    the `itur` Python package.

Important absolute-reference assumption
---------------------------------------
Without an external humidity reference or a complete telecom link budget,
the constant radio-link offset and absolute atmospheric attenuation cannot be
uniquely separated.

This script therefore implements an explicit INTERNAL anchor:

    B = low_quantile[ TL - L * gamma_P676(rho=0) ]

over safe-dry calibration samples.

The default low quantile is 0.01.  It means that the lower envelope of the
safe-dry calibration signal is anchored to zero water-vapour density.

This is NOT external humidity calibration, but it is a physical/reference
assumption and must be reported as such.  The output names the retrieval
`wvd_internal_anchor` for this reason.

The script also stores pre/post reference estimates so hardware/reference
stability can be assessed without using April-October humidity validation.

Output
------
A NetCDF with compact `series` dimension (one selected directional sublink
per series), `method = [graf, gauge]`, and April-October 2021 5-min time.

The script is streaming by series and does not load the full 3-D CML cube
into memory.

Dependencies
------------
numpy, pandas, xarray, netCDF4, itur

Example
-------
python 13_wvd_raw_60_100_independent.py \
  --cml-nc /home/jovyan/ESSD_publication_dataset/CZ_CML_RAW_1MIN_20200101_20211231.nc \
  --gauge-nc /home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_PRECIPITATION_2021.nc \
  --temperature-nc /home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_TEMPERATURE_2021.nc \
  --out-nc /home/jovyan/ESSD_publication_dataset/CZ_CML_WVD_60_100GHz_5MIN_APR_OCT_2021.nc \
  --diagnostics-csv /home/jovyan/ESSD_publication_dataset/CZ_CML_WVD_60_100GHz_DIAGNOSTICS.csv
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd
import xarray as xr

try:
    from netCDF4 import Dataset, date2num
except Exception as exc:
    raise ImportError("netCDF4 is required") from exc

try:
    from itur.models import itu676
except Exception as exc:
    raise ImportError(
        "itur is required. Install in your environment, e.g. `pip install itur==0.4.0`."
    ) from exc


SCRIPT_VERSION = "1.0"

METHODS = ("graf", "gauge")

EVAL_START = pd.Timestamp("2021-04-01 00:00:00")
EVAL_END_EXCLUSIVE = pd.Timestamp("2021-11-01 00:00:00")

TEMP_CAL_1_START = pd.Timestamp("2021-01-01 00:00:00")
TEMP_CAL_1_END = pd.Timestamp("2021-04-01 00:00:00")
TEMP_CAL_2_START = pd.Timestamp("2021-11-01 00:00:00")
TEMP_CAL_2_END = pd.Timestamp("2022-01-01 00:00:00")

# Publication/raw-signal sentinel rules already applied in the released raw
# file are re-applied defensively here. No plateau filtering, blackout filling,
# gap interpolation, rainfall detection-limit screening, or final rainfall QC.
TX_MIN_VALID_EXCLUSIVE = -99.0
RX_MIN_VALID_EXCLUSIVE = -99.9
RX_MAX_VALID_INCLUSIVE = -20.0
TX_FILL_CODE = 255.0


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

def require_file(path: str, label: str) -> Path:
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"{label} not found: {p}")
    return p


def decode_strings(values) -> np.ndarray:
    a = np.asarray(values)
    if a.dtype.kind == "S":
        out = []
        for x in a.reshape(-1):
            out.append(
                x.decode("utf-8", errors="ignore").replace("\x00", "").strip()
                if isinstance(x, (bytes, np.bytes_))
                else str(x).replace("\x00", "").strip()
            )
        return np.asarray(out, dtype=object)
    return np.asarray([str(x).replace("\x00", "").strip() for x in a.reshape(-1)], dtype=object)


def detect_name(ds: xr.Dataset, candidates: Iterable[str], label: str) -> str:
    for name in candidates:
        if name in ds.variables or name in ds.coords:
            return name
    raise KeyError(
        f"Could not find {label}. Tried {list(candidates)}. "
        f"Available variables: {list(ds.variables)}"
    )


def to_ghz(f) -> np.ndarray:
    x = np.asarray(f, dtype=np.float64)
    med = np.nanmedian(x)
    if med > 1e6:
        return x / 1e9
    if med > 1e3:
        return x / 1e3
    return x


def to_km(length) -> np.ndarray:
    x = np.asarray(length, dtype=np.float64)
    med = np.nanmedian(x)
    return x / 1000.0 if med > 100 else x


def point_to_segment_distances_km(
    station_lat: np.ndarray,
    station_lon: np.ndarray,
    lat_a: float,
    lon_a: float,
    lat_b: float,
    lon_b: float,
) -> np.ndarray:
    """
    Approximate point-to-geodesic-segment distance with a local equirectangular
    projection. Accurate enough for station assignment over Czech CML paths.
    """
    station_lat = np.asarray(station_lat, dtype=float)
    station_lon = np.asarray(station_lon, dtype=float)

    lat0 = math.radians((lat_a + lat_b) / 2.0)
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * math.cos(lat0)

    ax = lon_a * km_per_deg_lon
    ay = lat_a * km_per_deg_lat
    bx = lon_b * km_per_deg_lon
    by = lat_b * km_per_deg_lat

    px = station_lon * km_per_deg_lon
    py = station_lat * km_per_deg_lat

    vx = bx - ax
    vy = by - ay
    vv = vx * vx + vy * vy

    if not np.isfinite(vv) or vv <= 0:
        return np.sqrt((px - ax) ** 2 + (py - ay) ** 2)

    t = ((px - ax) * vx + (py - ay) * vy) / vv
    t = np.clip(t, 0.0, 1.0)
    qx = ax + t * vx
    qy = ay + t * vy
    return np.sqrt((px - qx) ** 2 + (py - qy) ** 2)


def nearest_station_to_path(
    station_lat: np.ndarray,
    station_lon: np.ndarray,
    lat_a: float,
    lon_a: float,
    lat_b: float,
    lon_b: float,
) -> Tuple[int, float]:
    d = point_to_segment_distances_km(
        station_lat, station_lon, lat_a, lon_a, lat_b, lon_b
    )
    finite = np.isfinite(d)
    if not finite.any():
        return -1, np.nan
    j = int(np.nanargmin(d))
    return j, float(d[j])


def rolling_std_centered(
    x: np.ndarray,
    window: int,
    min_periods: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Centered rolling sample standard deviation using cumulative sums.

    Window convention matches a 60-sample centered window:
        left = window//2
        right = window-left
        [i-left, i+right)
    """
    x = np.asarray(x, dtype=np.float64)
    n = x.size

    finite = np.isfinite(x)
    x0 = np.where(finite, x, 0.0)

    cs = np.concatenate(([0.0], np.cumsum(x0, dtype=np.float64)))
    cs2 = np.concatenate(([0.0], np.cumsum(x0 * x0, dtype=np.float64)))
    cc = np.concatenate(([0], np.cumsum(finite.astype(np.int64))))

    i = np.arange(n, dtype=np.int64)
    left_n = window // 2
    right_n = window - left_n
    lo = np.maximum(0, i - left_n)
    hi = np.minimum(n, i + right_n)

    count = cc[hi] - cc[lo]
    sx = cs[hi] - cs[lo]
    sx2 = cs2[hi] - cs2[lo]

    out = np.full(n, np.nan, dtype=np.float64)
    ok = count >= max(2, min_periods)

    c = count[ok].astype(np.float64)
    numer = sx2[ok] - (sx[ok] * sx[ok] / c)
    numer = np.maximum(numer, 0.0)
    out[ok] = np.sqrt(numer / np.maximum(c - 1.0, 1.0))
    return out, count


def dilate_boolean_centered(mask: np.ndarray, pad_samples: int) -> np.ndarray:
    """True if any True occurs within +/- pad_samples, inclusive."""
    m = np.asarray(mask, dtype=bool)
    n = m.size
    if pad_samples <= 0 or n == 0:
        return m.copy()

    c = np.concatenate(([0], np.cumsum(m.astype(np.int64))))
    i = np.arange(n, dtype=np.int64)
    lo = np.maximum(0, i - pad_samples)
    hi = np.minimum(n, i + pad_samples + 1)
    return (c[hi] - c[lo]) > 0


def graf_masks(
    tl: np.ndarray,
    window: int,
    quantile: float,
    factor: float,
    threshold_min_periods: int,
    detection_min_periods: int,
    wet_pad_minutes: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    """
    Return:
      wet_raw, detector_known, wet_buffered, q80, threshold
    """
    rstd_thresh, count_thresh = rolling_std_centered(
        tl, window=window, min_periods=threshold_min_periods
    )
    complete = count_thresh >= threshold_min_periods
    vals = rstd_thresh[complete & np.isfinite(rstd_thresh)]
    if vals.size == 0:
        nanb = np.zeros(tl.size, dtype=bool)
        return nanb, nanb, nanb, np.nan, np.nan

    q = float(np.nanquantile(vals, quantile))
    threshold = factor * q

    rstd_det, count_det = rolling_std_centered(
        tl, window=window, min_periods=detection_min_periods
    )
    known = (count_det >= detection_min_periods) & np.isfinite(rstd_det)
    wet = known & (rstd_det > threshold)
    wet_buffered = dilate_boolean_centered(wet, wet_pad_minutes)
    return wet, known, wet_buffered, q, threshold


def gauge_masks_on_cml_minutes(
    cml_times: pd.DatetimeIndex,
    gauge_times: pd.DatetimeIndex,
    gauge_amount_mm: np.ndarray,
    wet_pad_minutes: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    CHMI value stamped t represents [t-10 min, t).
    Convert to start-labelled minute intervals by assigning the 10-min amount
    to t-10, ..., t-1 minutes.

    Missing gauge values remain unknown; they are never interpreted as dry.
    """
    amount = np.asarray(gauge_amount_mm, dtype=np.float64)

    # Each amount belongs to interval starting 10 min before its timestamp.
    start_idx = gauge_times - pd.Timedelta(minutes=10)
    s = pd.Series(amount, index=start_idx)

    # Propagate a value only across the nine following one-minute starts.
    aligned = s.reindex(
        cml_times,
        method="ffill",
        tolerance=pd.Timedelta(minutes=9),
    ).to_numpy(dtype=np.float64)

    known = np.isfinite(aligned)
    wet = known & (aligned > 0.0)
    wet_buffered = dilate_boolean_centered(wet, wet_pad_minutes)
    return wet, known, wet_buffered


def aggregate_5min_complete(
    times: pd.DatetimeIndex,
    values: np.ndarray,
    safe: np.ndarray,
    eval_start: pd.Timestamp,
    eval_end_exclusive: pd.Timestamp,
    min_valid_minutes: int = 5,
) -> Tuple[pd.DatetimeIndex, np.ndarray, np.ndarray, np.ndarray]:
    """
    Aggregate 1-min values to exact left-labelled 5-min bins.

    A 5-min TL value is emitted only if at least `min_valid_minutes` samples
    are finite AND every minute in the bin is classified safe.
    """
    in_eval = (times >= eval_start) & (times < eval_end_exclusive)
    t = times[in_eval]
    x = np.asarray(values[in_eval], dtype=np.float64)
    s = np.asarray(safe[in_eval], dtype=bool)

    target = pd.date_range(
        eval_start,
        eval_end_exclusive - pd.Timedelta(minutes=5),
        freq="5min",
    )

    if t.size == 0:
        return (
            target,
            np.full(target.size, np.nan, dtype=float),
            np.zeros(target.size, dtype=bool),
            np.zeros(target.size, dtype=np.uint8),
        )

    # The raw publication grid is nominally regular 1 min. Reindex explicitly
    # so acquisition gaps remain missing rather than shifting bin membership.
    sx = pd.Series(x, index=t).reindex(
        pd.date_range(eval_start, eval_end_exclusive - pd.Timedelta(minutes=1), freq="1min")
    )
    ss = pd.Series(s.astype(np.uint8), index=t).reindex(sx.index)

    arr = sx.to_numpy(dtype=np.float64)
    safe_arr = ss.to_numpy(dtype=np.float64)

    nbin = target.size
    expected = nbin * 5
    if arr.size != expected:
        raise RuntimeError(
            f"Unexpected evaluation-grid length: {arr.size} != {expected}"
        )

    arr2 = arr.reshape(nbin, 5)
    safe2 = safe_arr.reshape(nbin, 5)

    finite_count = np.isfinite(arr2).sum(axis=1)
    safe_all = np.isfinite(safe2).all(axis=1) & (safe2 == 1).all(axis=1)
    valid = (finite_count >= min_valid_minutes) & safe_all

    mean = np.full(nbin, np.nan, dtype=np.float64)
    candidate = np.full(nbin, np.nan, dtype=np.float64)
    has_data = finite_count > 0
    candidate[has_data] = (
        np.nansum(arr2[has_data], axis=1) / finite_count[has_data]
    )
    mean[valid] = candidate[valid]

    # 0 unknown/unsafe, 1 complete safe-dry bin
    safe5 = valid.astype(np.uint8)
    return target, mean, valid, safe5


def as_float_quantity_value(q) -> np.ndarray:
    if hasattr(q, "value"):
        return np.asarray(q.value, dtype=np.float64)
    return np.asarray(q, dtype=np.float64)



_P676_TEMP_CACHE: Dict[Tuple[float, float, float, float, float, float, float], Tuple[np.ndarray, np.ndarray, np.ndarray, bool, np.ndarray]] = {}


def p676_temperature_lookup(
    frequency_ghz: float,
    pressure_hpa: float,
    temp_min_c: float,
    temp_max_c: float,
    temp_step_c: float,
    rho_max: float,
    rho_step: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, bool, np.ndarray]:
    """
    Build/cache gamma(T, rho) [dB/km] at fixed frequency and pressure.

    Temperature is evaluated on a regular grid. Retrieval uses interpolation
    between neighbouring temperature-grid rows, so measured station T(t)
    remains time varying while pressure is the only constant thermodynamic
    input.
    """
    key = (
        round(float(frequency_ghz), 6),
        round(float(pressure_hpa), 6),
        round(float(temp_min_c), 6),
        round(float(temp_max_c), 6),
        round(float(temp_step_c), 6),
        round(float(rho_max), 6),
        round(float(rho_step), 6),
    )
    cached = _P676_TEMP_CACHE.get(key)
    if cached is not None:
        return cached

    if temp_step_c <= 0 or rho_step <= 0:
        raise ValueError("Temperature/rho lookup steps must be > 0")
    if temp_max_c <= temp_min_c:
        raise ValueError("Temperature lookup max must exceed min")

    temp_grid = np.arange(
        temp_min_c,
        temp_max_c + 0.5 * temp_step_c,
        temp_step_c,
        dtype=np.float64,
    )
    rho_grid = np.arange(
        0.0,
        rho_max + 0.5 * rho_step,
        rho_step,
        dtype=np.float64,
    )

    gamma = np.empty((temp_grid.size, rho_grid.size), dtype=np.float32)
    sensitivity = np.full(temp_grid.size, np.nan, dtype=np.float64)

    rho_ref = 10.0
    kref = int(np.searchsorted(rho_grid, rho_ref))
    kref = max(1, min(kref, rho_grid.size - 2))

    monotonic = True
    for i, tc in enumerate(temp_grid):
        g = as_float_quantity_value(
            itu676.gamma_exact(
                float(frequency_ghz),
                float(pressure_hpa),
                rho_grid,
                float(tc) + 273.15,
            )
        ).reshape(-1)

        if (
            g.size != rho_grid.size
            or not np.isfinite(g).all()
            or np.any(np.diff(g) < -1e-10)
        ):
            monotonic = False

        gamma[i, :] = g.astype(np.float32)
        sensitivity[i] = float(
            (g[kref + 1] - g[kref - 1])
            / (rho_grid[kref + 1] - rho_grid[kref - 1])
        )

    result = (temp_grid, rho_grid, gamma, monotonic, sensitivity)
    _P676_TEMP_CACHE[key] = result
    return result


def gamma_rho0_for_temperature(
    temperature_c: np.ndarray,
    temp_grid: np.ndarray,
    gamma_table: np.ndarray,
) -> np.ndarray:
    """Interpolate P.676 gamma(rho=0) to measured T(t)."""
    t = np.asarray(temperature_c, dtype=np.float64)
    out = np.full(t.shape, np.nan, dtype=np.float64)
    ok = (
        np.isfinite(t)
        & (t >= temp_grid[0])
        & (t <= temp_grid[-1])
    )
    if ok.any():
        out[ok] = np.interp(
            t[ok],
            temp_grid,
            np.asarray(gamma_table[:, 0], dtype=np.float64),
        )
    return out


def invert_gamma_lookup_temperature(
    target_gamma: np.ndarray,
    temperature_c: np.ndarray,
    temp_grid: np.ndarray,
    rho_grid: np.ndarray,
    gamma_table: np.ndarray,
    monotonic: bool,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Invert gamma(T,rho) -> rho using measured time-varying temperature.

    The inverse solution is evaluated in the two neighbouring temperature
    lookup rows and linearly interpolated in temperature.

    status:
      0 unavailable / missing T or attenuation
      1 valid
      2 below P.676 rho lookup range
      3 above P.676 rho lookup range
      4 non-monotonic P.676 lookup
      5 temperature outside P.676 lookup grid
    """
    target = np.asarray(target_gamma, dtype=np.float64)
    temp = np.asarray(temperature_c, dtype=np.float64)

    if target.shape != temp.shape:
        raise ValueError("target_gamma and temperature_c must have same shape")

    out = np.full(target.shape, np.nan, dtype=np.float64)
    status = np.zeros(target.shape, dtype=np.uint8)

    finite = np.isfinite(target) & np.isfinite(temp)
    if not monotonic:
        status[finite] = 4
        return out, status

    outside_t = finite & (
        (temp < temp_grid[0]) | (temp > temp_grid[-1])
    )
    status[outside_t] = 5

    usable = finite & ~outside_t
    if not usable.any():
        return out, status

    # Lower T-grid row and interpolation weight.
    j = np.searchsorted(temp_grid, temp[usable], side="right") - 1
    j = np.clip(j, 0, len(temp_grid) - 2)
    t0 = temp_grid[j]
    t1 = temp_grid[j + 1]
    w = (temp[usable] - t0) / (t1 - t0)

    idx_usable = np.flatnonzero(usable)
    gmin0 = np.asarray(gamma_table[j, 0], dtype=np.float64)
    gmin1 = np.asarray(gamma_table[j + 1, 0], dtype=np.float64)
    gmax0 = np.asarray(gamma_table[j, -1], dtype=np.float64)
    gmax1 = np.asarray(gamma_table[j + 1, -1], dtype=np.float64)

    gmin = (1.0 - w) * gmin0 + w * gmin1
    gmax = (1.0 - w) * gmax0 + w * gmax1
    tgt = target[usable]

    below_local = tgt < gmin
    above_local = tgt > gmax
    inside_local = ~(below_local | above_local)

    status[idx_usable[below_local]] = 2
    status[idx_usable[above_local]] = 3

    if inside_local.any():
        inside_global = idx_usable[inside_local]
        j_inside = j[inside_local]
        w_inside = w[inside_local]
        tgt_inside = tgt[inside_local]

        # Group by lower temperature row so np.interp can work vectorially.
        for jj in np.unique(j_inside):
            g = j_inside == jj
            tg = tgt_inside[g]
            r0 = np.interp(
                tg,
                np.asarray(gamma_table[jj, :], dtype=np.float64),
                rho_grid,
            )
            r1 = np.interp(
                tg,
                np.asarray(gamma_table[jj + 1, :], dtype=np.float64),
                rho_grid,
            )
            ww = w_inside[g]
            out[inside_global[g]] = (1.0 - ww) * r0 + ww * r1
            status[inside_global[g]] = 1

    return out, status

def calibration_masks(times: pd.DatetimeIndex) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Internal CML reference uses only 2021, because measured CHMI temperature
    is available for 2021.

    Main: all available 2021 safe-dry observations, including Apr-Oct CML
          signal samples. No humidity/WVD observations enter this estimate.
    Pre : Jan-Mar 2021, retained as a reference-stability diagnostic.
    Post: Nov-Dec 2021, retained as a reference-stability diagnostic.
    """
    main = (
        (times >= pd.Timestamp("2021-01-01 00:00:00"))
        & (times < pd.Timestamp("2022-01-01 00:00:00"))
    )
    pre = (
        (times >= pd.Timestamp("2021-01-01 00:00:00"))
        & (times < EVAL_START)
    )
    post = (
        (times >= EVAL_END_EXCLUSIVE)
        & (times < pd.Timestamp("2022-01-01 00:00:00"))
    )
    return np.asarray(main), np.asarray(pre), np.asarray(post)


def robust_quantile(x: np.ndarray, q: float, min_n: int) -> Tuple[float, int]:
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    n = int(x.size)
    if n < min_n:
        return np.nan, n
    return float(np.nanquantile(x, q)), n


def interpolate_temperature_to_minutes(
    native_times: pd.DatetimeIndex,
    native_values_c: np.ndarray,
    target_times: pd.DatetimeIndex,
    max_gap_minutes: float,
) -> np.ndarray:
    """
    Linearly interpolate station temperature to target times.

    Interpolation is only allowed between two finite native observations whose
    separation is <= max_gap_minutes. There is no extrapolation across the
    beginning/end of the record and no bridging of longer station gaps.
    """
    values = np.asarray(native_values_c, dtype=np.float64)
    native_ns = native_times.view("int64")
    target_ns = target_times.view("int64")

    good = np.isfinite(values)
    if good.sum() < 2:
        return np.full(target_times.size, np.nan, dtype=np.float64)

    x = native_ns[good]
    y = values[good]

    pos = np.searchsorted(x, target_ns, side="left")
    out = np.full(target_ns.size, np.nan, dtype=np.float64)

    # Exact native observations.
    exact = (pos < x.size)
    exact_idx = np.flatnonzero(exact)
    exact_match = x[pos[exact]] == target_ns[exact]
    if exact_match.any():
        use = exact_idx[exact_match]
        out[use] = y[pos[use]]

    # Interpolation strictly between valid observations.
    between = (pos > 0) & (pos < x.size) & ~np.isfinite(out)
    if between.any():
        idx = np.flatnonzero(between)
        left = pos[idx] - 1
        right = pos[idx]

        gap_min = (x[right] - x[left]) / (60.0 * 1e9)
        ok_gap = gap_min <= float(max_gap_minutes)

        if ok_gap.any():
            ii = idx[ok_gap]
            ll = left[ok_gap]
            rr = right[ok_gap]
            frac = (
                (target_ns[ii] - x[ll]).astype(np.float64)
                / (x[rr] - x[ll]).astype(np.float64)
            )
            out[ii] = y[ll] + frac * (y[rr] - y[ll])

    return out


def aggregate_temperature_5min(
    times: pd.DatetimeIndex,
    temperature_1min_c: np.ndarray,
    eval_start: pd.Timestamp,
    eval_end_exclusive: pd.Timestamp,
    min_valid_minutes: int,
) -> Tuple[pd.DatetimeIndex, np.ndarray]:
    """Mean measured/interpolated station T over each left-labelled 5-min bin."""
    in_eval = (times >= eval_start) & (times < eval_end_exclusive)
    t = times[in_eval]
    x = np.asarray(temperature_1min_c[in_eval], dtype=np.float64)

    target = pd.date_range(
        eval_start,
        eval_end_exclusive - pd.Timedelta(minutes=5),
        freq="5min",
    )
    minute_grid = pd.date_range(
        eval_start,
        eval_end_exclusive - pd.Timedelta(minutes=1),
        freq="1min",
    )
    sx = pd.Series(x, index=t).reindex(minute_grid)
    arr = sx.to_numpy(dtype=np.float64).reshape(target.size, 5)

    count = np.isfinite(arr).sum(axis=1)
    out = np.full(target.size, np.nan, dtype=np.float64)
    ok = count >= int(min_valid_minutes)
    if ok.any():
        out[ok] = np.nansum(arr[ok], axis=1) / count[ok]
    return target, out

@dataclass
class SeriesMeta:
    cml_index: int
    sublink_index: int
    cml_id: str
    sublink_id: str
    frequency_ghz: float
    length_km: float
    polarization: str
    lat_a: float
    lon_a: float
    lat_b: float
    lon_b: float


# -----------------------------------------------------------------------------
# Output writer
# -----------------------------------------------------------------------------

def create_output(
    path: Path,
    eval_time: pd.DatetimeIndex,
    series: list[SeriesMeta],
    args,
) -> Dataset:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if not args.overwrite:
            raise FileExistsError(f"Output exists: {path} (use --overwrite)")
        path.unlink()

    nc = Dataset(path, "w", format="NETCDF4")
    nc.createDimension("method", len(METHODS))
    nc.createDimension("series", len(series))
    nc.createDimension("time", len(eval_time))

    vmethod = nc.createVariable("method", str, ("method",))
    vmethod[:] = np.asarray(METHODS, dtype=object)

    vtime = nc.createVariable("time", "f8", ("time",))
    vtime.units = "seconds since 1970-01-01 00:00:00 UTC"
    vtime.calendar = "proleptic_gregorian"
    vtime[:] = date2num(
        eval_time.to_pydatetime(),
        units=vtime.units,
        calendar=vtime.calendar,
    )

    vcml = nc.createVariable("cml_id", str, ("series",))
    vsub = nc.createVariable("sublink_id", str, ("series",))
    vfreq = nc.createVariable("frequency_ghz", "f4", ("series",), fill_value=np.nan)
    vlen = nc.createVariable("length_km", "f4", ("series",), fill_value=np.nan)
    vpol = nc.createVariable("polarization", str, ("series",))

    for name in (
        "site_a_latitude", "site_a_longitude",
        "site_b_latitude", "site_b_longitude",
    ):
        v = nc.createVariable(name, "f4", ("series",), fill_value=np.nan)
        v.units = "degrees_north" if "latitude" in name else "degrees_east"

    vcml[:] = np.asarray([s.cml_id for s in series], dtype=object)
    vsub[:] = np.asarray([s.sublink_id for s in series], dtype=object)
    vfreq[:] = np.asarray([s.frequency_ghz for s in series], dtype=np.float32)
    vfreq.units = "GHz"
    vlen[:] = np.asarray([s.length_km for s in series], dtype=np.float32)
    vlen.units = "km"
    vpol[:] = np.asarray([s.polarization for s in series], dtype=object)
    nc["site_a_latitude"][:] = np.asarray([s.lat_a for s in series], dtype=np.float32)
    nc["site_a_longitude"][:] = np.asarray([s.lon_a for s in series], dtype=np.float32)
    nc["site_b_latitude"][:] = np.asarray([s.lat_b for s in series], dtype=np.float32)
    nc["site_b_longitude"][:] = np.asarray([s.lon_b for s in series], dtype=np.float32)

    # Assignments / thermodynamics
    for name, units in (
        ("gauge_distance_km", "km"),
        ("temperature_distance_km", "km"),
        ("temperature_ref_c", "degree_Celsius"),
        ("p676_sensitivity_dbkm_per_gm3", "dB km-1 (g m-3)-1"),
    ):
        v = nc.createVariable(name, "f4", ("series",), fill_value=np.nan)
        v.units = units

    nc.createVariable("gauge_station_id", str, ("series",))
    nc.createVariable("temperature_station_id", str, ("series",))
    nc.createVariable("temperature_ref_n", "i4", ("series",), fill_value=-1)

    vt = nc.createVariable(
        "temperature_5min_c",
        "f4", ("series", "time"),
        fill_value=np.nan, zlib=True, complevel=2, shuffle=True,
        chunksizes=(1, min(2016, len(eval_time))),
    )
    vt.units = "degree_Celsius"
    vt.long_name = (
        "Nearest CHMI station air temperature averaged over the same "
        "left-labelled five-minute retrieval interval"
    )

    nc.createVariable("p676_monotonic", "u1", ("series",), fill_value=np.uint8(255))
    nc.createVariable("graf_q80_db", "f4", ("series",), fill_value=np.nan)
    nc.createVariable("graf_threshold_db", "f4", ("series",), fill_value=np.nan)

    # Calibration/reference variables
    for name in (
        "baseline_B_db",
        "baseline_pre_db",
        "baseline_post_db",
        "baseline_delta_post_minus_pre_db",
    ):
        v = nc.createVariable(
            name, "f4", ("method", "series"),
            fill_value=np.nan, zlib=True, complevel=1, shuffle=True
        )
        v.units = "dB"

    for name in ("baseline_n", "baseline_pre_n", "baseline_post_n"):
        nc.createVariable(
            name, "i8", ("method", "series"),
            fill_value=-1, zlib=True, complevel=1, shuffle=True
        )

    # Main time-series variables.
    chunks = (1, 1, min(2016, len(eval_time)))

    v = nc.createVariable(
        "wvd_internal_anchor",
        "f4", ("method", "series", "time"),
        fill_value=np.nan, zlib=True, complevel=2, shuffle=True,
        chunksizes=chunks,
    )
    v.units = "g m-3"
    v.long_name = (
        "Water vapour density retrieved with ITU-R P.676 and an internal "
        "zero-humidity lower-envelope reference"
    )

    v = nc.createVariable(
        "total_loss_5min",
        "f4", ("method", "series", "time"),
        fill_value=np.nan, zlib=True, complevel=2, shuffle=True,
        chunksizes=chunks,
    )
    v.units = "dB"
    v.long_name = "Five-minute mean TX minus RX for complete safe-dry bins"

    v = nc.createVariable(
        "safe_dry_5min",
        "u1", ("method", "series", "time"),
        fill_value=np.uint8(255), zlib=True, complevel=2, shuffle=True,
        chunksizes=chunks,
    )
    v.flag_values = np.asarray([0, 1], dtype=np.uint8)
    v.flag_meanings = "not_safe_or_incomplete safe_dry"

    v = nc.createVariable(
        "retrieval_status",
        "u1", ("method", "series", "time"),
        fill_value=np.uint8(255), zlib=True, complevel=2, shuffle=True,
        chunksizes=chunks,
    )
    v.flag_values = np.asarray([0, 1, 2, 3, 4, 5], dtype=np.uint8)
    v.flag_meanings = (
        "unavailable valid below_p676_range above_p676_range "
        "nonmonotonic_p676_lookup temperature_outside_lookup"
    )

    # Global provenance.
    nc.title = "Czech 60-100 GHz CML water-vapour retrieval from raw TX/RX"
    nc.script = Path(__file__).name
    nc.script_version = SCRIPT_VERSION
    nc.source_cml = str(Path(args.cml_nc).expanduser())
    nc.source_precipitation = str(Path(args.gauge_nc).expanduser())
    nc.source_temperature = str(Path(args.temperature_nc).expanduser())
    nc.frequency_selection_ghz = f"{args.freq_min_ghz:g} to {args.freq_max_ghz:g}"
    nc.length_selection_km = f"{args.length_min_km:g} to {args.length_max_km:g}"
    nc.evaluation_period = "2021-04-01 through 2021-10-31 UTC"
    nc.calibration_period = (
        "all available 2021 safe-dry CML observations; April-October CML "
        "signal may contribute to the internal radio reference; no RH/WVD "
        "reference used"
    )
    nc.humidity_used_for_calibration = "no"
    nc.humidity_used_for_retrieval = "no"
    nc.pressure_assumption_hpa = float(args.pressure_hpa)
    nc.temperature_reference = (
        "nearest CHMI temperature station; measured time-varying T interpolated "
        "to the CML minute grid without bridging station gaps longer than the "
        "configured maximum; five-minute means used in retrieval"
    )
    nc.wet_dry_methods = "graf,gauge"
    nc.graf_definition = (
        f"centered {args.graf_window_min}-min rolling sample STD; "
        f"Q{int(args.graf_quantile*100)}; threshold={args.graf_factor:g}*Q; "
        f"detection_min_periods={args.graf_detection_min_periods}"
    )
    nc.gauge_definition = (
        "nearest CHMI precipitation station to CML path; wet when native "
        "10-min accumulation > 0 mm; timestamp t represents [t-10min,t)"
    )
    nc.wet_buffer_hours = float(args.wet_buffer_hours)
    nc.baseline_definition = (
        "low quantile of TL - L*gamma_P676(rho=0) over safe-dry calibration samples"
    )
    nc.baseline_quantile = float(args.baseline_quantile)
    nc.absolute_reference_warning = (
        "No external humidity calibration or full telecom link budget is used. "
        "Absolute WVD is internally anchored by assigning the lower safe-dry "
        "calibration envelope to rho=0 g m-3."
    )
    nc.p676_pressure_hpa = float(args.pressure_hpa)
    nc.p676_rho_grid_max_gm3 = float(args.rho_max_gm3)
    nc.p676_rho_grid_step_gm3 = float(args.rho_step_gm3)
    nc.p676_temperature_grid_min_c = float(args.temperature_grid_min_c)
    nc.p676_temperature_grid_max_c = float(args.temperature_grid_max_c)
    nc.p676_temperature_grid_step_c = float(args.temperature_grid_step_c)
    nc.temperature_interpolation_max_gap_minutes = float(args.temperature_max_gap_min)

    return nc


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Raw 60-100 GHz CML WVD retrieval with GRAF and gauge wet/dry branches.",
    )

    p.add_argument(
        "--cml-nc",
        default="/home/jovyan/ESSD_publication_dataset/CZ_CML_RAW_1MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--gauge-nc",
        default="/home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_PRECIPITATION_2021.nc",
    )
    p.add_argument(
        "--temperature-nc",
        default="/home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_TEMPERATURE_2021.nc",
    )
    p.add_argument(
        "--out-nc",
        default="/home/jovyan/ESSD_publication_dataset/CZ_CML_WVD_60_100GHz_5MIN_APR_OCT_2021.nc",
    )
    p.add_argument(
        "--diagnostics-csv",
        default="/home/jovyan/ESSD_publication_dataset/CZ_CML_WVD_60_100GHz_DIAGNOSTICS.csv",
    )

    p.add_argument("--freq-min-ghz", type=float, default=60.0)
    p.add_argument("--freq-max-ghz", type=float, default=100.0)
    p.add_argument("--length-min-km", type=float, default=2.0)
    p.add_argument("--length-max-km", type=float, default=10.0)

    p.add_argument("--pressure-hpa", type=float, default=980.0)

    p.add_argument("--graf-window-min", type=int, default=60)
    p.add_argument("--graf-quantile", type=float, default=0.80)
    p.add_argument("--graf-factor", type=float, default=1.12)
    p.add_argument("--graf-threshold-min-periods", type=int, default=60)
    p.add_argument("--graf-detection-min-periods", type=int, default=45)

    p.add_argument("--wet-buffer-hours", type=float, default=6.0)
    p.add_argument("--max-gauge-distance-km", type=float, default=5.0)
    p.add_argument("--max-temperature-distance-km", type=float, default=20.0)
    p.add_argument(
        "--temperature-max-gap-min",
        type=float,
        default=30.0,
        help="Maximum gap between finite 10-min station T observations that may be interpolated.",
    )
    p.add_argument("--temperature-grid-min-c", type=float, default=-35.0)
    p.add_argument("--temperature-grid-max-c", type=float, default=45.0)
    p.add_argument("--temperature-grid-step-c", type=float, default=0.25)
    p.add_argument(
        "--min-temperature-minutes-5min",
        type=int,
        default=3,
        help="Minimum finite 1-min temperature values needed for a 5-min T mean.",
    )

    p.add_argument(
        "--baseline-quantile",
        type=float,
        default=0.01,
        help="Internal lower-envelope quantile used to anchor rho=0.",
    )
    p.add_argument(
        "--min-baseline-samples",
        type=int,
        default=1440,
        help="Minimum safe-dry 1-min samples for main baseline.",
    )
    p.add_argument(
        "--min-prepost-samples",
        type=int,
        default=300,
        help="Minimum samples for pre/post baseline diagnostics.",
    )

    p.add_argument("--rho-max-gm3", type=float, default=40.0)
    p.add_argument("--rho-step-gm3", type=float, default=0.02)
    p.add_argument("--min-valid-minutes-5min", type=int, default=5)

    p.add_argument(
        "--max-series",
        type=int,
        default=None,
        help="Smoke-test limit. Omit for the full selected population.",
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not (0.0 < args.graf_quantile < 1.0):
        raise ValueError("--graf-quantile must be in (0,1)")
    if not (0.0 <= args.baseline_quantile <= 0.5):
        raise ValueError("--baseline-quantile must be between 0 and 0.5")
    if args.min_valid_minutes_5min < 1 or args.min_valid_minutes_5min > 5:
        raise ValueError("--min-valid-minutes-5min must be 1..5")
    if args.min_temperature_minutes_5min < 1 or args.min_temperature_minutes_5min > 5:
        raise ValueError("--min-temperature-minutes-5min must be 1..5")
    if args.temperature_grid_step_c <= 0:
        raise ValueError("--temperature-grid-step-c must be > 0")
    if args.temperature_grid_max_c <= args.temperature_grid_min_c:
        raise ValueError("temperature grid max must exceed min")

    cml_path = require_file(args.cml_nc, "CML raw NetCDF")
    gauge_path = require_file(args.gauge_nc, "CHMI precipitation NetCDF")
    temp_path = require_file(args.temperature_nc, "CHMI temperature NetCDF")
    out_path = Path(args.out_nc).expanduser().resolve()
    diag_path = Path(args.diagnostics_csv).expanduser().resolve()

    print("=" * 100)
    print("RAW CML WVD — 60-100 GHz — no external humidity calibration")
    print("=" * 100)
    print(f"CML        : {cml_path}")
    print(f"Gauge      : {gauge_path}")
    print(f"Temperature: {temp_path}")
    print(f"Output     : {out_path}")
    print()

    # ------------------------------------------------------------------
    # Open input metadata
    # ------------------------------------------------------------------
    cml = xr.open_dataset(cml_path, cache=False)
    gauge = xr.open_dataset(gauge_path, cache=False)
    temp = xr.open_dataset(temp_path, cache=False)

    required = [
        "tx", "rx", "frequency", "length", "cml_id", "sublink_id",
        "site_a_latitude", "site_a_longitude",
        "site_b_latitude", "site_b_longitude",
    ]
    missing = [v for v in required if v not in cml.variables and v not in cml.coords]
    if missing:
        raise KeyError(f"Missing required raw CML variables: {missing}")

    cml_times = pd.DatetimeIndex(pd.to_datetime(cml["time"].values))
    if not cml_times.is_monotonic_increasing:
        raise ValueError("CML time coordinate must be monotonically increasing")

    cml_ids = decode_strings(cml["cml_id"].values)
    sub_ids = decode_strings(cml["sublink_id"].values)

    freq = to_ghz(cml["frequency"].values)
    length = to_km(cml["length"].values)

    # Canonical frequency shape expected is (sublink, cml).
    if freq.shape == (len(cml_ids), len(sub_ids)):
        freq = freq.T
    if freq.shape != (len(sub_ids), len(cml_ids)):
        raise ValueError(
            f"Unexpected frequency shape {freq.shape}; expected "
            f"({len(sub_ids)}, {len(cml_ids)})"
        )

    lat_a = np.asarray(cml["site_a_latitude"].values, dtype=float).reshape(-1)
    lon_a = np.asarray(cml["site_a_longitude"].values, dtype=float).reshape(-1)
    lat_b = np.asarray(cml["site_b_latitude"].values, dtype=float).reshape(-1)
    lon_b = np.asarray(cml["site_b_longitude"].values, dtype=float).reshape(-1)

    if "polarization" in cml:
        pol = np.asarray(cml["polarization"].values)
        if pol.shape == (len(cml_ids), len(sub_ids)):
            pol = pol.T
        pol_flat = np.empty(pol.shape, dtype=object)
        for j in range(pol.shape[0]):
            pol_flat[j, :] = decode_strings(pol[j, :])
        pol = pol_flat
    else:
        pol = np.full(freq.shape, "", dtype=object)

    selected: list[SeriesMeta] = []
    for ci, cid in enumerate(cml_ids):
        L = float(length[ci])
        if not (np.isfinite(L) and args.length_min_km <= L <= args.length_max_km):
            continue
        for si, sid in enumerate(sub_ids):
            f = float(freq[si, ci])
            if not (
                np.isfinite(f)
                and args.freq_min_ghz <= f <= args.freq_max_ghz
            ):
                continue
            selected.append(
                SeriesMeta(
                    cml_index=ci,
                    sublink_index=si,
                    cml_id=str(cid),
                    sublink_id=str(sid),
                    frequency_ghz=f,
                    length_km=L,
                    polarization=str(pol[si, ci]),
                    lat_a=float(lat_a[ci]),
                    lon_a=float(lon_a[ci]),
                    lat_b=float(lat_b[ci]),
                    lon_b=float(lon_b[ci]),
                )
            )

    if args.max_series is not None:
        selected = selected[: args.max_series]

    if not selected:
        raise RuntimeError("No directional sublinks satisfy the selection")

    n_phys = len(set(s.cml_id for s in selected))
    print(f"Selected directional sublinks: {len(selected):,}")
    print(f"Selected physical CMLs        : {n_phys:,}")
    print(
        f"Selection                     : {args.freq_min_ghz:g}-{args.freq_max_ghz:g} GHz, "
        f"{args.length_min_km:g}-{args.length_max_km:g} km"
    )

    # ------------------------------------------------------------------
    # Gauge metadata
    # ------------------------------------------------------------------
    gauge_var = detect_name(
        gauge, ("SRA10M", "rr"), "precipitation accumulation variable"
    )
    gauge_da = gauge[gauge_var]
    gauge_time_dim = next(d for d in gauge_da.dims if "time" in d.lower())
    gauge_station_dim = next(d for d in gauge_da.dims if d != gauge_time_dim)

    gauge_lat_name = detect_name(gauge, ("lat", "latitude"), "gauge latitude")
    gauge_lon_name = detect_name(gauge, ("lon", "longitude"), "gauge longitude")
    gauge_lats = np.asarray(gauge[gauge_lat_name].values, dtype=float).reshape(-1)
    gauge_lons = np.asarray(gauge[gauge_lon_name].values, dtype=float).reshape(-1)

    if gauge_station_dim in gauge.coords:
        gauge_ids = decode_strings(gauge[gauge_station_dim].values)
    else:
        gauge_id_name = detect_name(
            gauge, ("station", "WSI", "station_id", "GH_ID"), "gauge station ID"
        )
        gauge_ids = decode_strings(gauge[gauge_id_name].values)

    gauge_times = pd.DatetimeIndex(pd.to_datetime(gauge[gauge_time_dim].values))

    # ------------------------------------------------------------------
    # Temperature metadata
    # ------------------------------------------------------------------
    temp_var = detect_name(temp, ("T",), "standard air temperature variable")
    temp_da = temp[temp_var]
    temp_time_dim = next(d for d in temp_da.dims if "time" in d.lower())
    temp_station_dim = next(d for d in temp_da.dims if d != temp_time_dim)

    temp_lat_name = detect_name(temp, ("lat", "latitude"), "temperature latitude")
    temp_lon_name = detect_name(temp, ("lon", "longitude"), "temperature longitude")
    temp_lats = np.asarray(temp[temp_lat_name].values, dtype=float).reshape(-1)
    temp_lons = np.asarray(temp[temp_lon_name].values, dtype=float).reshape(-1)

    temp_times = pd.DatetimeIndex(pd.to_datetime(temp[temp_time_dim].values))

    if temp_station_dim in temp.coords:
        temp_ids = decode_strings(temp[temp_station_dim].values)
    else:
        temp_id_name = detect_name(
            temp, ("station", "WSI", "station_id", "GH_ID"), "temperature station ID"
        )
        temp_ids = decode_strings(temp[temp_id_name].values)

    # ------------------------------------------------------------------
    # Common masks / output time
    # ------------------------------------------------------------------
    cal_main, cal_pre, cal_post = calibration_masks(cml_times)
    eval_time = pd.date_range(
        EVAL_START,
        EVAL_END_EXCLUSIVE - pd.Timedelta(minutes=5),
        freq="5min",
    )

    wet_pad_minutes = int(round(args.wet_buffer_hours * 60.0))

    # Cache gauge masks by station index: gauge station often reused.
    gauge_mask_cache: Dict[int, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    # Only keep the most recently used interpolated temperature station series.
    # Directional sublinks of one physical CML are adjacent, so this avoids
    # duplicating a ~1 million element array without allowing an unbounded cache.
    last_temp_station_index: Optional[int] = None
    last_temp_1min_c: Optional[np.ndarray] = None

    nc = create_output(out_path, eval_time, selected, args)

    diagnostics = []
    t0_run = time.time()

    # ------------------------------------------------------------------
    # Process every directional sublink independently
    # ------------------------------------------------------------------
    for out_i, s in enumerate(selected):
        started = time.time()

        tx = np.asarray(
            cml["tx"].isel(
                sublink_id=s.sublink_index,
                cml_id=s.cml_index,
            ).values,
            dtype=np.float64,
        )
        rx = np.asarray(
            cml["rx"].isel(
                sublink_id=s.sublink_index,
                cml_id=s.cml_index,
            ).values,
            dtype=np.float64,
        )

        valid_signal = (
            np.isfinite(tx)
            & np.isfinite(rx)
            & (tx > TX_MIN_VALID_EXCLUSIVE)
            & (tx != TX_FILL_CODE)
            & (rx > RX_MIN_VALID_EXCLUSIVE)
            & (rx <= RX_MAX_VALID_INCLUSIVE)
        )
        tl = np.full(tx.shape, np.nan, dtype=np.float64)
        tl[valid_signal] = tx[valid_signal] - rx[valid_signal]

        # Station assignment is per physical path.
        gj, gd = nearest_station_to_path(
            gauge_lats, gauge_lons,
            s.lat_a, s.lon_a, s.lat_b, s.lon_b
        )
        tj, td = nearest_station_to_path(
            temp_lats, temp_lons,
            s.lat_a, s.lon_a, s.lat_b, s.lon_b
        )

        gauge_ok = gj >= 0 and np.isfinite(gd) and gd <= args.max_gauge_distance_km
        temp_ok = tj >= 0 and np.isfinite(td) and td <= args.max_temperature_distance_km

        # Measured time-varying station temperature.
        if temp_ok:
            if last_temp_station_index == tj and last_temp_1min_c is not None:
                temp_1min_c = last_temp_1min_c
            else:
                native_temp = np.asarray(
                    temp_da.isel({temp_station_dim: tj}).values,
                    dtype=np.float64,
                )
                temp_1min_c = interpolate_temperature_to_minutes(
                    temp_times,
                    native_temp,
                    cml_times,
                    max_gap_minutes=args.temperature_max_gap_min,
                )
                last_temp_station_index = tj
                last_temp_1min_c = temp_1min_c

            t_ref_n = int(np.isfinite(temp_1min_c).sum())
            t_ref_c = (
                float(np.nanmedian(temp_1min_c))
                if t_ref_n > 0
                else np.nan
            )
        else:
            temp_1min_c = np.full(cml_times.size, np.nan, dtype=np.float64)
            t_ref_c, t_ref_n = np.nan, 0

        nc["gauge_station_id"][out_i] = str(gauge_ids[gj]) if gj >= 0 else ""
        nc["gauge_distance_km"][out_i] = np.float32(gd)
        nc["temperature_station_id"][out_i] = str(temp_ids[tj]) if tj >= 0 else ""
        nc["temperature_distance_km"][out_i] = np.float32(td)
        # Scalar temperature_ref_c is diagnostic only; retrieval uses T(t).
        nc["temperature_ref_c"][out_i] = np.float32(t_ref_c)
        nc["temperature_ref_n"][out_i] = np.int32(t_ref_n)

        t5_time, temp_5min_c = aggregate_temperature_5min(
            cml_times,
            temp_1min_c,
            EVAL_START,
            EVAL_END_EXCLUSIVE,
            min_valid_minutes=args.min_temperature_minutes_5min,
        )
        if not t5_time.equals(eval_time):
            raise RuntimeError("Internal 5-min temperature time mismatch")
        nc["temperature_5min_c"][out_i, :] = temp_5min_c.astype(np.float32)

        # P.676 table: pressure fixed; temperature remains time varying.
        if np.isfinite(temp_1min_c).any():
            (
                temp_grid,
                rho_grid,
                gamma_table,
                monotonic,
                sensitivity_by_temp,
            ) = p676_temperature_lookup(
                s.frequency_ghz,
                args.pressure_hpa,
                args.temperature_grid_min_c,
                args.temperature_grid_max_c,
                args.temperature_grid_step_c,
                args.rho_max_gm3,
                args.rho_step_gm3,
            )

            gamma_rho0_t = gamma_rho0_for_temperature(
                temp_1min_c,
                temp_grid,
                gamma_table,
            )
            sensitivity = (
                float(np.interp(t_ref_c, temp_grid, sensitivity_by_temp))
                if np.isfinite(t_ref_c)
                and temp_grid[0] <= t_ref_c <= temp_grid[-1]
                else np.nan
            )
        else:
            temp_grid = np.asarray(
                [args.temperature_grid_min_c, args.temperature_grid_max_c],
                dtype=np.float64,
            )
            rho_grid = np.asarray([0.0, args.rho_max_gm3], dtype=np.float64)
            gamma_table = np.full((2, 2), np.nan, dtype=np.float32)
            monotonic = False
            sensitivity = np.nan
            gamma_rho0_t = np.full(cml_times.size, np.nan, dtype=np.float64)

        nc["p676_monotonic"][out_i] = np.uint8(1 if monotonic else 0)
        nc["p676_sensitivity_dbkm_per_gm3"][out_i] = np.float32(sensitivity)

        # GRAF masks
        graf_wet, graf_known, graf_buffered, graf_q80, graf_thr = graf_masks(
            tl,
            window=args.graf_window_min,
            quantile=args.graf_quantile,
            factor=args.graf_factor,
            threshold_min_periods=args.graf_threshold_min_periods,
            detection_min_periods=args.graf_detection_min_periods,
            wet_pad_minutes=wet_pad_minutes,
        )
        nc["graf_q80_db"][out_i] = np.float32(graf_q80)
        nc["graf_threshold_db"][out_i] = np.float32(graf_thr)

        # Gauge masks
        if gauge_ok:
            if gj not in gauge_mask_cache:
                gauge_amount = np.asarray(
                    gauge_da.isel({gauge_station_dim: gj}).values,
                    dtype=np.float64,
                )
                gauge_mask_cache[gj] = gauge_masks_on_cml_minutes(
                    cml_times,
                    gauge_times,
                    gauge_amount,
                    wet_pad_minutes,
                )
            gauge_wet, gauge_known, gauge_buffered = gauge_mask_cache[gj]
        else:
            gauge_wet = np.zeros(tl.size, dtype=bool)
            gauge_known = np.zeros(tl.size, dtype=bool)
            gauge_buffered = np.zeros(tl.size, dtype=bool)

        method_masks = {
            "graf": (graf_known, graf_buffered),
            "gauge": (gauge_known, gauge_buffered),
        }

        row = {
            "series_index": out_i,
            "cml_id": s.cml_id,
            "sublink_id": s.sublink_id,
            "frequency_ghz": s.frequency_ghz,
            "length_km": s.length_km,
            "polarization": s.polarization,
            "raw_finite_tl": int(np.isfinite(tl).sum()),
            "gauge_station_id": str(gauge_ids[gj]) if gj >= 0 else "",
            "gauge_distance_km": gd,
            "gauge_within_max_distance": bool(gauge_ok),
            "temperature_station_id": str(temp_ids[tj]) if tj >= 0 else "",
            "temperature_distance_km": td,
            "temperature_ref_c": t_ref_c,
            "temperature_ref_n": t_ref_n,
            "temperature_eval_min_c": (
                float(np.nanmin(temp_5min_c))
                if np.isfinite(temp_5min_c).any() else np.nan
            ),
            "temperature_eval_median_c": (
                float(np.nanmedian(temp_5min_c))
                if np.isfinite(temp_5min_c).any() else np.nan
            ),
            "temperature_eval_max_c": (
                float(np.nanmax(temp_5min_c))
                if np.isfinite(temp_5min_c).any() else np.nan
            ),
            "p676_monotonic": bool(monotonic),
            "p676_sensitivity_dbkm_per_gm3": sensitivity,
            "graf_q80_db": graf_q80,
            "graf_threshold_db": graf_thr,
        }

        for mi, method in enumerate(METHODS):
            known, wet_buffered = method_masks[method]
            safe = np.isfinite(tl) & known & ~wet_buffered
            baseline_safe = safe & np.isfinite(gamma_rho0_t)

            # Internal anchor samples with measured T(t):
            # B_sample(t) = TL(t) - L * gamma_P676(f,P,T(t),rho=0).
            b_sample = tl - (s.length_km * gamma_rho0_t)

            main_values = b_sample[baseline_safe & cal_main]
            pre_values = b_sample[baseline_safe & cal_pre]
            post_values = b_sample[baseline_safe & cal_post]

            B, nB = robust_quantile(
                main_values, args.baseline_quantile, args.min_baseline_samples
            )
            Bpre, npre = robust_quantile(
                pre_values, args.baseline_quantile, args.min_prepost_samples
            )
            Bpost, npost = robust_quantile(
                post_values, args.baseline_quantile, args.min_prepost_samples
            )
            delta = Bpost - Bpre if np.isfinite(Bpost) and np.isfinite(Bpre) else np.nan

            nc["baseline_B_db"][mi, out_i] = np.float32(B)
            nc["baseline_pre_db"][mi, out_i] = np.float32(Bpre)
            nc["baseline_post_db"][mi, out_i] = np.float32(Bpost)
            nc["baseline_delta_post_minus_pre_db"][mi, out_i] = np.float32(delta)
            nc["baseline_n"][mi, out_i] = np.int64(nB)
            nc["baseline_pre_n"][mi, out_i] = np.int64(npre)
            nc["baseline_post_n"][mi, out_i] = np.int64(npost)

            t5, tl5, valid5, safe5 = aggregate_5min_complete(
                cml_times,
                tl,
                safe,
                EVAL_START,
                EVAL_END_EXCLUSIVE,
                min_valid_minutes=args.min_valid_minutes_5min,
            )
            if not t5.equals(eval_time):
                raise RuntimeError("Internal 5-min output time mismatch")

            nc["safe_dry_5min"][mi, out_i, :] = safe5

            # Baseline must be valid before retrieval. Pressure is constant;
            # measured station T varies at every five-minute retrieval interval.
            if np.isfinite(B) and monotonic:
                gamma_obs = (tl5 - B) / s.length_km
                wvd, status = invert_gamma_lookup_temperature(
                    gamma_obs,
                    temp_5min_c,
                    temp_grid,
                    rho_grid,
                    gamma_table,
                    monotonic,
                )
            else:
                wvd = np.full(eval_time.size, np.nan, dtype=np.float64)
                status = np.zeros(eval_time.size, dtype=np.uint8)
                if not monotonic:
                    status[np.isfinite(tl5) & np.isfinite(temp_5min_c)] = 4

            nc["total_loss_5min"][mi, out_i, :] = tl5.astype(np.float32)
            nc["wvd_internal_anchor"][mi, out_i, :] = wvd.astype(np.float32)
            nc["retrieval_status"][mi, out_i, :] = status

            row[f"{method}_wet_raw_minutes"] = int(
                (graf_wet if method == "graf" else gauge_wet).sum()
            )
            row[f"{method}_wet_buffered_minutes"] = int(wet_buffered.sum())
            row[f"{method}_safe_minutes"] = int(safe.sum())
            row[f"{method}_baseline_safe_with_temperature_minutes"] = int(
                (baseline_safe & cal_main).sum()
            )
            row[f"{method}_baseline_B_db"] = B
            row[f"{method}_baseline_n"] = nB
            row[f"{method}_baseline_pre_db"] = Bpre
            row[f"{method}_baseline_pre_n"] = npre
            row[f"{method}_baseline_post_db"] = Bpost
            row[f"{method}_baseline_post_n"] = npost
            row[f"{method}_baseline_delta_db"] = delta
            row[f"{method}_eval_safe_5min"] = int(valid5.sum())
            row[f"{method}_eval_valid_wvd_5min"] = int(np.isfinite(wvd).sum())

        diagnostics.append(row)

        nc.sync()

        elapsed = time.time() - started
        if (
            out_i < 10
            or (out_i + 1) % 25 == 0
            or out_i + 1 == len(selected)
        ):
            print(
                f"[{out_i+1:4d}/{len(selected):4d}] "
                f"{s.cml_id} sub={s.sublink_id} "
                f"f={s.frequency_ghz:.3f} GHz L={s.length_km:.2f} km "
                f"Tmed={t_ref_c:.1f} C gauge={gd:.2f} km "
                f"({elapsed:.1f}s)"
            )

    nc.close()
    cml.close()
    gauge.close()
    temp.close()

    diag_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(diagnostics)
    df.to_csv(diag_path, index=False)

    # Compact JSON summary next to NetCDF.
    summary = {
        "script_version": SCRIPT_VERSION,
        "n_series": len(selected),
        "n_physical_cml": n_phys,
        "frequency_min_ghz": args.freq_min_ghz,
        "frequency_max_ghz": args.freq_max_ghz,
        "length_min_km": args.length_min_km,
        "length_max_km": args.length_max_km,
        "evaluation_start": str(EVAL_START),
        "evaluation_end_exclusive": str(EVAL_END_EXCLUSIVE),
        "methods": list(METHODS),
        "pressure_hpa": args.pressure_hpa,
        "temperature_input": "nearest CHMI station T(t), time varying",
        "temperature_max_gap_min": args.temperature_max_gap_min,
        "wet_buffer_hours": args.wet_buffer_hours,
        "baseline_quantile": args.baseline_quantile,
        "no_external_humidity_calibration": True,
        "output_netcdf": str(out_path),
        "diagnostics_csv": str(diag_path),
        "elapsed_seconds": time.time() - t0_run,
    }
    summary_path = out_path.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print()
    print("=" * 100)
    print("DONE")
    print("=" * 100)
    print(f"NetCDF     : {out_path}")
    print(f"Diagnostics: {diag_path}")
    print(f"Summary    : {summary_path}")
    print(f"Elapsed    : {time.time() - t0_run:.1f} s")
    print()
    print("IMPORTANT:")
    print(
        "The WVD field is internally zero-anchored. No RH data were used, "
        "but this is not equivalent to knowing the absolute telecom link budget."
    )


if __name__ == "__main__":
    main()
