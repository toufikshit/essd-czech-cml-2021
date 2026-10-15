#!/usr/bin/env python3
"""
19_wvd_eband_progressive_baseline_v6.py
===========================

Causal progressive-baseline E-band CML water-vapour-density (WVD) retrieval (V6).

This is the repair of the single-baseline retrieval after the raw TX/RX/TL
diagnostic showed persistent telecom level changes in a substantial fraction
of E-band directional sublinks.

Core design
-----------
* CML source: RAW 1-min TX/RX; TL = TX - RX.
* Population by default: Fencl-style E-band only
      71--76 GHz OR 81--86 GHz
  and path length 2--10 km.
* Every directional sublink is processed independently.
* Hardware-state boundaries are detected ONLY from raw CML TX/RX/TL.
  No RH/WVD validation data are read.
* State detection is intentionally conservative:
    - persistent before/after daily-median TL comparison,
    - large TL steps are accepted directly,
    - smaller TL steps require TX or strong RX corroboration,
    - long acquisition gaps are candidate boundaries only if the post-gap
      persistent TL level changes sufficiently,
    - Q05 is retained only as a diagnostic; it does NOT define a boundary.
* A transition buffer is excluded around accepted hardware boundaries.
* Hardware boundaries are RESET points only; no recursive Q05-only state
  splitting is used.
* The internal baseline is causal and progressive. For each hardware regime,
  processing moves forward through time and uses only COMPLETED PAST days:
    1) daily low statistic Q05[TL - A_dry],
    2) trailing 45-day history of those daily lows,
    3) Q10 of the trailing daily-low history as the baseline candidate,
    4) a conservative progressive update of B(t).
* By default B(t) may move downward by at most 1 dB/day but may NOT drift
  upward during ordinary tracking. This prevents short/medium humid periods
  from automatically ratcheting the radio baseline upward.
* V6 adds a separate persistent unsupported-floor reinitialization rule. If
  the recent causal lower envelope remains >3 dB above the active baseline
  for about 14 days (with enough valid daily lows and low recent spread), the
  baseline is reinitialized upward from those completed past days. This is
  logged as a BASELINE REINITIALIZATION, not as a hardware event.
* After every confirmed hardware reset there is a 14-day causal warm-up;
  WVD remains missing until enough post-reset history exists.
* The same B(t) is used by GRAF and gauge. No wet/dry method can alter it.
* Wet/dry branches:
    - GRAF/RSTD CML detector,
    - nearest CHMI precipitation gauge.
* +/-6 h rain/wet protection retained.
* Temperature: measured time-varying nearest-station T(t).
* Pressure: constant, default 980 hPa.
* P.676 inversion at actual frequency and path length.
* No humidity observation enters state detection, baseline estimation,
  retrieval, or population selection.

Important
---------
The internal absolute reference remains an assumption. V6 uses a
single constant B_{s,k} with a causal progressive B_s(t) learned from the
lower envelope of TL - L*gamma_P676(rho=0), with confirmed telecom changes
used only as reset points.

This does NOT solve the fundamental lack of a complete telecom link budget.
Slow radio changes and atmospheric attenuation remain partly confounded.

Outputs
-------
1) NetCDF progressive-baseline Apr-Oct 2021 5-min WVD product.
2) Per-series/per-method diagnostics CSV.
3) Hardware reset-regime / progressive-baseline diagnostics CSV.
4) Accepted hardware change-point table CSV.
5) Progressive-baseline reinitialization table CSV.

The NetCDF keeps `wvd_internal_anchor` for compatibility with the existing
validation scripts, but it now contains the V6 PROGRESSIVE internal-anchor
retrieval. `baseline_progressive_5min_db` records the actual causal B(t) used;
`baseline_state_5min_db` is retained as a compatibility alias.

Example
-------
PY=~/my-conda-envs/tf-env/bin/python
ROOT=/home/jovyan/ESSD_publication_dataset

$PY "$ROOT/19_wvd_eband_progressive_baseline_v6.py" \
  --helper-script "$ROOT/13_wvd_raw_60_100_independent.py" \
  --cml-nc "$ROOT/CZ_CML_RAW_1MIN_20200101_20211231.nc" \
  --gauge-nc /home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_PRECIPITATION_2021.nc \
  --temperature-nc /home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_TEMPERATURE_2021.nc \
  --out-nc "$ROOT/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_5MIN_APR_OCT_2021.nc" \
  --diagnostics-csv "$ROOT/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv" \
  --states-csv "$ROOT/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_REGIMES.csv" \
  --events-csv "$ROOT/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_EVENTS.csv" \
  --reinit-csv "$ROOT/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_REINITIALIZATIONS.csv" \
  --overwrite
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd
import xarray as xr

try:
    from netCDF4 import date2num
except Exception as exc:
    raise ImportError("netCDF4 is required") from exc


METHODS = ("graf", "gauge")
EVAL_START = pd.Timestamp("2021-04-01 00:00:00")
EVAL_END_EXCLUSIVE = pd.Timestamp("2021-11-01 00:00:00")
STATE_START = pd.Timestamp("2021-01-01 00:00:00")
STATE_END = pd.Timestamp("2022-01-01 00:00:00")


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------

def load_module(path: Path, name: str):
    """Dynamic import compatible with Python 3.12 @dataclass modules."""
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not import {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
    except Exception:
        sys.modules.pop(name, None)
        raise
    return mod


def clean_time_index(values) -> pd.DatetimeIndex:
    """
    Force decoded timestamps to nanosecond precision.

    This avoids pandas/xarray datetime64[us] ambiguity that previously
    corrupted a NetCDF time coordinate when raw integer datetime values
    were divided by 1e9.
    """
    a = np.asarray(values)
    t = pd.DatetimeIndex(pd.to_datetime(a))
    if t.tz is not None:
        t = t.tz_convert("UTC").tz_localize(None)
    # Explicitly force ns irrespective of incoming datetime64 unit.
    return pd.DatetimeIndex(t.to_numpy(dtype="datetime64[ns]"))


def interpolate_temperature_to_minutes_safe(
    native_times: pd.DatetimeIndex,
    native_values_c: np.ndarray,
    target_times: pd.DatetimeIndex,
    max_gap_minutes: float,
) -> np.ndarray:
    """
    Linear interpolation without extrapolation and without bridging long gaps.
    Integer arithmetic is explicitly nanoseconds.
    """
    native_times = clean_time_index(native_times)
    target_times = clean_time_index(target_times)
    values = np.asarray(native_values_c, dtype=np.float64).reshape(-1)

    if len(values) != len(native_times):
        raise ValueError("Temperature value/time length mismatch")

    good = np.isfinite(values)
    if good.sum() < 2:
        return np.full(len(target_times), np.nan, dtype=np.float64)

    x = native_times.to_numpy(dtype="datetime64[ns]").astype("int64")[good]
    y = values[good]
    q = target_times.to_numpy(dtype="datetime64[ns]").astype("int64")

    pos = np.searchsorted(x, q, side="left")
    out = np.full(q.size, np.nan, dtype=np.float64)

    exact = pos < x.size
    idx = np.flatnonzero(exact)
    if idx.size:
        eq = x[pos[idx]] == q[idx]
        ii = idx[eq]
        if ii.size:
            out[ii] = y[pos[ii]]

    between = (pos > 0) & (pos < x.size) & ~np.isfinite(out)
    ii = np.flatnonzero(between)
    if ii.size:
        l = pos[ii] - 1
        r = pos[ii]
        gap_minutes = (x[r] - x[l]).astype(np.float64) / (60.0 * 1e9)
        keep = gap_minutes <= float(max_gap_minutes)
        jj = ii[keep]
        if jj.size:
            ll = l[keep]
            rr = r[keep]
            frac = (
                (q[jj] - x[ll]).astype(np.float64)
                / (x[rr] - x[ll]).astype(np.float64)
            )
            out[jj] = y[ll] + frac * (y[rr] - y[ll])

    return out


def robust_quantile(x, q, min_n):
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    n = int(x.size)
    if n < int(min_n):
        return np.nan, n
    return float(np.quantile(x, q)), n


# ---------------------------------------------------------------------------
# Conservative hardware-state detector
# ---------------------------------------------------------------------------

def daily_radio_statistics(
    times: pd.DatetimeIndex,
    tx: np.ndarray,
    rx: np.ndarray,
    valid: np.ndarray,
    min_valid_minutes: int,
    lower_quantile: float,
) -> pd.DataFrame:
    """Daily medians drive state detection; Q_low is diagnostic only."""
    in_year = (times >= STATE_START) & (times < STATE_END)
    t = times[in_year]
    v = valid[in_year]
    txv = np.where(v, np.asarray(tx[in_year], dtype=float), np.nan)
    rxv = np.where(v, np.asarray(rx[in_year], dtype=float), np.nan)
    tlv = txv - rxv

    df = pd.DataFrame({"tx": txv, "rx": rxv, "tl": tlv}, index=t)
    g = df.resample("1D")

    out = pd.DataFrame(index=pd.date_range(STATE_START, STATE_END - pd.Timedelta(days=1), freq="1D"))
    out["n"] = g["tl"].count().reindex(out.index, fill_value=0)
    out["tx_med"] = g["tx"].median().reindex(out.index)
    out["rx_med"] = g["rx"].median().reindex(out.index)
    out["tl_med"] = g["tl"].median().reindex(out.index)
    out["tl_low"] = g["tl"].quantile(float(lower_quantile)).reindex(out.index)

    bad = out["n"] < int(min_valid_minutes)
    out.loc[bad, ["tx_med", "rx_med", "tl_med", "tl_low"]] = np.nan
    return out


def side_medians(
    daily: pd.DataFrame,
    i: int,
    window_days: int,
    min_days_each_side: int,
) -> Optional[dict]:
    """Robust median radio levels in fixed windows before and after day i."""
    lo0 = max(0, i - int(window_days))
    hi1 = min(len(daily), i + int(window_days))

    pre = daily.iloc[lo0:i]
    post = daily.iloc[i:hi1]

    # Require TL support on both sides. Component medians then use same days.
    pre_valid = pre["tl_med"].dropna()
    post_valid = post["tl_med"].dropna()
    if (
        len(pre_valid) < int(min_days_each_side)
        or len(post_valid) < int(min_days_each_side)
    ):
        return None

    out = {}
    for col in ("tl_med", "tx_med", "rx_med", "tl_low"):
        a = pre[col].dropna()
        b = post[col].dropna()
        if len(a) == 0 or len(b) == 0:
            out[f"pre_{col}"] = np.nan
            out[f"post_{col}"] = np.nan
            out[f"delta_{col}"] = np.nan
        else:
            av = float(a.median())
            bv = float(b.median())
            out[f"pre_{col}"] = av
            out[f"post_{col}"] = bv
            out[f"delta_{col}"] = bv - av

    out["n_pre_days"] = int(len(pre_valid))
    out["n_post_days"] = int(len(post_valid))
    return out


def long_gap_boundaries(
    times: pd.DatetimeIndex,
    valid: np.ndarray,
    gap_hours: float,
) -> list[pd.Timestamp]:
    """First finite sample after every >= gap_hours interruption in 2021."""
    in_year = (times >= STATE_START) & (times < STATE_END)
    idx = np.flatnonzero(in_year & valid)
    if len(idx) < 2:
        return []

    t_ns = times.to_numpy(dtype="datetime64[ns]").astype("int64")
    hours = (t_ns[idx[1:]] - t_ns[idx[:-1]]) / (3600.0 * 1e9)
    p = np.flatnonzero(hours >= float(gap_hours))
    return [pd.Timestamp(times[idx[j + 1]]) for j in p]


def near_gap(boundary: pd.Timestamp, gaps: list[pd.Timestamp], hours: float = 36.0) -> bool:
    if not gaps:
        return False
    b = pd.Timestamp(boundary)
    tol = pd.Timedelta(hours=float(hours))
    return any(abs(g - b) <= tol for g in gaps)


def classify_state_event(dtl, dtx, drx):
    if not np.isfinite(dtl):
        return "unknown"
    if np.isfinite(dtx) and abs(dtx) >= 0.5:
        if np.isfinite(drx) and abs(drx) >= 0.5:
            return "TX+RX"
        return "TX"
    if np.isfinite(drx) and abs(drx) >= 1.5:
        return "RX"
    if abs(dtl) >= 3.0:
        return "large_TL"
    return "other"


def detect_hardware_events(
    times: pd.DatetimeIndex,
    tx: np.ndarray,
    rx: np.ndarray,
    valid: np.ndarray,
    args,
) -> Tuple[list[dict], pd.DataFrame]:
    """
    Accepted boundary rules.

    A) large persistent TL:
         |Δ median TL| >= --state-large-tl-step-db

    B) TX-confirmed:
         |Δ median TL| >= --state-min-tl-with-tx-db
         AND |Δ median TX| >= --state-tx-step-db

    C) strong RX-confirmed:
         |Δ median TL| >= --state-rx-tl-step-db
         AND |Δ median RX| >= --state-rx-step-db

    D) long-gap confirmed:
         raw finite-data gap >= --state-gap-hours
         AND |Δ median TL| >= --state-gap-tl-step-db

    Daily Q_low is saved only for diagnostics and never accepts a boundary.
    """
    daily = daily_radio_statistics(
        times, tx, rx, valid,
        min_valid_minutes=args.state_min_valid_minutes_per_day,
        lower_quantile=args.state_lower_quantile,
    )
    gaps = long_gap_boundaries(times, valid, args.state_gap_hours)

    candidates = []
    for i in range(1, len(daily)):
        z = side_medians(
            daily, i,
            args.state_window_days,
            args.state_min_days_each_side,
        )
        if z is None:
            continue

        dtl = z["delta_tl_med"]
        dtx = z["delta_tx_med"]
        drx = z["delta_rx_med"]
        if not np.isfinite(dtl):
            continue

        boundary = pd.Timestamp(daily.index[i])
        is_gap = near_gap(boundary, gaps)

        reasons = []
        if abs(dtl) >= float(args.state_large_tl_step_db):
            reasons.append("large_TL")
        if (
            np.isfinite(dtx)
            and abs(dtx) >= float(args.state_tx_step_db)
            and abs(dtl) >= float(args.state_min_tl_with_tx_db)
        ):
            reasons.append("TX_confirmed")
        if (
            np.isfinite(drx)
            and abs(drx) >= float(args.state_rx_step_db)
            and abs(dtl) >= float(args.state_rx_tl_step_db)
        ):
            reasons.append("RX_confirmed")
        if is_gap and abs(dtl) >= float(args.state_gap_tl_step_db):
            reasons.append("gap_confirmed")

        if not reasons:
            continue

        candidates.append(
            {
                "boundary_time": boundary,
                "reason": "+".join(reasons),
                "near_long_gap": int(is_gap),
                "origin": classify_state_event(dtl, dtx, drx),
                **z,
            }
        )

    # Add exact post-gap timestamps for qualifying gaps; this improves timing
    # when a long acquisition interruption brackets the state change.
    for g in gaps:
        # Compare around the day containing the first post-gap observation.
        day = g.floor("D")
        i = int(np.searchsorted(daily.index.values, day.to_datetime64()))
        i = min(max(i, 1), len(daily) - 1)
        z = side_medians(
            daily, i,
            args.state_window_days,
            args.state_min_days_each_side,
        )
        if z is None or not np.isfinite(z["delta_tl_med"]):
            continue
        if abs(z["delta_tl_med"]) < float(args.state_gap_tl_step_db):
            continue
        candidates.append(
            {
                "boundary_time": pd.Timestamp(g),
                "reason": "gap_confirmed_exact",
                "near_long_gap": 1,
                "origin": classify_state_event(
                    z["delta_tl_med"], z["delta_tx_med"], z["delta_rx_med"]
                ),
                **z,
            }
        )

    if not candidates:
        return [], daily

    # Merge multiple daily detections around the same physical transition.
    candidates = sorted(candidates, key=lambda r: pd.Timestamp(r["boundary_time"]))
    groups = [[candidates[0]]]
    merge_dt = pd.Timedelta(days=float(args.state_merge_days))
    for r in candidates[1:]:
        if (
            pd.Timestamp(r["boundary_time"])
            - pd.Timestamp(groups[-1][-1]["boundary_time"])
        ) <= merge_dt:
            groups[-1].append(r)
        else:
            groups.append([r])

    accepted = []
    for group in groups:
        # Prefer an exact gap timestamp when its persistent median shift is
        # essentially as strong as the strongest member.
        strongest = max(group, key=lambda r: abs(r["delta_tl_med"]))
        gap_exact = [r for r in group if r["reason"] == "gap_confirmed_exact"]
        if gap_exact:
            gx = max(gap_exact, key=lambda r: abs(r["delta_tl_med"]))
            if abs(gx["delta_tl_med"]) >= 0.8 * abs(strongest["delta_tl_med"]):
                best = gx
            else:
                best = strongest
        else:
            best = strongest

        best = dict(best)
        best["merged_candidate_count"] = len(group)
        best["all_reasons"] = "+".join(
            sorted(set(x["reason"] for x in group))
        )
        accepted.append(best)

    return accepted, daily


def assign_state_minutes(
    times: pd.DatetimeIndex,
    events: list[dict],
    transition_buffer_hours: float,
) -> Tuple[np.ndarray, np.ndarray, list[pd.Timestamp]]:
    """State ID for every raw minute; -1 outside 2021."""
    boundaries = sorted(
        pd.Timestamp(r["boundary_time"])
        for r in events
        if STATE_START < pd.Timestamp(r["boundary_time"]) < STATE_END
    )

    state = np.full(len(times), -1, dtype=np.int16)
    in_year = (times >= STATE_START) & (times < STATE_END)
    if in_year.any():
        tyear = times[in_year]
        b_ns = np.asarray(
            [b.to_datetime64().astype("datetime64[ns]").astype("int64") for b in boundaries],
            dtype=np.int64,
        )
        q_ns = tyear.to_numpy(dtype="datetime64[ns]").astype("int64")
        sid = np.searchsorted(b_ns, q_ns, side="right").astype(np.int16)
        state[in_year] = sid

    transition = np.zeros(len(times), dtype=bool)
    pad = pd.Timedelta(hours=float(transition_buffer_hours))
    for b in boundaries:
        transition |= (times >= b - pad) & (times <= b + pad)

    return state, transition, boundaries


def aggregate_state_to_5min(
    times: pd.DatetimeIndex,
    state_min: np.ndarray,
    transition_min: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return exact state ID and transition flag on the Apr-Oct 5-min grid.

    A state is emitted only if all 5 constituent minutes have the same
    non-negative state ID and none is in the transition buffer.
    """
    full = pd.date_range(
        EVAL_START,
        EVAL_END_EXCLUSIVE - pd.Timedelta(minutes=1),
        freq="1min",
    )
    s = pd.Series(np.asarray(state_min, dtype=np.int16), index=times).reindex(full)
    tr = pd.Series(
        np.asarray(transition_min, dtype=np.uint8), index=times
    ).reindex(full)

    a = s.to_numpy(dtype=float).reshape(-1, 5)
    b = tr.to_numpy(dtype=float).reshape(-1, 5)

    same = (
        np.isfinite(a).all(axis=1)
        & (a[:, 0, None] == a).all(axis=1)
        & (a[:, 0] >= 0)
    )
    trans = np.isfinite(b).all(axis=1) & (b > 0).any(axis=1)

    out = np.full(a.shape[0], -1, dtype=np.int16)
    out[same & ~trans] = a[same & ~trans, 0].astype(np.int16)
    return out, trans


def aggregate_baseline_to_5min(
    times: pd.DatetimeIndex,
    baseline_min: np.ndarray,
) -> np.ndarray:
    full = pd.date_range(
        EVAL_START,
        EVAL_END_EXCLUSIVE - pd.Timedelta(minutes=1),
        freq="1min",
    )
    s = pd.Series(np.asarray(baseline_min, dtype=float), index=times).reindex(full)
    a = s.to_numpy(dtype=float).reshape(-1, 5)
    good = np.isfinite(a).all(axis=1) & np.all(
        np.isclose(a, a[:, [0]], atol=1e-6, rtol=0.0) | ~np.isfinite(a),
        axis=1,
    )
    out = np.full(a.shape[0], np.nan, dtype=float)
    out[good] = np.mean(a[good], axis=1)
    return out


def state_intervals(boundaries: list[pd.Timestamp]) -> list[Tuple[int, pd.Timestamp, pd.Timestamp]]:
    edges = [STATE_START] + list(boundaries) + [STATE_END]
    return [(k, edges[k], edges[k + 1]) for k in range(len(edges) - 1)]



# ---------------------------------------------------------------------------
# Causal progressive CML-only baseline
# ---------------------------------------------------------------------------

def progressive_causal_baseline(
    times: pd.DatetimeIndex,
    b_sample: np.ndarray,
    state_min: np.ndarray,
    transition_min: np.ndarray,
    intervals: list[Tuple[int, pd.Timestamp, pd.Timestamp]],
    args,
) -> Tuple[np.ndarray, list[dict], list[dict], np.ndarray]:
    """
    Build the V6 causal progressive radio baseline B(t).

    b_sample = TL - L*gamma_P676(rho=0).

    For each CONFIRMED telecom reset regime, only completed PAST calendar
    days are eligible to affect the baseline used for the current day.

    Daily statistic
    ---------------
        D_d = Q_daily_low(b_sample on day d)

    Ordinary progressive candidate
    ------------------------------
        C_d = Q_history(D_j),  d-W <= j < d

    Ordinary update
    ---------------
      * first valid C_d initializes B;
      * downward movement is capped by --progressive-max-down-db-per-day;
      * ordinary upward movement is capped by
        --progressive-max-up-db-per-day (default 0).

    V6 persistent unsupported-floor reinitialization
    ------------------------------------------------
    A second, deliberately conservative causal test asks whether the RECENT
    daily lower envelope no longer supports the active baseline. For current
    day d, completed daily lows from the previous
    --progressive-reinit-persistence-days are examined. If:

      1) enough valid days are present,
      2) at least --progressive-reinit-min-support-fraction of those daily
         lows are > --progressive-reinit-threshold-db above active B,
      3) their recent lower quantile is also above that threshold,
      4) that recent floor is >= --progressive-reinit-min-floor-step-db above
         the immediately preceding lower-floor window,
      5) their IQR is <= --progressive-reinit-max-iqr-db, and
      6) cooldown / maximum-reinitialization constraints are satisfied,

    then B is reinitialized upward to that recent lower-floor quantile.

    This is NOT labelled a hardware event. It means only that the old radio
    floor is persistently unsupported by recent CML observations. Every such
    event is written to a separate diagnostic table.

    After reinitialization, the ordinary long-window history is restarted
    from the beginning of the recent support window so obsolete low-floor
    observations cannot immediately pull the new baseline back down.

    No GRAF/gauge/RH/WVD information enters this estimator.
    """
    times = clean_time_index(times)
    x = np.asarray(b_sample, dtype=np.float64)
    state_min = np.asarray(state_min, dtype=np.int16)
    transition_min = np.asarray(transition_min, dtype=bool)

    out = np.full(len(times), np.nan, dtype=np.float64)
    reinit_flag_min = np.zeros(len(times), dtype=bool)
    regime_rows: list[dict] = []
    reinit_rows: list[dict] = []

    one_day = pd.Timedelta(days=1)
    window = pd.Timedelta(days=float(args.progressive_window_days))
    warmup = pd.Timedelta(days=float(args.progressive_warmup_days))
    reinit_span = pd.Timedelta(days=float(args.progressive_reinit_persistence_days))
    cooldown = pd.Timedelta(days=float(args.progressive_reinit_cooldown_days))

    for sid, a0, b0 in intervals:
        a = pd.Timestamp(a0)
        b = pd.Timestamp(b0)
        if b <= a:
            continue

        regime = (
            (state_min == int(sid))
            & ~transition_min
            & np.isfinite(x)
            & (times >= a)
            & (times < b)
        )

        first_day = a.floor('D')
        last_day = (b - pd.Timedelta(nanoseconds=1)).floor('D')
        days = pd.date_range(first_day, last_day, freq='D')

        daily = []
        for day in days:
            lo = max(pd.Timestamp(day), a)
            hi = min(pd.Timestamp(day) + one_day, b)
            m = regime & (times >= lo) & (times < hi)
            vals = x[m]
            vals = vals[np.isfinite(vals)]
            n = int(vals.size)
            if n >= int(args.progressive_min_valid_minutes_per_day):
                q = float(np.quantile(vals, float(args.progressive_daily_low_quantile)))
                daily.append({'day': pd.Timestamp(day), 'low_db': q, 'n_minutes': n})

        current_B = np.nan
        first_valid_time = pd.NaT
        n_updates = 0
        n_down_updates = 0
        n_up_updates = 0
        n_reinitializations = 0
        candidate_values = []
        reinit_jumps = []
        first_reinit_time = pd.NaT
        last_reinit_time = pd.NaT

        # The normal 45-day candidate ignores obsolete pre-reinitialization
        # history after an unsupported-floor reset.
        history_start = a

        for day in days:
            day = pd.Timestamp(day)

            # Strictly causal: history ends before the day being retrieved.
            eligible = [
                r for r in daily
                if (r['day'] < day)
                and (r['day'] >= day - window)
                and (r['day'] >= history_start)
            ]

            enough_calendar_warmup = (day - a) >= warmup
            enough_days = len(eligible) >= int(args.progressive_min_history_days)
            n_hist_minutes = int(sum(r['n_minutes'] for r in eligible))
            enough_minutes = n_hist_minutes >= int(args.progressive_min_history_minutes)

            candidate = np.nan
            if enough_calendar_warmup and enough_days and enough_minutes:
                lows = np.asarray([r['low_db'] for r in eligible], dtype=float)
                lows = lows[np.isfinite(lows)]
                if lows.size >= int(args.progressive_min_history_days):
                    candidate = float(
                        np.quantile(lows, float(args.progressive_history_quantile))
                    )

            if np.isfinite(candidate):
                candidate_values.append(candidate)
                if not np.isfinite(current_B):
                    current_B = candidate
                    first_valid_time = day
                    n_updates += 1
                else:
                    delta = candidate - current_B
                    if delta < 0:
                        step = max(delta, -float(args.progressive_max_down_db_per_day))
                        if step != 0:
                            current_B += step
                            n_updates += 1
                            n_down_updates += 1
                    elif delta > 0:
                        # Ordinary upward drift remains disabled by default.
                        step = min(delta, float(args.progressive_max_up_db_per_day))
                        if step != 0:
                            current_B += step
                            n_updates += 1
                            n_up_updates += 1

            # --------------------------------------------------------------
            # V6 persistent unsupported-floor reinitialization.
            # Uses COMPLETED PAST daily lows only; never current-day data.
            # --------------------------------------------------------------
            if (
                np.isfinite(current_B)
                and bool(args.progressive_reinit_enabled)
                and n_reinitializations < int(args.progressive_reinit_max_per_regime)
            ):
                cooldown_ok = (
                    pd.isna(last_reinit_time)
                    or (day - pd.Timestamp(last_reinit_time)) >= cooldown
                )

                recent_start = day - reinit_span
                recent = [
                    r for r in daily
                    if (r['day'] < day)
                    and (r['day'] >= recent_start)
                    and (r['day'] >= history_start)
                ]

                if cooldown_ok and len(recent) >= int(args.progressive_reinit_min_valid_days):
                    recent_lows = np.asarray([r['low_db'] for r in recent], dtype=float)
                    recent_lows = recent_lows[np.isfinite(recent_lows)]

                    # Require the evidence to span most of the requested
                    # persistence interval, not merely be clustered in a few
                    # adjacent available days.
                    recent_days = [pd.Timestamp(r['day']) for r in recent]
                    recent_calendar_span = (
                        (max(recent_days) - min(recent_days)) / one_day + 1.0
                        if recent_days else 0.0
                    )
                    required_span = max(
                        1.0,
                        float(args.progressive_reinit_persistence_days)
                        - float(args.progressive_reinit_max_missing_days),
                    )

                    if (
                        recent_lows.size >= int(args.progressive_reinit_min_valid_days)
                        and recent_calendar_span >= required_span
                    ):
                        recent_floor = float(
                            np.quantile(
                                recent_lows,
                                float(args.progressive_reinit_recent_quantile),
                            )
                        )

                        previous_start = day - 2 * reinit_span
                        previous_end = day - reinit_span
                        previous = [
                            r for r in daily
                            if (r['day'] < previous_end)
                            and (r['day'] >= previous_start)
                            and (r['day'] >= history_start)
                        ]
                        previous_lows = np.asarray(
                            [r['low_db'] for r in previous], dtype=float
                        )
                        previous_lows = previous_lows[np.isfinite(previous_lows)]
                        previous_floor = (
                            float(np.quantile(
                                previous_lows,
                                float(args.progressive_reinit_recent_quantile),
                            ))
                            if previous_lows.size >= int(args.progressive_reinit_min_previous_days)
                            else np.nan
                        )
                        recent_floor_step = (
                            recent_floor - previous_floor
                            if np.isfinite(previous_floor) else np.nan
                        )

                        q25, q75 = np.quantile(recent_lows, [0.25, 0.75])
                        recent_iqr = float(q75 - q25)
                        threshold_level = current_B + float(args.progressive_reinit_threshold_db)
                        recent_support_fraction = float(
                            np.mean(recent_lows >= threshold_level)
                        )
                        excess = recent_floor - current_B

                        if (
                            recent_support_fraction >= float(args.progressive_reinit_min_support_fraction)
                            and excess >= float(args.progressive_reinit_threshold_db)
                            and np.isfinite(recent_floor_step)
                            and recent_floor_step >= float(args.progressive_reinit_min_floor_step_db)
                            and recent_iqr <= float(args.progressive_reinit_max_iqr_db)
                        ):
                            old_B = float(current_B)
                            new_B = float(recent_floor)
                            jump = new_B - old_B

                            current_B = new_B
                            n_reinitializations += 1
                            n_updates += 1
                            reinit_jumps.append(jump)
                            if pd.isna(first_reinit_time):
                                first_reinit_time = day
                            last_reinit_time = day

                            # Keep only the recent supporting floor history for
                            # subsequent ordinary candidate updates.
                            history_start = min(recent_days)

                            lo = max(day, a)
                            hi = min(day + one_day, b)
                            mre = (
                                (state_min == int(sid))
                                & (times >= lo)
                                & (times < hi)
                            )
                            reinit_flag_min[mre] = True

                            reinit_rows.append(
                                {
                                    'state_id': int(sid),
                                    'reinit_index_within_state': int(n_reinitializations - 1),
                                    'reinit_time': day,
                                    'reason': 'persistent_unsupported_lower_floor',
                                    'old_B_db': old_B,
                                    'new_B_db': new_B,
                                    'delta_B_db': jump,
                                    'recent_floor_db': recent_floor,
                                    'previous_floor_db': previous_floor,
                                    'recent_floor_step_db': recent_floor_step,
                                    'recent_iqr_db': recent_iqr,
                                    'previous_n_valid_days': int(previous_lows.size),
                                    'recent_n_valid_days': int(recent_lows.size),
                                    'recent_support_fraction': recent_support_fraction,
                                    'recent_calendar_span_days': float(recent_calendar_span),
                                    'recent_quantile': float(args.progressive_reinit_recent_quantile),
                                    'threshold_db': float(args.progressive_reinit_threshold_db),
                                    'persistence_days': float(args.progressive_reinit_persistence_days),
                                    'history_restart_day': history_start,
                                    'hardware_event': 0,
                                }
                            )

            # Baseline is constant within the current calendar day. Every
            # value used here was determined only from completed past days.
            if np.isfinite(current_B):
                lo = max(day, a)
                hi = min(day + one_day, b)
                mday = (
                    (state_min == int(sid))
                    & ~transition_min
                    & (times >= lo)
                    & (times < hi)
                )
                out[mday] = current_B

        m_state_out = (state_min == int(sid)) & np.isfinite(out)
        bvals = out[m_state_out]
        bvals = bvals[np.isfinite(bvals)]
        daily_lows = np.asarray([r['low_db'] for r in daily], dtype=float)
        total_candidate_minutes = int(sum(r['n_minutes'] for r in daily))
        valid_days = int(
            len(pd.DatetimeIndex(times[m_state_out]).floor('D').unique())
        ) if m_state_out.any() else 0

        regime_rows.append(
            {
                'state_id': int(sid),
                'state_start': a,
                'state_end_exclusive': b,
                'state_duration_days': float((b - a) / one_day),
                'baseline_valid': int(bvals.size > 0),
                'baseline_first_valid_time': first_valid_time,
                'baseline_first_db': float(bvals[0]) if bvals.size else np.nan,
                'baseline_last_db': float(bvals[-1]) if bvals.size else np.nan,
                'baseline_median_db': float(np.median(bvals)) if bvals.size else np.nan,
                'baseline_min_db': float(np.min(bvals)) if bvals.size else np.nan,
                'baseline_max_db': float(np.max(bvals)) if bvals.size else np.nan,
                'baseline_valid_days': valid_days,
                'baseline_n_daily_lows': int(daily_lows.size),
                'baseline_n_candidate_minutes': total_candidate_minutes,
                'baseline_n_updates': int(n_updates),
                'baseline_n_down_updates': int(n_down_updates),
                'baseline_n_up_updates': int(n_up_updates),
                'baseline_n_reinitializations': int(n_reinitializations),
                'baseline_first_reinit_time': first_reinit_time,
                'baseline_last_reinit_time': last_reinit_time,
                'baseline_reinit_total_up_db': float(np.sum(reinit_jumps)) if reinit_jumps else 0.0,
                'baseline_reinit_max_up_db': float(np.max(reinit_jumps)) if reinit_jumps else 0.0,
                'candidate_min_db': float(np.min(candidate_values)) if candidate_values else np.nan,
                'candidate_max_db': float(np.max(candidate_values)) if candidate_values else np.nan,
            }
        )

    out[transition_min] = np.nan
    reinit_flag_min[transition_min] = False
    return out, regime_rows, reinit_rows, reinit_flag_min

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="E-band causal progressive-baseline raw-CML WVD retrieval with no humidity calibration.",
    )
    root = "/home/jovyan/ESSD_publication_dataset"

    p.add_argument(
        "--helper-script",
        default=f"{root}/13_wvd_raw_60_100_independent.py",
        help="Existing time-varying-T retrieval script; its physical/GRAF/gauge helpers are reused.",
    )
    p.add_argument(
        "--cml-nc",
        default=f"{root}/CZ_CML_RAW_1MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--gauge-nc",
        default="/home/jovyan/opendata.chmi.cz/processed/"
                "CHMI_10MIN_PRECIPITATION_2021.nc",
    )
    p.add_argument(
        "--temperature-nc",
        default="/home/jovyan/opendata.chmi.cz/processed/"
                "CHMI_10MIN_TEMPERATURE_2021.nc",
    )
    p.add_argument(
        "--out-nc",
        default=f"{root}/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_5MIN_APR_OCT_2021.nc",
    )
    p.add_argument(
        "--diagnostics-csv",
        default=f"{root}/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv",
    )
    p.add_argument(
        "--states-csv",
        default=f"{root}/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_REGIMES.csv",
    )
    p.add_argument(
        "--events-csv",
        default=f"{root}/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_EVENTS.csv",
    )
    p.add_argument(
        "--reinit-csv",
        default=f"{root}/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_REINITIALIZATIONS.csv",
        help="CML-only progressive-baseline reinitialization diagnostics.",
    )

    # Population.
    p.add_argument("--freq-min-ghz", type=float, default=71.0)
    p.add_argument("--freq-max-ghz", type=float, default=86.0)
    p.add_argument("--band1-min-ghz", type=float, default=71.0)
    p.add_argument("--band1-max-ghz", type=float, default=76.0)
    p.add_argument("--band2-min-ghz", type=float, default=81.0)
    p.add_argument("--band2-max-ghz", type=float, default=86.0)
    p.add_argument(
        "--include-gap-76-81",
        action="store_true",
        help="Normally false: only 71-76 and 81-86 GHz are selected.",
    )
    p.add_argument("--length-min-km", type=float, default=2.0)
    p.add_argument("--length-max-km", type=float, default=10.0)

    p.add_argument("--cml-id", action="append", default=[],
                   help="Restrict to exact cml_id; repeat as needed.")
    p.add_argument("--sublink-id", action="append", default=[],
                   help="Optional exact sublink_id; repeat as needed.")

    # Wet/dry.
    p.add_argument("--graf-window-min", type=int, default=60)
    p.add_argument("--graf-quantile", type=float, default=0.80)
    p.add_argument("--graf-factor", type=float, default=1.12)
    p.add_argument("--graf-threshold-min-periods", type=int, default=60)
    p.add_argument("--graf-detection-min-periods", type=int, default=45)
    p.add_argument("--wet-buffer-hours", type=float, default=6.0)
    p.add_argument("--max-gauge-distance-km", type=float, default=5.0)

    # T/P/P.676.
    p.add_argument("--pressure-hpa", type=float, default=980.0)
    p.add_argument("--max-temperature-distance-km", type=float, default=20.0)
    p.add_argument("--temperature-max-gap-min", type=float, default=30.0)
    p.add_argument("--temperature-grid-min-c", type=float, default=-35.0)
    p.add_argument("--temperature-grid-max-c", type=float, default=45.0)
    p.add_argument("--temperature-grid-step-c", type=float, default=0.25)
    p.add_argument("--min-temperature-minutes-5min", type=int, default=3)
    p.add_argument("--rho-max-gm3", type=float, default=60.0)
    p.add_argument("--rho-step-gm3", type=float, default=0.02)

    # Progressive causal baseline.
    # The legacy baseline arguments are retained because helper.create_output
    # expects them, but the V6 retrieval below uses the progressive settings.
    p.add_argument("--baseline-quantile", type=float, default=0.01,
                   help="Compatibility attribute only in V6; not used to estimate B(t).")
    p.add_argument("--min-state-baseline-samples", type=int, default=2880,
                   help="Compatibility attribute only in V6.")
    p.add_argument("--min-baseline-samples", type=int, default=1440)
    p.add_argument("--min-prepost-samples", type=int, default=300)
    p.add_argument("--min-state-duration-days", type=float, default=3.0,
                   help="Compatibility attribute only in V6.")
    p.add_argument("--min-valid-minutes-5min", type=int, default=5)

    p.add_argument("--progressive-window-days", type=float, default=45.0)
    p.add_argument("--progressive-warmup-days", type=float, default=14.0)
    p.add_argument("--progressive-daily-low-quantile", type=float, default=0.05)
    p.add_argument("--progressive-history-quantile", type=float, default=0.10)
    p.add_argument("--progressive-min-valid-minutes-per-day", type=int, default=120)
    p.add_argument("--progressive-min-history-days", type=int, default=7)
    p.add_argument("--progressive-min-history-minutes", type=int, default=1440)
    p.add_argument(
        "--progressive-max-down-db-per-day",
        type=float,
        default=1.0,
        help="Largest normal downward B(t) update in one day.",
    )
    p.add_argument(
        "--progressive-max-up-db-per-day",
        type=float,
        default=0.0,
        help="Largest normal upward B(t) update in one day; default 0 prevents humid-season ratcheting.",
    )

    # V6 persistent unsupported-floor reinitialization. This is deliberately
    # separate from the strong TX/RX/TL hardware-event detector.
    p.add_argument(
        "--progressive-reinit-enabled",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Allow causal upward baseline reinitialization after persistent unsupported-floor evidence.",
    )
    p.add_argument("--progressive-reinit-threshold-db", type=float, default=3.0)
    p.add_argument("--progressive-reinit-persistence-days", type=float, default=14.0)
    p.add_argument("--progressive-reinit-min-valid-days", type=int, default=12)
    p.add_argument(
        "--progressive-reinit-min-support-fraction",
        type=float,
        default=0.90,
        help="Fraction of recent valid daily lows that must individually exceed active B by the reinit threshold.",
    )
    p.add_argument(
        "--progressive-reinit-min-floor-step-db",
        type=float,
        default=3.0,
        help="Minimum increase of recent lower floor over the immediately preceding floor window.",
    )
    p.add_argument(
        "--progressive-reinit-min-previous-days",
        type=int,
        default=7,
        help="Minimum valid daily lows required in the preceding comparison window.",
    )
    p.add_argument("--progressive-reinit-max-missing-days", type=float, default=3.0)
    p.add_argument("--progressive-reinit-recent-quantile", type=float, default=0.10)
    p.add_argument(
        "--progressive-reinit-max-iqr-db",
        type=float,
        default=3.0,
        help="Maximum IQR of recent daily lower-envelope values allowed for a reinitialization.",
    )
    p.add_argument("--progressive-reinit-cooldown-days", type=float, default=14.0)
    p.add_argument("--progressive-reinit-max-per-regime", type=int, default=3)

    # Conservative state detector.
    p.add_argument("--state-window-days", type=int, default=7)
    p.add_argument("--state-min-days-each-side", type=int, default=4)
    p.add_argument("--state-min-valid-minutes-per-day", type=int, default=120)
    p.add_argument(
        "--state-lower-quantile",
        type=float,
        default=0.05,
        help="Diagnostic only; state acceptance is based on persistent medians.",
    )
    p.add_argument(
        "--state-large-tl-step-db",
        type=float,
        default=3.0,
        help="Persistent TL median step accepted without further corroboration.",
    )
    p.add_argument("--state-tx-step-db", type=float, default=2.0)
    p.add_argument("--state-min-tl-with-tx-db", type=float, default=2.0)
    p.add_argument("--state-rx-step-db", type=float, default=3.0)
    p.add_argument("--state-rx-tl-step-db", type=float, default=3.0)
    p.add_argument("--state-gap-hours", type=float, default=12.0)
    p.add_argument("--state-gap-tl-step-db", type=float, default=1.5)
    p.add_argument("--state-merge-days", type=float, default=5.0)
    p.add_argument("--hardware-transition-buffer-hours", type=float, default=12.0)

    p.add_argument("--max-series", type=int, default=None)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not (0.0 < args.baseline_quantile <= 0.5):
        raise ValueError("--baseline-quantile must be in (0,0.5]")
    if not (0.0 < args.progressive_daily_low_quantile <= 0.5):
        raise ValueError("--progressive-daily-low-quantile must be in (0,0.5]")
    if not (0.0 < args.progressive_history_quantile <= 0.5):
        raise ValueError("--progressive-history-quantile must be in (0,0.5]")
    if args.progressive_window_days <= 0 or args.progressive_warmup_days < 0:
        raise ValueError("Invalid progressive window/warm-up settings")
    if args.progressive_min_history_days < 1 or args.progressive_min_history_minutes < 1:
        raise ValueError("Progressive history requirements must be positive")
    if args.progressive_max_down_db_per_day < 0 or args.progressive_max_up_db_per_day < 0:
        raise ValueError("Progressive B(t) update caps must be non-negative")
    if args.progressive_reinit_threshold_db <= 0:
        raise ValueError("--progressive-reinit-threshold-db must be >0")
    if args.progressive_reinit_persistence_days <= 0:
        raise ValueError("--progressive-reinit-persistence-days must be >0")
    if args.progressive_reinit_min_valid_days < 1:
        raise ValueError("--progressive-reinit-min-valid-days must be >=1")
    if not (0.0 < args.progressive_reinit_min_support_fraction <= 1.0):
        raise ValueError("--progressive-reinit-min-support-fraction must be in (0,1]")
    if args.progressive_reinit_min_floor_step_db <= 0:
        raise ValueError("--progressive-reinit-min-floor-step-db must be >0")
    if args.progressive_reinit_min_previous_days < 1:
        raise ValueError("--progressive-reinit-min-previous-days must be >=1")
    if not (0.0 < args.progressive_reinit_recent_quantile <= 0.5):
        raise ValueError("--progressive-reinit-recent-quantile must be in (0,0.5]")
    if args.progressive_reinit_max_iqr_db < 0:
        raise ValueError("--progressive-reinit-max-iqr-db must be >=0")
    if args.progressive_reinit_cooldown_days < 0:
        raise ValueError("--progressive-reinit-cooldown-days must be >=0")
    if args.progressive_reinit_max_per_regime < 0:
        raise ValueError("--progressive-reinit-max-per-regime must be >=0")
    if args.rho_max_gm3 <= 0 or args.rho_step_gm3 <= 0:
        raise ValueError("Invalid rho lookup settings")

    helper = load_module(
        Path(args.helper_script), "wvd_stateaware_physics_helper"
    )

    cml_path = helper.require_file(args.cml_nc, "CML raw NetCDF")
    gauge_path = helper.require_file(args.gauge_nc, "CHMI precipitation NetCDF")
    temp_path = helper.require_file(args.temperature_nc, "CHMI temperature NetCDF")

    out_path = Path(args.out_nc).expanduser().resolve()
    diag_path = Path(args.diagnostics_csv).expanduser().resolve()
    states_path = Path(args.states_csv).expanduser().resolve()
    events_path = Path(args.events_csv).expanduser().resolve()
    reinit_path = Path(args.reinit_csv).expanduser().resolve()

    for p in (out_path, diag_path, states_path, events_path, reinit_path):
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}; use --overwrite")

    print("=" * 110)
    print("PROGRESSIVE-BASELINE E-BAND CML WVD — raw TX/RX — no external humidity calibration")
    print("=" * 110)
    print("CML        :", cml_path)
    print("Gauge      :", gauge_path)
    print("Temperature:", temp_path)
    print("Output     :", out_path)
    print(
        "State rule : persistent median levels; "
        f"large TL>={args.state_large_tl_step_db:g} dB, "
        f"TX>={args.state_tx_step_db:g} dB with TL>={args.state_min_tl_with_tx_db:g} dB, "
        f"RX/TL>={args.state_rx_tl_step_db:g} dB, "
        f"gap>={args.state_gap_hours:g} h with TL>={args.state_gap_tl_step_db:g} dB"
    )
    print(
        "Baseline   : causal daily lower-envelope; "
        f"warm-up={args.progressive_warmup_days:g} d, "
        f"window={args.progressive_window_days:g} d, "
        f"Qday={args.progressive_daily_low_quantile:g}, "
        f"Qhist={args.progressive_history_quantile:g}, "
        f"up<= {args.progressive_max_up_db_per_day:g} dB/d, "
        f"down<= {args.progressive_max_down_db_per_day:g} dB/d"
    )
    print(
        "Reinit     : persistent unsupported floor; "
        f">={args.progressive_reinit_threshold_db:g} dB, "
        f"~{args.progressive_reinit_persistence_days:g} d, "
        f"min_days={args.progressive_reinit_min_valid_days}, "
        f"support>={args.progressive_reinit_min_support_fraction:g}, "
        f"floor_step>={args.progressive_reinit_min_floor_step_db:g} dB, "
        f"Qrecent={args.progressive_reinit_recent_quantile:g}, "
        f"IQR<={args.progressive_reinit_max_iqr_db:g} dB, "
        f"enabled={bool(args.progressive_reinit_enabled)}"
    )
    print("Humidity   : NEVER READ")
    print()

    cml = xr.open_dataset(cml_path, cache=False)
    gauge = xr.open_dataset(gauge_path, cache=False)
    temp = xr.open_dataset(temp_path, cache=False)

    try:
        cml_times = clean_time_index(cml["time"].values)
        if not cml_times.is_monotonic_increasing:
            raise ValueError("CML time must be monotonic")

        cml_ids = helper.decode_strings(cml["cml_id"].values)
        sub_ids = helper.decode_strings(cml["sublink_id"].values)
        freq = helper.to_ghz(cml["frequency"].values)
        length = helper.to_km(cml["length"].values)

        if freq.shape == (len(cml_ids), len(sub_ids)):
            freq = freq.T
        if freq.shape != (len(sub_ids), len(cml_ids)):
            raise ValueError(f"Unexpected frequency shape {freq.shape}")

        lat_a = np.asarray(cml["site_a_latitude"].values, dtype=float).reshape(-1)
        lon_a = np.asarray(cml["site_a_longitude"].values, dtype=float).reshape(-1)
        lat_b = np.asarray(cml["site_b_latitude"].values, dtype=float).reshape(-1)
        lon_b = np.asarray(cml["site_b_longitude"].values, dtype=float).reshape(-1)

        if "polarization" in cml:
            pol = np.asarray(cml["polarization"].values)
            if pol.shape == (len(cml_ids), len(sub_ids)):
                pol = pol.T
            pol2 = np.empty(pol.shape, dtype=object)
            for j in range(pol.shape[0]):
                pol2[j, :] = helper.decode_strings(pol[j, :])
            pol = pol2
        else:
            pol = np.full(freq.shape, "", dtype=object)

        selected = []
        for ci, cid in enumerate(cml_ids):
            L = float(length[ci])
            if not (
                np.isfinite(L)
                and args.length_min_km <= L <= args.length_max_km
            ):
                continue
            for si, sid in enumerate(sub_ids):
                f = float(freq[si, ci])
                if not np.isfinite(f):
                    continue
                if args.include_gap_76_81:
                    f_ok = args.freq_min_ghz <= f <= args.freq_max_ghz
                else:
                    f_ok = (
                        args.band1_min_ghz <= f <= args.band1_max_ghz
                        or args.band2_min_ghz <= f <= args.band2_max_ghz
                    )
                if not f_ok:
                    continue
                selected.append(
                    helper.SeriesMeta(
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

        wanted_cml = set(str(x) for x in args.cml_id)
        wanted_sub = set(str(x) for x in args.sublink_id)
        if wanted_cml:
            selected = [x for x in selected if str(x.cml_id) in wanted_cml]
        if wanted_sub:
            selected = [x for x in selected if str(x.sublink_id) in wanted_sub]

        if args.max_series is not None:
            selected = selected[: int(args.max_series)]
        if not selected:
            raise RuntimeError("No E-band directional sublinks selected")

        print(f"Selected directional sublinks: {len(selected):,}")
        print(f"Selected physical CMLs        : {len(set(s.cml_id for s in selected)):,}")
        print("Population                     : 71-76 OR 81-86 GHz, 2-10 km")
        if wanted_cml:
            print("Target CML IDs                 :", ", ".join(sorted(wanted_cml)))
        if wanted_sub:
            print("Target sublink IDs             :", ", ".join(sorted(wanted_sub)))

        # Gauge metadata.
        gauge_var = helper.detect_name(
            gauge, ("SRA10M", "rr"), "precipitation accumulation variable"
        )
        gauge_da = gauge[gauge_var]
        gauge_time_dim = next(d for d in gauge_da.dims if "time" in d.lower())
        gauge_station_dim = next(d for d in gauge_da.dims if d != gauge_time_dim)
        gauge_lat_name = helper.detect_name(gauge, ("lat", "latitude"), "gauge latitude")
        gauge_lon_name = helper.detect_name(gauge, ("lon", "longitude"), "gauge longitude")
        gauge_lats = np.asarray(gauge[gauge_lat_name].values, dtype=float).reshape(-1)
        gauge_lons = np.asarray(gauge[gauge_lon_name].values, dtype=float).reshape(-1)
        gauge_times = clean_time_index(gauge[gauge_time_dim].values)
        if gauge_station_dim in gauge.coords:
            gauge_ids = helper.decode_strings(gauge[gauge_station_dim].values)
        else:
            gn = helper.detect_name(
                gauge, ("station", "WSI", "station_id", "GH_ID"), "gauge station ID"
            )
            gauge_ids = helper.decode_strings(gauge[gn].values)

        # Temperature metadata.
        temp_var = helper.detect_name(temp, ("T",), "standard air temperature variable")
        temp_da = temp[temp_var]
        temp_time_dim = next(d for d in temp_da.dims if "time" in d.lower())
        temp_station_dim = next(d for d in temp_da.dims if d != temp_time_dim)
        temp_lat_name = helper.detect_name(temp, ("lat", "latitude"), "temperature latitude")
        temp_lon_name = helper.detect_name(temp, ("lon", "longitude"), "temperature longitude")
        temp_lats = np.asarray(temp[temp_lat_name].values, dtype=float).reshape(-1)
        temp_lons = np.asarray(temp[temp_lon_name].values, dtype=float).reshape(-1)
        temp_times = clean_time_index(temp[temp_time_dim].values)
        if temp_station_dim in temp.coords:
            temp_ids = helper.decode_strings(temp[temp_station_dim].values)
        else:
            tn = helper.detect_name(
                temp, ("station", "WSI", "station_id", "GH_ID"), "temperature station ID"
            )
            temp_ids = helper.decode_strings(temp[tn].values)

        eval_time = pd.date_range(
            EVAL_START,
            EVAL_END_EXCLUSIVE - pd.Timedelta(minutes=5),
            freq="5min",
        )
        cal_2021 = (cml_times >= STATE_START) & (cml_times < STATE_END)
        wet_pad_minutes = int(round(args.wet_buffer_hours * 60.0))

        # Use existing writer, then correct/extend it.
        nc = helper.create_output(out_path, eval_time, selected, args)

        # Always overwrite time using netCDF date2num; never raw datetime integers.
        nc["time"][:] = date2num(
            eval_time.to_pydatetime(),
            units=nc["time"].units,
            calendar=nc["time"].calendar,
        )

        # Hardware-reset and progressive-baseline variables.
        chunks_st = (1, min(2016, len(eval_time)))
        chunks_m = (1, 1, min(2016, len(eval_time)))

        v = nc.createVariable(
            "hardware_state_id",
            "i2", ("series", "time"),
            fill_value=np.int16(-32768),
            zlib=True, complevel=2, shuffle=True,
            chunksizes=chunks_st,
        )
        v.long_name = "Conservative CML-only hardware state identifier"
        v.comment = "-1 means no stable state assigned to this 5-min bin"

        v = nc.createVariable(
            "hardware_transition_buffer_5min",
            "u1", ("series", "time"),
            fill_value=np.uint8(255),
            zlib=True, complevel=2, shuffle=True,
            chunksizes=chunks_st,
        )
        v.flag_values = np.asarray([0, 1], dtype=np.uint8)
        v.flag_meanings = "outside_transition_buffer inside_transition_buffer"

        v = nc.createVariable(
            "baseline_state_5min_db",
            "f4", ("method", "series", "time"),
            fill_value=np.nan,
            zlib=True, complevel=2, shuffle=True,
            chunksizes=chunks_m,
        )
        v.units = "dB"
        v.long_name = "Compatibility alias of causal progressive CML-only baseline actually used for retrieval"

        v = nc.createVariable(
            "baseline_progressive_5min_db",
            "f4", ("method", "series", "time"),
            fill_value=np.nan,
            zlib=True, complevel=2, shuffle=True,
            chunksizes=chunks_m,
        )
        v.units = "dB"
        v.long_name = "Causal progressive CML-only internal baseline B(t) actually used for retrieval"
        v.comment = "Identical for GRAF and gauge; computed only from completed past CML days within each confirmed telecom regime."

        v = nc.createVariable(
            "baseline_reinit_5min",
            "u1", ("series", "time"),
            fill_value=np.uint8(255),
            zlib=True, complevel=2, shuffle=True,
            chunksizes=chunks_st,
        )
        v.flag_values = np.asarray([0, 1], dtype=np.uint8)
        v.flag_meanings = "no_reinitialization baseline_reinitialized_at_start_of_day"
        v.long_name = "V6 persistent unsupported-floor baseline reinitialization flag"

        nc.createVariable(
            "baseline_n_reinitializations", "i2", ("series",), fill_value=np.int16(-1)
        )

        nc.createVariable(
            "hardware_n_states", "i2", ("series",), fill_value=np.int16(-1)
        )
        nc.createVariable(
            "hardware_n_change_points", "i2", ("series",), fill_value=np.int16(-1)
        )
        vm = nc.createVariable(
            "hardware_max_abs_step_db", "f4", ("series",), fill_value=np.nan
        )
        vm.units = "dB"

        # Expand status semantics.
        nc["retrieval_status"].flag_values = np.asarray(
            [0, 1, 2, 3, 4, 5, 6, 7], dtype=np.uint8
        )
        nc["retrieval_status"].flag_meanings = (
            "unavailable valid below_p676_range above_p676_range "
            "nonmonotonic_p676_lookup temperature_outside_lookup "
            "no_valid_progressive_baseline hardware_transition_buffer"
        )

        nc.title = "Czech E-band CML V6 causal progressive-baseline water-vapour retrieval from raw TX/RX"
        nc.frequency_selection_ghz = "71-76 and 81-86 GHz"
        nc.calibration_period = (
            "Causal 2021 CML-only baseline history. Each current day uses only completed "
            "past days in the same confirmed telecom regime; no RH/WVD reference used."
        )
        nc.baseline_definition = (
            "V6 B(t): daily Q05 of TL-L*gamma_P676(rho=0), trailing 45-day Q10 of completed "
            "past daily lows by default; 14-day warm-up after strong telecom reset; ordinary "
            "upward drift disabled, but a separate persistent unsupported-floor rule may "
            "causally reinitialize B upward; same B(t) for GRAF and gauge"
        )
        nc.hardware_state_detection = (
            "CML-only persistent daily-median TX/RX/TL level changes; "
            "Q05 retained only as diagnostic and never sufficient by itself"
        )
        nc.hardware_state_window_days = int(args.state_window_days)
        nc.hardware_transition_buffer_hours = float(args.hardware_transition_buffer_hours)
        nc.progressive_baseline_window_days = float(args.progressive_window_days)
        nc.progressive_baseline_warmup_days = float(args.progressive_warmup_days)
        nc.progressive_daily_low_quantile = float(args.progressive_daily_low_quantile)
        nc.progressive_history_quantile = float(args.progressive_history_quantile)
        nc.progressive_max_down_db_per_day = float(args.progressive_max_down_db_per_day)
        nc.progressive_max_up_db_per_day = float(args.progressive_max_up_db_per_day)
        nc.progressive_reinit_enabled = str(bool(args.progressive_reinit_enabled)).lower()
        nc.progressive_reinit_threshold_db = float(args.progressive_reinit_threshold_db)
        nc.progressive_reinit_persistence_days = float(args.progressive_reinit_persistence_days)
        nc.progressive_reinit_min_valid_days = int(args.progressive_reinit_min_valid_days)
        nc.progressive_reinit_min_support_fraction = float(args.progressive_reinit_min_support_fraction)
        nc.progressive_reinit_min_floor_step_db = float(args.progressive_reinit_min_floor_step_db)
        nc.progressive_reinit_min_previous_days = int(args.progressive_reinit_min_previous_days)
        nc.progressive_reinit_recent_quantile = float(args.progressive_reinit_recent_quantile)
        nc.progressive_reinit_max_iqr_db = float(args.progressive_reinit_max_iqr_db)
        nc.progressive_reinit_cooldown_days = float(args.progressive_reinit_cooldown_days)
        nc.progressive_baseline_causal = "yes; current-day baseline uses completed past days only"
        nc.humidity_used_for_state_detection = "no"
        nc.humidity_used_for_calibration = "no"
        nc.humidity_used_for_retrieval = "no"
        nc.state_aware_warning = (
            "Confirmed telecom states are reset regimes only. After each reset, retrieval remains "
            "missing until the causal progressive baseline has enough post-reset history. "
            "Baseline reinitializations are logged separately and are not labelled hardware events."
        )

        gauge_mask_cache: Dict[int, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        last_temp_station_index: Optional[int] = None
        last_temp_1min_c: Optional[np.ndarray] = None

        diag_rows = []
        state_rows = []
        event_rows = []
        reinit_rows = []
        run_t0 = time.time()

        for out_i, s in enumerate(selected):
            t_start = time.time()

            tx = np.asarray(
                cml["tx"].isel(
                    sublink_id=s.sublink_index, cml_id=s.cml_index
                ).values,
                dtype=np.float64,
            ).reshape(-1)
            rx = np.asarray(
                cml["rx"].isel(
                    sublink_id=s.sublink_index, cml_id=s.cml_index
                ).values,
                dtype=np.float64,
            ).reshape(-1)

            valid_signal = (
                np.isfinite(tx)
                & np.isfinite(rx)
                & (tx > helper.TX_MIN_VALID_EXCLUSIVE)
                & (tx != helper.TX_FILL_CODE)
                & (rx > helper.RX_MIN_VALID_EXCLUSIVE)
                & (rx <= helper.RX_MAX_VALID_INCLUSIVE)
            )
            tl = np.full(tx.shape, np.nan, dtype=np.float64)
            tl[valid_signal] = tx[valid_signal] - rx[valid_signal]

            # ---- Hardware states: RAW CML ONLY.
            events, daily_state_diag = detect_hardware_events(
                cml_times, tx, rx, valid_signal, args
            )
            state_min, transition_min, boundaries = assign_state_minutes(
                cml_times, events, args.hardware_transition_buffer_hours
            )
            intervals = state_intervals(boundaries)

            max_abs_step = (
                max(abs(float(e["delta_tl_med"])) for e in events)
                if events else 0.0
            )
            nc["hardware_n_states"][out_i] = np.int16(len(intervals))
            nc["hardware_n_change_points"][out_i] = np.int16(len(events))
            nc["hardware_max_abs_step_db"][out_i] = np.float32(max_abs_step)

            state5, transition5 = aggregate_state_to_5min(
                cml_times, state_min, transition_min
            )
            nc["hardware_state_id"][out_i, :] = state5
            nc["hardware_transition_buffer_5min"][out_i, :] = transition5.astype(np.uint8)

            for ei, e in enumerate(events):
                event_rows.append(
                    {
                        "series_index": out_i,
                        "cml_id": s.cml_id,
                        "sublink_id": s.sublink_id,
                        "frequency_ghz": s.frequency_ghz,
                        "length_km": s.length_km,
                        "event_index": ei,
                        **e,
                    }
                )

            # ---- Station assignments.
            gj, gd = helper.nearest_station_to_path(
                gauge_lats, gauge_lons,
                s.lat_a, s.lon_a, s.lat_b, s.lon_b,
            )
            tj, td = helper.nearest_station_to_path(
                temp_lats, temp_lons,
                s.lat_a, s.lon_a, s.lat_b, s.lon_b,
            )
            gauge_ok = (
                gj >= 0 and np.isfinite(gd)
                and gd <= args.max_gauge_distance_km
            )
            temp_ok = (
                tj >= 0 and np.isfinite(td)
                and td <= args.max_temperature_distance_km
            )

            if temp_ok:
                if last_temp_station_index == tj and last_temp_1min_c is not None:
                    temp_1min_c = last_temp_1min_c
                else:
                    native_temp = np.asarray(
                        temp_da.isel({temp_station_dim: tj}).values,
                        dtype=np.float64,
                    ).reshape(-1)
                    temp_1min_c = interpolate_temperature_to_minutes_safe(
                        temp_times, native_temp, cml_times,
                        args.temperature_max_gap_min,
                    )
                    last_temp_station_index = tj
                    last_temp_1min_c = temp_1min_c
                t_ref_n = int(np.isfinite(temp_1min_c).sum())
                t_ref_c = float(np.nanmedian(temp_1min_c)) if t_ref_n else np.nan
            else:
                temp_1min_c = np.full(len(cml_times), np.nan, dtype=float)
                t_ref_n = 0
                t_ref_c = np.nan

            nc["gauge_station_id"][out_i] = str(gauge_ids[gj]) if gj >= 0 else ""
            nc["gauge_distance_km"][out_i] = np.float32(gd)
            nc["temperature_station_id"][out_i] = str(temp_ids[tj]) if tj >= 0 else ""
            nc["temperature_distance_km"][out_i] = np.float32(td)
            nc["temperature_ref_c"][out_i] = np.float32(t_ref_c)
            nc["temperature_ref_n"][out_i] = np.int32(t_ref_n)

            t5_time, temp_5min_c = helper.aggregate_temperature_5min(
                cml_times,
                temp_1min_c,
                EVAL_START,
                EVAL_END_EXCLUSIVE,
                min_valid_minutes=args.min_temperature_minutes_5min,
            )
            if not clean_time_index(t5_time).equals(clean_time_index(eval_time)):
                raise RuntimeError("5-min temperature time mismatch")
            nc["temperature_5min_c"][out_i, :] = temp_5min_c.astype(np.float32)

            if np.isfinite(temp_1min_c).any():
                (
                    temp_grid,
                    rho_grid,
                    gamma_table,
                    monotonic,
                    sensitivity_by_temp,
                ) = helper.p676_temperature_lookup(
                    s.frequency_ghz,
                    args.pressure_hpa,
                    args.temperature_grid_min_c,
                    args.temperature_grid_max_c,
                    args.temperature_grid_step_c,
                    args.rho_max_gm3,
                    args.rho_step_gm3,
                )
                gamma0_t = helper.gamma_rho0_for_temperature(
                    temp_1min_c, temp_grid, gamma_table
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
                    dtype=float,
                )
                rho_grid = np.asarray([0.0, args.rho_max_gm3], dtype=float)
                gamma_table = np.full((2, 2), np.nan, dtype=np.float32)
                monotonic = False
                sensitivity = np.nan
                gamma0_t = np.full(len(cml_times), np.nan, dtype=float)

            nc["p676_monotonic"][out_i] = np.uint8(1 if monotonic else 0)
            nc["p676_sensitivity_dbkm_per_gm3"][out_i] = np.float32(sensitivity)

            # ---- Wet/dry masks.
            (
                graf_wet,
                graf_known,
                graf_buffered,
                graf_q80,
                graf_thr,
            ) = helper.graf_masks(
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

            if gauge_ok:
                if gj not in gauge_mask_cache:
                    gauge_amount = np.asarray(
                        gauge_da.isel({gauge_station_dim: gj}).values,
                        dtype=np.float64,
                    ).reshape(-1)
                    gauge_mask_cache[gj] = helper.gauge_masks_on_cml_minutes(
                        cml_times, gauge_times, gauge_amount, wet_pad_minutes
                    )
                gauge_wet, gauge_known, gauge_buffered = gauge_mask_cache[gj]
            else:
                gauge_wet = np.zeros(len(tl), dtype=bool)
                gauge_known = np.zeros(len(tl), dtype=bool)
                gauge_buffered = np.zeros(len(tl), dtype=bool)

            method_masks = {
                "graf": (graf_known, graf_buffered, graf_wet),
                "gauge": (gauge_known, gauge_buffered, gauge_wet),
            }

            b_sample = tl - (s.length_km * gamma0_t)

            diag = {
                "series_index": out_i,
                "cml_id": s.cml_id,
                "sublink_id": s.sublink_id,
                "frequency_ghz": s.frequency_ghz,
                "length_km": s.length_km,
                "raw_finite_tl": int(np.isfinite(tl).sum()),
                "gauge_station_id": str(gauge_ids[gj]) if gj >= 0 else "",
                "gauge_distance_km": gd,
                "gauge_within_max_distance": bool(gauge_ok),
                "temperature_station_id": str(temp_ids[tj]) if tj >= 0 else "",
                "temperature_distance_km": td,
                "temperature_ref_c": t_ref_c,
                "p676_monotonic": bool(monotonic),
                "p676_sensitivity_dbkm_per_gm3": sensitivity,
                "graf_q80_db": graf_q80,
                "graf_threshold_db": graf_thr,
                "hardware_n_states": len(intervals),
                "hardware_n_change_points": len(events),
                "hardware_max_abs_step_db": max_abs_step,
                "hardware_transition_minutes": int(transition_min.sum()),
            }

            # ------------------------------------------------------------------
            # ONE common CAUSAL PROGRESSIVE baseline, independent of GRAF/gauge.
            #
            # Strong CML-only hardware events define reset regimes. Within each
            # regime the baseline is updated chronologically from completed PAST
            # daily lower-envelope statistics only. No humidity or wet/dry branch
            # enters B(t).
            # ------------------------------------------------------------------
            (
                baseline_min,
                progressive_rows,
                progressive_reinit_rows,
                reinit_flag_min,
            ) = progressive_causal_baseline(
                cml_times,
                b_sample,
                state_min,
                transition_min,
                intervals,
                args,
            )

            for rr in progressive_reinit_rows:
                reinit_rows.append(
                    {
                        "series_index": out_i,
                        "cml_id": s.cml_id,
                        "sublink_id": s.sublink_id,
                        "frequency_ghz": s.frequency_ghz,
                        "length_km": s.length_km,
                        **rr,
                    }
                )

            for r in progressive_rows:
                state_rows.append(
                    {
                        "series_index": out_i,
                        "baseline_source": "progressive_causal_cml_lower_envelope_v6",
                        "cml_id": s.cml_id,
                        "sublink_id": s.sublink_id,
                        "frequency_ghz": s.frequency_ghz,
                        "length_km": s.length_km,
                        "progressive_window_days": args.progressive_window_days,
                        "progressive_warmup_days": args.progressive_warmup_days,
                        "progressive_daily_low_quantile": args.progressive_daily_low_quantile,
                        "progressive_history_quantile": args.progressive_history_quantile,
                        "progressive_max_down_db_per_day": args.progressive_max_down_db_per_day,
                        "progressive_max_up_db_per_day": args.progressive_max_up_db_per_day,
                        "progressive_reinit_enabled": bool(args.progressive_reinit_enabled),
                        "progressive_reinit_threshold_db": args.progressive_reinit_threshold_db,
                        "progressive_reinit_persistence_days": args.progressive_reinit_persistence_days,
                        "progressive_reinit_min_support_fraction": args.progressive_reinit_min_support_fraction,
            "progressive_reinit_min_floor_step_db": args.progressive_reinit_min_floor_step_db,
            "progressive_reinit_min_previous_days": args.progressive_reinit_min_previous_days,
                        "progressive_reinit_min_floor_step_db": args.progressive_reinit_min_floor_step_db,
                        "progressive_reinit_min_previous_days": args.progressive_reinit_min_previous_days,
                        "progressive_reinit_recent_quantile": args.progressive_reinit_recent_quantile,
                        "progressive_reinit_max_iqr_db": args.progressive_reinit_max_iqr_db,
                        **r,
                    }
                )

            B5_common = aggregate_baseline_to_5min(cml_times, baseline_min)
            # A reinitialization is applied at the start of a UTC day. Mark all
            # 5-min bins on that day so downstream diagnostics can identify it.
            full_eval_min = pd.date_range(
                EVAL_START,
                EVAL_END_EXCLUSIVE - pd.Timedelta(minutes=1),
                freq="1min",
            )
            re_ser = pd.Series(
                np.asarray(reinit_flag_min, dtype=np.uint8), index=cml_times
            ).reindex(full_eval_min, fill_value=0)
            re_a = re_ser.to_numpy(dtype=np.uint8).reshape(-1, 5)
            reinit5 = (re_a > 0).any(axis=1).astype(np.uint8)
            nc["baseline_reinit_5min"][out_i, :] = reinit5
            nc["baseline_n_reinitializations"][out_i] = np.int16(
                len(progressive_reinit_rows)
            )

            finite_B = baseline_min[np.isfinite(baseline_min)]
            Brep = float(np.median(finite_B)) if finite_B.size else np.nan
            total_baseline_n = int(
                sum(int(r["baseline_n_candidate_minutes"]) for r in progressive_rows)
            )
            n_valid_regimes = int(sum(int(r["baseline_valid"]) for r in progressive_rows))
            first_valid_candidates = [
                pd.Timestamp(r["baseline_first_valid_time"])
                for r in progressive_rows
                if not pd.isna(r["baseline_first_valid_time"])
            ]
            first_progressive_B = (
                min(first_valid_candidates) if first_valid_candidates else pd.NaT
            )

            for mi, method in enumerate(METHODS):
                known, wet_buffered, wet_raw = method_masks[method]

                # Wet/dry only controls which retrieval samples are accepted.
                # It does NOT change the progressive CML-only baseline.
                safe = (
                    np.isfinite(tl)
                    & known
                    & ~wet_buffered
                    & ~transition_min
                    & (state_min >= 0)
                )

                B5 = B5_common
                nc["baseline_state_5min_db"][mi, out_i, :] = B5.astype(np.float32)
                nc["baseline_progressive_5min_db"][mi, out_i, :] = B5.astype(np.float32)

                # Compatibility variables: identical across methods by design.
                nc["baseline_B_db"][mi, out_i] = np.float32(Brep)
                nc["baseline_n"][mi, out_i] = np.int64(total_baseline_n)
                nc["baseline_pre_db"][mi, out_i] = np.float32(np.nan)
                nc["baseline_post_db"][mi, out_i] = np.float32(np.nan)
                nc["baseline_delta_post_minus_pre_db"][mi, out_i] = np.float32(np.nan)
                nc["baseline_pre_n"][mi, out_i] = np.int64(0)
                nc["baseline_post_n"][mi, out_i] = np.int64(0)

                t5, tl5, valid5, safe5 = helper.aggregate_5min_complete(
                    cml_times,
                    tl,
                    safe,
                    EVAL_START,
                    EVAL_END_EXCLUSIVE,
                    min_valid_minutes=args.min_valid_minutes_5min,
                )
                if not clean_time_index(t5).equals(clean_time_index(eval_time)):
                    raise RuntimeError("Internal 5-min output time mismatch")

                nc["safe_dry_5min"][mi, out_i, :] = safe5
                nc["total_loss_5min"][mi, out_i, :] = tl5.astype(np.float32)

                if monotonic:
                    gamma_obs = (tl5 - B5) / s.length_km
                    wvd, status = helper.invert_gamma_lookup_temperature(
                        gamma_obs,
                        temp_5min_c,
                        temp_grid,
                        rho_grid,
                        gamma_table,
                        monotonic,
                    )

                    # Distinguish unavailable progressive baseline from generic unavailable.
                    no_B = (
                        np.isfinite(tl5)
                        & np.isfinite(temp_5min_c)
                        & ~np.isfinite(B5)
                    )
                    status[no_B] = 6

                    # Transition status, even though TL is deliberately
                    # withheld by the safe-dry aggregator.
                    status[
                        transition5 & np.isfinite(temp_5min_c)
                    ] = 7
                else:
                    wvd = np.full(len(eval_time), np.nan, dtype=float)
                    status = np.zeros(len(eval_time), dtype=np.uint8)
                    q = np.isfinite(tl5) & np.isfinite(temp_5min_c)
                    status[q] = 4

                nc["wvd_internal_anchor"][mi, out_i, :] = wvd.astype(np.float32)
                nc["retrieval_status"][mi, out_i, :] = status

                n_eval_states = int(
                    len(set(int(x) for x in state5[state5 >= 0]))
                )
                diag[f"{method}_wet_raw_minutes"] = int(wet_raw.sum())
                diag[f"{method}_wet_buffered_minutes"] = int(wet_buffered.sum())
                diag[f"{method}_safe_minutes"] = int(safe.sum())
                # Compatibility name retained; value now means reset regimes
                # with at least one valid progressive baseline day.
                diag[f"{method}_valid_state_baselines"] = n_valid_regimes
                diag[f"{method}_valid_progressive_regimes"] = n_valid_regimes
                diag[f"{method}_eval_states"] = n_eval_states
                diag[f"{method}_representative_B_db"] = Brep
                diag[f"{method}_baseline_source"] = "progressive_causal_cml_lower_envelope_v6"
                diag[f"{method}_progressive_first_valid_time"] = first_progressive_B
                diag[f"{method}_progressive_valid_5min"] = int(np.isfinite(B5).sum())
                diag[f"{method}_progressive_n_reinitializations"] = int(
                    len(progressive_reinit_rows)
                )
                diag[f"{method}_eval_safe_5min"] = int(valid5.sum())
                diag[f"{method}_eval_valid_wvd_5min"] = int(np.isfinite(wvd).sum())
                diag[f"{method}_eval_no_baseline_5min"] = int(
                    np.sum(status == 6)
                )
                diag[f"{method}_eval_transition_5min"] = int(
                    np.sum(status == 7)
                )

            diag_rows.append(diag)
            nc.sync()

            if (
                out_i < 10
                or (out_i + 1) % 25 == 0
                or out_i + 1 == len(selected)
            ):
                elapsed = time.time() - t_start
                print(
                    f"[{out_i+1:4d}/{len(selected):4d}] "
                    f"{s.cml_id} sub={s.sublink_id} "
                    f"f={s.frequency_ghz:.3f} GHz L={s.length_km:.2f} km "
                    f"reset_regimes={len(intervals)} changes={len(events)} "
                    f"reinit={len(progressive_reinit_rows)} "
                    f"B5valid={int(np.isfinite(B5_common).sum())} "
                    f"max|dTL|={max_abs_step:.2f} dB "
                    f"({elapsed:.1f}s)"
                )

        # Finish files.
        nc.history = (
            f"Created {pd.Timestamp.now("UTC").isoformat()} by "
            f"{Path(__file__).name}"
        )
        nc.close()

        pd.DataFrame(diag_rows).to_csv(diag_path, index=False)
        pd.DataFrame(state_rows).to_csv(states_path, index=False)
        event_columns = [
            "series_index", "cml_id", "sublink_id", "frequency_ghz", "length_km",
            "event_index", "boundary_time", "reason", "near_long_gap", "origin",
            "pre_tl_med", "post_tl_med", "delta_tl_med",
            "pre_tx_med", "post_tx_med", "delta_tx_med",
            "pre_rx_med", "post_rx_med", "delta_rx_med",
            "pre_tl_low", "post_tl_low", "delta_tl_low",
            "n_pre_days", "n_post_days",
            "merged_candidate_count", "all_reasons",
        ]
        pd.DataFrame(event_rows, columns=event_columns).to_csv(events_path, index=False)

        reinit_columns = [
            "series_index", "cml_id", "sublink_id", "frequency_ghz", "length_km",
            "state_id", "reinit_index_within_state", "reinit_time", "reason",
            "old_B_db", "new_B_db", "delta_B_db", "recent_floor_db",
            "previous_floor_db", "recent_floor_step_db", "recent_iqr_db",
            "previous_n_valid_days", "recent_n_valid_days", "recent_support_fraction", "recent_calendar_span_days",
            "recent_quantile", "threshold_db", "persistence_days",
            "history_restart_day", "hardware_event",
        ]
        pd.DataFrame(reinit_rows, columns=reinit_columns).to_csv(reinit_path, index=False)

        summary_path = out_path.with_suffix(".summary.json")
        summary = {
            "n_series": len(selected),
            "n_physical_cml": len(set(s.cml_id for s in selected)),
            "n_hardware_events": len(event_rows),
            "n_baseline_reinitializations": len(reinit_rows),
            "n_state_rows_common": len(state_rows),
            "frequency_population": "71-76 and 81-86 GHz",
            "length_population_km": [args.length_min_km, args.length_max_km],
            "pressure_hpa": args.pressure_hpa,
            "baseline_mode": "progressive_causal_cml_lower_envelope_v6",
            "progressive_window_days": args.progressive_window_days,
            "progressive_warmup_days": args.progressive_warmup_days,
            "progressive_daily_low_quantile": args.progressive_daily_low_quantile,
            "progressive_history_quantile": args.progressive_history_quantile,
            "progressive_min_valid_minutes_per_day": args.progressive_min_valid_minutes_per_day,
            "progressive_min_history_days": args.progressive_min_history_days,
            "progressive_min_history_minutes": args.progressive_min_history_minutes,
            "progressive_max_down_db_per_day": args.progressive_max_down_db_per_day,
            "progressive_max_up_db_per_day": args.progressive_max_up_db_per_day,
            "progressive_reinit_enabled": bool(args.progressive_reinit_enabled),
            "progressive_reinit_threshold_db": args.progressive_reinit_threshold_db,
            "progressive_reinit_persistence_days": args.progressive_reinit_persistence_days,
            "progressive_reinit_min_valid_days": args.progressive_reinit_min_valid_days,
            "progressive_reinit_min_support_fraction": args.progressive_reinit_min_support_fraction,
            "progressive_reinit_recent_quantile": args.progressive_reinit_recent_quantile,
            "progressive_reinit_max_iqr_db": args.progressive_reinit_max_iqr_db,
            "progressive_reinit_cooldown_days": args.progressive_reinit_cooldown_days,
            "hardware_transition_buffer_hours": args.hardware_transition_buffer_hours,
            "humidity_used": False,
            "elapsed_seconds": time.time() - run_t0,
        }
        summary_path.write_text(json.dumps(summary, indent=2))

        print("\n" + "=" * 110)
        print("DONE")
        print("=" * 110)
        print("NetCDF     :", out_path)
        print("Diagnostics:", diag_path)
        print("States     :", states_path)
        print("Events     :", events_path)
        print("Reinit     :", reinit_path)
        print("Summary    :", summary_path)
        print(f"Elapsed    : {time.time() - run_t0:.1f} s")
        print()
        print("IMPORTANT:")
        print("  * Strong telecom reset events were detected from RAW CML TX/RX/TL only.")
        print("  * No RH/WVD validation data were read.")
        print("  * B(t) is causal: each day uses completed PAST CML days only.")
        print("  * Persistent unsupported-floor reinitializations are CML-only and logged separately from hardware events.")
        print("  * GRAF/gauge use the SAME progressive B(t); only their wet/dry masks differ.")
        print("  * After a hardware reset, WVD remains missing until the progressive warm-up/history criteria are met.")

    finally:
        # nc may already be closed.
        try:
            if "nc" in locals() and nc.isopen():
                nc.close()
        except Exception:
            pass
        cml.close()
        gauge.close()
        temp.close()


if __name__ == "__main__":
    main()
