#!/usr/bin/env python3
"""
cml_qc_filter_standalone_v3.py
================================
Standalone, streamed quality control for commercial microwave-link NetCDF data.

The program preserves the input dimensions, coordinates, metadata, and all
non-signal variables. The output TX/RX variables retain their original names
but contain the quality-controlled signal series. The original input file is
never modified.

QC steps
--------
1. Mask missing/fill/physically invalid TX and RX values.
2. Mask low, stable RSL plateaus and two neighbouring samples on each side.
3. Fill blackout gaps of at most 60 samples when adjacent RSL is below -65 dBm.
4. Linearly interpolate total-loss gaps of at most 5 samples.
5. Flag sublinks with excessive 5-hour or 1-hour total-loss fluctuations.
6. Flag sublinks whose legacy ITU detection limit is >= 2 mm h-1.

The output contains separate flags for every QC step, a sample-level final-valid
flag, final sublink/link flags, and the QC total loss. A compact bit-mask can be
added with --write-bitmask when needed.
Steps 5 and 6 are sublink-level flags; the final link flag passes when at least
one directional sublink passes, matching the original link aggregation logic.

Processing is streamed in CML batches. Expensive per-sublink calculations are
compiled and parallelized with Numba; NetCDF reads/writes remain in the main
process for HDF5 safety and efficient sequential I/O.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    from numba import njit, prange, set_num_threads, get_num_threads
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "This standalone optimized script requires numba. Install it with "
        "`python -m pip install numba` in the active environment."
    ) from exc

try:
    from netCDF4 import Dataset, chartostring
except Exception as exc:  # pragma: no cover
    Dataset = None
    _NETCDF4_IMPORT_ERROR = exc
else:
    _NETCDF4_IMPORT_ERROR = None


# -----------------------------------------------------------------------------
# QC constants: kept explicit so the output records the exact configuration.
# -----------------------------------------------------------------------------
TX_MIN_DBM = -99.0
TX_SENTINEL = 255.0
RX_MIN_DBM = -99.9
RX_MAX_DBM = -20.0
HUGE_FILL_ABS = 1.0e35

PLATEAU_WINDOW = 3
PLATEAU_MAX_DBM = -85.0
PLATEAU_STD_DB = 0.5
PLATEAU_PAD = 2

BLACKOUT_LOW_DBM = -65.0
BLACKOUT_MAX_GAP = 60
SHORT_GAP_MAX = 5

LONG_STD_WINDOW = 300       # 5 hours for 1-minute data
LONG_STD_THRESHOLD = 2.0
LONG_STD_MAX_FRACTION = 0.10
SHORT_STD_WINDOW = 60       # 1 hour for 1-minute data
SHORT_STD_THRESHOLD = 0.8
SHORT_STD_MAX_FRACTION = 0.33
ROLLING_MIN_PERIODS = 60

A_QUANTIZATION_DB = 0.33
DETECTION_LIMIT_MAX_MMPH = 2.0

# Per-sample compact bitmask. Bit 4/5 represent sublink-level failures broadcast
# over the samples; bit 6 means the sample is not usable in the final QC product.
BIT_STEP1_INVALID = np.uint16(1)
BIT_STEP2_PLATEAU = np.uint16(2)
BIT_STEP3_BLACKOUT_FILLED = np.uint16(4)
BIT_STEP4_INTERPOLATED = np.uint16(8)
BIT_STEP5_FLUCTUATION_FAIL = np.uint16(16)
BIT_STEP6_DETECTION_FAIL = np.uint16(32)
BIT_FINAL_INVALID = np.uint16(64)


# -----------------------------------------------------------------------------
# ITU-R P.838-3 k coefficient, implemented locally to avoid pycomlink dependency.
# For horizontal/vertical terrestrial links at zero elevation, k is k_H or k_V.
# -----------------------------------------------------------------------------
_KH_A = np.array([-5.33980, -0.35351, -0.23789, -0.94158], dtype=np.float64)
_KH_B = np.array([-0.10008, 1.26970, 0.86036, 0.64552], dtype=np.float64)
_KH_C = np.array([1.13098, 0.45400, 0.15354, 0.16817], dtype=np.float64)
_KH_M = -0.18961
_KH_C0 = 0.71147

_KV_A = np.array([-3.80595, -3.44965, -0.39902, 0.50167], dtype=np.float64)
_KV_B = np.array([0.56934, -0.22911, 0.73042, 1.07319], dtype=np.float64)
_KV_C = np.array([0.81061, 0.51059, 0.11899, 0.27195], dtype=np.float64)
_KV_M = -0.16398
_KV_C0 = 0.63297


def _itu_k(f_ghz: float, pol_code: int) -> float:
    """Return ITU-R P.838-3 k_H (pol_code=0) or k_V (pol_code=1)."""
    if not np.isfinite(f_ghz) or f_ghz <= 0:
        return np.nan
    x = math.log10(float(f_ghz))
    if int(pol_code) == 1:
        aa, bb, cc, mm, c0 = _KV_A, _KV_B, _KV_C, _KV_M, _KV_C0
    else:
        aa, bb, cc, mm, c0 = _KH_A, _KH_B, _KH_C, _KH_M, _KH_C0
    y = float(np.sum(aa * np.exp(-((x - bb) / cc) ** 2)) + mm * x + c0)
    return 10.0 ** y


def normalize_frequency_to_ghz(values: np.ndarray) -> Tuple[np.ndarray, str]:
    """Normalize a frequency array to GHz using a conservative magnitude rule."""
    arr = np.asarray(values, dtype=np.float64)
    finite = np.abs(arr[np.isfinite(arr)])
    if finite.size == 0:
        return arr, "unknown"
    med = float(np.nanmedian(finite))
    if med > 1.0e8:
        return arr / 1.0e9, "Hz"
    if med > 1.0e5:
        return arr / 1.0e6, "kHz"
    if med > 1.0e3:
        return arr / 1.0e3, "MHz"
    return arr, "GHz"


def pol_to_HV(pol_raw) -> str:
    """Normalize one polarization label to H or V.

    This intentionally retains the legacy fallback to H for missing, mixed, or
    unknown labels.
    """
    if isinstance(pol_raw, (bytes, np.bytes_)):
        pol_raw = pol_raw.decode("utf-8", errors="ignore")
    s = str(pol_raw).replace("\x00", "").strip().lower()
    if "v" in s and "h" not in s:
        return "V"
    if "h" in s and "v" not in s:
        return "H"
    return "H"


def polarization_codes(values: np.ndarray) -> np.ndarray:
    """Map common polarization encodings to 0=H and 1=V."""
    arr = np.asarray(values)
    out = np.zeros(arr.shape, dtype=np.int8)
    if arr.dtype.kind in "iufb":
        # Common numeric convention: 0/H, 1/V. Any other value defaults to H.
        out[np.asarray(arr, dtype=np.float64) == 1] = 1
        return out
    flat = arr.reshape(-1)
    of = out.reshape(-1)
    for i, raw in enumerate(flat):
        of[i] = 1 if pol_to_HV(raw) == "V" else 0
    return out


def detection_limits_legacy(
    frequency_ghz: np.ndarray,
    polarization_code: np.ndarray,
    length_km: np.ndarray,
) -> np.ndarray:
    """
    Reproduce the legacy detection-limit formula used by the original script:
        detection_limit = 0.33 / (k_ITU * link_length_km)
    """
    f = np.asarray(frequency_ghz, dtype=np.float64)
    p = np.asarray(polarization_code, dtype=np.int8)
    L = np.asarray(length_km, dtype=np.float64)
    f, p, L = np.broadcast_arrays(f, p, L)
    out = np.full(f.shape, np.nan, dtype=np.float32)
    it = np.nditer([f, p, L, out], op_flags=[["readonly"], ["readonly"], ["readonly"], ["writeonly"]])
    for fi, pi, li, oi in it:
        k = _itu_k(float(fi), int(pi))
        sens = k * float(li)
        oi[...] = A_QUANTIZATION_DB / sens if np.isfinite(sens) and sens > 0 else np.nan
    return out


# -----------------------------------------------------------------------------
# Fast Numba kernels. Arrays use canonical shape (series, time), where series
# is flattened (sublink, cml) for the current link batch.
# -----------------------------------------------------------------------------
@njit(cache=True, inline="always")
def _is_bad_tx(v: float) -> bool:
    return (not math.isfinite(v)) or v <= TX_MIN_DBM or abs(v - TX_SENTINEL) < 1.0e-6 or abs(v) >= HUGE_FILL_ABS


@njit(cache=True, inline="always")
def _is_bad_rx(v: float) -> bool:
    return (not math.isfinite(v)) or v <= RX_MIN_DBM or v > RX_MAX_DBM or abs(v) >= HUGE_FILL_ABS


@njit(cache=True)
def _rolling_std_exceed_fraction(x: np.ndarray, window: int, threshold: float, min_periods: int) -> Tuple[float, int]:
    """NaN-aware centered rolling std (ddof=0), matching pandas window bounds."""
    n = x.size
    cs = np.zeros(n + 1, dtype=np.float64)
    cs2 = np.zeros(n + 1, dtype=np.float64)
    cc = np.zeros(n + 1, dtype=np.int64)
    for i in range(n):
        v = x[i]
        cs[i + 1] = cs[i]
        cs2[i + 1] = cs2[i]
        cc[i + 1] = cc[i]
        if math.isfinite(v):
            cs[i + 1] += v
            cs2[i + 1] += v * v
            cc[i + 1] += 1

    left = window // 2
    right = window - left
    n_good = 0
    n_exceed = 0
    for i in range(n):
        a = i - left
        if a < 0:
            a = 0
        b = i + right
        if b > n:
            b = n
        cnt = cc[b] - cc[a]
        if cnt >= min_periods:
            sm = cs[b] - cs[a]
            sm2 = cs2[b] - cs2[a]
            mean = sm / cnt
            var = sm2 / cnt - mean * mean
            if var < 0.0 and var > -1.0e-10:
                var = 0.0
            std = math.sqrt(var) if var >= 0.0 else 0.0
            n_good += 1
            if std > threshold:
                n_exceed += 1
    if n_good == 0:
        return math.nan, 0
    return n_exceed / n_good, n_good


@njit(cache=True, parallel=True)
def qc_kernel(
    tx_in: np.ndarray,
    rx_in: np.ndarray,
    proc_code: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Process canonical arrays of shape (n_series, n_time).

    proc_code: 0=no_filter, 1=graf_2020, 2=full.
    """
    n_series, n_time = tx_in.shape
    tx = tx_in.copy()
    rx = rx_in.copy()
    trsl = np.empty((n_series, n_time), dtype=tx_in.dtype)

    step1 = np.zeros((n_series, n_time), dtype=np.uint8)
    step2 = np.zeros((n_series, n_time), dtype=np.uint8)
    step3 = np.zeros((n_series, n_time), dtype=np.uint8)
    step4 = np.zeros((n_series, n_time), dtype=np.uint8)
    step5_fail = np.zeros(n_series, dtype=np.uint8)
    long_frac = np.full(n_series, np.nan, dtype=np.float32)
    short_frac = np.full(n_series, np.nan, dtype=np.float32)

    for s in prange(n_series):
        # no_filter records raw invalidity but intentionally leaves TX/RX unchanged.
        for i in range(n_time):
            bad_tx = _is_bad_tx(float(tx[s, i]))
            bad_rx = _is_bad_rx(float(rx[s, i]))
            if bad_tx or bad_rx:
                step1[s, i] = 1
            if proc_code >= 1:
                if bad_tx:
                    tx[s, i] = np.nan
                if bad_rx:
                    rx[s, i] = np.nan

        if proc_code == 2:
            # Base plateau mask using a centered 3-sample rolling max/std.
            base = np.zeros(n_time, dtype=np.uint8)
            for i in range(n_time):
                count = 0
                sm = 0.0
                sm2 = 0.0
                vmax = -math.inf
                a = i - 1
                if a < 0:
                    a = 0
                b = i + 2
                if b > n_time:
                    b = n_time
                for j in range(a, b):
                    v = float(rx[s, j])
                    if math.isfinite(v):
                        count += 1
                        sm += v
                        sm2 += v * v
                        if v > vmax:
                            vmax = v
                if count > 0:
                    mean = sm / count
                    var = sm2 / count - mean * mean
                    if var < 0.0 and var > -1.0e-10:
                        var = 0.0
                    std = math.sqrt(var) if var >= 0.0 else 0.0
                    if vmax < PLATEAU_MAX_DBM and std < PLATEAU_STD_DB:
                        base[i] = 1

            # Two neighbours on each side, equivalent to centered rolling max size 5.
            for i in range(n_time):
                a = i - PLATEAU_PAD
                if a < 0:
                    a = 0
                b = i + PLATEAU_PAD + 1
                if b > n_time:
                    b = n_time
                bad = 0
                for j in range(a, b):
                    if base[j] == 1:
                        bad = 1
                        break
                if bad == 1:
                    step2[s, i] = 1
                    tx[s, i] = np.nan
                    rx[s, i] = np.nan

            # Whole-series extrema after plateau filtering.
            tx_max = -math.inf
            rx_min = math.inf
            have_tx = False
            have_rx = False
            for i in range(n_time):
                tv = float(tx[s, i])
                rv = float(rx[s, i])
                if math.isfinite(tv):
                    have_tx = True
                    if tv > tx_max:
                        tx_max = tv
                if math.isfinite(rv):
                    have_rx = True
                    if rv < rx_min:
                        rx_min = rv

            # Fill complete NaN runs in RSL, if <=60 and adjacent to low signal.
            if have_tx and have_rx:
                i = 0
                while i < n_time:
                    if math.isfinite(float(rx[s, i])):
                        i += 1
                        continue
                    j = i + 1
                    while j < n_time and not math.isfinite(float(rx[s, j])):
                        j += 1
                    gap = j - i
                    if gap <= BLACKOUT_MAX_GAP:
                        left_low = i > 0 and math.isfinite(float(rx[s, i - 1])) and float(rx[s, i - 1]) < BLACKOUT_LOW_DBM
                        right_low = j < n_time and math.isfinite(float(rx[s, j])) and float(rx[s, j]) < BLACKOUT_LOW_DBM
                        if left_low or right_low:
                            for q in range(i, j):
                                tx[s, q] = tx_max
                                rx[s, q] = rx_min
                                step3[s, q] = 1
                    i = j

        # Total loss after steps 1-3, or raw total loss for no_filter.
        for i in range(n_time):
            tv = float(tx[s, i])
            rv = float(rx[s, i])
            trsl[s, i] = tv - rv if math.isfinite(tv) and math.isfinite(rv) else np.nan

        # Step 4 for graf/full: bounded linear interpolation of gaps <=5 samples.
        if proc_code >= 1:
            i = 0
            while i < n_time:
                if math.isfinite(float(trsl[s, i])):
                    i += 1
                    continue
                j = i + 1
                while j < n_time and not math.isfinite(float(trsl[s, j])):
                    j += 1
                gap = j - i
                if gap <= SHORT_GAP_MAX and i > 0 and j < n_time:
                    left = float(trsl[s, i - 1])
                    right = float(trsl[s, j])
                    if math.isfinite(left) and math.isfinite(right):
                        denom = gap + 1.0
                        for q in range(gap):
                            trsl[s, i + q] = left + (right - left) * ((q + 1.0) / denom)
                            step4[s, i + q] = 1
                i = j

            # Step 5: Blettner fluctuation tests.
            frac_long, ng_long = _rolling_std_exceed_fraction(
                trsl[s], LONG_STD_WINDOW, LONG_STD_THRESHOLD, ROLLING_MIN_PERIODS
            )
            frac_short, ng_short = _rolling_std_exceed_fraction(
                trsl[s], SHORT_STD_WINDOW, SHORT_STD_THRESHOLD, ROLLING_MIN_PERIODS
            )
            long_frac[s] = frac_long
            short_frac[s] = frac_short
            pass_long = ng_long > 0 and math.isfinite(frac_long) and frac_long <= LONG_STD_MAX_FRACTION
            pass_short = ng_short > 0 and math.isfinite(frac_short) and frac_short <= SHORT_STD_MAX_FRACTION
            if not (pass_long and pass_short):
                step5_fail[s] = 1

    return tx, rx, trsl, step1, step2, step3, step4, step5_fail, long_frac, short_frac


# -----------------------------------------------------------------------------
# NetCDF helpers
# -----------------------------------------------------------------------------
def _require_netcdf4() -> None:
    if Dataset is None:
        raise SystemExit(
            "The production streamed writer requires netCDF4. Install it in the "
            "active environment with `python -m pip install netCDF4`. Original "
            f"import error: {_NETCDF4_IMPORT_ERROR}"
        )


def _safe_attr_dict(var) -> Dict[str, object]:
    return {name: var.getncattr(name) for name in var.ncattrs() if name != "_FillValue"}


def _fill_value_of(var):
    return var.getncattr("_FillValue") if "_FillValue" in var.ncattrs() else None


def _create_variable_like(dst, name: str, src_var, *, compression: bool, compression_level: int, chunksizes=None):
    dtype = src_var.datatype
    fill = _fill_value_of(src_var)
    kwargs = {}
    if fill is not None:
        kwargs["fill_value"] = fill
    if compression and np.dtype(src_var.dtype).kind not in "OSU":
        kwargs.update(zlib=True, complevel=int(compression_level), shuffle=np.dtype(src_var.dtype).kind in "fiu")
    if chunksizes is not None and len(chunksizes) == len(src_var.dimensions):
        kwargs["chunksizes"] = tuple(int(x) for x in chunksizes)
    try:
        out = dst.createVariable(name, dtype, src_var.dimensions, **kwargs)
    except Exception:
        # String/VLEN variables do not accept compression/chunk options consistently.
        kwargs.pop("zlib", None)
        kwargs.pop("complevel", None)
        kwargs.pop("shuffle", None)
        kwargs.pop("chunksizes", None)
        out = dst.createVariable(name, dtype, src_var.dimensions, **kwargs)
    out.setncatts(_safe_attr_dict(src_var))
    return out


def _copy_data_streamed(src_var, dst_var, cml_dim: str, time_dim: str, link_batch: int, time_batch: int = 10080) -> None:
    """Copy a non-signal variable without loading large arrays unnecessarily."""
    dims = src_var.dimensions
    if len(dims) == 0:
        dst_var.assignValue(src_var.getValue())
        return
    if cml_dim in dims:
        ax = dims.index(cml_dim)
        n = src_var.shape[ax]
        for i0 in range(0, n, link_batch):
            i1 = min(i0 + link_batch, n)
            sl = [slice(None)] * len(dims)
            sl[ax] = slice(i0, i1)
            dst_var[tuple(sl)] = src_var[tuple(sl)]
        return
    if time_dim in dims and src_var.size > 5_000_000:
        ax = dims.index(time_dim)
        n = src_var.shape[ax]
        for i0 in range(0, n, time_batch):
            i1 = min(i0 + time_batch, n)
            sl = [slice(None)] * len(dims)
            sl[ax] = slice(i0, i1)
            dst_var[tuple(sl)] = src_var[tuple(sl)]
        return
    dst_var[:] = src_var[:]


def _masked_to_ndarray(data, dtype=None) -> np.ndarray:
    if np.ma.isMaskedArray(data):
        # Character metadata cannot be filled with NaN. Use an empty character
        # so fixed-length NetCDF strings remain decodable.
        kind = np.asarray(data.data).dtype.kind
        if kind == "S":
            data = data.filled(b"")
        elif kind == "U":
            data = data.filled("")
        else:
            data = data.filled(np.nan)
    arr = np.asarray(data)
    if dtype is not None:
        arr = arr.astype(dtype, copy=False)
    return arr


def _guess_name(available: Iterable[str], requested: Optional[str], candidates: Sequence[str], role: str) -> str:
    names = list(available)
    if requested:
        if requested not in names:
            raise KeyError(f"Requested {role} variable {requested!r} not found. Available: {names}")
        return requested
    lower = {n.lower(): n for n in names}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    raise KeyError(f"Could not auto-detect {role}. Supply its command-line option. Available: {names}")


def _guess_dim(dims: Sequence[str], requested: Optional[str], candidates: Sequence[str], role: str) -> str:
    if requested:
        if requested not in dims:
            raise KeyError(f"Requested {role} dimension {requested!r} not in signal dimensions {tuple(dims)}")
        return requested
    lower = {d.lower(): d for d in dims}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    for d in dims:
        dl = d.lower()
        if any(token in dl for token in candidates):
            return d
    raise KeyError(f"Could not auto-detect {role} dimension from {tuple(dims)}")


def _canonical_read(var, time_dim: str, sublink_dim: str, cml_dim: str, c0: int, c1: int) -> np.ndarray:
    dims = var.dimensions
    sl = [slice(None)] * len(dims)
    sl[dims.index(cml_dim)] = slice(c0, c1)
    arr = _masked_to_ndarray(var[tuple(sl)])
    axes = [dims.index(time_dim), dims.index(sublink_dim), dims.index(cml_dim)]
    arr = np.transpose(arr, axes)
    # (time, sublink, cml) -> flattened series x time.
    arr = np.transpose(arr, (1, 2, 0)).reshape(-1, arr.shape[0])
    if arr.dtype.kind not in "fc":
        arr = arr.astype(np.float32)
    return np.ascontiguousarray(arr)


def _canonical_to_native(arr_series_time: np.ndarray, native_dims: Sequence[str], time_dim: str, sublink_dim: str, cml_dim: str, n_sublink: int, n_links: int) -> np.ndarray:
    n_time = arr_series_time.shape[1]
    arr = arr_series_time.reshape(n_sublink, n_links, n_time).transpose(2, 0, 1)  # time, sublink, cml
    canonical_dims = [time_dim, sublink_dim, cml_dim]
    perm = [canonical_dims.index(d) for d in native_dims]
    return np.transpose(arr, perm)


def _batch_slice(native_dims: Sequence[str], cml_dim: str, c0: int, c1: int) -> Tuple[slice, ...]:
    sl = [slice(None)] * len(native_dims)
    sl[native_dims.index(cml_dim)] = slice(c0, c1)
    return tuple(sl)


def _chunks_for_dims(dims: Sequence[str], sizes: Dict[str, int], time_dim: str, sublink_dim: str, cml_dim: str, time_chunk: int, link_chunk: int) -> Tuple[int, ...]:
    out = []
    for d in dims:
        if d == time_dim:
            out.append(min(sizes[d], int(time_chunk)))
        elif d == sublink_dim:
            out.append(sizes[d])
        elif d == cml_dim:
            out.append(min(sizes[d], int(link_chunk)))
        else:
            out.append(min(sizes[d], 1))
    return tuple(max(1, int(v)) for v in out)


def _new_numeric_var(dst, name: str, dtype: str, dims: Sequence[str], attrs: Dict[str, object], *, compression: bool, level: int, chunksizes: Optional[Tuple[int, ...]] = None, fill_value=None):
    kwargs = {}
    if fill_value is not None:
        kwargs["fill_value"] = fill_value
    if compression:
        kwargs.update(zlib=True, complevel=int(level), shuffle=True)
    if chunksizes:
        kwargs["chunksizes"] = chunksizes
    v = dst.createVariable(name, dtype, tuple(dims), **kwargs)
    v.setncatts(attrs)
    return v


def _append_ancillary(var, names: Sequence[str]) -> None:
    existing = str(getattr(var, "ancillary_variables", "")).strip()
    merged = " ".join([x for x in [existing, *names] if x])
    var.setncattr("ancillary_variables", merged)


def _decode_fixed_length_strings(arr: np.ndarray, dims: Tuple[str, ...], sublink_dim: str, cml_dim: str) -> Tuple[np.ndarray, Tuple[str, ...]]:
    """Collapse a trailing NetCDF character dimension into full strings.

    Example:
        (sublink_id, cml_id, string10) S1
    becomes:
        (sublink_id, cml_id) string
    """
    if (
        arr.ndim >= 2
        and len(dims) == arr.ndim
        and dims[-1] not in (sublink_dim, cml_dim)
        and arr.dtype.kind in "SU"
        and arr.dtype.itemsize in (1, 4)
    ):
        if arr.dtype.kind == "S":
            try:
                decoded = np.asarray(chartostring(arr), dtype=object)
            except Exception:
                shape = arr.shape[:-1]
                rows = arr.reshape(-1, arr.shape[-1])
                decoded = np.array(
                    [b"".join(np.asarray(row, dtype="S1").tolist()).decode("utf-8", "ignore").rstrip("\x00") for row in rows],
                    dtype=object,
                ).reshape(shape)
        else:
            shape = arr.shape[:-1]
            rows = arr.reshape(-1, arr.shape[-1])
            decoded = np.array(
                ["".join(str(x) for x in row).rstrip("\x00") for row in rows],
                dtype=object,
            ).reshape(shape)
        return decoded, tuple(dims[:-1])
    return arr, tuple(dims)


def _metadata_array(var, sublink_dim: str, cml_dim: str, n_sublink: int, n_cml: int) -> np.ndarray:
    arr = _masked_to_ndarray(var[:])
    dims = tuple(var.dimensions)
    arr, dims = _decode_fixed_length_strings(arr, dims, sublink_dim, cml_dim)

    if dims == (sublink_dim, cml_dim):
        return arr
    if dims == (cml_dim, sublink_dim):
        return arr.T
    if dims == (cml_dim,):
        return np.broadcast_to(arr[None, :], (n_sublink, n_cml))
    if dims == (sublink_dim,):
        return np.broadcast_to(arr[:, None], (n_sublink, n_cml))
    if arr.ndim == 0:
        return np.full((n_sublink, n_cml), arr.item())
    raise ValueError(
        f"Metadata variable {getattr(var, 'name', '?')} has unsupported dimensions {dims}; "
        f"expected ({sublink_dim},{cml_dim}), ({cml_dim},{sublink_dim}), or one-dimensional."
    )


def _worker_candidates(max_workers: int) -> List[int]:
    """Return up to five deterministic worker-count trials."""
    maximum = max(1, int(max_workers))
    if maximum == 1:
        return [1]
    base = [1, 4, 8, 16, maximum]
    candidates = sorted({min(maximum, max(1, x)) for x in base})
    # Keep a maximum of five trials and always retain the requested maximum.
    if len(candidates) > 5:
        candidates = candidates[:4] + [maximum]
    return candidates


def _run_benchmark(
    tx: np.ndarray,
    rx: np.ndarray,
    proc_code: int,
    requested_threads: int,
    estimate_links: int,
    n_sublink: int,
    repeats: int,
    auto_tune: bool,
) -> Tuple[Dict[str, object], int]:
    """Benchmark 1-5 thread-count iterations and return the fastest setting."""
    warm_n = min(tx.shape[0], max(1, n_sublink))
    # Compile kernels and initialize Numba before any timing.
    set_num_threads(1)
    qc_kernel(tx[:warm_n, : min(tx.shape[1], 512)], rx[:warm_n, : min(rx.shape[1], 512)], proc_code)

    candidates = _worker_candidates(requested_threads) if auto_tune else [max(1, requested_threads)]
    trials: List[Dict[str, object]] = []
    repeats = max(1, int(repeats))
    for iteration, workers in enumerate(candidates, start=1):
        set_num_threads(int(workers))
        # Untimed call initializes/reconfigures the parallel runtime.
        qc_kernel(tx, rx, proc_code)
        timings = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            qc_kernel(tx, rx, proc_code)
            timings.append(float(time.perf_counter() - t0))
        median_s = float(np.median(np.asarray(timings, dtype=np.float64)))
        trials.append({
            "iteration": int(iteration),
            "threads": int(workers),
            "median_seconds_compute_only": median_s,
            "runs_seconds": timings,
        })

    best = min(trials, key=lambda row: float(row["median_seconds_compute_only"]))
    selected_threads = int(best["threads"])
    best_seconds = float(best["median_seconds_compute_only"])
    links_tested = max(1, tx.shape[0] // max(1, n_sublink))
    estimate = best_seconds * float(estimate_links) / float(links_tested) if estimate_links > 0 else np.nan
    report: Dict[str, object] = {
        "benchmark_links": int(links_tested),
        "benchmark_samples": int(tx.size),
        "worker_trials": trials,
        "selected_threads": selected_threads,
        "best_seconds_compute_only": best_seconds,
        "best_seconds_per_link_compute_only": float(best_seconds / links_tested),
        "estimated_seconds_for_requested_links_compute_only": float(estimate),
        "important": "Estimate excludes reading, compression, copying other variables, and NetCDF writing.",
    }
    set_num_threads(selected_threads)
    return report, selected_threads


def process_file(args: argparse.Namespace) -> None:
    _require_netcdf4()
    src_path = Path(args.cml_nc).expanduser().resolve()
    out_path = Path(args.out_nc).expanduser().resolve() if args.out_nc else src_path.with_name(src_path.stem + "_QC.nc")
    if not src_path.exists():
        raise FileNotFoundError(src_path)
    if out_path.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists: {out_path}. Add --overwrite to replace it.")
    if src_path == out_path:
        raise ValueError("Input and output paths must differ; the raw input is never modified in place.")

    proc_codes = {"no_filter": 0, "graf_2020": 1, "full": 2}
    proc_code = proc_codes[args.proc_line]
    n_workers = max(1, min(int(args.n_workers), os.cpu_count() or 1))
    set_num_threads(n_workers)

    temp_path = out_path.with_name(out_path.name + ".partial")
    if temp_path.exists():
        temp_path.unlink()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    t_all = time.perf_counter()
    with Dataset(src_path, "r") as src:
        src.set_auto_mask(True)
        tx_name = _guess_name(src.variables, args.tx_var, ["tx", "tsl", "tx_signal", "transmitted_signal_level"], "TX")
        rx_name = _guess_name(src.variables, args.rx_var, ["rx", "rsl", "rx_signal", "received_signal_level"], "RX")
        tx_src = src.variables[tx_name]
        rx_src = src.variables[rx_name]
        if tx_src.dimensions != rx_src.dimensions:
            raise ValueError(f"TX/RX dimensions differ: {tx_src.dimensions} vs {rx_src.dimensions}")
        if len(tx_src.dimensions) != 3:
            raise ValueError(f"TX/RX must be 3-D. Found dimensions {tx_src.dimensions}")

        native_dims = tx_src.dimensions
        time_dim = _guess_dim(native_dims, args.time_dim, ["time"], "time")
        sublink_dim = _guess_dim(native_dims, args.sublink_dim, ["sublink_id", "sublink", "channel"], "sublink")
        cml_dim = _guess_dim(native_dims, args.cml_dim, ["cml_id", "cml", "link"], "CML")
        if len({time_dim, sublink_dim, cml_dim}) != 3:
            raise ValueError("The time, sublink, and CML dimensions must be distinct.")

        sizes = {d: len(src.dimensions[d]) for d in src.dimensions}
        n_time = sizes[time_dim]
        n_sublink = sizes[sublink_dim]
        n_cml = sizes[cml_dim]

        freq_name = _guess_name(src.variables, args.frequency_var, ["frequency", "freq", "frequency_GHz", "frequency_MHz"], "frequency")
        pol_name = _guess_name(src.variables, args.polarization_var, ["polarization", "pol"], "polarization")
        length_name = _guess_name(src.variables, args.length_var, ["length", "length_km", "path_length"], "link length")

        freq_raw = _metadata_array(src.variables[freq_name], sublink_dim, cml_dim, n_sublink, n_cml)
        freq_ghz, inferred_freq_unit = normalize_frequency_to_ghz(freq_raw)
        pol = polarization_codes(_metadata_array(src.variables[pol_name], sublink_dim, cml_dim, n_sublink, n_cml))
        length = _metadata_array(src.variables[length_name], sublink_dim, cml_dim, n_sublink, n_cml).astype(np.float64)
        detection = detection_limits_legacy(freq_ghz, pol, length)
        step6_fail_all = (~np.isfinite(detection) | (detection >= DETECTION_LIMIT_MAX_MMPH)).astype(np.uint8)
        if proc_code < 2:
            step6_applied_all = np.zeros_like(step6_fail_all)
        else:
            step6_applied_all = step6_fail_all

        print("=" * 88)
        print("Standalone streamed CML quality control")
        print(f"Input:             {src_path}")
        print(f"Output:            {out_path}")
        print(f"Signals:           TX={tx_name}, RX={rx_name}, dims={native_dims}")
        print(f"Shape:             time={n_time:,}, sublink={n_sublink}, cml={n_cml:,}")
        print(f"Processing line:   {args.proc_line}")
        print(f"Numba threads:     {get_num_threads()}")
        print(f"Link batch:        {args.link_batch}")
        print(f"Frequency inferred:{inferred_freq_unit}; GHz range={np.nanmin(freq_ghz):.3f}–{np.nanmax(freq_ghz):.3f}")
        print("=" * 88)

        if args.benchmark_links > 0:
            b1 = min(n_cml, int(args.benchmark_links))
            tx_b = _canonical_read(tx_src, time_dim, sublink_dim, cml_dim, 0, b1)
            rx_b = _canonical_read(rx_src, time_dim, sublink_dim, cml_dim, 0, b1)
            bench, selected_workers = _run_benchmark(
                tx_b, rx_b, proc_code, n_workers, int(args.estimate_total_links),
                n_sublink, int(args.benchmark_repeats), not args.no_auto_tune_workers,
            )
            n_workers = int(selected_workers)
            set_num_threads(n_workers)
            print("Benchmark and worker tuning (compute only; excludes NetCDF I/O):")
            print(json.dumps(bench, indent=2))
            print(f"Selected {n_workers} compute threads for the full run.")
            del tx_b, rx_b
            if args.benchmark_only:
                return

        compression = not args.no_compression
        chunks_native = _chunks_for_dims(
            native_dims, sizes, time_dim, sublink_dim, cml_dim,
            args.time_chunk, min(args.link_batch, args.link_chunk)
        )

        with Dataset(temp_path, "w", format=args.output_format) as dst:
            dst.set_fill_off()
            # Dimensions
            for name, dim in src.dimensions.items():
                dst.createDimension(name, None if dim.isunlimited() else len(dim))
            # Global metadata
            dst.setncatts({name: src.getncattr(name) for name in src.ncattrs()})
            dst.setncattr("quality_control_applied", args.proc_line)
            dst.setncattr("quality_control_program", "cml_qc_filter_optimized.py standalone streamed implementation")
            dst.setncattr("quality_control_created_utc", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            dst.setncattr("quality_control_frequency_input_unit_inferred", inferred_freq_unit)
            dst.setncattr("quality_control_link_aggregation", "final link passes when any directional sublink passes")
            dst.setncattr("quality_control_thresholds_json", json.dumps({
                "tx_min_dbm": TX_MIN_DBM,
                "tx_sentinel": TX_SENTINEL,
                "rx_min_dbm": RX_MIN_DBM,
                "rx_max_dbm": RX_MAX_DBM,
                "plateau_max_dbm": PLATEAU_MAX_DBM,
                "plateau_std_db": PLATEAU_STD_DB,
                "plateau_pad_samples": PLATEAU_PAD,
                "blackout_low_dbm": BLACKOUT_LOW_DBM,
                "blackout_max_gap_samples": BLACKOUT_MAX_GAP,
                "short_gap_max_samples": SHORT_GAP_MAX,
                "long_std_window_samples": LONG_STD_WINDOW,
                "long_std_threshold_db": LONG_STD_THRESHOLD,
                "long_std_max_fraction": LONG_STD_MAX_FRACTION,
                "short_std_window_samples": SHORT_STD_WINDOW,
                "short_std_threshold_db": SHORT_STD_THRESHOLD,
                "short_std_max_fraction": SHORT_STD_MAX_FRACTION,
                "detection_limit_max_mmph": DETECTION_LIMIT_MAX_MMPH,
                "detection_limit_formula": "0.33/(ITU_P838_3_k*length_km), legacy-compatible",
            }, sort_keys=True))

            # Create and copy all original variables. TX/RX are created now but written
            # batchwise after QC. This preserves coordinate and metadata functionality.
            out_vars = {}
            for name, var in src.variables.items():
                use_chunks = chunks_native if name in (tx_name, rx_name) else None
                out_vars[name] = _create_variable_like(
                    dst, name, var,
                    compression=compression and name in (tx_name, rx_name),
                    compression_level=args.compression_level,
                    chunksizes=use_chunks,
                )

            for name, var in src.variables.items():
                if name in (tx_name, rx_name):
                    continue
                _copy_data_streamed(var, out_vars[name], cml_dim, time_dim, args.link_batch)

            if args.keep_original_signals:
                tx_raw = _new_numeric_var(
                    dst, tx_name + "_raw", str(tx_src.dtype), native_dims,
                    {"long_name": f"Original unmodified {tx_name}", "source_variable": tx_name},
                    compression=compression, level=args.compression_level, chunksizes=chunks_native,
                    fill_value=_fill_value_of(tx_src),
                )
                rx_raw = _new_numeric_var(
                    dst, rx_name + "_raw", str(rx_src.dtype), native_dims,
                    {"long_name": f"Original unmodified {rx_name}", "source_variable": rx_name},
                    compression=compression, level=args.compression_level, chunksizes=chunks_native,
                    fill_value=_fill_value_of(rx_src),
                )
            else:
                tx_raw = rx_raw = None

            # New sample-level variables follow the original TX/RX dimension order.
            flag_attrs = {
                "flag_values": np.array([0, 1], dtype=np.uint8),
                "flag_meanings": "not_affected affected",
                "valid_min": np.uint8(0),
                "valid_max": np.uint8(1),
            }
            trsl_var = _new_numeric_var(
                dst, "trsl_qc", "f4", native_dims,
                {"long_name": "Quality-controlled total signal loss (TX minus RX)", "units": getattr(tx_src, "units", "dB")},
                compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.float32(np.nan),
            )
            step1_var = _new_numeric_var(dst, "qc_step1_invalid_signal", "u1", native_dims,
                {**flag_attrs, "long_name": "QC step 1: missing, fill, or physically invalid TX/RX sample"},
                compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))
            step2_var = _new_numeric_var(dst, "qc_step2_plateau", "u1", native_dims,
                {**flag_attrs, "long_name": "QC step 2: low stable RSL plateau sample including two-sample padding"},
                compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))
            step3_var = _new_numeric_var(dst, "qc_step3_blackout_filled", "u1", native_dims,
                {**flag_attrs, "long_name": "QC step 3: sample reconstructed as part of a short blackout gap"},
                compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))
            step4_var = _new_numeric_var(dst, "qc_step4_short_gap_interpolated", "u1", native_dims,
                {**flag_attrs, "long_name": "QC step 4: total-loss sample linearly interpolated across a gap of at most five samples"},
                compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))
            final_var = _new_numeric_var(dst, "qc_final_valid", "u1", native_dims,
                {"flag_values": np.array([0, 1], dtype=np.uint8), "flag_meanings": "invalid valid",
                 "long_name": "Final sample-level QC validity after point and sublink filters"},
                compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))
            if args.write_bitmask:
                bitmask_var = _new_numeric_var(dst, "qc_flag", "u2", native_dims,
                    {"long_name": "Compact QC bitmask", "flag_masks": np.array([1,2,4,8,16,32,64], dtype=np.uint16),
                     "flag_meanings": "step1_invalid_signal step2_plateau step3_blackout_filled step4_short_gap_interpolated step5_fluctuation_fail step6_detection_limit_fail final_invalid"},
                    compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint16(65535))
            else:
                bitmask_var = None

            small_dims = (sublink_dim, cml_dim)
            small_chunks = (n_sublink, min(n_cml, max(1, args.link_chunk)))
            step5_var = _new_numeric_var(dst, "qc_step5_fluctuation_fail", "u1", small_dims,
                {**flag_attrs, "long_name": "QC step 5: sublink failed one or both total-loss fluctuation tests"},
                compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.uint8(255))
            step6_var = _new_numeric_var(dst, "qc_step6_detection_limit_fail", "u1", small_dims,
                {**flag_attrs, "long_name": "QC step 6: sublink detection limit is missing or at least 2 mm h-1"},
                compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.uint8(255))
            det_var = _new_numeric_var(dst, "qc_detection_limit_mmph", "f4", small_dims,
                {"long_name": "Legacy-compatible rainfall detection limit", "units": "mm h-1"},
                compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.float32(np.nan))
            long_frac_var = _new_numeric_var(dst, "qc_step5_fraction_long_std_exceed", "f4", small_dims,
                {"long_name": "Fraction of valid centered 5-hour windows with total-loss standard deviation above 2 dB", "units": "1"},
                compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.float32(np.nan))
            short_frac_var = _new_numeric_var(dst, "qc_step5_fraction_short_std_exceed", "f4", small_dims,
                {"long_name": "Fraction of valid centered 1-hour windows with total-loss standard deviation above 0.8 dB", "units": "1"},
                compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.float32(np.nan))
            final_sl_var = _new_numeric_var(dst, "qc_final_sublink_pass", "u1", small_dims,
                {"flag_values": np.array([0,1],dtype=np.uint8), "flag_meanings": "fail pass",
                 "long_name": "Final sublink-level QC result"},
                compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.uint8(255))
            final_link_var = _new_numeric_var(dst, "qc_final_link_pass", "u1", (cml_dim,),
                {"flag_values": np.array([0,1],dtype=np.uint8), "flag_meanings": "fail pass",
                 "long_name": "Final CML-level QC result; pass when at least one directional sublink passes"},
                compression=compression, level=args.compression_level, chunksizes=(min(n_cml, max(1,args.link_chunk)),), fill_value=np.uint8(255))

            ancillary = ["qc_final_valid", "qc_final_sublink_pass", "qc_final_link_pass"]
            if bitmask_var is not None:
                ancillary.insert(0, "qc_flag")
            _append_ancillary(out_vars[tx_name], ancillary)
            _append_ancillary(out_vars[rx_name], ancillary)
            out_vars[tx_name].setncattr("quality_controlled", "true")
            out_vars[rx_name].setncattr("quality_controlled", "true")

            # Small variables known before the batch loop.
            det_var[:, :] = detection.astype(np.float32)
            step6_var[:, :] = step6_applied_all.astype(np.uint8)

            counts = {
                "step1_invalid_samples": 0,
                "step2_plateau_samples": 0,
                "step3_blackout_filled_samples": 0,
                "step4_interpolated_samples": 0,
                "step5_failed_sublinks": 0,
                "step6_failed_sublinks": int(step6_applied_all.sum()),
                "final_valid_samples": 0,
                "final_pass_sublinks": 0,
                "final_pass_links": 0,
            }

            t_compute = 0.0
            t_write = 0.0
            for c0 in range(0, n_cml, int(args.link_batch)):
                c1 = min(c0 + int(args.link_batch), n_cml)
                nb = c1 - c0
                t0 = time.perf_counter()
                tx_in = _canonical_read(tx_src, time_dim, sublink_dim, cml_dim, c0, c1)
                rx_in = _canonical_read(rx_src, time_dim, sublink_dim, cml_dim, c0, c1)
                tx_qc, rx_qc, trsl, f1, f2, f3, f4, f5_flat, long_frac, short_frac = qc_kernel(tx_in, rx_in, proc_code)
                t_compute += time.perf_counter() - t0

                f6 = step6_applied_all[:, c0:c1].reshape(-1).astype(np.uint8)
                if proc_code == 0:
                    sl_pass = np.ones_like(f5_flat, dtype=np.uint8)
                elif proc_code == 1:
                    sl_pass = (f5_flat == 0).astype(np.uint8)
                else:
                    sl_pass = ((f5_flat == 0) & (f6 == 0)).astype(np.uint8)
                valid = (np.isfinite(trsl) & (sl_pass[:, None] == 1)).astype(np.uint8)
                link_pass = (sl_pass.reshape(n_sublink, nb).any(axis=0)).astype(np.uint8)

                if bitmask_var is not None:
                    bitmask = (
                        f1.astype(np.uint16) * BIT_STEP1_INVALID
                        | f2.astype(np.uint16) * BIT_STEP2_PLATEAU
                        | f3.astype(np.uint16) * BIT_STEP3_BLACKOUT_FILLED
                        | f4.astype(np.uint16) * BIT_STEP4_INTERPOLATED
                        | f5_flat[:, None].astype(np.uint16) * BIT_STEP5_FLUCTUATION_FAIL
                        | f6[:, None].astype(np.uint16) * BIT_STEP6_DETECTION_FAIL
                        | (1 - valid).astype(np.uint16) * BIT_FINAL_INVALID
                    )
                else:
                    bitmask = None

                t0 = time.perf_counter()
                native_slice = _batch_slice(native_dims, cml_dim, c0, c1)
                out_vars[tx_name][native_slice] = _canonical_to_native(tx_qc, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                out_vars[rx_name][native_slice] = _canonical_to_native(rx_qc, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                trsl_var[native_slice] = _canonical_to_native(trsl.astype(np.float32), native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                step1_var[native_slice] = _canonical_to_native(f1, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                step2_var[native_slice] = _canonical_to_native(f2, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                step3_var[native_slice] = _canonical_to_native(f3, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                step4_var[native_slice] = _canonical_to_native(f4, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                final_var[native_slice] = _canonical_to_native(valid, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                if bitmask_var is not None:
                    bitmask_var[native_slice] = _canonical_to_native(bitmask, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                if tx_raw is not None:
                    tx_raw[native_slice] = _canonical_to_native(tx_in, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                    rx_raw[native_slice] = _canonical_to_native(rx_in, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)

                step5_var[:, c0:c1] = f5_flat.reshape(n_sublink, nb)
                long_frac_var[:, c0:c1] = long_frac.reshape(n_sublink, nb)
                short_frac_var[:, c0:c1] = short_frac.reshape(n_sublink, nb)
                final_sl_var[:, c0:c1] = sl_pass.reshape(n_sublink, nb)
                final_link_var[c0:c1] = link_pass
                t_write += time.perf_counter() - t0

                counts["step1_invalid_samples"] += int(f1.sum())
                counts["step2_plateau_samples"] += int(f2.sum())
                counts["step3_blackout_filled_samples"] += int(f3.sum())
                counts["step4_interpolated_samples"] += int(f4.sum())
                counts["step5_failed_sublinks"] += int(f5_flat.sum())
                counts["final_valid_samples"] += int(valid.sum())
                counts["final_pass_sublinks"] += int(sl_pass.sum())
                counts["final_pass_links"] += int(link_pass.sum())

                elapsed = time.perf_counter() - t_all
                rate = c1 / elapsed if elapsed > 0 else np.nan
                eta = (n_cml - c1) / rate if rate > 0 else np.nan
                print(f"Processed links {c0:5d}:{c1:5d}/{n_cml} | elapsed={elapsed:8.1f}s | ETA≈{eta:8.1f}s", flush=True)

            dst.setncattr("quality_control_summary_json", json.dumps(counts, sort_keys=True))
            dst.setncattr("quality_control_compute_seconds", float(t_compute))
            dst.setncattr("quality_control_write_seconds", float(t_write))
            dst.sync()

    # Atomic completion: a partial file is never mistaken for a finished product.
    if out_path.exists():
        out_path.unlink()
    os.replace(temp_path, out_path)
    total_elapsed = time.perf_counter() - t_all
    print("=" * 88)
    print(f"Saved quality-controlled NetCDF: {out_path}")
    print(f"Total elapsed: {total_elapsed:.1f} s")
    print(f"Output size: {out_path.stat().st_size / (1024**3):.3f} GiB")
    print("=" * 88)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Standalone streamed CML QC with per-step and final NetCDF flags.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--cml-nc", required=True, help="Input CML NetCDF.")
    ap.add_argument("--out-nc", default=None, help="Output QC NetCDF; default is INPUT_STEM_QC.nc.")
    ap.add_argument("--proc-line", choices=["no_filter", "graf_2020", "full"], default="full")
    ap.add_argument("--n-workers", type=int, default=max(1, min(16, os.cpu_count() or 1)), help="Numba compute threads; NetCDF I/O remains single-process.")
    ap.add_argument("--link-batch", type=int, default=64, help="Number of CML links read, processed, and written per batch.")
    ap.add_argument("--link-chunk", type=int, default=32, help="CML chunk length in the output NetCDF.")
    ap.add_argument("--time-chunk", type=int, default=1440, help="Time chunk length in the output NetCDF.")
    ap.add_argument("--compression-level", type=int, choices=range(0, 10), default=1, help="zlib compression level; level 1 is usually the best speed/size compromise.")
    ap.add_argument("--no-compression", action="store_true", help="Disable compression for fastest writing and larger files.")
    ap.add_argument("--keep-original-signals", action="store_true", help="Also store raw TX/RX copies; substantially increases output size and write time.")
    ap.add_argument("--write-bitmask", action="store_true", help="Also write the redundant uint16 compact qc_flag bitmask. Separate step flags are always written.")
    ap.add_argument("--benchmark-links", type=int, default=0, help="Benchmark the first N links before the full run.")
    ap.add_argument("--benchmark-repeats", type=int, default=3, help="Timed repeats for each worker-count tuning iteration.")
    ap.add_argument("--estimate-total-links", type=int, default=2000, help="Scale compute-only benchmark to this link count.")
    ap.add_argument("--no-auto-tune-workers", action="store_true", help="Benchmark only the requested --n-workers value instead of testing up to five values.")
    ap.add_argument("--benchmark-only", action="store_true", help="Run benchmark and exit without creating an output file.")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--output-format", choices=["NETCDF4"], default="NETCDF4")

    ap.add_argument("--tx-var", default=None)
    ap.add_argument("--rx-var", default=None)
    ap.add_argument("--frequency-var", default=None)
    ap.add_argument("--polarization-var", default=None)
    ap.add_argument("--length-var", default=None)
    ap.add_argument("--time-dim", default=None)
    ap.add_argument("--sublink-dim", default=None)
    ap.add_argument("--cml-dim", default=None)
    return ap


def main() -> None:
    args = build_parser().parse_args()
    try:
        process_file(args)
    except KeyboardInterrupt:
        print("Interrupted. Any .partial output can be removed safely.", file=sys.stderr)
        raise SystemExit(130)


if __name__ == "__main__":
    main()
