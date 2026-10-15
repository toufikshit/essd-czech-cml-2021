#!/usr/bin/env python3
"""
cml_qc_filter_standalone_v5.py
================================
Complete standalone quality control and strict compaction of a raw commercial
microwave-link NetCDF dataset. This V5 program does not require pre-existing QC
flags and therefore uses --cml-nc rather than --qc-nc.

The program computes the full per-sublink QC directly from TX/RX, physically
removes physical links for which no directional sublink passes, preserves the
original time x sublink_id x cml_id structure, and stores failed sublinks and
invalid observations as missing values. This is required because different
sublinks fail at different timestamps and therefore cannot share a physically
shortened time dimension.

QC sequence
-----------
1. Invalid TX/RX masking.
2. Low stable RSL plateau removal.
3. Short blackout reconstruction.
4. Short TRSL gap interpolation.
5. Per-sublink long/short fluctuation screening.
6. Per-sublink legacy-compatible ITU detection-limit screening.
7. Per-sublink valid-data availability screening.

A physical CML is retained when at least one directional sublink passes all
sublink-level filters. Only observations with finite processed TRSL on a passed
sublink are valid in the strict output.
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



def _canonical_read_indices(var, time_dim: str, sublink_dim: str, cml_dim: str, indices: np.ndarray) -> np.ndarray:
    """Read arbitrary source-link indices and return canonical (series, time)."""
    dims = var.dimensions
    sl = [slice(None)] * len(dims)
    sl[dims.index(cml_dim)] = np.asarray(indices, dtype=np.int64)
    arr = _masked_to_ndarray(var[tuple(sl)])
    axes = [dims.index(time_dim), dims.index(sublink_dim), dims.index(cml_dim)]
    arr = np.transpose(arr, axes)
    arr = np.transpose(arr, (1, 2, 0)).reshape(-1, arr.shape[0])
    if arr.dtype.kind not in "fc":
        arr = arr.astype(np.float32)
    return np.ascontiguousarray(arr)


def _copy_selected_variable(
    src_var,
    dst_var,
    *,
    cml_dim: str,
    time_dim: str,
    source_indices: np.ndarray,
    link_batch: int,
    time_batch: int = 10080,
) -> None:
    """Copy a variable while compacting its CML dimension."""
    dims = src_var.dimensions
    if len(dims) == 0:
        dst_var.assignValue(src_var.getValue())
        return

    if cml_dim in dims:
        ax = dims.index(cml_dim)
        n_out = int(source_indices.size)
        for o0 in range(0, n_out, int(link_batch)):
            o1 = min(o0 + int(link_batch), n_out)
            src_sel = np.asarray(source_indices[o0:o1], dtype=np.int64)
            src_sl = [slice(None)] * len(dims)
            dst_sl = [slice(None)] * len(dims)
            src_sl[ax] = src_sel
            dst_sl[ax] = slice(o0, o1)
            dst_var[tuple(dst_sl)] = src_var[tuple(src_sl)]
        return

    if time_dim in dims and src_var.size > 5_000_000:
        ax = dims.index(time_dim)
        n = src_var.shape[ax]
        for t0 in range(0, n, int(time_batch)):
            t1 = min(t0 + int(time_batch), n)
            sl = [slice(None)] * len(dims)
            sl[ax] = slice(t0, t1)
            dst_var[tuple(sl)] = src_var[tuple(sl)]
        return

    dst_var[:] = src_var[:]


def _write_summary_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "stage",
            "retained_links",
            "dropped_links_at_stage",
            "retained_sublinks",
            "dropped_sublinks_at_stage",
        ])
        writer.writeheader()
        writer.writerows(rows)


def process_file(args: argparse.Namespace) -> None:
    """Compute full QC from the raw file and write a compact strict-QC file."""
    _require_netcdf4()
    src_path = Path(args.cml_nc).expanduser().resolve()
    out_path = Path(args.out_nc).expanduser().resolve()
    summary_json = Path(args.summary_json).expanduser().resolve() if args.summary_json else out_path.with_name(out_path.stem + "_summary.json")
    summary_csv = Path(args.summary_csv).expanduser().resolve() if args.summary_csv else out_path.with_name(out_path.stem + "_filtering_steps.csv")

    if not src_path.exists():
        raise FileNotFoundError(src_path)
    if src_path == out_path:
        raise ValueError("Input and output paths must differ; the original file is never modified.")
    if out_path.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists: {out_path}. Add --overwrite to replace it.")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = out_path.with_name(out_path.name + ".partial")
    if temp_path.exists():
        temp_path.unlink()

    n_workers = max(1, min(int(args.n_workers), os.cpu_count() or 1))
    set_num_threads(n_workers)
    proc_code = 2  # complete six-step QC
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
            raise ValueError(f"TX/RX must be 3-D. Found {tx_src.dimensions}")
        if np.dtype(tx_src.dtype).kind not in "fc" or np.dtype(rx_src.dtype).kind not in "fc":
            raise TypeError("TX and RX must be floating-point variables so rejected observations can be stored as NaN.")

        native_dims = tx_src.dimensions
        time_dim = _guess_dim(native_dims, args.time_dim, ["time"], "time")
        sublink_dim = _guess_dim(native_dims, args.sublink_dim, ["sublink_id", "sublink", "channel"], "sublink")
        cml_dim = _guess_dim(native_dims, args.cml_dim, ["cml_id", "cml", "link"], "CML")
        if len({time_dim, sublink_dim, cml_dim}) != 3:
            raise ValueError("time, sublink and CML dimensions must be distinct")

        sizes_src = {d: len(src.dimensions[d]) for d in src.dimensions}
        n_time = sizes_src[time_dim]
        n_sublink = sizes_src[sublink_dim]
        n_cml = sizes_src[cml_dim]

        freq_name = _guess_name(src.variables, args.frequency_var, ["frequency", "freq", "frequency_GHz", "frequency_MHz"], "frequency")
        pol_name = _guess_name(src.variables, args.polarization_var, ["polarization", "pol"], "polarization")
        length_name = _guess_name(src.variables, args.length_var, ["length", "length_km", "path_length"], "link length")

        freq_raw = _metadata_array(src.variables[freq_name], sublink_dim, cml_dim, n_sublink, n_cml)
        freq_ghz, inferred_freq_unit = normalize_frequency_to_ghz(freq_raw)
        pol = polarization_codes(_metadata_array(src.variables[pol_name], sublink_dim, cml_dim, n_sublink, n_cml))
        length = _metadata_array(src.variables[length_name], sublink_dim, cml_dim, n_sublink, n_cml).astype(np.float64)
        detection = detection_limits_legacy(freq_ghz, pol, length)
        step6_fail_all = (~np.isfinite(detection) | (detection >= DETECTION_LIMIT_MAX_MMPH)).astype(np.uint8)

        print("=" * 92)
        print("Running cml_qc_filter_standalone_v5.py version 5.0")
        print("Complete raw-data QC followed by strict compact output")
        print(f"Input                 : {src_path}")
        print(f"Output                : {out_path}")
        print(f"Signals               : TX={tx_name}, RX={rx_name}")
        print(f"Input dimensions      : time={n_time:,}, sublinks={n_sublink}, links={n_cml:,}")
        print(f"Minimum valid fraction: {args.min_valid_fraction:.3f}")
        print(f"Minimum valid samples : {args.min_valid_observations:,}")
        print(f"Numba threads         : {get_num_threads()}")
        print(f"Link batch            : {args.link_batch}")
        print(f"Frequency input       : inferred {inferred_freq_unit}")
        print("=" * 92)

        # ------------------------------------------------------------------
        # PASS 1: determine final sublink and physical-link population.
        # ------------------------------------------------------------------
        fluct_pass = np.zeros((n_sublink, n_cml), dtype=np.uint8)
        det_pass = (step6_fail_all == 0).astype(np.uint8)
        availability_pass = np.zeros((n_sublink, n_cml), dtype=np.uint8)
        availability_fraction = np.full((n_sublink, n_cml), np.nan, dtype=np.float32)
        availability_count = np.zeros((n_sublink, n_cml), dtype=np.int64)
        long_fraction_all = np.full((n_sublink, n_cml), np.nan, dtype=np.float32)
        short_fraction_all = np.full((n_sublink, n_cml), np.nan, dtype=np.float32)

        print("\nPASS 1/2: computing final link and sublink selection", flush=True)
        t_pass1 = time.perf_counter()
        for c0 in range(0, n_cml, int(args.link_batch)):
            c1 = min(c0 + int(args.link_batch), n_cml)
            nb = c1 - c0
            tx_in = _canonical_read(tx_src, time_dim, sublink_dim, cml_dim, c0, c1)
            rx_in = _canonical_read(rx_src, time_dim, sublink_dim, cml_dim, c0, c1)
            _, _, trsl, _, _, _, _, f5_flat, long_frac, short_frac = qc_kernel(tx_in, rx_in, proc_code)

            finite_count = np.isfinite(trsl).sum(axis=1).astype(np.int64)
            finite_fraction = finite_count.astype(np.float64) / float(n_time)
            avail_flat = (
                (finite_fraction >= float(args.min_valid_fraction))
                & (finite_count >= int(args.min_valid_observations))
            ).astype(np.uint8)

            fluct_pass[:, c0:c1] = (f5_flat.reshape(n_sublink, nb) == 0).astype(np.uint8)
            availability_pass[:, c0:c1] = avail_flat.reshape(n_sublink, nb)
            availability_fraction[:, c0:c1] = finite_fraction.reshape(n_sublink, nb).astype(np.float32)
            availability_count[:, c0:c1] = finite_count.reshape(n_sublink, nb)
            long_fraction_all[:, c0:c1] = long_frac.reshape(n_sublink, nb)
            short_fraction_all[:, c0:c1] = short_frac.reshape(n_sublink, nb)

            elapsed = time.perf_counter() - t_pass1
            rate = c1 / elapsed if elapsed > 0 else np.nan
            eta = (n_cml - c1) / rate if rate > 0 else np.nan
            print(f"  selection links {c0:5d}:{c1:5d}/{n_cml} | ETA≈{eta:8.1f}s", flush=True)

        pass_after_fluct = fluct_pass.astype(bool)
        pass_after_detection = pass_after_fluct & det_pass.astype(bool)
        final_sublink_pass = pass_after_detection & availability_pass.astype(bool)

        link_after_fluct = pass_after_fluct.any(axis=0)
        link_after_detection = pass_after_detection.any(axis=0)
        final_link_pass = final_sublink_pass.any(axis=0)
        retained_source_indices = np.flatnonzero(final_link_pass).astype(np.int64)
        n_retained = int(retained_source_indices.size)
        if n_retained == 0:
            raise RuntimeError("No physical links passed the complete QC. Check units, metadata and thresholds.")

        stage_rows = [
            {
                "stage": "Original population",
                "retained_links": int(n_cml),
                "dropped_links_at_stage": 0,
                "retained_sublinks": int(n_sublink * n_cml),
                "dropped_sublinks_at_stage": 0,
            },
            {
                "stage": "After fluctuation screening",
                "retained_links": int(link_after_fluct.sum()),
                "dropped_links_at_stage": int(n_cml - link_after_fluct.sum()),
                "retained_sublinks": int(pass_after_fluct.sum()),
                "dropped_sublinks_at_stage": int(n_sublink * n_cml - pass_after_fluct.sum()),
            },
            {
                "stage": "After detection-limit screening",
                "retained_links": int(link_after_detection.sum()),
                "dropped_links_at_stage": int(link_after_fluct.sum() - link_after_detection.sum()),
                "retained_sublinks": int(pass_after_detection.sum()),
                "dropped_sublinks_at_stage": int(pass_after_fluct.sum() - pass_after_detection.sum()),
            },
            {
                "stage": "After availability screening (final)",
                "retained_links": int(final_link_pass.sum()),
                "dropped_links_at_stage": int(link_after_detection.sum() - final_link_pass.sum()),
                "retained_sublinks": int(final_sublink_pass.sum()),
                "dropped_sublinks_at_stage": int(pass_after_detection.sum() - final_sublink_pass.sum()),
            },
        ]

        print("\n=== LINK/SUBLINK FILTERING ===")
        for row in stage_rows:
            print(
                f"{row['stage']:<46s} "
                f"links={row['retained_links']:5d} "
                f"drop_links={row['dropped_links_at_stage']:5d} | "
                f"sublinks={row['retained_sublinks']:5d} "
                f"drop_sublinks={row['dropped_sublinks_at_stage']:5d}"
            )

        # ------------------------------------------------------------------
        # Create compact output with retained physical links only.
        # ------------------------------------------------------------------
        compression = not args.no_compression
        sizes_out = dict(sizes_src)
        sizes_out[cml_dim] = n_retained
        chunks_native = _chunks_for_dims(
            native_dims,
            sizes_out,
            time_dim,
            sublink_dim,
            cml_dim,
            args.time_chunk,
            min(args.link_batch, args.link_chunk),
        )

        total_valid_observations = 0
        total_invalid_masked = 0
        failed_sublinks_inside_retained = int(
            (~final_sublink_pass[:, retained_source_indices]).sum()
        )

        with Dataset(temp_path, "w", format="NETCDF4") as dst:
            dst.set_fill_off()
            for name, dim in src.dimensions.items():
                if name == cml_dim:
                    dst.createDimension(name, n_retained)
                else:
                    dst.createDimension(name, None if dim.isunlimited() else len(dim))

            dst.setncatts({name: src.getncattr(name) for name in src.ncattrs()})
            dst.setncattr("quality_control_version", "5.0")
            dst.setncattr("quality_control_program", "cml_qc_filter_standalone_v5.py")
            dst.setncattr("quality_control_created_utc", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            dst.setncattr("quality_control_source_file", str(src_path))
            dst.setncattr("quality_control_processing", "Complete six-step per-sublink QC followed by availability screening and strict masking")
            dst.setncattr("quality_control_link_rule", "A physical link is retained when at least one directional sublink passes fluctuation, detection-limit, and availability screening.")
            dst.setncattr("quality_control_sublink_rule", "A failed directional sublink remains on the shared sublink dimension but all of its TX, RX, and TRSL observations are missing.")
            dst.setncattr("quality_control_observation_rule", "An observation is valid only when processed TRSL is finite and its directional sublink passes all sublink-level filters. Invalid values are represented by missing values because different sublinks fail at different times.")
            dst.setncattr("quality_control_original_structure", "Preserved as time x sublink_id x cml_id; the cml_id dimension contains only retained physical links.")
            dst.setncattr("quality_control_frequency_input_unit_inferred", inferred_freq_unit)
            dst.setncattr("quality_control_min_valid_fraction", float(args.min_valid_fraction))
            dst.setncattr("quality_control_min_valid_observations", int(args.min_valid_observations))
            dst.setncattr("quality_control_thresholds_json", json.dumps({
                "tx_min_dbm": TX_MIN_DBM,
                "tx_sentinel": TX_SENTINEL,
                "rx_min_dbm": RX_MIN_DBM,
                "rx_max_dbm": RX_MAX_DBM,
                "plateau_window_samples": PLATEAU_WINDOW,
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
                "rolling_min_periods": ROLLING_MIN_PERIODS,
                "attenuation_quantization_db": A_QUANTIZATION_DB,
                "detection_limit_max_mmph": DETECTION_LIMIT_MAX_MMPH,
                "detection_limit_formula": "0.33/(ITU_R_P838_3_k*length_km), legacy-compatible",
            }, sort_keys=True))

            out_vars = {}
            for name, var in src.variables.items():
                use_chunks = chunks_native if name in (tx_name, rx_name) else None
                out_vars[name] = _create_variable_like(
                    dst,
                    name,
                    var,
                    compression=compression and np.dtype(var.dtype).kind not in "OSU",
                    compression_level=args.compression_level,
                    chunksizes=use_chunks,
                )

            # Copy everything except signal variables; cml-dependent data are compacted.
            for name, var in src.variables.items():
                if name in (tx_name, rx_name):
                    continue
                _copy_selected_variable(
                    var,
                    out_vars[name],
                    cml_dim=cml_dim,
                    time_dim=time_dim,
                    source_indices=retained_source_indices,
                    link_batch=args.link_batch,
                    time_batch=args.time_chunk,
                )

            if args.keep_original_signals:
                tx_raw_var = _new_numeric_var(
                    dst,
                    tx_name + "_raw",
                    str(tx_src.dtype),
                    native_dims,
                    {"long_name": f"Original unmodified {tx_name} for retained links", "source_variable": tx_name},
                    compression=compression,
                    level=args.compression_level,
                    chunksizes=chunks_native,
                    fill_value=_fill_value_of(tx_src),
                )
                rx_raw_var = _new_numeric_var(
                    dst,
                    rx_name + "_raw",
                    str(rx_src.dtype),
                    native_dims,
                    {"long_name": f"Original unmodified {rx_name} for retained links", "source_variable": rx_name},
                    compression=compression,
                    level=args.compression_level,
                    chunksizes=chunks_native,
                    fill_value=_fill_value_of(rx_src),
                )
            else:
                tx_raw_var = rx_raw_var = None

            flag_attrs = {
                "flag_values": np.array([0, 1], dtype=np.uint8),
                "flag_meanings": "not_affected affected",
                "valid_min": np.uint8(0),
                "valid_max": np.uint8(1),
            }
            trsl_var = _new_numeric_var(
                dst,
                "trsl_qc",
                "f4",
                native_dims,
                {"long_name": "Strict quality-controlled total signal loss (TX minus RX)", "units": getattr(tx_src, "units", "dB")},
                compression=compression,
                level=args.compression_level,
                chunksizes=chunks_native,
                fill_value=np.float32(np.nan),
            )
            step1_var = _new_numeric_var(dst, "qc_step1_invalid_signal", "u1", native_dims, {**flag_attrs, "long_name": "QC step 1: invalid TX or RX sample"}, compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))
            step2_var = _new_numeric_var(dst, "qc_step2_plateau", "u1", native_dims, {**flag_attrs, "long_name": "QC step 2: low stable RSL plateau including padding"}, compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))
            step3_var = _new_numeric_var(dst, "qc_step3_blackout_filled", "u1", native_dims, {**flag_attrs, "long_name": "QC step 3: reconstructed short blackout sample"}, compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))
            step4_var = _new_numeric_var(dst, "qc_step4_short_gap_interpolated", "u1", native_dims, {**flag_attrs, "long_name": "QC step 4: linearly interpolated short TRSL gap"}, compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))
            final_valid_var = _new_numeric_var(dst, "qc_final_valid", "u1", native_dims, {"flag_values": np.array([0, 1], dtype=np.uint8), "flag_meanings": "invalid valid", "long_name": "Final observation-level validity after all QC"}, compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint8(255))

            small_dims = (sublink_dim, cml_dim)
            small_chunks = (n_sublink, min(n_retained, max(1, args.link_chunk)))
            step5_var = _new_numeric_var(dst, "qc_step5_fluctuation_fail", "u1", small_dims, {**flag_attrs, "long_name": "QC step 5: sublink failed fluctuation screening"}, compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.uint8(255))
            step6_var = _new_numeric_var(dst, "qc_step6_detection_limit_fail", "u1", small_dims, {**flag_attrs, "long_name": "QC step 6: sublink detection limit missing or >= 2 mm h-1"}, compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.uint8(255))
            step7_var = _new_numeric_var(dst, "qc_step7_availability_fail", "u1", small_dims, {**flag_attrs, "long_name": "QC step 7: sublink failed minimum valid-data availability"}, compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.uint8(255))
            det_var = _new_numeric_var(dst, "qc_detection_limit_mmph", "f4", small_dims, {"long_name": "Legacy-compatible rainfall detection limit", "units": "mm h-1"}, compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.float32(np.nan))
            avail_frac_var = _new_numeric_var(dst, "qc_availability_fraction", "f4", small_dims, {"long_name": "Fraction of timestamps with finite processed total loss before sublink-level rejection", "units": "1"}, compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.float32(np.nan))
            avail_count_var = _new_numeric_var(dst, "qc_availability_count", "i8", small_dims, {"long_name": "Number of timestamps with finite processed total loss before sublink-level rejection", "units": "samples"}, compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.int64(-1))
            long_frac_var = _new_numeric_var(dst, "qc_step5_fraction_long_std_exceed", "f4", small_dims, {"long_name": "Fraction of valid centered 5-hour windows exceeding 2 dB standard deviation", "units": "1"}, compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.float32(np.nan))
            short_frac_var = _new_numeric_var(dst, "qc_step5_fraction_short_std_exceed", "f4", small_dims, {"long_name": "Fraction of valid centered 1-hour windows exceeding 0.8 dB standard deviation", "units": "1"}, compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.float32(np.nan))
            final_sl_var = _new_numeric_var(dst, "qc_final_sublink_pass", "u1", small_dims, {"flag_values": np.array([0, 1], dtype=np.uint8), "flag_meanings": "fail pass", "long_name": "Final directional-sublink QC result"}, compression=compression, level=args.compression_level, chunksizes=small_chunks, fill_value=np.uint8(255))
            final_link_var = _new_numeric_var(dst, "qc_final_link_pass", "u1", (cml_dim,), {"flag_values": np.array([0, 1], dtype=np.uint8), "flag_meanings": "fail pass", "long_name": "Final physical-link QC result; all links in this compact file pass"}, compression=compression, level=args.compression_level, chunksizes=(min(n_retained, max(1, args.link_chunk)),), fill_value=np.uint8(255))
            source_index_var = _new_numeric_var(dst, "source_cml_index", "i4", (cml_dim,), {"long_name": "Zero-based CML index in the original source dataset"}, compression=compression, level=args.compression_level, chunksizes=(min(n_retained, max(1, args.link_chunk)),), fill_value=np.int32(-1))

            if args.write_bitmask:
                bitmask_var = _new_numeric_var(dst, "qc_flag", "u2", native_dims, {"long_name": "Compact QC bitmask", "flag_masks": np.array([1, 2, 4, 8, 16, 32, 64], dtype=np.uint16), "flag_meanings": "step1_invalid_signal step2_plateau step3_blackout_filled step4_short_gap_interpolated step5_fluctuation_fail step6_detection_limit_fail final_invalid"}, compression=compression, level=args.compression_level, chunksizes=chunks_native, fill_value=np.uint16(65535))
            else:
                bitmask_var = None

            _append_ancillary(out_vars[tx_name], ["qc_final_valid", "qc_final_sublink_pass", "qc_final_link_pass"])
            _append_ancillary(out_vars[rx_name], ["qc_final_valid", "qc_final_sublink_pass", "qc_final_link_pass"])
            out_vars[tx_name].setncattr("quality_controlled", "true")
            out_vars[rx_name].setncattr("quality_controlled", "true")
            out_vars[tx_name].setncattr("invalid_observations", "stored as missing values")
            out_vars[rx_name].setncattr("invalid_observations", "stored as missing values")

            # Write precomputed small variables for retained links.
            ri = retained_source_indices
            step5_var[:, :] = (~pass_after_fluct[:, ri]).astype(np.uint8)
            step6_var[:, :] = step6_fail_all[:, ri].astype(np.uint8)
            step7_var[:, :] = (~availability_pass[:, ri].astype(bool)).astype(np.uint8)
            det_var[:, :] = detection[:, ri].astype(np.float32)
            avail_frac_var[:, :] = availability_fraction[:, ri].astype(np.float32)
            avail_count_var[:, :] = availability_count[:, ri].astype(np.int64)
            long_frac_var[:, :] = long_fraction_all[:, ri].astype(np.float32)
            short_frac_var[:, :] = short_fraction_all[:, ri].astype(np.float32)
            final_sl_var[:, :] = final_sublink_pass[:, ri].astype(np.uint8)
            final_link_var[:] = np.ones(n_retained, dtype=np.uint8)
            source_index_var[:] = ri.astype(np.int32)

            print("\nPASS 2/2: writing compact strict-QC dataset", flush=True)
            t_pass2 = time.perf_counter()
            for o0 in range(0, n_retained, int(args.link_batch)):
                o1 = min(o0 + int(args.link_batch), n_retained)
                src_indices = retained_source_indices[o0:o1]
                nb = o1 - o0
                tx_in = _canonical_read_indices(tx_src, time_dim, sublink_dim, cml_dim, src_indices)
                rx_in = _canonical_read_indices(rx_src, time_dim, sublink_dim, cml_dim, src_indices)
                tx_qc, rx_qc, trsl, f1, f2, f3, f4, f5_flat, _, _ = qc_kernel(tx_in, rx_in, proc_code)

                sl_pass_flat = final_sublink_pass[:, src_indices].reshape(-1)
                valid = (np.isfinite(trsl) & (sl_pass_flat[:, None] == 1)).astype(np.uint8)

                # Strict dataset: no scientific value remains at an invalid sample.
                tx_strict = tx_qc.copy()
                rx_strict = rx_qc.copy()
                trsl_strict = trsl.astype(np.float32, copy=True)
                invalid = valid == 0
                tx_strict[invalid] = np.nan
                rx_strict[invalid] = np.nan
                trsl_strict[invalid] = np.nan

                total_valid_observations += int(valid.sum())
                total_invalid_masked += int(valid.size - valid.sum())

                if bitmask_var is not None:
                    f6_flat = step6_fail_all[:, src_indices].reshape(-1).astype(np.uint8)
                    bitmask = (
                        f1.astype(np.uint16) * BIT_STEP1_INVALID
                        | f2.astype(np.uint16) * BIT_STEP2_PLATEAU
                        | f3.astype(np.uint16) * BIT_STEP3_BLACKOUT_FILLED
                        | f4.astype(np.uint16) * BIT_STEP4_INTERPOLATED
                        | f5_flat[:, None].astype(np.uint16) * BIT_STEP5_FLUCTUATION_FAIL
                        | f6_flat[:, None].astype(np.uint16) * BIT_STEP6_DETECTION_FAIL
                        | (1 - valid).astype(np.uint16) * BIT_FINAL_INVALID
                    )
                else:
                    bitmask = None

                out_slice = _batch_slice(native_dims, cml_dim, o0, o1)
                out_vars[tx_name][out_slice] = _canonical_to_native(tx_strict, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                out_vars[rx_name][out_slice] = _canonical_to_native(rx_strict, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                trsl_var[out_slice] = _canonical_to_native(trsl_strict, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                step1_var[out_slice] = _canonical_to_native(f1, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                step2_var[out_slice] = _canonical_to_native(f2, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                step3_var[out_slice] = _canonical_to_native(f3, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                step4_var[out_slice] = _canonical_to_native(f4, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                final_valid_var[out_slice] = _canonical_to_native(valid, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                if bitmask_var is not None:
                    bitmask_var[out_slice] = _canonical_to_native(bitmask, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                if tx_raw_var is not None:
                    tx_raw_var[out_slice] = _canonical_to_native(tx_in, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)
                    rx_raw_var[out_slice] = _canonical_to_native(rx_in, native_dims, time_dim, sublink_dim, cml_dim, n_sublink, nb)

                elapsed = time.perf_counter() - t_pass2
                rate = o1 / elapsed if elapsed > 0 else np.nan
                eta = (n_retained - o1) / rate if rate > 0 else np.nan
                print(f"  writing links {o0:5d}:{o1:5d}/{n_retained} | ETA≈{eta:8.1f}s", flush=True)

            total_possible_output = int(n_time * n_sublink * n_retained)
            removed_with_dropped_links = int(n_time * n_sublink * (n_cml - n_retained))
            summary = {
                "version": "5.0",
                "source_file": str(src_path),
                "output_file": str(out_path),
                "original_links": int(n_cml),
                "retained_links": int(n_retained),
                "dropped_links": int(n_cml - n_retained),
                "original_directional_sublinks": int(n_sublink * n_cml),
                "retained_directional_sublinks": int(final_sublink_pass[:, ri].sum()),
                "failed_sublinks_masked_inside_retained_links": failed_sublinks_inside_retained,
                "observations_removed_with_dropped_links": removed_with_dropped_links,
                "possible_observations_in_output": total_possible_output,
                "valid_observations_retained": int(total_valid_observations),
                "invalid_observations_masked_in_output": int(total_invalid_masked),
                "min_valid_fraction": float(args.min_valid_fraction),
                "min_valid_observations": int(args.min_valid_observations),
                "filtering_stages": stage_rows,
            }
            dst.setncattr("quality_control_summary_json", json.dumps(summary, sort_keys=True))
            dst.sync()

    if out_path.exists():
        out_path.unlink()
    os.replace(temp_path, out_path)
    _write_summary_csv(summary_csv, stage_rows)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    with summary_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    total_elapsed = time.perf_counter() - t_all
    print("\n" + "=" * 92)
    print("STRICT QC DATASET COMPLETE")
    print(f"Retained links                  : {summary['retained_links']:,}/{summary['original_links']:,}")
    print(f"Dropped links                   : {summary['dropped_links']:,}")
    print(f"Retained directional sublinks   : {summary['retained_directional_sublinks']:,}")
    print(f"Failed sublinks masked in output: {summary['failed_sublinks_masked_inside_retained_links']:,}")
    print(f"Valid observations retained     : {summary['valid_observations_retained']:,}")
    print(f"Invalid observations masked     : {summary['invalid_observations_masked_in_output']:,}")
    print(f"Saved NetCDF                    : {out_path}")
    print(f"Saved summary JSON              : {summary_json}")
    print(f"Saved filtering CSV             : {summary_csv}")
    print(f"Elapsed                         : {total_elapsed:.1f} s")
    print("=" * 92)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="cml_qc_filter_standalone_v5.py",
        description=(
            "Compute complete per-sublink CML QC directly from a raw NetCDF, "
            "physically remove failed physical links, preserve the original "
            "time x sublink x cml structure, and mask failed sublinks and invalid observations."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--version", action="version", version="%(prog)s 5.0")
    ap.add_argument("--cml-nc", required=True, help="Original/raw CML NetCDF containing TX, RX and metadata.")
    ap.add_argument("--out-nc", required=True, help="Strictly filtered output NetCDF.")
    ap.add_argument("--summary-json", default=None)
    ap.add_argument("--summary-csv", default=None)

    ap.add_argument("--tx-var", default=None)
    ap.add_argument("--rx-var", default=None)
    ap.add_argument("--frequency-var", default=None)
    ap.add_argument("--polarization-var", default=None)
    ap.add_argument("--length-var", default=None)
    ap.add_argument("--time-dim", default=None)
    ap.add_argument("--sublink-dim", default=None)
    ap.add_argument("--cml-dim", default=None)

    ap.add_argument("--min-valid-fraction", type=float, default=0.5, help="Minimum fraction of finite processed TRSL samples required per directional sublink.")
    ap.add_argument("--min-valid-observations", type=int, default=1, help="Minimum number of finite processed TRSL observations required per directional sublink.")
    ap.add_argument("--n-workers", type=int, default=max(1, min(16, os.cpu_count() or 1)), help="Numba compute threads.")
    ap.add_argument("--link-batch", type=int, default=16)
    ap.add_argument("--link-chunk", type=int, default=8)
    ap.add_argument("--time-chunk", type=int, default=10080)
    ap.add_argument("--compression-level", type=int, choices=range(0, 10), default=1)
    ap.add_argument("--no-compression", action="store_true")
    ap.add_argument("--keep-original-signals", action="store_true", help="Also store raw TX/RX for retained links; increases file size substantially.")
    ap.add_argument("--write-bitmask", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    if not 0.0 <= args.min_valid_fraction <= 1.0:
        raise SystemExit("--min-valid-fraction must be between 0 and 1")
    if args.min_valid_observations < 0:
        raise SystemExit("--min-valid-observations must be >= 0")
    try:
        process_file(args)
    except KeyboardInterrupt:
        print("Interrupted. Any .partial output can be removed safely.", file=sys.stderr)
        raise SystemExit(130)


if __name__ == "__main__":
    main()
