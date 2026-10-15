#!/usr/bin/env python3
"""
End-to-end publishable GRAF/RSTD CML rainfall product and evaluation pipeline.

Purpose
-------
1. Rebuild a GRAF-conditioned CML rain-rate product for every CML link that
   passes the original Blettner link-level filter. Links are NOT selected by
   radar availability, so links with no radar reference remain available for
   rain-gauge comparison.
2. Keep the two directional sublinks as independent observations. A failed
   sublink is masked; sublinks are never averaged in the validation.
3. Mask every native one-minute observation that does not pass the original
   QC validity flag. The common time axis is preserved because a rectangular
   NetCDF cannot physically delete a different set of timestamps per link.
4. Evaluate GRAF wet/dry detection against path-associated radar at 5 minutes
   and the closest rain gauge at 10 minutes using contingency tables.
5. Evaluate GRAF-conditioned rain rates with pooled-sublink hexbin plots.

Reference units
---------------
* Radar ``rr`` is assumed to be a 5-minute amount in mm and is multiplied by
  12 to obtain mm h-1. Use ``--radar-scale 1`` if it is already a rate.
* Gauge ``rr`` is assumed to be a 10-minute amount in mm and is multiplied by
  6 to obtain mm h-1. Use ``--gauge-scale 1`` if it is already a rate.
* CML rain rates returned by the k-R inversion are already in mm h-1. They are
  averaged, not multiplied by 60, when moved to 5- or 10-minute resolution.

Fair-comparison rules
---------------------
* The filtered V2 CML product is independent of reference availability.
* Radar and gauge populations are formed independently.
* Only exact, complete native CML intervals are used by default: all five
  one-minute samples for radar and all ten for gauges must be valid.
* Missing CML or reference values are excluded, never recoded as dry.
* GRAF interval wetness defaults to ``any``: an interval is CML wet when at
  least one valid native minute is wet. This is configurable.
* Contingency tables use all paired wet and dry intervals.
* Hexbins default to reference-wet intervals (reference >= 0.1 mm h-1), while
  retaining CML zeros so missed rainfall affects the metrics.
* Hexbin metrics use all finite paired values, including values above the
  displayed 150 mm h-1 limit.

V2 NetCDF representation
------------------------
The output contains only links with ``qc_final_link_pass == 1``. Because the
same sublink dimension is shared by all links, a failed sublink cannot be
removed independently without a ragged representation. It is therefore kept
as a coordinate position but all of its time-dependent observations are
masked. The same applies to failed one-minute observations: the common time
axis is retained and failed values are stored as missing with explicit flags.

Default scientific configuration
--------------------------------
* Wet/dry source: graf_rstd_wet_1min and graf_rstd_valid_1min
* Native QC: qc_final_link_pass, qc_final_sublink_pass, qc_final_valid
* Baseline: pycomlink baseline_constant, one most recent dry value
* WAA: Pastorek et al. (2021), A_max=14, zeta=0.55, d=0.1
* k-R: ITU_2005 (ITU_2003 is also accepted)

The script is streamed by link batches. Rainy paired values are written to
small temporary binary files and memory-mapped for plotting, avoiding a large
in-memory pair table.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import tempfile
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import xarray as xr

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    def tqdm(x: Iterable, **_: object):
        return x

SCRIPT_VERSION = "2.0"
ALLOWED_KR_MODELS = ("ITU_2003", "ITU_2005")


# =============================================================================
# Generic helpers
# =============================================================================
def _decode_scalar(value: object) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return bytes(value).decode("utf-8", errors="ignore").replace("\x00", "").strip()
    return str(value).replace("\x00", "").strip()


def _decode_ids(values: np.ndarray) -> np.ndarray:
    a = np.asarray(values)
    if a.ndim == 2 and a.dtype.kind in "SU":
        # Fixed-width character array. Join the final axis.
        out = []
        for row in a:
            out.append("".join(_decode_scalar(v) for v in row).strip())
        return np.asarray(out, dtype=object)
    return np.asarray([_decode_scalar(v) for v in a.reshape(-1)], dtype=object)


def _get_ids(ds: xr.Dataset, dim: str) -> np.ndarray:
    candidates = (dim, "cml_id", "link_id", "station_id", "WSI", "id")
    for name in candidates:
        if name in ds.variables or name in ds.coords:
            da = ds[name]
            if dim in da.dims:
                arr = np.asarray(da.values)
                ids = _decode_ids(arr)
                if len(ids) == ds.sizes[dim]:
                    return ids
    return np.asarray([str(i) for i in range(ds.sizes[dim])], dtype=object)


def _find_dim(dims: Iterable[str], requested: Optional[str], candidates: Sequence[str], label: str) -> str:
    dims = tuple(dims)
    if requested:
        if requested not in dims:
            raise KeyError(f"Requested {label} dimension {requested!r} not found in {dims}")
        return requested
    for c in candidates:
        if c in dims:
            return c
    raise KeyError(f"Could not infer {label} dimension from {dims}; tried {candidates}")


def _find_var(ds: xr.Dataset, requested: Optional[str], candidates: Sequence[str], label: str) -> str:
    if requested:
        if requested not in ds.variables:
            raise KeyError(f"Requested {label} variable {requested!r} not found")
        return requested
    for c in candidates:
        if c in ds.variables:
            return c
    raise KeyError(f"Could not infer {label}; tried {candidates}")


def _normalize_times(values: np.ndarray, label: str) -> pd.DatetimeIndex:
    idx = pd.DatetimeIndex(values)
    if idx.hasnans:
        raise ValueError(f"{label} contains invalid timestamps")
    if not idx.is_monotonic_increasing:
        raise ValueError(f"{label} must be monotonically increasing")
    if not idx.is_unique:
        raise ValueError(f"{label} contains duplicated timestamps")
    return idx


def _time_slice_indices(times: pd.DatetimeIndex, tstart: Optional[str], tend: Optional[str]) -> tuple[int, int]:
    ns = times.asi8
    left = 0 if tstart is None else int(np.searchsorted(ns, pd.Timestamp(tstart).value, side="left"))
    right = len(times) if tend is None else int(np.searchsorted(ns, pd.Timestamp(tend).value, side="right"))
    if left >= right:
        raise ValueError("Requested time period contains no CML samples")
    return left, right


def _drop_singleton_extra_dims(da: xr.DataArray, keep: set[str], label: str) -> xr.DataArray:
    work = da
    for d in list(work.dims):
        if d in keep:
            continue
        if work.sizes[d] != 1:
            raise ValueError(f"{label} has unsupported extra dimension {d!r} of size {work.sizes[d]}")
        work = work.isel({d: 0}, drop=True)
    return work


def pol_to_HV(pol_raw: object) -> str:
    if isinstance(pol_raw, (bytes, np.bytes_)):
        pol_raw = bytes(pol_raw).decode("utf-8", errors="ignore")
    s = str(pol_raw).strip().lower()
    if "v" in s and "h" not in s:
        return "V"
    if "h" in s and "v" not in s:
        return "H"
    return "H"


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
    if float(np.nanmedian(finite)) > 100.0:
        return x / 1000.0, "m"
    return x, "km"


def _numeric_metadata_2d(
    da: xr.DataArray,
    sublink_dim: str,
    cml_dim: str,
    cml_positions: np.ndarray,
    nsub: int,
) -> np.ndarray:
    work = da
    if cml_dim not in work.dims:
        raise ValueError(f"{da.name}.dims={da.dims}; expected {cml_dim!r}")
    work = work.isel({cml_dim: cml_positions})
    for d in list(work.dims):
        if d not in {sublink_dim, cml_dim}:
            if work.sizes[d] != 1:
                raise ValueError(f"Numeric metadata {da.name!r} has unsupported dimension {d!r}")
            work = work.isel({d: 0}, drop=True)
    if sublink_dim not in work.dims:
        arr = np.asarray(work.transpose(cml_dim).values, dtype=np.float64)[None, :]
        return np.broadcast_to(arr, (nsub, len(cml_positions))).copy()
    return np.asarray(work.transpose(sublink_dim, cml_dim).values, dtype=np.float64)


def _string_metadata_2d(
    da: xr.DataArray,
    sublink_dim: str,
    cml_dim: str,
    cml_positions: np.ndarray,
    nsub: int,
) -> np.ndarray:
    work = da.isel({cml_dim: cml_positions}) if cml_dim in da.dims else da
    stored_dims = list(work.dims)
    arr = np.asarray(work.values)

    # Character axis may be an extra dimension. Arrange signal dimensions first.
    signal_dims = []
    if sublink_dim in stored_dims:
        signal_dims.append(sublink_dim)
    if cml_dim in stored_dims:
        signal_dims.append(cml_dim)
    extra_dims = [d for d in stored_dims if d not in signal_dims]
    desired = signal_dims + extra_dims
    axes = [stored_dims.index(d) for d in desired]
    if axes != list(range(len(axes))):
        arr = np.transpose(arr, axes)

    if sublink_dim not in signal_dims:
        # cml x optional-char
        if extra_dims:
            rows = arr.reshape(len(cml_positions), -1)
            vals = np.asarray(["".join(_decode_scalar(v) for v in row) for row in rows], dtype=object)
        else:
            vals = np.asarray([_decode_scalar(v) for v in arr.reshape(-1)], dtype=object)
        return np.broadcast_to(vals[None, :], (nsub, len(cml_positions))).copy()

    # sublink x cml x optional-char
    if extra_dims:
        rows = arr.reshape(nsub * len(cml_positions), -1)
        vals = np.asarray(["".join(_decode_scalar(v) for v in row) for row in rows], dtype=object)
        return vals.reshape(nsub, len(cml_positions))
    return np.asarray([_decode_scalar(v) for v in arr.reshape(-1)], dtype=object).reshape(nsub, len(cml_positions))


def _read_static_link_flag(ds: xr.Dataset, name: str, cml_dim: str) -> np.ndarray:
    if name not in ds.variables:
        raise KeyError(f"Required link flag {name!r} not found")
    da = ds[name]
    if cml_dim not in da.dims:
        raise ValueError(f"{name}.dims={da.dims}; expected {cml_dim!r}")
    work = da
    for d in list(work.dims):
        if d == cml_dim:
            continue
        work = work.max(dim=d, skipna=True)
    arr = np.asarray(work.transpose(cml_dim).values, dtype=np.float64)
    return np.isfinite(arr) & (arr > 0.5)


def _read_static_sublink_flag(
    ds: xr.Dataset, name: str, sublink_dim: str, cml_dim: str, positions: np.ndarray
) -> np.ndarray:
    if name not in ds.variables:
        raise KeyError(f"Required sublink flag {name!r} not found")
    da = ds[name].isel({cml_dim: positions})
    for d in list(da.dims):
        if d not in {sublink_dim, cml_dim}:
            da = da.max(dim=d, skipna=True)
    if sublink_dim not in da.dims:
        arr = np.asarray(da.transpose(cml_dim).values, dtype=np.float64)[None, :]
        arr = np.broadcast_to(arr, (ds.sizes[sublink_dim], len(positions)))
    else:
        arr = np.asarray(da.transpose(sublink_dim, cml_dim).values, dtype=np.float64)
    return np.isfinite(arr) & (arr > 0.5)


def _read_signal_block(
    da: xr.DataArray,
    time_dim: str,
    sublink_dim: str,
    cml_dim: str,
    time_slice: slice,
    positions: np.ndarray,
    dtype: np.dtype,
) -> np.ndarray:
    work = da.isel({time_dim: time_slice, cml_dim: positions})
    work = _drop_singleton_extra_dims(work, {time_dim, sublink_dim, cml_dim}, da.name or "signal")
    if sublink_dim not in work.dims:
        raise ValueError(f"{da.name}.dims={da.dims}; expected sublink dimension {sublink_dim!r}")
    return np.asarray(work.transpose(time_dim, sublink_dim, cml_dim).values, dtype=dtype)


# =============================================================================
# pycomlink compatibility and GRAF rain retrieval
# =============================================================================
def _call_with_kr_model(fn, model: str, **kwargs):
    errors = []
    for keyword in ("a_b_approximation", "approx_type"):
        try:
            return fn(**kwargs, **{keyword: model})
        except (TypeError, ValueError) as exc:
            errors.append(f"{keyword}: {type(exc).__name__}: {exc}")
    raise RuntimeError(f"pycomlink rejected {model}: {' | '.join(errors)}")


def _calc_rain_from_a(pycml, A: xr.DataArray, L_km: float, f_GHz: float, pol: str, model: str):
    return _call_with_kr_model(
        pycml.processing.k_R_relation.calc_R_from_A,
        model,
        A=A,
        L_km=float(L_km),
        f_GHz=float(f_GHz),
        pol=pol_to_HV(pol),
        R_min=0,
    )


def _estimate_one_series(
    pycml,
    trsl: np.ndarray,
    wet: np.ndarray,
    valid: np.ndarray,
    times: np.ndarray,
    length_km: float,
    frequency_hz: float,
    polarization: str,
    kr_model: str,
    n_average_last_dry: int,
    A_max: float,
    zeta: float,
    d: float,
) -> np.ndarray:
    out = np.full(trsl.shape, np.nan, dtype=np.float32)
    valid = np.asarray(valid, dtype=bool) & np.isfinite(trsl)
    if not valid.any():
        return out
    wet = np.asarray(wet, dtype=bool) & valid
    dry = valid & ~wet
    out[dry] = 0.0
    if not wet.any():
        return out

    pol = pol_to_HV(polarization)
    xtrsl = xr.DataArray(
        np.where(valid, trsl, np.nan),
        dims=("time",),
        coords={"time": np.asarray(times)},
    )
    xwet = xr.DataArray(wet, dims=("time",), coords={"time": np.asarray(times)})
    baseline = pycml.processing.baseline.baseline_constant(
        trsl=xtrsl,
        wet=xwet,
        n_average_last_dry=int(n_average_last_dry),
    )
    A_obs = xr.where((xtrsl - baseline) > 0, xtrsl - baseline, 0.0)
    waa = pycml.processing.wet_antenna.waa_pastorek_2021_from_A_obs(
        A_obs=A_obs,
        f_Hz=float(frequency_hz),
        pol=pol,
        L_km=float(length_km),
        A_max=float(A_max),
        zeta=float(zeta),
        d=float(d),
    )
    # Match the user's own implementation: physically constrain WAA.
    waa = xr.where(waa < 0, 0.0, waa)
    waa = xr.where(waa > A_obs, A_obs, waa)
    A_rain = xr.where((A_obs - waa) > 0, A_obs - waa, 0.0)
    R = _calc_rain_from_a(
        pycml,
        A=A_rain,
        L_km=float(length_km),
        f_GHz=float(frequency_hz) / 1e9,
        pol=pol,
        model=kr_model,
    )
    values = np.asarray(getattr(R, "values", R), dtype=np.float64)
    good = wet & np.isfinite(values) & (values >= 0)
    out[good] = values[good].astype(np.float32)
    return out


# =============================================================================
# V2 NetCDF writer
# =============================================================================
def _compression(level: int, chunks: tuple[int, ...]) -> dict[str, object]:
    return {
        "zlib": int(level) > 0,
        "complevel": int(level),
        "shuffle": True,
        "chunksizes": tuple(max(1, int(v)) for v in chunks),
    }


@dataclass
class V2Handles:
    nc: object
    R: object
    wet: object
    wet_valid: object
    qc_valid: object
    trsl: Optional[object]


def _create_v2_output(
    output: Path,
    args: argparse.Namespace,
    times: pd.DatetimeIndex,
    selected_ids: np.ndarray,
    selected_positions: np.ndarray,
    sublink_values: np.ndarray,
    sublink_pass: np.ndarray,
    length_km: np.ndarray,
    frequency_input: np.ndarray,
    frequency_hz: np.ndarray,
    polarization: np.ndarray,
    lon0: np.ndarray,
    lat0: np.ndarray,
    lon1: np.ndarray,
    lat1: np.ndarray,
    source_cml_file: Path,
    input_link_count: int,
    link_drop_count: int,
) -> tuple[V2Handles, Path]:
    try:
        from netCDF4 import Dataset, date2num
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("netCDF4 is required for streamed V2 output") from exc

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists: {output}; add --overwrite")
    partial = output.with_suffix(output.suffix + ".partial")
    if partial.exists():
        partial.unlink()

    nt = len(times)
    ns = len(sublink_values)
    nl = len(selected_ids)
    nc = Dataset(partial, "w", format="NETCDF4")
    nc.createDimension("time", nt)
    nc.createDimension("sublink_id", ns)
    nc.createDimension("cml_id", nl)

    vt = nc.createVariable("time", "i8", ("time",))
    vt.units = "seconds since 1970-01-01 00:00:00 UTC"
    vt.calendar = "standard"
    vt[:] = np.asarray(date2num(times.to_pydatetime(), vt.units, vt.calendar), dtype=np.int64)

    vc = nc.createVariable("cml_id", str, ("cml_id",))
    vc[:] = np.asarray(selected_ids, dtype=object)
    vs = nc.createVariable("sublink_id", str, ("sublink_id",))
    vs[:] = np.asarray([_decode_scalar(v) for v in sublink_values], dtype=object)
    vsrc = nc.createVariable("source_cml_index", "i4", ("cml_id",))
    vsrc[:] = selected_positions.astype(np.int32)

    chunks3 = (min(nt, int(args.time_chunk_1min)), 1, min(nl, int(args.output_link_chunk)))
    chunks2 = (1, min(nl, int(args.output_link_chunk)))
    chunk_link = (min(nl, int(args.output_link_chunk)),)

    R = nc.createVariable(
        "R_cml", "f4", ("time", "sublink_id", "cml_id"),
        fill_value=np.float32(np.nan), **_compression(args.compression_level, chunks3)
    )
    R.units = "mm h-1"
    R.long_name = "GRAF-conditioned CML rain rate after complete native QC masking"

    wet = nc.createVariable(
        "graf_wet", "u1", ("time", "sublink_id", "cml_id"),
        fill_value=np.uint8(255), **_compression(args.compression_level, chunks3)
    )
    wet.flag_values = np.asarray([0, 1], dtype=np.uint8)
    wet.flag_meanings = "dry wet"
    wet.long_name = "GRAF/RSTD one-minute wet-dry classification"

    wet_valid = nc.createVariable(
        "graf_wet_valid", "u1", ("time", "sublink_id", "cml_id"),
        **_compression(args.compression_level, chunks3)
    )
    wet_valid.flag_values = np.asarray([0, 1], dtype=np.uint8)
    wet_valid.flag_meanings = "invalid valid"

    qv = nc.createVariable(
        "qc_final_valid", "u1", ("time", "sublink_id", "cml_id"),
        **_compression(args.compression_level, chunks3)
    )
    qv.flag_values = np.asarray([0, 1], dtype=np.uint8)
    qv.flag_meanings = "invalid valid"
    qv.long_name = "Original one-minute final QC validity after link and sublink masking"

    tr = None
    if args.write_trsl:
        tr = nc.createVariable(
            "trsl_qc", "f4", ("time", "sublink_id", "cml_id"),
            fill_value=np.float32(np.nan), **_compression(args.compression_level, chunks3)
        )
        tr.units = "dB"
        tr.long_name = "Quality-controlled total received signal loss"

    vlp = nc.createVariable("qc_final_link_pass", "u1", ("cml_id",), **_compression(args.compression_level, chunk_link))
    vlp[:] = np.ones(nl, dtype=np.uint8)
    vlp.flag_values = np.asarray([0, 1], dtype=np.uint8)
    vlp.flag_meanings = "fail pass"

    vsp = nc.createVariable(
        "qc_final_sublink_pass", "u1", ("sublink_id", "cml_id"),
        **_compression(args.compression_level, chunks2)
    )
    vsp[:] = sublink_pass.astype(np.uint8)
    vsp.flag_values = np.asarray([0, 1], dtype=np.uint8)
    vsp.flag_meanings = "fail pass"

    vlen = nc.createVariable("length_km", "f4", ("sublink_id", "cml_id"), **_compression(args.compression_level, chunks2))
    vlen[:] = length_km.astype(np.float32)
    vlen.units = "km"
    vfin = nc.createVariable("frequency_input", "f8", ("sublink_id", "cml_id"), **_compression(args.compression_level, chunks2))
    vfin[:] = frequency_input.astype(np.float64)
    vfghz = nc.createVariable("frequency_GHz", "f4", ("sublink_id", "cml_id"), **_compression(args.compression_level, chunks2))
    vfghz[:] = (frequency_hz / 1e9).astype(np.float32)
    vfghz.units = "GHz"
    vpol = nc.createVariable("polarization", str, ("sublink_id", "cml_id"))
    vpol[:] = polarization.astype(object)

    for name, vals in (
        ("site_a_longitude", lon0),
        ("site_a_latitude", lat0),
        ("site_b_longitude", lon1),
        ("site_b_latitude", lat1),
    ):
        v = nc.createVariable(name, "f8", ("cml_id",), **_compression(args.compression_level, chunk_link))
        v[:] = np.asarray(vals, dtype=np.float64)
        v.units = "degrees_east" if "longitude" in name else "degrees_north"

    nc.setncattr("title", "Blettner-filtered GRAF-conditioned CML rainfall product")
    nc.setncattr("script", Path(__file__).name)
    nc.setncattr("script_version", SCRIPT_VERSION)
    nc.setncattr("created_utc", pd.Timestamp.now(tz="UTC").isoformat())
    nc.setncattr("source_cml_file", str(source_cml_file))
    nc.setncattr("reference_availability_used_for_link_selection", "false")
    nc.setncattr("link_selection", f"Only links with {args.link_pass_var} == 1")
    nc.setncattr("sublink_selection", f"Failed {args.sublink_pass_var} entries retained only as masked coordinate positions")
    nc.setncattr("observation_selection", f"Only samples with {args.qc_valid_var} == 1, {args.graf_valid_var} == 1, and finite {args.trsl_var} are usable")
    nc.setncattr("rectangular_dataset_note", "A shared time and sublink grid is retained; rejected sublinks and observations are represented by missing data and explicit validity flags")
    nc.setncattr("input_link_count", int(input_link_count))
    nc.setncattr("retained_link_count", int(nl))
    nc.setncattr("dropped_link_count", int(link_drop_count))
    nc.setncattr("wet_dry_method", "GRAF/RSTD flags supplied by source dataset")
    nc.setncattr("baseline_method", "pycomlink baseline_constant")
    nc.setncattr("n_average_last_dry", int(args.n_average_last_dry))
    nc.setncattr("waa_model", "Pastorek_2021")
    nc.setncattr("pastorek_A_max", float(args.pastorek_a_max))
    nc.setncattr("pastorek_zeta", float(args.pastorek_zeta))
    nc.setncattr("pastorek_d", float(args.pastorek_d))
    nc.setncattr("kr_model", args.kr_model)
    nc.setncattr("rain_rate_units", "mm h-1")
    nc.setncattr("rain_rate_formula", "R=(max(A_obs-WAA,0)/(a*L_km))**(1/b)")
    nc.setncattr("temporal_resolution", "native one-minute CML resolution")
    nc.setncattr("temporal_resampling", "none in product; reference-specific averaging occurs only during evaluation")
    nc.setncattr("dry_value_policy", "Valid GRAF-dry samples are 0 mm h-1; invalid samples are missing")

    return V2Handles(nc, R, wet, wet_valid, qv, tr), partial


def _finalize_v2(handles: V2Handles, partial: Path, output: Path) -> None:
    handles.nc.sync()
    handles.nc.close()
    if output.exists():
        output.unlink()
    partial.replace(output)


# =============================================================================
# Exact temporal mapping and aggregation
# =============================================================================
@dataclass
class IntervalMap:
    reference_times: pd.DatetimeIndex
    starts: pd.DatetimeIndex
    indices: np.ndarray
    complete: np.ndarray
    expected_count: int
    minutes: int


def _reference_interval_starts(times: pd.DatetimeIndex, minutes: int, label: str) -> pd.DatetimeIndex:
    delta = pd.Timedelta(minutes=minutes)
    if label == "start":
        return times
    if label == "end":
        return times - delta
    return times - delta / 2


def _build_interval_map(
    cml_times: pd.DatetimeIndex,
    reference_times: pd.DatetimeIndex,
    minutes: int,
    timestamp_label: str,
) -> IntervalMap:
    starts = _reference_interval_starts(reference_times, minutes, timestamp_label)
    offsets = np.arange(minutes, dtype=np.int64) * int(60e9)
    wanted = starts.asi8[:, None] + offsets[None, :]
    flat = cml_times.get_indexer(pd.DatetimeIndex(wanted.reshape(-1)))
    indices = flat.reshape(len(starts), minutes).astype(np.int64)
    complete = np.all(indices >= 0, axis=1)
    return IntervalMap(reference_times, starts, indices, complete, minutes, minutes)


def _aggregate_native(
    rain_1min: np.ndarray,
    wet_1min: np.ndarray,
    valid_1min: np.ndarray,
    mapping: IntervalMap,
    wet_aggregation: str,
    min_valid_fraction: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Aggregate time x sublink x link arrays to exact reference intervals."""
    nref = len(mapping.reference_times)
    _, ns, nl = rain_1min.shape
    safe = np.where(mapping.indices >= 0, mapping.indices, 0)
    values = rain_1min[safe, :, :]          # nref x expected x sub x link
    wets = wet_1min[safe, :, :]
    valids = valid_1min[safe, :, :]
    valids &= mapping.complete[:, None, None, None]
    valids &= np.isfinite(values)
    count = np.sum(valids, axis=1)
    required = int(math.ceil(mapping.expected_count * float(min_valid_fraction) - 1e-12))
    interval_valid = count >= required

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        masked = np.where(valids, values, np.nan)
        rain = np.nanmean(masked, axis=1).astype(np.float32)

    wet_count = np.sum(wets & valids, axis=1)
    if wet_aggregation == "any":
        wet = wet_count > 0
    elif wet_aggregation == "majority":
        wet = wet_count > (count / 2.0)
    elif wet_aggregation == "all":
        wet = wet_count == count
    else:  # pragma: no cover
        raise ValueError(wet_aggregation)
    wet &= interval_valid
    rain[~interval_valid] = np.nan
    return rain, wet, interval_valid


# =============================================================================
# Gauge geometry and reference loading
# =============================================================================
def _infer_station_dim(da: xr.DataArray, time_dim: str, explicit: Optional[str]) -> str:
    if explicit:
        if explicit not in da.dims:
            raise KeyError(f"Gauge station dimension {explicit!r} not found in {da.dims}")
        return explicit
    candidates = [d for d in da.dims if d != time_dim]
    for preferred in ("WSI", "station", "station_id", "gauge", "site", "id"):
        if preferred in candidates:
            return preferred
    if len(candidates) != 1:
        raise ValueError(f"Cannot infer gauge station dimension from {da.dims}")
    return candidates[0]


def _read_numeric_for_dim(ds: xr.Dataset, name: str, dim: str) -> np.ndarray:
    if name not in ds.variables and name not in ds.coords:
        raise KeyError(f"Coordinate {name!r} not found")
    da = ds[name]
    if dim not in da.dims:
        raise ValueError(f"{name}.dims={da.dims}; expected {dim!r}")
    work = da
    for d in list(work.dims):
        if d != dim:
            work = work.mean(dim=d, skipna=True)
    return np.asarray(work.transpose(dim).values, dtype=np.float64).reshape(-1)


def _assign_closest_gauges(
    cml_ids: np.ndarray,
    lon0: np.ndarray,
    lat0: np.ndarray,
    lon1: np.ndarray,
    lat1: np.ndarray,
    gauge_ids: np.ndarray,
    gauge_lon: np.ndarray,
    gauge_lat: np.ndarray,
    projection: str,
    max_distance_m: float,
    batch_size: int,
) -> pd.DataFrame:
    try:
        from pyproj import Transformer
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("pyproj is required for closest-gauge assignment") from exc

    valid_g = np.isfinite(gauge_lon) & np.isfinite(gauge_lat)
    if not valid_g.any():
        raise RuntimeError("No gauge has finite coordinates")
    gauge_positions = np.flatnonzero(valid_g)
    transformer = Transformer.from_crs("EPSG:4326", projection, always_xy=True)
    gx, gy = transformer.transform(gauge_lon[valid_g], gauge_lat[valid_g])
    gx = np.asarray(gx, dtype=np.float64)
    gy = np.asarray(gy, dtype=np.float64)

    valid_c = np.isfinite(lon0) & np.isfinite(lat0) & np.isfinite(lon1) & np.isfinite(lat1)
    x0 = np.full(len(cml_ids), np.nan)
    y0 = np.full(len(cml_ids), np.nan)
    x1 = np.full(len(cml_ids), np.nan)
    y1 = np.full(len(cml_ids), np.nan)
    if valid_c.any():
        x0[valid_c], y0[valid_c] = transformer.transform(lon0[valid_c], lat0[valid_c])
        x1[valid_c], y1[valid_c] = transformer.transform(lon1[valid_c], lat1[valid_c])

    nearest = np.full(len(cml_ids), -1, dtype=np.int64)
    dist_out = np.full(len(cml_ids), np.nan, dtype=np.float64)
    positions = np.flatnonzero(valid_c)
    for b0 in range(0, len(positions), int(batch_size)):
        idx = positions[b0:b0 + int(batch_size)]
        ax, ay = x0[idx][:, None], y0[idx][:, None]
        bx, by = x1[idx][:, None], y1[idx][:, None]
        dx, dy = bx - ax, by - ay
        denom = dx * dx + dy * dy
        px, py = gx[None, :] - ax, gy[None, :] - ay
        t = np.divide(px * dx + py * dy, denom, out=np.zeros_like(px), where=denom > 0)
        t = np.clip(t, 0.0, 1.0)
        qx, qy = ax + t * dx, ay + t * dy
        distances = np.hypot(gx[None, :] - qx, gy[None, :] - qy)
        local = np.argmin(distances, axis=1)
        nearest[idx] = gauge_positions[local]
        dist_out[idx] = distances[np.arange(len(idx)), local]

    accepted = (nearest >= 0) & np.isfinite(dist_out) & (dist_out <= float(max_distance_m))
    rows = []
    for i, cid in enumerate(cml_ids):
        gp = int(nearest[i]) if accepted[i] else -1
        rows.append({
            "cml_id": str(cid),
            "cml_position": i,
            "gauge_id": str(gauge_ids[gp]) if gp >= 0 else "",
            "gauge_position": gp,
            "distance_m": float(dist_out[i]) if np.isfinite(dist_out[i]) else np.nan,
            "within_max_distance": int(gp >= 0),
        })
    return pd.DataFrame(rows)


# =============================================================================
# Streaming metrics, contingency, and pair storage
# =============================================================================
@dataclass
class ContinuousAccumulator:
    n: int = 0
    sum_error: float = 0.0
    sum_sq_error: float = 0.0
    sum_x: float = 0.0
    sum_y: float = 0.0
    sum_x2: float = 0.0
    sum_y2: float = 0.0
    sum_xy: float = 0.0

    def update(self, reference: np.ndarray, cml: np.ndarray) -> None:
        x = np.asarray(reference, dtype=np.float64).reshape(-1)
        y = np.asarray(cml, dtype=np.float64).reshape(-1)
        m = np.isfinite(x) & np.isfinite(y) & (x >= 0) & (y >= 0)
        if not m.any():
            return
        x, y = x[m], y[m]
        e = y - x
        self.n += int(len(x))
        self.sum_error += float(np.sum(e, dtype=np.float64))
        self.sum_sq_error += float(np.sum(e * e, dtype=np.float64))
        self.sum_x += float(np.sum(x, dtype=np.float64))
        self.sum_y += float(np.sum(y, dtype=np.float64))
        self.sum_x2 += float(np.sum(x * x, dtype=np.float64))
        self.sum_y2 += float(np.sum(y * y, dtype=np.float64))
        self.sum_xy += float(np.sum(x * y, dtype=np.float64))

    def result(self) -> dict[str, float | int]:
        if self.n == 0:
            return {"n": 0, "rmse_mmph": np.nan, "mbe_mmph": np.nan, "pcc": np.nan}
        n = float(self.n)
        rmse = math.sqrt(self.sum_sq_error / n)
        mbe = self.sum_error / n
        cov = self.sum_xy - self.sum_x * self.sum_y / n
        vx = self.sum_x2 - self.sum_x * self.sum_x / n
        vy = self.sum_y2 - self.sum_y * self.sum_y / n
        pcc = cov / math.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan
        return {"n": self.n, "rmse_mmph": rmse, "mbe_mmph": mbe, "pcc": pcc}


@dataclass
class ContingencyAccumulator:
    hit: int = 0
    false_alarm: int = 0
    miss: int = 0
    correct_negative: int = 0

    def update(self, predicted_wet: np.ndarray, reference_wet: np.ndarray, valid: np.ndarray) -> None:
        p = np.asarray(predicted_wet, dtype=bool)
        r = np.asarray(reference_wet, dtype=bool)
        v = np.asarray(valid, dtype=bool)
        self.hit += int(np.count_nonzero(v & p & r))
        self.false_alarm += int(np.count_nonzero(v & p & ~r))
        self.miss += int(np.count_nonzero(v & ~p & r))
        self.correct_negative += int(np.count_nonzero(v & ~p & ~r))

    def counts(self) -> dict[str, int]:
        return {
            "hit": self.hit,
            "false_alarm": self.false_alarm,
            "miss": self.miss,
            "correct_negative": self.correct_negative,
        }

    def metrics(self) -> dict[str, float | int]:
        h, f, m, n = self.hit, self.false_alarm, self.miss, self.correct_negative
        total = h + f + m + n
        pod = h / (h + m) if h + m else np.nan
        far = f / (h + f) if h + f else np.nan
        denom = math.sqrt((h + f) * (h + m) * (n + f) * (n + m))
        mcc = (h * n - f * m) / denom if denom else np.nan
        return {"N": total, "POD": pod, "FAR": far, "MCC": mcc}


class PairStore:
    def __init__(self, directory: Path, stem: str):
        directory.mkdir(parents=True, exist_ok=True)
        self.ref_path = directory / f"{stem}_reference.float32.bin"
        self.cml_path = directory / f"{stem}_cml.float32.bin"
        self.ref_f = open(self.ref_path, "wb")
        self.cml_f = open(self.cml_path, "wb")
        self.n = 0

    def append(self, ref: np.ndarray, cml: np.ndarray) -> None:
        x = np.asarray(ref, dtype=np.float32).reshape(-1)
        y = np.asarray(cml, dtype=np.float32).reshape(-1)
        if len(x) != len(y):
            raise ValueError("Pair lengths differ")
        if len(x) == 0:
            return
        x.tofile(self.ref_f)
        y.tofile(self.cml_f)
        self.n += int(len(x))

    def close(self) -> None:
        if not self.ref_f.closed:
            self.ref_f.close()
        if not self.cml_f.closed:
            self.cml_f.close()

    def memmaps(self) -> tuple[np.memmap, np.memmap]:
        self.close()
        return (
            np.memmap(self.ref_path, dtype=np.float32, mode="r", shape=(self.n,)),
            np.memmap(self.cml_path, dtype=np.float32, mode="r", shape=(self.n,)),
        )


# =============================================================================
# Plotting
# =============================================================================
def _plot_contingency(counts: dict[str, int], reference_name: str, output: Path, dpi: int, save_pdf: bool) -> None:
    h = int(counts["hit"])
    f = int(counts["false_alarm"])
    m = int(counts["miss"])
    n = int(counts["correct_negative"])
    ref_wet = h + m
    ref_dry = f + n
    values = np.asarray([[h, f], [m, n]], dtype=np.float64)
    labels = [["Hit", "False alarm"], ["Miss", "Correct negative"]]
    percentages = [
        [100 * h / ref_wet if ref_wet else np.nan, 100 * f / ref_dry if ref_dry else np.nan],
        [100 * m / ref_wet if ref_wet else np.nan, 100 * n / ref_dry if ref_dry else np.nan],
    ]

    fig, ax = plt.subplots(figsize=(8.2, 8.2))
    ax.imshow(values, cmap="Blues", aspect="equal")
    ax.set_xticks([0, 1], ["Wet", "Dry"], fontsize=21)
    ax.set_yticks([0, 1], ["Wet", "Dry"], fontsize=21)
    ax.set_xlabel(f"{reference_name} reference", fontsize=23)
    ax.set_ylabel("CML GRAF detection", fontsize=23)
    for i in range(2):
        for j in range(2):
            rgba = plt.cm.Blues((values[i, j] - np.nanmin(values)) / (np.nanmax(values) - np.nanmin(values) + 1e-12))
            color = "white" if rgba[0] < 0.55 else "#202020"
            ax.text(j, i - 0.19, labels[i][j], ha="center", va="center", fontsize=22, fontweight="bold", color=color)
            ax.text(j, i + 0.05, f"{int(values[i,j]):,}", ha="center", va="center", fontsize=22, color=color)
            denom_label = f"{reference_name.lower()} {'wet' if j == 0 else 'dry'}"
            pct = percentages[i][j]
            ax.text(j, i + 0.29, f"{pct:.1f}% of {denom_label}" if np.isfinite(pct) else "N/A", ha="center", va="center", fontsize=16, color=color)
    ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 2, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=4)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def _metric_annotation(metrics: dict[str, float | int]) -> str:
    pcc = float(metrics["pcc"])
    pcc_text = f"{pcc:.2f}" if np.isfinite(pcc) else "NaN"
    return (
        f"N={int(metrics['n']):,}\n"
        f"RMSE={float(metrics['rmse_mmph']):.2f}\n"
        f"MBE={float(metrics['mbe_mmph']):.2f}\n"
        f"PCC={pcc_text}"
    )


def _plot_hexbin(
    store: PairStore,
    metrics: dict[str, float | int],
    reference_symbol: str,
    output: Path,
    plot_max: float,
    gridsize: int,
    dpi: int,
    save_pdf: bool,
) -> None:
    if store.n == 0:
        raise RuntimeError(f"No pairs available for {output.name}")
    ref, cml = store.memmaps()
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    hb = ax.hexbin(
        ref,
        cml,
        gridsize=int(gridsize),
        extent=(0, float(plot_max), 0, float(plot_max)),
        mincnt=1,
        bins="log",
        cmap="jet",
        linewidths=0,
    )
    ax.plot([0, plot_max], [0, plot_max], "--", color="0.25", linewidth=1.5)
    ax.set_xlim(0, plot_max)
    ax.set_ylim(0, plot_max)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(rf"$RR_{{{reference_symbol}}}$ [mm h$^{{-1}}$]", fontsize=21)
    ax.set_ylabel(r"$RR_{CML}$ [mm h$^{-1}$]", fontsize=21)
    ax.tick_params(labelsize=16)
    ax.text(
        0.98, 0.04, _metric_annotation(metrics),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=14,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "0.6", "alpha": 0.9},
    )
    cbar = fig.colorbar(hb, ax=ax, pad=0.04)
    cbar.set_label("Count", fontsize=19)
    cbar.ax.tick_params(labelsize=14)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    del ref, cml


# =============================================================================
# Stage 1: build filtered all-link GRAF rainfall product
# =============================================================================
def build_v2_product(args: argparse.Namespace) -> dict[str, object]:
    try:
        import pycomlink as pycml
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("pycomlink is required") from exc

    source = Path(args.cml_qc_nc).resolve()
    output = Path(args.output_v2_nc).resolve()
    ds = xr.open_dataset(source, decode_times=True, cache=False)
    handles: Optional[V2Handles] = None
    partial: Optional[Path] = None
    started = time.perf_counter()
    try:
        time_dim = _find_dim(ds.dims, args.time_dim, ("time",), "CML time")
        cml_dim = _find_dim(ds.dims, args.cml_dim, ("cml_id", "cml", "link"), "CML")
        sub_dim = _find_dim(ds.dims, args.sublink_dim, ("sublink_id", "sublink"), "sublink")
        for name in (
            args.trsl_var, args.graf_wet_var, args.graf_valid_var,
            args.link_pass_var, args.sublink_pass_var, args.qc_valid_var,
            args.frequency_var, args.length_var, args.polarization_var,
            args.site_a_lon, args.site_a_lat, args.site_b_lon, args.site_b_lat,
        ):
            if name not in ds.variables and name not in ds.coords:
                raise KeyError(f"Required CML variable {name!r} not found")

        all_times = _normalize_times(ds[time_dim].values, "CML time")
        out_left, out_right = _time_slice_indices(all_times, args.tstart, args.tend)
        output_times = all_times[out_left:out_right]
        proc_start_time = output_times[0] - pd.Timedelta(hours=float(args.baseline_pad_hours))
        proc_left = int(np.searchsorted(all_times.asi8, proc_start_time.value, side="left"))
        proc_right = out_right
        proc_times = all_times[proc_left:proc_right]
        output_slice = slice(out_left - proc_left, out_right - proc_left)

        all_ids = _get_ids(ds, cml_dim)
        link_pass_all = _read_static_link_flag(ds, args.link_pass_var, cml_dim)
        selected_positions = np.flatnonzero(link_pass_all)
        if args.max_links is not None:
            selected_positions = selected_positions[: int(args.max_links)]
        if len(selected_positions) == 0:
            raise RuntimeError("No link passes the requested link-level QC flag")
        selected_ids = all_ids[selected_positions]
        nsub = ds.sizes[sub_dim]
        sublink_values = np.asarray(ds[sub_dim].values if sub_dim in ds.variables else np.arange(nsub))
        sub_pass = _read_static_sublink_flag(ds, args.sublink_pass_var, sub_dim, cml_dim, selected_positions)

        frequency_input = _numeric_metadata_2d(ds[args.frequency_var], sub_dim, cml_dim, selected_positions, nsub)
        frequency_hz, inferred_frequency = _frequency_to_hz(frequency_input, args.frequency_unit)
        length_input = _numeric_metadata_2d(ds[args.length_var], sub_dim, cml_dim, selected_positions, nsub)
        length_km, inferred_length = _length_to_km(length_input, args.length_unit)
        polarization_raw = _string_metadata_2d(ds[args.polarization_var], sub_dim, cml_dim, selected_positions, nsub)
        polarization = np.asarray([[pol_to_HV(v) for v in row] for row in polarization_raw], dtype=object)

        def endpoint(name: str) -> np.ndarray:
            da = ds[name].isel({cml_dim: selected_positions})
            for d in list(da.dims):
                if d != cml_dim:
                    da = da.mean(dim=d, skipna=True)
            return np.asarray(da.transpose(cml_dim).values, dtype=np.float64)

        lon0, lat0 = endpoint(args.site_a_lon), endpoint(args.site_a_lat)
        lon1, lat1 = endpoint(args.site_b_lon), endpoint(args.site_b_lat)

        handles, partial = _create_v2_output(
            output, args, output_times, selected_ids, selected_positions,
            sublink_values, sub_pass, length_km, frequency_input, frequency_hz,
            polarization, lon0, lat0, lon1, lat1, source,
            input_link_count=len(all_ids),
            link_drop_count=int(len(all_ids) - np.count_nonzero(link_pass_all)),
        )

        valid_obs_total = 0
        possible_obs_total = int(len(output_times) * nsub * len(selected_positions))
        failed_sublinks = int(sub_pass.size - np.count_nonzero(sub_pass))
        log_rows = []
        batch_starts = range(0, len(selected_positions), int(args.link_batch))
        pbar = tqdm(batch_starts, total=math.ceil(len(selected_positions) / int(args.link_batch)), desc="GRAF V2 rainfall", unit="batch")
        for batch_no, b0 in enumerate(pbar):
            b1 = min(b0 + int(args.link_batch), len(selected_positions))
            src_pos = selected_positions[b0:b1]
            nb = b1 - b0
            read_slice = slice(proc_left, proc_right)
            trsl = _read_signal_block(ds[args.trsl_var], time_dim, sub_dim, cml_dim, read_slice, src_pos, np.float32)
            wet_raw = _read_signal_block(ds[args.graf_wet_var], time_dim, sub_dim, cml_dim, read_slice, src_pos, np.float32)
            gv_raw = _read_signal_block(ds[args.graf_valid_var], time_dim, sub_dim, cml_dim, read_slice, src_pos, np.float32)
            qv_raw = _read_signal_block(ds[args.qc_valid_var], time_dim, sub_dim, cml_dim, read_slice, src_pos, np.float32)

            static_pass = sub_pass[:, b0:b1][None, :, :]
            wet_binary = np.isfinite(wet_raw) & (wet_raw == 1)
            graf_valid = np.isfinite(gv_raw) & (gv_raw == 1)
            qc_valid = np.isfinite(qv_raw) & (qv_raw == 1)
            usable = static_pass & graf_valid & qc_valid & np.isfinite(trsl)

            rain = np.full(trsl.shape, np.nan, dtype=np.float32)
            tasks = []
            for local_c in range(nb):
                for si in range(nsub):
                    if not sub_pass[si, b0 + local_c]:
                        continue
                    tasks.append((si, local_c))

            def task(si: int, local_c: int):
                return si, local_c, _estimate_one_series(
                    pycml,
                    trsl[:, si, local_c],
                    wet_binary[:, si, local_c],
                    usable[:, si, local_c],
                    proc_times.values,
                    float(length_km[si, b0 + local_c]),
                    float(frequency_hz[si, b0 + local_c]),
                    str(polarization[si, b0 + local_c]),
                    args.kr_model,
                    args.n_average_last_dry,
                    args.pastorek_a_max,
                    args.pastorek_zeta,
                    args.pastorek_d,
                )

            if int(args.n_workers) <= 1:
                results = [task(si, lc) for si, lc in tasks]
            else:
                results = []
                with ThreadPoolExecutor(max_workers=int(args.n_workers)) as ex:
                    futures = [ex.submit(task, si, lc) for si, lc in tasks]
                    for fut in as_completed(futures):
                        results.append(fut.result())
            for si, lc, values in results:
                rain[:, si, lc] = values

            rain_out = rain[output_slice]
            usable_out = usable[output_slice]
            wet_out_bool = wet_binary[output_slice] & usable_out
            trsl_out = trsl[output_slice]
            trsl_out = np.where(usable_out, trsl_out, np.nan).astype(np.float32)
            rain_out[~usable_out] = np.nan
            valid_obs_total += int(np.count_nonzero(usable_out))

            handles.R[:, :, b0:b1] = rain_out
            wet_store = np.full(wet_out_bool.shape, 255, dtype=np.uint8)
            wet_store[usable_out] = wet_out_bool[usable_out].astype(np.uint8)
            handles.wet[:, :, b0:b1] = wet_store
            handles.wet_valid[:, :, b0:b1] = usable_out.astype(np.uint8)
            handles.qc_valid[:, :, b0:b1] = (qc_valid[output_slice] & static_pass).astype(np.uint8)
            if handles.trsl is not None:
                handles.trsl[:, :, b0:b1] = trsl_out

            for lc in range(nb):
                for si in range(nsub):
                    log_rows.append({
                        "cml_id": str(selected_ids[b0 + lc]),
                        "source_cml_index": int(src_pos[lc]),
                        "sublink_position": int(si),
                        "qc_final_sublink_pass": int(sub_pass[si, b0 + lc]),
                        "usable_one_minute_observations": int(np.count_nonzero(usable_out[:, si, lc])),
                        "finite_rain_rate_observations": int(np.count_nonzero(np.isfinite(rain_out[:, si, lc]))),
                    })

            if (batch_no + 1) % int(args.sync_every_batches) == 0:
                handles.nc.sync()
            pbar.set_postfix_str(f"links {b0}:{b1} | valid {np.count_nonzero(usable_out):,}")

        handles.nc.setncattr("possible_one_minute_sublink_observations", possible_obs_total)
        handles.nc.setncattr("retained_one_minute_sublink_observations", int(valid_obs_total))
        handles.nc.setncattr("masked_one_minute_sublink_observations", int(possible_obs_total - valid_obs_total))
        handles.nc.setncattr("retained_static_sublinks", int(np.count_nonzero(sub_pass)))
        handles.nc.setncattr("masked_static_sublinks", failed_sublinks)
        handles.nc.setncattr("frequency_input_unit", inferred_frequency)
        handles.nc.setncattr("length_input_unit", inferred_length)
        handles.nc.setncattr("elapsed_generation_seconds", float(time.perf_counter() - started))
        _finalize_v2(handles, partial, output)
        handles = None

        pd.DataFrame(log_rows).to_csv(Path(args.out_dir) / "v2_product_per_sublink_log.csv", index=False)
        summary = {
            "input_links": int(len(all_ids)),
            "links_passing_qc_final_link_pass": int(len(selected_positions)),
            "links_dropped_by_qc_final_link_pass": int(len(all_ids) - np.count_nonzero(link_pass_all)),
            "static_sublinks_total_in_retained_links": int(sub_pass.size),
            "static_sublinks_passing": int(np.count_nonzero(sub_pass)),
            "static_sublinks_dropped": failed_sublinks,
            "possible_one_minute_sublink_observations": possible_obs_total,
            "retained_one_minute_sublink_observations": int(valid_obs_total),
            "masked_one_minute_sublink_observations": int(possible_obs_total - valid_obs_total),
            "reference_availability_used_for_selection": False,
            "output_v2_nc": str(output),
            "elapsed_seconds": float(time.perf_counter() - started),
        }
        return summary
    except Exception:
        if handles is not None:
            try:
                handles.nc.close()
            except Exception:
                pass
        if partial is not None and partial.exists():
            partial.unlink()
        raise
    finally:
        ds.close()


# =============================================================================
# Stage 2: radar/gauge detection and rain-rate evaluation
# =============================================================================
def evaluate_v2(args: argparse.Namespace, generation_summary: dict[str, object]) -> dict[str, object]:
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    radar_dir = out / "radar"
    gauge_dir = out / "gauge"
    radar_dir.mkdir(exist_ok=True)
    gauge_dir.mkdir(exist_ok=True)
    tmp_dir = Path(tempfile.mkdtemp(prefix="graf_pairs_", dir=out))
    started = time.perf_counter()

    cml = xr.open_dataset(args.output_v2_nc, decode_times=True, cache=False)
    radar = xr.open_dataset(args.radar_nc, decode_times=True, cache=False)
    gauge = xr.open_dataset(args.gauge_nc, decode_times=True, cache=False)
    radar_store = PairStore(tmp_dir, "radar")
    gauge_store = PairStore(tmp_dir, "gauge")
    try:
        cml_times = _normalize_times(cml["time"].values, "V2 CML time")
        cml_ids = _get_ids(cml, "cml_id")
        nlink = len(cml_ids)
        nsub = cml.sizes["sublink_id"]

        # Radar structure and exact map.
        radar_var_name = _find_var(radar, args.radar_var, ("rr", "radar_mmph", "rain_rate"), "radar rain")
        radar_da = radar[radar_var_name]
        radar_time_dim = _find_dim(radar_da.dims, args.radar_time_dim, ("time", "time_5min"), "radar time")
        radar_cml_dim = _find_dim(radar_da.dims, args.radar_cml_dim, ("cml_id", "cml", "link"), "radar CML")
        radar_da = _drop_singleton_extra_dims(radar_da, {radar_time_dim, radar_cml_dim}, radar_var_name).transpose(radar_time_dim, radar_cml_dim)
        radar_times = _normalize_times(radar[radar_time_dim].values, "radar time")
        radar_ids = _get_ids(radar, radar_cml_dim)
        radar_id_pos = {str(cid): i for i, cid in enumerate(radar_ids)}
        radar_positions = np.asarray([radar_id_pos.get(str(cid), -1) for cid in cml_ids], dtype=np.int64)
        radar_map = _build_interval_map(cml_times, radar_times, args.radar_minutes, args.radar_timestamp_label)

        # Gauge structure and exact map.
        gauge_var_name = _find_var(gauge, args.gauge_var, ("rr", "rainfall", "precipitation"), "gauge rain")
        gauge_da = gauge[gauge_var_name]
        gauge_time_dim = _find_dim(gauge_da.dims, args.gauge_time_dim, ("time", "time_gauge"), "gauge time")
        gauge_station_dim = _infer_station_dim(gauge_da, gauge_time_dim, args.gauge_station_dim)
        gauge_da = _drop_singleton_extra_dims(gauge_da, {gauge_time_dim, gauge_station_dim}, gauge_var_name).transpose(gauge_time_dim, gauge_station_dim)
        gauge_times = _normalize_times(gauge[gauge_time_dim].values, "gauge time")
        gauge_ids = _get_ids(gauge, gauge_station_dim)
        gauge_lon_name = _find_var(gauge, args.gauge_lon, ("lon", "longitude", "station_lon", "station_longitude"), "gauge longitude")
        gauge_lat_name = _find_var(gauge, args.gauge_lat, ("lat", "latitude", "station_lat", "station_latitude"), "gauge latitude")
        gauge_lon = _read_numeric_for_dim(gauge, gauge_lon_name, gauge_station_dim)
        gauge_lat = _read_numeric_for_dim(gauge, gauge_lat_name, gauge_station_dim)
        gauge_map = _build_interval_map(cml_times, gauge_times, args.gauge_minutes, args.gauge_timestamp_label)

        lon0 = np.asarray(cml["site_a_longitude"].values, dtype=np.float64)
        lat0 = np.asarray(cml["site_a_latitude"].values, dtype=np.float64)
        lon1 = np.asarray(cml["site_b_longitude"].values, dtype=np.float64)
        lat1 = np.asarray(cml["site_b_latitude"].values, dtype=np.float64)
        gauge_assignment = _assign_closest_gauges(
            cml_ids, lon0, lat0, lon1, lat1,
            gauge_ids, gauge_lon, gauge_lat,
            args.projection, args.max_gauge_distance_m, args.spatial_batch,
        )
        gauge_assignment.to_csv(gauge_dir / "closest_gauge_assignment.csv", index=False)
        gauge_positions = gauge_assignment["gauge_position"].to_numpy(dtype=np.int64)
        used_gauges = np.asarray(sorted(set(int(v) for v in gauge_positions if int(v) >= 0)), dtype=np.int64)
        gauge_col = {int(pos): i for i, pos in enumerate(used_gauges)}
        if len(used_gauges):
            gauge_values_used = np.asarray(
                gauge_da.isel({gauge_station_dim: used_gauges}).values,
                dtype=np.float32,
            ) * float(args.gauge_scale)
            gauge_values_used[~np.isfinite(gauge_values_used) | (gauge_values_used < 0)] = np.nan
        else:
            gauge_values_used = np.empty((len(gauge_times), 0), dtype=np.float32)

        radar_cont = ContingencyAccumulator()
        gauge_cont = ContingencyAccumulator()
        radar_metrics = ContinuousAccumulator()
        gauge_metrics = ContinuousAccumulator()
        radar_links_cont = np.zeros(nlink, dtype=bool)
        radar_links_hex = np.zeros(nlink, dtype=bool)
        gauge_links_cont = np.zeros(nlink, dtype=bool)
        gauge_links_hex = np.zeros(nlink, dtype=bool)

        batch_starts = range(0, nlink, int(args.eval_link_batch))
        pbar = tqdm(batch_starts, total=math.ceil(nlink / int(args.eval_link_batch)), desc="reference evaluation", unit="batch")
        for b0 in pbar:
            b1 = min(b0 + int(args.eval_link_batch), nlink)
            nb = b1 - b0
            rain1 = np.asarray(cml[args.cml_rain_var].isel(cml_id=slice(b0, b1)).transpose("time", "sublink_id", "cml_id").values, dtype=np.float32)
            wet_raw = np.asarray(cml[args.cml_wet_var].isel(cml_id=slice(b0, b1)).transpose("time", "sublink_id", "cml_id").values)
            wet1 = wet_raw == 1
            valid1 = np.asarray(cml[args.cml_valid_var].isel(cml_id=slice(b0, b1)).transpose("time", "sublink_id", "cml_id").values) == 1
            valid1 &= np.isfinite(rain1)

            rain5, wet5, valid5 = _aggregate_native(
                rain1, wet1, valid1, radar_map,
                args.wet_aggregation, args.min_cml_valid_fraction,
            )
            rain10, wet10, valid10 = _aggregate_native(
                rain1, wet1, valid1, gauge_map,
                args.wet_aggregation, args.min_cml_valid_fraction,
            )

            # Radar values, independently mapped by CML ID.
            rr = np.full((len(radar_times), nb), np.nan, dtype=np.float32)
            local_positions = radar_positions[b0:b1]
            good_local = np.flatnonzero(local_positions >= 0)
            if len(good_local):
                vals = np.asarray(
                    radar_da.isel({radar_cml_dim: local_positions[good_local]}).values,
                    dtype=np.float32,
                ) * float(args.radar_scale)
                vals[~np.isfinite(vals) | (vals < 0)] = np.nan
                rr[:, good_local] = vals
            rr_valid = np.isfinite(rr) & (rr >= 0)
            rr_wet = rr_valid & (rr >= float(args.wet_threshold))
            rr3 = rr[:, None, :]
            rr_valid3 = rr_valid[:, None, :]
            rr_wet3 = rr_wet[:, None, :]
            det_valid = valid5 & rr_valid3
            radar_cont.update(wet5, rr_wet3, det_valid)
            radar_links_cont[b0:b1] |= np.any(det_valid, axis=(0, 1))

            rain_valid = valid5 & rr_valid3 & np.isfinite(rain5) & (rain5 >= 0)
            if args.hexbin_scope == "reference_wet":
                rain_valid &= rr_wet3
            ref_rep = np.broadcast_to(rr3, rain5.shape)
            radar_metrics.update(ref_rep[rain_valid], rain5[rain_valid])
            radar_store.append(ref_rep[rain_valid], rain5[rain_valid])
            radar_links_hex[b0:b1] |= np.any(rain_valid, axis=(0, 1))

            # Gauge values, independently assigned; radar availability is irrelevant.
            gg = np.full((len(gauge_times), nb), np.nan, dtype=np.float32)
            for lc, gp in enumerate(gauge_positions[b0:b1]):
                if int(gp) >= 0:
                    gg[:, lc] = gauge_values_used[:, gauge_col[int(gp)]]
            gg_valid = np.isfinite(gg) & (gg >= 0)
            gg_wet = gg_valid & (gg >= float(args.wet_threshold))
            gg3 = gg[:, None, :]
            gg_valid3 = gg_valid[:, None, :]
            gg_wet3 = gg_wet[:, None, :]
            gdet_valid = valid10 & gg_valid3
            gauge_cont.update(wet10, gg_wet3, gdet_valid)
            gauge_links_cont[b0:b1] |= np.any(gdet_valid, axis=(0, 1))

            grain_valid = valid10 & gg_valid3 & np.isfinite(rain10) & (rain10 >= 0)
            if args.hexbin_scope == "reference_wet":
                grain_valid &= gg_wet3
            gref_rep = np.broadcast_to(gg3, rain10.shape)
            gauge_metrics.update(gref_rep[grain_valid], rain10[grain_valid])
            gauge_store.append(gref_rep[grain_valid], rain10[grain_valid])
            gauge_links_hex[b0:b1] |= np.any(grain_valid, axis=(0, 1))

        radar_store.close()
        gauge_store.close()
        radar_counts = radar_cont.counts()
        gauge_counts = gauge_cont.counts()
        radar_det_metrics = radar_cont.metrics()
        gauge_det_metrics = gauge_cont.metrics()
        radar_rain_metrics = radar_metrics.result()
        gauge_rain_metrics = gauge_metrics.result()

        # Save tables and plots.
        pd.DataFrame([{**radar_counts, **radar_det_metrics}]).to_csv(radar_dir / "graf_vs_radar_contingency_metrics.csv", index=False)
        pd.DataFrame([{**gauge_counts, **gauge_det_metrics}]).to_csv(gauge_dir / "graf_vs_gauge_contingency_metrics.csv", index=False)
        pd.DataFrame([{**radar_rain_metrics, "scope": args.hexbin_scope, "resolution_minutes": args.radar_minutes}]).to_csv(radar_dir / "graf_rain_vs_radar_metrics.csv", index=False)
        pd.DataFrame([{**gauge_rain_metrics, "scope": args.hexbin_scope, "resolution_minutes": args.gauge_minutes}]).to_csv(gauge_dir / "graf_rain_vs_gauge_metrics.csv", index=False)

        _plot_contingency(radar_counts, "Radar", radar_dir / "graf_vs_radar_contingency.png", args.dpi, args.save_pdf)
        _plot_contingency(gauge_counts, "Rain gauge", gauge_dir / "graf_vs_gauge_contingency.png", args.dpi, args.save_pdf)
        _plot_hexbin(radar_store, radar_rain_metrics, "Radar", radar_dir / "graf_rain_vs_radar_hexbin_pooled_sublinks.png", args.plot_max, args.gridsize, args.dpi, args.save_pdf)
        _plot_hexbin(gauge_store, gauge_rain_metrics, "Gauge", gauge_dir / "graf_rain_vs_gauge_hexbin_pooled_sublinks.png", args.plot_max, args.gridsize, args.dpi, args.save_pdf)

        reference_summary = pd.DataFrame([
            {
                "reference": "radar",
                "reference_resolution_minutes": args.radar_minutes,
                "links_in_v2": nlink,
                "links_with_reference_id_mapping": int(np.count_nonzero(radar_positions >= 0)),
                "links_with_paired_contingency_observations": int(np.count_nonzero(radar_links_cont)),
                "links_with_hexbin_observations": int(np.count_nonzero(radar_links_hex)),
                "unique_reference_stations": np.nan,
                "pooled_sublink_contingency_N": int(radar_det_metrics["N"]),
                "pooled_sublink_hexbin_N": int(radar_rain_metrics["n"]),
            },
            {
                "reference": "rain_gauge",
                "reference_resolution_minutes": args.gauge_minutes,
                "links_in_v2": nlink,
                "links_with_reference_id_mapping": int(np.count_nonzero(gauge_positions >= 0)),
                "links_with_paired_contingency_observations": int(np.count_nonzero(gauge_links_cont)),
                "links_with_hexbin_observations": int(np.count_nonzero(gauge_links_hex)),
                "unique_reference_stations": int(len(used_gauges)),
                "pooled_sublink_contingency_N": int(gauge_det_metrics["N"]),
                "pooled_sublink_hexbin_N": int(gauge_rain_metrics["n"]),
            },
        ])
        reference_summary.to_csv(out / "reference_population_summary.csv", index=False)

        filtering_rows = [
            {"stage": "Input CML links", "retained": int(generation_summary["input_links"]), "dropped_at_stage": 0, "unit": "links"},
            {"stage": "After qc_final_link_pass", "retained": int(generation_summary["links_passing_qc_final_link_pass"]), "dropped_at_stage": int(generation_summary["links_dropped_by_qc_final_link_pass"]), "unit": "links"},
            {"stage": "Static sublinks after qc_final_sublink_pass", "retained": int(generation_summary["static_sublinks_passing"]), "dropped_at_stage": int(generation_summary["static_sublinks_dropped"]), "unit": "sublinks"},
            {"stage": "Native observations after qc_final_valid + GRAF validity + finite TRSL", "retained": int(generation_summary["retained_one_minute_sublink_observations"]), "dropped_at_stage": int(generation_summary["masked_one_minute_sublink_observations"]), "unit": "one-minute sublink observations"},
            {"stage": "Radar links with any paired contingency observation", "retained": int(np.count_nonzero(radar_links_cont)), "dropped_at_stage": int(nlink - np.count_nonzero(radar_links_cont)), "unit": "links relative to V2"},
            {"stage": "Gauge links with any paired contingency observation", "retained": int(np.count_nonzero(gauge_links_cont)), "dropped_at_stage": int(nlink - np.count_nonzero(gauge_links_cont)), "unit": "links relative to V2"},
        ]
        pd.DataFrame(filtering_rows).to_csv(out / "filtering_steps.csv", index=False)

        summary = {
            "script": Path(__file__).name,
            "script_version": SCRIPT_VERSION,
            "v2_product": str(Path(args.output_v2_nc).resolve()),
            "filtering": generation_summary,
            "reference_selection_independent": True,
            "radar": {
                "raw_scale_to_mmph": args.radar_scale,
                "resolution_minutes": args.radar_minutes,
                "links_with_id_mapping": int(np.count_nonzero(radar_positions >= 0)),
                "links_with_contingency_pairs": int(np.count_nonzero(radar_links_cont)),
                "links_with_hexbin_pairs": int(np.count_nonzero(radar_links_hex)),
                "contingency_counts": radar_counts,
                "contingency_metrics": radar_det_metrics,
                "rain_metrics": radar_rain_metrics,
            },
            "gauge": {
                "raw_scale_to_mmph": args.gauge_scale,
                "resolution_minutes": args.gauge_minutes,
                "maximum_distance_m": args.max_gauge_distance_m,
                "links_assigned": int(np.count_nonzero(gauge_positions >= 0)),
                "unique_gauges_used": int(len(used_gauges)),
                "links_with_contingency_pairs": int(np.count_nonzero(gauge_links_cont)),
                "links_with_hexbin_pairs": int(np.count_nonzero(gauge_links_hex)),
                "contingency_counts": gauge_counts,
                "contingency_metrics": gauge_det_metrics,
                "rain_metrics": gauge_rain_metrics,
            },
            "pooled_sublinks": True,
            "sublink_averaging": False,
            "wet_threshold_mmph": args.wet_threshold,
            "wet_aggregation": args.wet_aggregation,
            "minimum_native_valid_fraction": args.min_cml_valid_fraction,
            "hexbin_scope": args.hexbin_scope,
            "plot_limit_mmph": args.plot_max,
            "metrics_clipped_to_plot_limit": False,
            "elapsed_evaluation_seconds": float(time.perf_counter() - started),
        }
        with open(out / "run_summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2, allow_nan=True)

        print("\n=== BLETTNER / V2 FILTERING ===")
        print(f"Input CML links                         : {generation_summary['input_links']:,}")
        print(f"Links passing qc_final_link_pass        : {generation_summary['links_passing_qc_final_link_pass']:,}")
        print(f"Links dropped by qc_final_link_pass     : {generation_summary['links_dropped_by_qc_final_link_pass']:,}")
        print(f"Static sublinks passing                 : {generation_summary['static_sublinks_passing']:,}")
        print(f"Static sublinks dropped                 : {generation_summary['static_sublinks_dropped']:,}")
        print(f"Native observations retained            : {generation_summary['retained_one_minute_sublink_observations']:,}")
        print(f"Native observations masked              : {generation_summary['masked_one_minute_sublink_observations']:,}")
        print("\n=== RADAR POPULATION ===")
        print(f"V2 links with radar ID mapping          : {np.count_nonzero(radar_positions >= 0):,}")
        print(f"Links used in radar contingency         : {np.count_nonzero(radar_links_cont):,}")
        print(f"Links used in radar hexbin              : {np.count_nonzero(radar_links_hex):,}")
        print(f"Radar contingency pooled N              : {radar_det_metrics['N']:,}")
        print(f"Radar hexbin pooled N                   : {radar_rain_metrics['n']:,}")
        print("\n=== GAUGE POPULATION ===")
        print(f"Gauge stations in source                : {len(gauge_ids):,}")
        print(f"Unique gauges used within distance      : {len(used_gauges):,}")
        print(f"V2 links assigned to a gauge            : {np.count_nonzero(gauge_positions >= 0):,}")
        print(f"Links used in gauge contingency         : {np.count_nonzero(gauge_links_cont):,}")
        print(f"Links used in gauge hexbin              : {np.count_nonzero(gauge_links_hex):,}")
        print(f"Gauge contingency pooled N              : {gauge_det_metrics['N']:,}")
        print(f"Gauge hexbin pooled N                   : {gauge_rain_metrics['n']:,}")
        return summary
    finally:
        radar_store.close()
        gauge_store.close()
        cml.close()
        radar.close()
        gauge.close()
        if not args.keep_pair_files:
            shutil.rmtree(tmp_dir, ignore_errors=True)


# =============================================================================
# CLI
# =============================================================================
def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build a Blettner-filtered all-link GRAF rainfall V2 product and publishable radar/gauge contingency and hexbin figures."
    )
    p.add_argument("--cml-qc-nc", required=True)
    p.add_argument("--radar-nc", required=True)
    p.add_argument("--gauge-nc", required=True)
    p.add_argument("--output-v2-nc", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--skip-v2-generation", action="store_true", help="Use an already generated V2 product; filtering summary is read from its attributes.")

    p.add_argument("--time-dim", default="time")
    p.add_argument("--cml-dim", default="cml_id")
    p.add_argument("--sublink-dim", default="sublink_id")
    p.add_argument("--trsl-var", default="trsl_qc")
    p.add_argument("--graf-wet-var", default="graf_rstd_wet_1min")
    p.add_argument("--graf-valid-var", default="graf_rstd_valid_1min")
    p.add_argument("--link-pass-var", default="qc_final_link_pass")
    p.add_argument("--sublink-pass-var", default="qc_final_sublink_pass")
    p.add_argument("--qc-valid-var", default="qc_final_valid")
    p.add_argument("--frequency-var", default="frequency")
    p.add_argument("--length-var", default="length")
    p.add_argument("--polarization-var", default="polarization")
    p.add_argument("--site-a-lon", default="site_a_longitude")
    p.add_argument("--site-a-lat", default="site_a_latitude")
    p.add_argument("--site-b-lon", default="site_b_longitude")
    p.add_argument("--site-b-lat", default="site_b_latitude")
    p.add_argument("--frequency-unit", choices=("Hz", "MHz", "GHz", "auto"), default="MHz")
    p.add_argument("--length-unit", choices=("m", "km", "auto"), default="km")
    p.add_argument("--kr-model", choices=ALLOWED_KR_MODELS, default="ITU_2005")
    p.add_argument("--n-average-last-dry", type=int, default=1)
    p.add_argument("--pastorek-a-max", type=float, default=14.0)
    p.add_argument("--pastorek-zeta", type=float, default=0.55)
    p.add_argument("--pastorek-d", type=float, default=0.1)
    p.add_argument("--baseline-pad-hours", type=float, default=48.0)
    p.add_argument("--tstart", default=None)
    p.add_argument("--tend", default=None)
    p.add_argument("--max-links", type=int, default=None, help="Testing only; omit for every QC-passed link.")
    p.add_argument("--n-workers", type=int, default=min(16, max(1, os.cpu_count() or 1)))
    p.add_argument("--link-batch", type=int, default=8)
    p.add_argument("--compression-level", type=int, default=1)
    p.add_argument("--time-chunk-1min", type=int, default=10080)
    p.add_argument("--output-link-chunk", type=int, default=8)
    p.add_argument("--sync-every-batches", type=int, default=4)
    p.add_argument("--write-trsl", action="store_true")

    p.add_argument("--cml-rain-var", default="R_cml")
    p.add_argument("--cml-wet-var", default="graf_wet")
    p.add_argument("--cml-valid-var", default="graf_wet_valid")
    p.add_argument("--eval-link-batch", type=int, default=8)
    p.add_argument("--wet-aggregation", choices=("any", "majority", "all"), default="any")
    p.add_argument("--min-cml-valid-fraction", type=float, default=1.0)
    p.add_argument("--wet-threshold", type=float, default=0.1)
    p.add_argument("--hexbin-scope", choices=("reference_wet", "all_paired"), default="reference_wet")

    p.add_argument("--radar-var", default="rr")
    p.add_argument("--radar-time-dim", default=None)
    p.add_argument("--radar-cml-dim", default=None)
    p.add_argument("--radar-scale", type=float, default=12.0)
    p.add_argument("--radar-minutes", type=int, default=5)
    p.add_argument("--radar-timestamp-label", choices=("start", "end", "center"), default="start")

    p.add_argument("--gauge-var", default="rr")
    p.add_argument("--gauge-time-dim", default=None)
    p.add_argument("--gauge-station-dim", default=None)
    p.add_argument("--gauge-lon", default=None)
    p.add_argument("--gauge-lat", default=None)
    p.add_argument("--gauge-scale", type=float, default=6.0)
    p.add_argument("--gauge-minutes", type=int, default=10)
    p.add_argument("--gauge-timestamp-label", choices=("start", "end", "center"), default="start")
    p.add_argument("--projection", default="EPSG:25833")
    p.add_argument("--max-gauge-distance-m", type=float, default=1000.0)
    p.add_argument("--spatial-batch", type=int, default=128)

    p.add_argument("--plot-max", type=float, default=150.0)
    p.add_argument("--gridsize", type=int, default=55)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--save-pdf", action="store_true")
    p.add_argument("--keep-pair-files", action="store_true")

    args = p.parse_args(argv)
    if args.kr_model not in ALLOWED_KR_MODELS:
        p.error(f"--kr-model must be one of {ALLOWED_KR_MODELS}")
    if args.n_workers < 1 or args.link_batch < 1 or args.eval_link_batch < 1:
        p.error("worker and batch sizes must be positive")
    if not (0 < args.min_cml_valid_fraction <= 1):
        p.error("--min-cml-valid-fraction must be in (0,1]")
    if args.radar_minutes < 1 or args.gauge_minutes < 1:
        p.error("reference resolutions must be positive")
    if args.max_gauge_distance_m <= 0:
        p.error("--max-gauge-distance-m must be positive")
    if args.plot_max <= 0 or args.gridsize < 5:
        p.error("invalid plot settings")
    return args


def _generation_summary_from_v2(path: str) -> dict[str, object]:
    with xr.open_dataset(path, decode_times=False) as ds:
        attrs = ds.attrs
        retained_links = int(ds.sizes["cml_id"])
        return {
            "input_links": int(attrs.get("input_link_count", retained_links)),
            "links_passing_qc_final_link_pass": retained_links,
            "links_dropped_by_qc_final_link_pass": int(attrs.get("dropped_link_count", 0)),
            "static_sublinks_total_in_retained_links": int(ds.sizes["sublink_id"] * retained_links),
            "static_sublinks_passing": int(attrs.get("retained_static_sublinks", np.asarray(ds["qc_final_sublink_pass"].values).sum())),
            "static_sublinks_dropped": int(attrs.get("masked_static_sublinks", 0)),
            "possible_one_minute_sublink_observations": int(attrs.get("possible_one_minute_sublink_observations", 0)),
            "retained_one_minute_sublink_observations": int(attrs.get("retained_one_minute_sublink_observations", 0)),
            "masked_one_minute_sublink_observations": int(attrs.get("masked_one_minute_sublink_observations", 0)),
            "reference_availability_used_for_selection": False,
            "output_v2_nc": str(Path(path).resolve()),
        }


def main() -> None:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    args = parse_args()
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    print(f"Running {Path(__file__).name} version {SCRIPT_VERSION}")
    if args.skip_v2_generation:
        if not Path(args.output_v2_nc).exists():
            raise FileNotFoundError(args.output_v2_nc)
        generation_summary = _generation_summary_from_v2(args.output_v2_nc)
    else:
        generation_summary = build_v2_product(args)
    with open(Path(args.out_dir) / "v2_filtering_summary.json", "w", encoding="utf-8") as f:
        json.dump(generation_summary, f, indent=2, allow_nan=True)
    evaluate_v2(args, generation_summary)
    print(f"\nSaved V2 product: {Path(args.output_v2_nc).resolve()}")
    print(f"Saved results   : {Path(args.out_dir).resolve()}")


if __name__ == "__main__":
    main()
