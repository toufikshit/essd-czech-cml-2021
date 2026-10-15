#!/usr/bin/env python3
"""
Blettner/Graf-style CML processing + rainfall estimation + event IDW maps.

This script applies the processing logic from the Example_Processing notebook,
but in a streamed way for the full CML dataset:

  1) period plateau filter for the full processing line
  2) blackout gap filling for the full processing line
  3) TRSL = TSL - RSL
  4) 5 min linear interpolation of short TRSL gaps for graf_2020 and full
  5) CML sanity labels: detection limit, long std, short std, availability
  6) derive_rain() using the same pycomlink baseline/WAA/k-R logic
  7) resample 1 min rain rates to 5 min
  8) interpolate CML, radar path averages, and gauges to event rain maps
  9) write NetCDF outputs and an optional GIF

The goal is speed and memory safety without changing the scientific logic.
The script avoids building long lists of single-link xarray datasets and avoids
large xr.concat operations during processing.
"""
from __future__ import annotations

import argparse
import gc
import time
from pathlib import Path
from typing import Optional, Sequence, Tuple, Dict, List

import numpy as np
import pandas as pd
import xarray as xr

try:
    import numba
except Exception:
    numba = None

try:
    from tqdm import tqdm
except Exception:
    def tqdm(x, **kwargs):
        return x

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import imageio.v2 as imageio
except Exception:
    imageio = None


# ============================================================
# Generic helpers
# ============================================================
def find_dim(ds: xr.Dataset, candidates: Sequence[str]) -> Optional[str]:
    for c in candidates:
        if c in ds.dims:
            return c
    return None


def find_var(ds: xr.Dataset, candidates: Sequence[str]) -> Optional[str]:
    for c in candidates:
        if c in ds.data_vars or c in ds.coords:
            return c
    return None


def time_slice_indices(times: pd.DatetimeIndex, tstart: pd.Timestamp, tend: pd.Timestamp) -> Tuple[int, int]:
    vals = times.values
    i0 = int(np.searchsorted(vals, np.datetime64(tstart), "left"))
    i1 = int(np.searchsorted(vals, np.datetime64(tend), "right"))
    return i0, i1


def as_str_array(a):
    return np.asarray([str(x) for x in np.asarray(a).ravel()])


def safe_float(x) -> float:
    try:
        return float(np.asarray(x).squeeze())
    except Exception:
        return float("nan")


def pol_to_str(pol_raw) -> str:
    if isinstance(pol_raw, (bytes, np.bytes_)):
        pol_raw = pol_raw.decode("utf-8", errors="ignore")
    s = str(np.asarray(pol_raw).squeeze()).strip().lower()
    if "v" in s and "h" not in s:
        return "V"
    if "h" in s and "v" not in s:
        return "H"
    return "H"


def frequency_to_hz(f_raw: float) -> float:
    """Accept Hz, MHz, or GHz and return Hz.

    Dataset convention is often MHz. Notebook functions expect Hz.
    """
    f = float(f_raw)
    if not np.isfinite(f):
        return np.nan
    if f > 1e8:      # already Hz
        return f
    if f > 1e3:      # MHz
        return f * 1e6
    return f * 1e9   # GHz


def select_scalar_or_index(da: xr.DataArray, cml_dim: str, cidx: int, sublink_dim: Optional[str], sl: int):
    sel = {}
    if cml_dim in da.dims:
        sel[cml_dim] = cidx
    if sublink_dim is not None and sublink_dim in da.dims:
        sel[sublink_dim] = sl
    return da.isel(sel).values


def extract_signal_1d(da: xr.DataArray, time_dim: str, cml_dim: str, cidx: int,
                      sublink_dim: Optional[str], sl: int, i0: int, i1: int) -> np.ndarray:
    sel = {time_dim: slice(i0, i1)}
    if cml_dim in da.dims:
        sel[cml_dim] = cidx
    if sublink_dim is not None and sublink_dim in da.dims:
        sel[sublink_dim] = sl
    out = da.isel(sel)
    if time_dim in out.dims:
        out = out.transpose(time_dim)
    return np.asarray(out.values, dtype="float64").reshape(-1)


def resample_1min_rate_to_5min_rate(values_1min: np.ndarray, time_1min: pd.DatetimeIndex,
                                    target_5min: pd.DatetimeIndex) -> np.ndarray:
    """Convert 1 min rain rate to 5 min rain rate by conserving accumulation."""
    s = pd.Series(np.asarray(values_1min, dtype="float64"), index=time_1min)
    out = (s / 60.0).resample("5min", label="left", closed="left").sum(min_count=1) * 12.0
    return out.reindex(target_5min).to_numpy(dtype="float32")


# ============================================================
# Exact processing functions from the notebook logic
# ============================================================
def detect_plateaus(series, time_span=3, max_thld=-85, std_thld=0.5, pad=5):
    low = series.rolling(time=time_span, center=True).max("time") < max_thld
    steady = series.rolling(time=time_span, center=True).std("time") < std_thld
    isSane = ~(low & steady)
    isSane = ~(isSane.rolling(time=pad, center=True).min("time") == 0)
    return isSane


def label_timeseries(cml):
    cml["isSane_plateaus_rsl"] = (["time"], detect_plateaus(cml.rsl).values)
    return cml


def get_blackout_start_and_end_v1(rsl, rsl_threshold):
    with np.errstate(invalid="ignore"):
        rsl_below_threshold = rsl < rsl_threshold
    rsl_nan = np.isnan(rsl)
    gap_start = np.roll(rsl_nan, -1) & rsl_below_threshold
    gap_end = np.roll(rsl_nan, 1) & rsl_below_threshold
    return gap_start, gap_end


if numba is not None:
    @numba.jit(nopython=True)
    def created_blackout_gap_mask_from_start_end_markers(rsl, gap_start, gap_end, max_gap_length=60):
        mask = np.zeros(rsl.shape, dtype=numba.boolean)
        in_blackout_gap = False
        cnt_gap_length = 0
        for i in range(len(rsl)):
            if gap_start[i] == True:
                in_blackout_gap = True
                cnt_gap_length = 0
            if (gap_start[i] == False) and ~np.isnan(rsl[i]):
                in_blackout_gap = False
            if gap_end[i] == True:
                in_blackout_gap = False
            if in_blackout_gap and np.isnan(rsl[i]) and (cnt_gap_length < max_gap_length):
                mask[i] = True
                cnt_gap_length += 1
                if cnt_gap_length > max_gap_length:
                    mask[(i-max_gap_length):(i+1)] = False
                    in_blackout_gap = False
        return mask
else:
    def created_blackout_gap_mask_from_start_end_markers(rsl, gap_start, gap_end, max_gap_length=60):
        mask = np.zeros(rsl.shape, dtype=bool)
        in_blackout_gap = False
        cnt_gap_length = 0
        for i in range(len(rsl)):
            if gap_start[i] == True:
                in_blackout_gap = True
                cnt_gap_length = 0
            if (gap_start[i] == False) and ~np.isnan(rsl[i]):
                in_blackout_gap = False
            if gap_end[i] == True:
                in_blackout_gap = False
            if in_blackout_gap and np.isnan(rsl[i]) and (cnt_gap_length < max_gap_length):
                mask[i] = True
                cnt_gap_length += 1
                if cnt_gap_length > max_gap_length:
                    mask[(i-max_gap_length):(i+1)] = False
                    in_blackout_gap = False
        return mask


def blackout_filling_one(tsl_in: np.ndarray, rsl_in: np.ndarray, tsl_max: float, rsl_min: float):
    tsl_out = np.copy(tsl_in)
    rsl_out = np.copy(rsl_in)

    gap_start, gap_end = get_blackout_start_and_end_v1(rsl=rsl_out, rsl_threshold=-65)
    mask_forward = created_blackout_gap_mask_from_start_end_markers(rsl_out, gap_start, gap_end)
    mask_reverse = created_blackout_gap_mask_from_start_end_markers(rsl_out[::-1], gap_end[::-1], gap_start[::-1])
    mask = mask_forward | mask_reverse[::-1]

    tsl_out[mask] = tsl_max
    rsl_out[mask] = rsl_min
    return tsl_out, rsl_out, mask


def detection_limit_from_kR(L_km, freq_Hz, pol, A_quantization_dB=0.33):
    import pycomlink as pycml
    try:
        a, b = pycml.processing.k_R_relation.a_b(
            f_GHz=freq_Hz / 1e9, pol=pol, approx_type="ITU_2003"
        )
    except TypeError:
        a, b = pycml.processing.k_R_relation.a_b(
            f_GHz=freq_Hz / 1e9, pol=pol, a_b_approximation="ITU_2003"
        )
    sensitivity = a * L_km
    detection_limit = A_quantization_dB / sensitivity
    return detection_limit


def temporal_sanity_check(trsl, tspan=300, thld=2, perc=0.1):
    roll_std = trsl.rolling(time=tspan, center=True).std()
    high_std = roll_std > thld
    low_std = roll_std <= thld
    isSane = (
        high_std.sum(dim="time") / (high_std.sum(dim="time") + low_std.sum(dim="time"))
        <= perc
    )
    return isSane


def availability_check(trsl, perc_thld=0.5):
    isSane = (trsl.isnull().sum("time") / len(trsl.time)) < perc_thld
    return isSane


def derive_rain(cml, waa_model="Leijnse"):
    """Same rainfall derivation logic as the Example_Processing notebook."""
    import pycomlink as pycml

    wet_threshold = 1.12 * cml.trsl.rolling(
        time=60, center=True, min_periods=60
    ).std().quantile(0.8, dim="time")

    wet = cml.trsl_raw.rolling(
        time=60, center=True, min_periods=45
    ).std(dim="time") > wet_threshold

    if "flag_blackout_fill" in cml.data_vars or "flag_blackout_fill" in cml.coords:
        wet = wet.where(~cml.flag_blackout_fill, True)

    baseline = pycml.processing.baseline.baseline_constant(
        trsl=cml.trsl,
        wet=wet,
        n_average_last_dry=1,
    )

    if waa_model == "Pastorek":
        A_obs = cml.trsl - baseline
        A_obs = A_obs.where(~(A_obs < 0), 0)
        waa = pycml.processing.wet_antenna.waa_pastorek_2021_from_A_obs(
            A_obs=A_obs,
            f_Hz=cml.frequency,
            L_km=cml.length,
            A_max=14,
            zeta=0.55,
            d=0.1,
        )

    elif waa_model == "Schleiss":
        waa = pycml.processing.wet_antenna.waa_schleiss_2013(
            rsl=cml.trsl,
            baseline=baseline,
            wet=wet,
            waa_max=2.5,
            delta_t=1,
            tau=15,
        )

    elif waa_model == "Leijnse":
        A_obs = cml.trsl - baseline
        A_obs = A_obs.where(~(A_obs < 0), 0)
        waa = pycml.processing.wet_antenna.waa_leijnse_2008_from_A_obs(
            A_obs=A_obs,
            f_Hz=cml.frequency,
            L_km=cml.length,
            T_K=293.0,
            gamma=1.47e-5,
            delta=0.36,
            n_antenna=complex(1.73, 0.014),
            l_antenna=0.041,
        )
    else:
        raise ValueError("waa_model must be Leijnse, Pastorek, or Schleiss")

    A = cml.trsl - baseline - waa
    if "quantile" in A.coords:
        A = A.drop_vars("quantile")

    cml["R_cml"] = pycml.processing.k_R_relation.calc_R_from_A(
        A=A,
        L_km=float(cml.length),
        f_GHz=cml.frequency / 1e9,
        pol=cml.polarization,
    )
    return cml


# ============================================================
# Per-series pipeline
# ============================================================
def process_one_sublink(
    tsl: np.ndarray,
    rsl: np.ndarray,
    proc_time: pd.DatetimeIndex,
    out_rel0: int,
    out_rel1: int,
    length_km: float,
    frequency_hz: float,
    polarization: str,
    proc_line: str,
    waa_model: str,
) -> Tuple[np.ndarray, Dict[str, float]]:
    """Apply notebook processing to one sublink time series."""
    cml = xr.Dataset(
        data_vars={
            "tsl": ("time", np.asarray(tsl, dtype="float64")),
            "rsl": ("time", np.asarray(rsl, dtype="float64")),
        },
        coords={
            "time": proc_time.values,
            "length": float(length_km),
            "frequency": float(frequency_hz),
            "polarization": polarization,
        },
    )

    info: Dict[str, float] = {}
    info["valid_input"] = int(np.isfinite(tsl).sum() + np.isfinite(rsl).sum())

    # Period filter and blackout gap filling only for full, as in notebook.
    if proc_line in ["full"]:
        cml = label_timeseries(cml)
        cml["tsl"] = cml.tsl.where(cml["isSane_plateaus_rsl"])
        cml["rsl"] = cml.rsl.where(cml["isSane_plateaus_rsl"])

        tsl_max = safe_float(cml.tsl.max("time").values)
        rsl_min = safe_float(cml.rsl.min("time").values)
        tsl_b, rsl_b, mask_b = blackout_filling_one(cml.tsl.values, cml.rsl.values, tsl_max, rsl_min)
        cml["flag_blackout_fill"] = ("time", cml.rsl.isnull().values & np.isfinite(rsl_b))
        cml["tsl"] = ("time", tsl_b)
        cml["rsl"] = ("time", rsl_b)
    else:
        cml["isSane_plateaus_rsl"] = ("time", np.ones(len(proc_time), dtype=bool))
        cml["flag_blackout_fill"] = ("time", np.zeros(len(proc_time), dtype=bool))

    cml["trsl"] = cml["tsl"] - cml["rsl"]
    cml["trsl_raw"] = cml.trsl.copy()

    if proc_line in ["graf_2020", "full"]:
        cml["trsl"] = cml.trsl.interpolate_na(dim="time", method="linear", max_gap="5min")

    # General CML labels.
    det_limit = detection_limit_from_kR(
        L_km=float(length_km),
        freq_Hz=float(frequency_hz),
        pol=polarization,
    )
    cml["detection_limit"] = det_limit
    cml["isSane_detection_limit"] = bool(det_limit < 2)
    cml["isSane_std_long"] = temporal_sanity_check(cml.trsl, tspan=300, thld=2, perc=0.1)
    cml["isSane_std_short"] = temporal_sanity_check(cml.trsl, tspan=60, thld=0.8, perc=0.33)
    cml["isSane_available"] = availability_check(cml.trsl, perc_thld=0.5)

    if proc_line in ["graf_2020"]:
        is_sane = bool(cml.isSane_std_long.values) and bool(cml.isSane_std_short.values)
    elif proc_line in ["full"]:
        is_sane = (
            bool(cml.isSane_std_long.values)
            and bool(cml.isSane_std_short.values)
            and bool(cml.isSane_detection_limit.values)
            and bool(cml.isSane_available.values)
        )
    else:
        is_sane = True

    info["detection_limit"] = float(det_limit)
    info["isSane"] = int(is_sane)
    info["isSane_std_long"] = int(bool(cml.isSane_std_long.values))
    info["isSane_std_short"] = int(bool(cml.isSane_std_short.values))
    info["isSane_detection_limit"] = int(bool(cml.isSane_detection_limit.values))
    info["isSane_available"] = int(bool(cml.isSane_available.values))
    info["blackout_filled"] = int(np.asarray(cml.flag_blackout_fill.values).sum())

    if not is_sane:
        return np.full(out_rel1 - out_rel0, np.nan, dtype="float32"), info

    try:
        cml = derive_rain(cml, waa_model=waa_model)
        R = np.asarray(cml.R_cml.values, dtype="float64")
        R = np.where(np.isfinite(R) & (R >= 0), R, np.nan)
        return R[out_rel0:out_rel1].astype("float32"), info
    except Exception as exc:
        info["error"] = str(exc)
        return np.full(out_rel1 - out_rel0, np.nan, dtype="float32"), info


# ============================================================
# IDW and plotting
# ============================================================
def build_grid(lon: np.ndarray, lat: np.ndarray, resolution: float, pad: float):
    xmin, xmax = np.nanmin(lon) - pad, np.nanmax(lon) + pad
    ymin, ymax = np.nanmin(lat) - pad, np.nanmax(lat) + pad
    x = np.arange(xmin, xmax + resolution, resolution)
    y = np.arange(ymin, ymax + resolution, resolution)
    return x, y, np.meshgrid(x, y)


def idw_grid(lon, lat, z, xg, yg, nnear=15, p=2, max_dist=0.3):
    from scipy.spatial import cKDTree
    lon = np.asarray(lon, dtype="float64")
    lat = np.asarray(lat, dtype="float64")
    z = np.asarray(z, dtype="float64")
    valid = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(z)
    out = np.full(xg.shape, np.nan, dtype="float32")
    if valid.sum() == 0:
        return out
    xy = np.c_[lon[valid], lat[valid]]
    tree = cKDTree(xy)
    q = np.c_[xg.ravel(), yg.ravel()]
    k = min(int(nnear), xy.shape[0])
    d, ix = tree.query(q, k=k)
    if k == 1:
        d = d[:, None]
        ix = ix[:, None]
    far = d[:, 0] > max_dist if max_dist is not None else np.zeros(d.shape[0], dtype=bool)
    exact = d[:, 0] < 1e-12
    with np.errstate(divide="ignore", invalid="ignore"):
        w = 1.0 / (d ** p)
        w = w / np.sum(w, axis=1, keepdims=True)
    vals = np.sum(z[valid][ix] * w, axis=1)
    vals[exact] = z[valid][ix[exact, 0]]
    vals[far] = np.nan
    return vals.reshape(xg.shape).astype("float32")


def make_gif(ds_grid: xr.Dataset, out_gif: str, tstart: str, tend: str, frame_freq: str, vmax: Optional[float]):
    if imageio is None:
        print("imageio not installed; skipping GIF")
        return
    out_gif = Path(out_gif)
    frames_dir = out_gif.with_suffix("")
    frames_dir.mkdir(parents=True, exist_ok=True)

    ds_sel = ds_grid.sel(time=slice(pd.Timestamp(tstart), pd.Timestamp(tend)))
    target_times = pd.date_range(pd.Timestamp(tstart), pd.Timestamp(tend), freq=frame_freq)
    available = pd.DatetimeIndex(ds_sel.time.values)
    times = [t for t in target_times if t in available]
    if not times:
        times = list(available)

    vars_to_plot = []
    for v, title in [("R_cml", "CML"), ("R_radar", "Radar path avg"), ("R_gauge", "Gauge")]:
        if v in ds_sel.data_vars:
            vars_to_plot.append((v, title))

    if vmax is None:
        vals = []
        for v, _ in vars_to_plot:
            a = ds_sel[v].values
            a = a[np.isfinite(a)]
            if a.size:
                vals.append(a)
        vmax = max(5.0, float(np.nanpercentile(np.concatenate(vals), 99))) if vals else 10.0

    pngs = []
    for t in times:
        fig, axes = plt.subplots(1, len(vars_to_plot), figsize=(5.2 * len(vars_to_plot), 4.6), squeeze=False)
        axes = axes[0]
        for ax, (v, title) in zip(axes, vars_to_plot):
            data = ds_sel[v].sel(time=t).values.astype("float64")
            data[data < 0.1] = np.nan
            im = ax.pcolormesh(ds_sel.x.values, ds_sel.y.values, data, shading="auto", cmap="turbo", vmin=0.1, vmax=vmax)
            ax.set_title(title)
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.set_aspect("equal", adjustable="box")
            fig.colorbar(im, ax=ax, label="mm h$^{-1}$", shrink=0.78)
        fig.suptitle(pd.Timestamp(t).strftime("%Y-%m-%d %H:%M UTC"))
        fig.tight_layout()
        png = frames_dir / f"frame_{pd.Timestamp(t).strftime('%Y%m%d_%H%M')}.png"
        fig.savefig(png, dpi=130, bbox_inches="tight")
        plt.close(fig)
        pngs.append(png)

    imageio.mimsave(out_gif, [imageio.imread(p) for p in pngs], duration=0.35)
    print(f"Saved GIF: {out_gif}")


# ============================================================
# Main pipeline
# ============================================================
def run(args):
    t0 = time.perf_counter()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    est_nc = out_dir / args.estimation_nc
    grid_nc = out_dir / args.gridded_nc
    log_csv = out_dir / args.log_csv
    gif_out = out_dir / args.gif_name if args.gif_name else None

    print("Opening datasets...")
    ds_c = xr.open_dataset(args.cml_nc)
    ds_r = xr.open_dataset(args.radar_nc) if args.radar_nc else None
    ds_g = xr.open_dataset(args.gauge_nc) if args.gauge_nc else None

    try:
        cml_dim = find_dim(ds_c, ["cml_id", "cml", "link", "link_id"])
        time_dim = find_dim(ds_c, ["time", "t"])
        sublink_dim = find_dim(ds_c, ["sublink_id", "sublink", "channel"])
        if cml_dim is None or time_dim is None:
            raise ValueError(f"Cannot find CML/time dims. dims={list(ds_c.dims)}")

        tx_name = find_var(ds_c, ["tsl", "tx", "TSL"])
        rx_name = find_var(ds_c, ["rsl", "rx", "RSL"])
        if tx_name is None or rx_name is None:
            raise ValueError(f"Cannot find tx/rx or tsl/rsl variables. data_vars={list(ds_c.data_vars)}")

        len_name = find_var(ds_c, ["length", "path_length_km", "L_km"])
        freq_name = find_var(ds_c, ["frequency", "freq", "frequency_GHz", "frequency_MHz"])
        pol_name = find_var(ds_c, ["polarization", "polarisation", "pol"])
        lat_a_name = find_var(ds_c, ["site_a_latitude", "lat_a", "latitude_a"])
        lon_a_name = find_var(ds_c, ["site_a_longitude", "lon_a", "longitude_a"])
        lat_b_name = find_var(ds_c, ["site_b_latitude", "lat_b", "latitude_b"])
        lon_b_name = find_var(ds_c, ["site_b_longitude", "lon_b", "longitude_b"])

        if len_name is None or freq_name is None or pol_name is None:
            raise ValueError("Missing length/frequency/polarization variables in CML dataset.")
        if None in [lat_a_name, lon_a_name, lat_b_name, lon_b_name]:
            raise ValueError("Missing CML endpoint coordinates.")

        cml_ids = as_str_array(ds_c[cml_dim].values)
        if sublink_dim is not None:
            sublink_ids = np.asarray(ds_c[sublink_dim].values)
            n_sublink = len(sublink_ids)
        else:
            sublink_ids = np.array([0])
            n_sublink = 1

        full_time = pd.DatetimeIndex(ds_c[time_dim].values)
        tstart = pd.Timestamp(args.tstart)
        tend = pd.Timestamp(args.tend)
        gtstart = pd.Timestamp(args.grid_tstart)
        gtend = pd.Timestamp(args.grid_tend)

        # Extra padding for centered rolling windows and baseline access.
        p0 = tstart - pd.Timedelta(minutes=args.pad_minutes)
        p1 = tend + pd.Timedelta(minutes=args.pad_minutes)
        proc0, proc1 = time_slice_indices(full_time, p0, p1)
        out0, out1 = time_slice_indices(full_time, tstart, tend)
        if proc1 <= proc0 or out1 <= out0:
            raise ValueError("Empty time window.")
        proc_time = full_time[proc0:proc1]
        out_time = full_time[out0:out1]
        out_rel0 = out0 - proc0
        out_rel1 = out1 - proc0
        time_5min = pd.date_range(tstart.floor("5min"), tend.floor("5min"), freq="5min")
        time_5min = time_5min[(time_5min >= tstart) & (time_5min <= tend)]
        grid_times = time_5min[(time_5min >= gtstart) & (time_5min <= gtend)]

        n_cml = len(cml_ids)
        if args.max_links is not None:
            n_cml = min(n_cml, int(args.max_links))
            cml_ids = cml_ids[:n_cml]

        print("Configuration")
        print(f"  CML file     : {args.cml_nc}")
        print(f"  tx/rx vars   : {tx_name}, {rx_name}")
        print(f"  proc_line    : {args.proc_line}")
        print(f"  WAA model    : {args.waa_model}")
        print(f"  CMLs         : {n_cml}")
        print(f"  sublinks     : {n_sublink}")
        print(f"  process time : {proc_time[0]} -> {proc_time[-1]}")
        print(f"  output time  : {out_time[0]} -> {out_time[-1]}")
        print(f"  grid time    : {grid_times[0]} -> {grid_times[-1]}")

        R_sub_5 = np.full((n_cml, n_sublink, len(time_5min)), np.nan, dtype="float32")
        lengths = np.full(n_cml, np.nan, dtype="float32")
        freqs_hz = np.full((n_cml, n_sublink), np.nan, dtype="float64")
        pols = np.empty((n_cml, n_sublink), dtype=object)
        mid_lon = np.full(n_cml, np.nan, dtype="float64")
        mid_lat = np.full(n_cml, np.nan, dtype="float64")
        logs: List[Dict] = []

        for i in tqdm(range(n_cml), desc="CML links"):
            cid = cml_ids[i]
            length_km = safe_float(select_scalar_or_index(ds_c[len_name], cml_dim, i, sublink_dim, 0))
            lengths[i] = length_km
            lon_a = safe_float(select_scalar_or_index(ds_c[lon_a_name], cml_dim, i, None, 0))
            lon_b = safe_float(select_scalar_or_index(ds_c[lon_b_name], cml_dim, i, None, 0))
            lat_a = safe_float(select_scalar_or_index(ds_c[lat_a_name], cml_dim, i, None, 0))
            lat_b = safe_float(select_scalar_or_index(ds_c[lat_b_name], cml_dim, i, None, 0))
            mid_lon[i] = 0.5 * (lon_a + lon_b)
            mid_lat[i] = 0.5 * (lat_a + lat_b)

            if not np.isfinite(length_km) or length_km < args.min_link_km:
                logs.append({"cml_id": cid, "status": "skip_length", "length_km": length_km})
                continue

            for sl in range(n_sublink):
                tsl = extract_signal_1d(ds_c[tx_name], time_dim, cml_dim, i, sublink_dim, sl, proc0, proc1)
                rsl = extract_signal_1d(ds_c[rx_name], time_dim, cml_dim, i, sublink_dim, sl, proc0, proc1)
                f_raw = safe_float(select_scalar_or_index(ds_c[freq_name], cml_dim, i, sublink_dim, sl))
                f_hz = frequency_to_hz(f_raw)
                pol = pol_to_str(select_scalar_or_index(ds_c[pol_name], cml_dim, i, sublink_dim, sl))
                freqs_hz[i, sl] = f_hz
                pols[i, sl] = pol

                if np.isfinite(tsl).sum() < args.min_valid_samples or np.isfinite(rsl).sum() < args.min_valid_samples:
                    logs.append({"cml_id": cid, "sublink_id": int(sl), "status": "skip_no_signal"})
                    continue

                R1, info = process_one_sublink(
                    tsl=tsl,
                    rsl=rsl,
                    proc_time=proc_time,
                    out_rel0=out_rel0,
                    out_rel1=out_rel1,
                    length_km=length_km,
                    frequency_hz=f_hz,
                    polarization=pol,
                    proc_line=args.proc_line,
                    waa_model=args.waa_model,
                )
                R_sub_5[i, sl, :] = resample_1min_rate_to_5min_rate(R1, out_time, time_5min)

                row = {"cml_id": cid, "sublink_id": int(sl), "status": "processed", "length_km": length_km,
                       "frequency_hz": f_hz, "polarization": pol,
                       "rr_finite_1min": int(np.isfinite(R1).sum()),
                       "rr_positive_1min": int(np.sum(np.isfinite(R1) & (R1 > 0)))}
                row.update(info)
                logs.append(row)

            if i % 50 == 0:
                gc.collect()

        with np.errstate(invalid="ignore"):
            R_mean_5 = np.nanmean(R_sub_5, axis=1).astype("float32")

        # Radar path average aligned to CML IDs and 5 min time.
        R_radar = np.full((n_cml, len(time_5min)), np.nan, dtype="float32")
        if ds_r is not None:
            r_cml_dim = find_dim(ds_r, ["cml_id", "cml", "link", "link_id"])
            r_time_dim = find_dim(ds_r, ["time", "t"])
            r_var = find_var(ds_r, ["rr", "R_radar", "rain_rate", "rainfall"])
            if r_cml_dim is None or r_time_dim is None or r_var is None:
                print("WARNING: radar dims/var not found; skipping radar extraction")
            else:
                ds_r2 = ds_r.assign_coords({r_cml_dim: as_str_array(ds_r[r_cml_dim].values)})
                radar_time = pd.DatetimeIndex(ds_r2[r_time_dim].values)
                for i, cid in enumerate(tqdm(cml_ids, desc="Radar align")):
                    try:
                        da = ds_r2[r_var].sel({r_cml_dim: cid})
                        s = pd.Series(np.asarray(da.values, dtype="float64") * args.radar_factor, index=radar_time)
                        R_radar[i, :] = s.reindex(time_5min).to_numpy(dtype="float32")
                    except Exception:
                        pass

        ds_est = xr.Dataset(
            data_vars={
                "R_sublink_5min": (("cml_id", "sublink_id", "time_5min"), R_sub_5),
                "R_mean_5min": (("cml_id", "time_5min"), R_mean_5),
                "R_radar": (("cml_id", "time_5min"), R_radar),
                "length": (("cml_id",), lengths),
                "frequency_hz": (("cml_id", "sublink_id"), freqs_hz.astype("float32")),
                "mid_lon": (("cml_id",), mid_lon.astype("float32")),
                "mid_lat": (("cml_id",), mid_lat.astype("float32")),
            },
            coords={"cml_id": cml_ids, "sublink_id": sublink_ids, "time_5min": time_5min.values},
            attrs={
                "processing": "Example_Processing notebook logic: plateau filter, blackout filling, TRSL interpolation, sanity labels, derive_rain",
                "proc_line": args.proc_line,
                "waa_model": args.waa_model,
                "radar_factor": args.radar_factor,
            },
        )
        ds_est.to_netcdf(est_nc)
        pd.DataFrame(logs).to_csv(log_csv, index=False)
        print(f"Saved estimation NetCDF: {est_nc}")
        print(f"Saved log CSV: {log_csv}")

        # Gridding for event time window.
        grid_idx = np.where((time_5min >= gtstart) & (time_5min <= gtend))[0]
        if len(grid_idx) == 0:
            raise ValueError("No 5 min times inside grid window")

        good_cml_xy = np.isfinite(mid_lon) & np.isfinite(mid_lat)
        all_lon = [mid_lon[good_cml_xy]]
        all_lat = [mid_lat[good_cml_xy]]

        gauge_meta = None
        if ds_g is not None:
            g_time_dim = find_dim(ds_g, ["time", "t"])
            g_var = find_var(ds_g, ["rr", "rain", "rainfall", "precip", "R"])
            g_lat_name = find_var(ds_g, ["lat", "latitude"])
            g_lon_name = find_var(ds_g, ["lon", "longitude"])
            if g_time_dim and g_var and g_lat_name and g_lon_name:
                g_da = ds_g[g_var]
                g_station_dim = [d for d in g_da.dims if d != g_time_dim][0]
                glon = np.asarray(ds_g[g_lon_name].values, dtype="float64")
                glat = np.asarray(ds_g[g_lat_name].values, dtype="float64")
                if glon.ndim > 1: glon = glon.reshape(-1)
                if glat.ndim > 1: glat = glat.reshape(-1)
                gauge_meta = (g_da, g_time_dim, g_station_dim, glon, glat)
                all_lon.append(glon[np.isfinite(glon) & np.isfinite(glat)])
                all_lat.append(glat[np.isfinite(glon) & np.isfinite(glat)])
            else:
                print("WARNING: gauge variables not found; skipping gauges")

        x, y, (xg, yg) = build_grid(np.concatenate(all_lon), np.concatenate(all_lat), args.resolution, args.grid_pad)
        R_cml_grid = np.full((len(grid_idx), len(y), len(x)), np.nan, dtype="float32")
        R_radar_grid = np.full_like(R_cml_grid, np.nan)
        R_gauge_grid = np.full_like(R_cml_grid, np.nan) if gauge_meta is not None else None

        for jj, ti in enumerate(tqdm(grid_idx, desc="IDW maps")):
            R_cml_grid[jj] = idw_grid(mid_lon, mid_lat, R_mean_5[:, ti], xg, yg, args.nnear, args.p, args.max_dist)
            R_radar_grid[jj] = idw_grid(mid_lon, mid_lat, R_radar[:, ti], xg, yg, args.nnear, args.p, args.max_dist)

            if gauge_meta is not None:
                g_da, g_time_dim, g_station_dim, glon, glat = gauge_meta
                gt = pd.Timestamp(time_5min[ti])
                try:
                    # Gauge is 10 min accumulation in mm. Convert to mm h-1 with factor 6.
                    gg = g_da.sel({g_time_dim: gt}, method="nearest", tolerance=pd.Timedelta(args.gauge_tolerance)).values
                    gg = np.asarray(gg, dtype="float64").reshape(-1) * args.gauge_factor
                    R_gauge_grid[jj] = idw_grid(glon, glat, gg, xg, yg, args.nnear, args.p, args.max_dist)
                except Exception:
                    pass

        data_vars = {
            "R_cml": (("time", "y", "x"), R_cml_grid, {"units": "mm h-1"}),
            "R_radar": (("time", "y", "x"), R_radar_grid, {"units": "mm h-1"}),
        }
        if R_gauge_grid is not None:
            data_vars["R_gauge"] = (("time", "y", "x"), R_gauge_grid, {"units": "mm h-1"})

        ds_grid = xr.Dataset(
            data_vars=data_vars,
            coords={"time": time_5min[grid_idx].values, "x": x, "y": y},
            attrs={"interpolation": f"IDW nnear={args.nnear}, p={args.p}, max_dist={args.max_dist}"},
        )
        ds_grid.to_netcdf(grid_nc)
        print(f"Saved gridded NetCDF: {grid_nc}")

        if gif_out is not None:
            make_gif(ds_grid, str(gif_out), args.grid_tstart, args.grid_tend, args.frame_freq, args.vmax)

        print(f"Done in {(time.perf_counter() - t0)/60:.1f} min")

    finally:
        ds_c.close()
        if ds_r is not None:
            ds_r.close()
        if ds_g is not None:
            ds_g.close()


def parse_args():
    ap = argparse.ArgumentParser(description="Blettner/Graf notebook-style CML rainfall event maps")
    ap.add_argument("--cml-nc", required=True)
    ap.add_argument("--radar-nc", default=None)
    ap.add_argument("--gauge-nc", default=None)
    ap.add_argument("--out-dir", default="outputs/blettner_event_20210621")
    ap.add_argument("--estimation-nc", default="cml_rainrate_blettner_full_event.nc")
    ap.add_argument("--gridded-nc", default="rainfall_maps_blettner_full_event.nc")
    ap.add_argument("--log-csv", default="cml_rainrate_blettner_full_event_log.csv")
    ap.add_argument("--gif-name", default="rainfall_maps_blettner_full_event.gif")

    ap.add_argument("--tstart", default="2021-06-21 12:00:00")
    ap.add_argument("--tend", default="2021-06-22 05:30:00")
    ap.add_argument("--grid-tstart", default="2021-06-21 21:50:00")
    ap.add_argument("--grid-tend", default="2021-06-22 04:50:00")
    ap.add_argument("--pad-minutes", type=int, default=360)

    ap.add_argument("--proc-line", choices=["no_filter", "graf_2020", "full"], default="full")
    ap.add_argument("--waa-model", choices=["Leijnse", "Pastorek", "Schleiss"], default="Leijnse")
    ap.add_argument("--min-link-km", type=float, default=0.5)
    ap.add_argument("--min-valid-samples", type=int, default=100)

    ap.add_argument("--radar-factor", type=float, default=12.0, help="Multiply radar rr by this factor to get mm h-1")
    ap.add_argument("--gauge-factor", type=float, default=6.0, help="10 min accumulation mm to mm h-1")
    ap.add_argument("--gauge-tolerance", default="5min")

    ap.add_argument("--resolution", type=float, default=0.01)
    ap.add_argument("--grid-pad", type=float, default=0.05)
    ap.add_argument("--nnear", type=int, default=15)
    ap.add_argument("--p", type=int, default=2)
    ap.add_argument("--max-dist", type=float, default=0.3)
    ap.add_argument("--frame-freq", default="10min")
    ap.add_argument("--vmax", type=float, default=None)
    ap.add_argument("--max-links", type=int, default=None)
    return ap.parse_args()


if __name__ == "__main__":
    run(parse_args())
