#!/usr/bin/env python3
"""
Merged CML rainfall estimation + streamed parallel gridding for CML, radar, and rain gauges.

This file combines the user's tested rainfall-estimation pipeline with the streamed
128-GB-safe IDW gridding script, and adds rain-gauge IDW interpolation on the same
regular grid. The scientific estimation methods are copied from the original script;
the new code only orchestrates stages and streams/parallelizes gridding chunks.
"""
from __future__ import annotations

"""
CML rainfall-estimation pipeline with:
  * multiple baseline methods
  * multiple k-R / ITU models
  * multiple wet-antenna attenuation corrections

Main new options
----------------
baseline_method:
  pycml_constant : pycomlink baseline_constant
  last_dry       : custom optimized last-dry baseline

kr_model:
  ITU_2003
  ITU_2005

waa_method:
  none
  leijnse_2008
  pastorek_2021
  schleiss_2013

Output dimensions, recommended 5-min product:
  R_sublink_5min(cml_id, sublink_id, baseline_method, kr_model, waa_method, time_5min)
  R_mean_5min(cml_id, baseline_method, kr_model, waa_method, time_5min)
  R_radar(cml_id, time_5min)

Optional 1-min output:
  Add --write-1min to also write:
  R_sublink(cml_id, sublink_id, baseline_method, kr_model, waa_method, time)
  R_mean(cml_id, baseline_method, kr_model, waa_method, time)

Important
---------
The custom last_dry baseline does NOT average over the 48 hour pad.
The 48 hour pad only gives the method access to dry samples before the output
period starts. The actual last_dry baseline at time t is the closest previous
finite dry TRSL value. Rain attenuation is then:

    A_obs = max(TRSL_current - baseline_last_dry, 0)

because this workflow uses TRSL = TSL - RSL, which increases during attenuation.
"""


import argparse
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple, List

import numpy as np
import pandas as pd
import xarray as xr

try:
    from tqdm import tqdm
except Exception:
    def tqdm(x, **kwargs):
        return x

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ============================================================
# Defaults
# ============================================================
N_AVERAGE_LAST_DRY_DEFAULT = 1
RADAR_WET_THR_5MIN_DEFAULT = 0.00833
MIN_LINK_KM_DEFAULT = 0.5
MIN_VALID_SAMPLES_DEFAULT = 100
BASELINE_PAD_HOURS_DEFAULT = 48

BASELINE_METHODS_ALLOWED = ("pycml_constant", "last_dry")
BASELINE_METHODS_DEFAULT = BASELINE_METHODS_ALLOWED

KR_MODELS_ALLOWED = ("ITU_2003", "ITU_2005")
KR_MODELS_DEFAULT = KR_MODELS_ALLOWED

WAA_METHODS_ALLOWED = ("none", "leijnse_2008", "pastorek_2021", "schleiss_2013")
WAA_METHODS_DEFAULT = WAA_METHODS_ALLOWED

# Leijnse 2008 defaults from your tested workflow.
WAA_GAMMA_DEFAULT = (2.06e-05) / 1.4
WAA_L_ANT_DEFAULT = 0.0043
WAA_DELTA_DEFAULT = 0.24 * 1.5

# Pastorek 2021 defaults.
PASTOREK_A_MAX_DEFAULT = 14.0
PASTOREK_ZETA_DEFAULT = 0.55
PASTOREK_D_DEFAULT = 0.1

# Schleiss 2013 WAA defaults.
SCHLEISS_WAA_MAX_DEFAULT = 2.0
SCHLEISS_DELTA_T_DEFAULT = 1.0
SCHLEISS_TAU_DEFAULT = 15.0


# ============================================================
# Small helpers
# ============================================================
def find_dim(ds: xr.Dataset, candidates: Sequence[str]) -> Optional[str]:
    for c in candidates:
        if c in ds.dims:
            return c
    return None


def variable_exists(ds: xr.Dataset, name: str) -> bool:
    return name in ds.variables or name in ds.data_vars


def pol_to_HV(pol_raw) -> str:
    if isinstance(pol_raw, (bytes, np.bytes_)):
        pol_raw = pol_raw.decode("utf-8", errors="ignore")
    s = str(pol_raw).strip().lower()
    if "v" in s and "h" not in s:
        return "V"
    if "h" in s and "v" not in s:
        return "H"
    return "H"


def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def time_slice_indices(times: pd.DatetimeIndex, tstart: pd.Timestamp, tend: pd.Timestamp) -> Tuple[int, int]:
    vals = times.values
    i0 = int(np.searchsorted(vals, np.datetime64(tstart), "left"))
    i1 = int(np.searchsorted(vals, np.datetime64(tend), "right"))
    return i0, i1


def radar_values_to_mm_h(raw: np.ndarray) -> np.ndarray:
    """Exact convention from your working event script: radar_mmph = radar_raw * 12."""
    return np.asarray(raw, dtype="float64") * 12.0


def forward_fill_radar_wet_to_cml_time(
    radar_time: pd.DatetimeIndex,
    cml_time: pd.DatetimeIndex,
    wet_5min: np.ndarray,
) -> np.ndarray:
    idx = np.searchsorted(radar_time.values, cml_time.values, "right") - 1
    pre_first = idx < 0
    idx = np.clip(idx, 0, None)
    out = wet_5min[idx].astype(bool)
    out[pre_first] = False
    return out


def resample_1min_rate_to_5min_rate(
    values_1min: np.ndarray,
    time_1min: pd.DatetimeIndex,
    target_5min: pd.DatetimeIndex,
) -> np.ndarray:
    """Rate -> amount per minute -> 5-min sum -> rate."""
    s = pd.Series(values_1min, index=time_1min)
    out = (s / 60.0).resample("5min", label="left", closed="left").sum(min_count=1) * 12.0
    return out.reindex(target_5min).to_numpy(dtype="float32")


def normalize_2d_time_sublink(da: xr.DataArray, sublink_dim: str) -> np.ndarray:
    if sublink_dim in da.dims:
        return da.transpose("time", sublink_dim).values
    return da.transpose("time").values[:, None]


def finite_quantiles(a: np.ndarray) -> Dict[str, float]:
    a = np.asarray(a, dtype="float64")
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {"min": np.nan, "p50": np.nan, "p95": np.nan, "p99": np.nan, "max": np.nan}
    q = np.nanpercentile(a, [50, 95, 99])
    return {
        "min": float(np.nanmin(a)),
        "p50": float(q[0]),
        "p95": float(q[1]),
        "p99": float(q[2]),
        "max": float(np.nanmax(a)),
    }


# ============================================================
# Custom baseline
# ============================================================
def baseline_last_dry_optimized(trsl: np.ndarray, wet: np.ndarray) -> np.ndarray:
    """Closest previous finite dry TRSL baseline.

    This is the exact custom logic:
      * During dry finite samples, the baseline is the current TRSL.
      * During wet samples, the baseline is the closest previous finite dry TRSL.
      * If no previous dry value exists, baseline is NaN.
      * Missing dry TRSL values do not update the last dry reference.

    It is vectorized using np.maximum.accumulate, so it is fast for long 1-min
    time series.

    Because this workflow uses TRSL = TSL - RSL, attenuation is:
        A_obs = max(TRSL - baseline, 0)
    """
    trsl = np.asarray(trsl, dtype="float64")
    wet = np.asarray(wet, dtype=bool)
    n = trsl.size

    dry_valid = (~wet) & np.isfinite(trsl)
    dry_idx = np.where(dry_valid, np.arange(n, dtype=np.int64), -1)
    last_dry_idx = np.maximum.accumulate(dry_idx)

    baseline = np.full(n, np.nan, dtype="float64")
    ok = last_dry_idx >= 0
    baseline[ok] = trsl[last_dry_idx[ok]]
    return baseline


# ============================================================
# Backend
# ============================================================
@dataclass
class EstimationResult:
    R: Dict[Tuple[str, str, str], np.ndarray]          # (baseline, kr, waa) -> R
    waa: Dict[Tuple[str, str], np.ndarray]             # (baseline, waa) -> waa
    A_obs: Dict[str, np.ndarray]                       # baseline -> A_obs
    A_rain: Dict[Tuple[str, str], np.ndarray]          # (baseline, waa) -> A_rain


class PycomlinkBackend:
    def __init__(
        self,
        baseline_methods: Sequence[str],
        kr_models: Sequence[str],
        waa_methods: Sequence[str],
        n_average_last_dry: int,
        waa_gamma: float,
        waa_l_antenna: float,
        waa_delta: float,
        pastorek_a_max: float,
        pastorek_zeta: float,
        pastorek_d: float,
        schleiss_waa_max: float,
        schleiss_delta_t: float,
        schleiss_tau: float,
    ):
        import pycomlink as pycml
        import pycomlink.processing.baseline
        import pycomlink.processing.wet_antenna
        import pycomlink.processing.k_R_relation

        self.pycml = pycml
        self.baseline_methods = tuple(baseline_methods)
        self.kr_models = tuple(kr_models)
        self.waa_methods = tuple(waa_methods)
        self.n_average_last_dry = int(n_average_last_dry)

        self.waa_gamma = float(waa_gamma)
        self.waa_l_antenna = float(waa_l_antenna)
        self.waa_delta = float(waa_delta)

        self.pastorek_a_max = float(pastorek_a_max)
        self.pastorek_zeta = float(pastorek_zeta)
        self.pastorek_d = float(pastorek_d)

        self.schleiss_waa_max = float(schleiss_waa_max)
        self.schleiss_delta_t = float(schleiss_delta_t)
        self.schleiss_tau = float(schleiss_tau)

    def _baseline(self, method: str, xr_trsl: xr.DataArray, xr_wet: xr.DataArray) -> xr.DataArray:
        if method == "pycml_constant":
            return self.pycml.processing.baseline.baseline_constant(
                trsl=xr_trsl,
                wet=xr_wet,
                n_average_last_dry=self.n_average_last_dry,
            )
        if method == "last_dry":
            b = baseline_last_dry_optimized(xr_trsl.values, xr_wet.values)
            return xr.DataArray(b, coords=xr_trsl.coords, dims=xr_trsl.dims)
        raise ValueError(f"Unknown baseline method: {method}")

    def _calc_R(self, A_rain: xr.DataArray, kr_model: str, f_GHz: float, pol: str, length_km: float) -> xr.DataArray:
        if kr_model not in KR_MODELS_ALLOWED:
            raise ValueError(f"Unknown k-R model: {kr_model}")

        # pycomlink versions differ in keyword name.
        try:
            return self.pycml.processing.k_R_relation.calc_R_from_A(
                A=A_rain,
                L_km=length_km,
                f_GHz=f_GHz,
                pol=pol,
                a_b_approximation=kr_model,
            )
        except TypeError:
            return self.pycml.processing.k_R_relation.calc_R_from_A(
                A=A_rain,
                L_km=length_km,
                f_GHz=f_GHz,
                pol=pol,
                approx_type=kr_model,
            )

    @staticmethod
    def _clip_waa(waa: xr.DataArray, A_obs: xr.DataArray) -> xr.DataArray:
        waa = xr.where(waa < 0, 0.0, waa)
        waa = xr.where(waa > A_obs, A_obs, waa)
        return waa

    def _waa(self, method: str, A_obs: xr.DataArray, xr_trsl: xr.DataArray,
             baseline: xr.DataArray, xr_wet: xr.DataArray, f_GHz: float,
             pol: str, length_km: float) -> xr.DataArray:
        if method == "none":
            waa = xr.zeros_like(A_obs)
        elif method == "leijnse_2008":
            waa = self.pycml.processing.wet_antenna.waa_leijnse_2008_from_A_obs(
                A_obs=A_obs,
                f_Hz=f_GHz * 1e9,
                pol=pol,
                L_km=length_km,
                gamma=self.waa_gamma,
                l_antenna=self.waa_l_antenna,
                delta=self.waa_delta,
            )
        elif method == "pastorek_2021":
            waa = self.pycml.processing.wet_antenna.waa_pastorek_2021_from_A_obs(
                A_obs=A_obs,
                f_Hz=f_GHz * 1e9,
                pol=pol,
                L_km=length_km,
                A_max=self.pastorek_a_max,
                zeta=self.pastorek_zeta,
                d=self.pastorek_d,
            )
        elif method == "schleiss_2013":
            waa = self.pycml.processing.wet_antenna.waa_schleiss_2013(
                rsl=xr_trsl,
                baseline=baseline,
                wet=xr_wet,
                waa_max=self.schleiss_waa_max,
                delta_t=self.schleiss_delta_t,
                tau=self.schleiss_tau,
            )
        else:
            raise ValueError(f"Unknown WAA method: {method}")

        return self._clip_waa(waa, A_obs)

    def estimate(
        self,
        trsl: np.ndarray,
        wet: np.ndarray,
        time_index: pd.DatetimeIndex,
        f_GHz: float,
        pol: str,
        length_km: float,
    ) -> EstimationResult:
        xr_trsl = xr.DataArray(trsl.astype("float64"), coords={"time": time_index.values}, dims=["time"])
        xr_wet = xr.DataArray(wet.astype(bool), coords={"time": time_index.values}, dims=["time"])

        out_R: Dict[Tuple[str, str, str], np.ndarray] = {}
        out_waa: Dict[Tuple[str, str], np.ndarray] = {}
        out_A_obs: Dict[str, np.ndarray] = {}
        out_A_rain: Dict[Tuple[str, str], np.ndarray] = {}

        for b_method in self.baseline_methods:
            baseline = self._baseline(b_method, xr_trsl, xr_wet)
            A_obs = xr.where((xr_trsl - baseline) < 0, 0.0, xr_trsl - baseline)
            out_A_obs[b_method] = np.asarray(A_obs.values, dtype="float32")

            for w_method in self.waa_methods:
                waa = self._waa(w_method, A_obs, xr_trsl, baseline, xr_wet, f_GHz, pol, length_km)
                A_rain = xr.where((A_obs - waa) < 0, 0.0, A_obs - waa)

                out_waa[(b_method, w_method)] = np.asarray(waa.values, dtype="float32")
                out_A_rain[(b_method, w_method)] = np.asarray(A_rain.values, dtype="float32")

                for kr_model in self.kr_models:
                    R = self._calc_R(A_rain=A_rain, kr_model=kr_model, f_GHz=f_GHz, pol=pol, length_km=length_km)
                    out_R[(b_method, kr_model, w_method)] = np.asarray(R.values, dtype="float32")

        return EstimationResult(R=out_R, waa=out_waa, A_obs=out_A_obs, A_rain=out_A_rain)


# ============================================================
# Streaming writer
# ============================================================
class StreamingNetCDFWriter:
    def __init__(
        self,
        cml_ids,
        sublink_ids,
        baseline_methods,
        kr_models,
        waa_methods,
        time_out,
        time_5min,
        out_path,
        attrs,
        write_1min=False,
        overwrite=False,
        complevel=4,
        chunk_time=100_000,
    ):
        try:
            import netCDF4
        except Exception as exc:
            raise RuntimeError("Streaming output requires the netCDF4 Python package.") from exc

        self.out_path = Path(out_path)
        self.out_path.parent.mkdir(parents=True, exist_ok=True)
        if self.out_path.exists() and overwrite:
            self.out_path.unlink()
        if self.out_path.exists():
            raise FileExistsError(f"Output exists: {self.out_path}. Use --overwrite.")

        self.write_1min = bool(write_1min)

        self.nc = netCDF4.Dataset(self.out_path, "w", format="NETCDF4")
        C = len(cml_ids)
        S = len(sublink_ids)
        B = len(baseline_methods)
        K = len(kr_models)
        W = len(waa_methods)
        T = len(time_out)
        T5 = len(time_5min)

        self.C, self.S, self.B, self.K, self.W, self.T, self.T5 = C, S, B, K, W, T, T5

        self.nc.createDimension("cml_id", C)
        self.nc.createDimension("sublink_id", S)
        self.nc.createDimension("baseline_method", B)
        self.nc.createDimension("kr_model", K)
        self.nc.createDimension("waa_method", W)
        self.nc.createDimension("time_5min", T5)
        if self.write_1min:
            self.nc.createDimension("time", T)

        self.nc.createVariable("cml_id", str, ("cml_id",))[:] = np.asarray(cml_ids, dtype=object)
        self.nc.createVariable("sublink_id", "i4", ("sublink_id",))[:] = np.asarray(sublink_ids, dtype="int32")
        self.nc.createVariable("baseline_method", str, ("baseline_method",))[:] = np.asarray(baseline_methods, dtype=object)
        self.nc.createVariable("kr_model", str, ("kr_model",))[:] = np.asarray(kr_models, dtype=object)
        self.nc.createVariable("waa_method", str, ("waa_method",))[:] = np.asarray(waa_methods, dtype=object)

        units = "seconds since 1970-01-01 00:00:00 UTC"
        v_t5 = self.nc.createVariable("time_5min", "i8", ("time_5min",))
        v_t5.units = units
        v_t5.calendar = "standard"
        v_t5[:] = (pd.DatetimeIndex(time_5min).view("int64") // 10**9).astype("int64")
        if self.write_1min:
            v_t = self.nc.createVariable("time", "i8", ("time",))
            v_t.units = units
            v_t.calendar = "standard"
            v_t[:] = (pd.DatetimeIndex(time_out).view("int64") // 10**9).astype("int64")

        zlib = True
        ct = min(int(chunk_time), max(T, 1))
        ct5 = min(max(1, int(chunk_time // 5)), max(T5, 1))

        def f4var(name, dims, chunksizes):
            return self.nc.createVariable(
                name, "f4", dims, zlib=zlib, complevel=int(complevel),
                chunksizes=chunksizes, fill_value=np.float32(np.nan)
            )

        if self.write_1min:
            self.v_Rs = f4var(
                "R_sublink",
                ("cml_id", "sublink_id", "baseline_method", "kr_model", "waa_method", "time"),
                (1, 1, 1, 1, 1, ct),
            )
            self.v_Rs.units = "mm h-1"
            self.v_Rm = f4var(
                "R_mean",
                ("cml_id", "baseline_method", "kr_model", "waa_method", "time"),
                (1, 1, 1, 1, ct),
            )
            self.v_Rm.units = "mm h-1"
        else:
            self.v_Rs = None
            self.v_Rm = None

        self.v_Rs5 = f4var(
            "R_sublink_5min",
            ("cml_id", "sublink_id", "baseline_method", "kr_model", "waa_method", "time_5min"),
            (1, 1, 1, 1, 1, ct5),
        )
        self.v_Rs5.units = "mm h-1"

        self.v_Rm5 = f4var(
            "R_mean_5min",
            ("cml_id", "baseline_method", "kr_model", "waa_method", "time_5min"),
            (1, 1, 1, 1, ct5),
        )
        self.v_Rm5.units = "mm h-1"

        self.v_rad = f4var("R_radar", ("cml_id", "time_5min"), (1, ct5))
        self.v_rad.units = "mm h-1"

        self.v_waa = f4var(
            "waa_mean",
            ("cml_id", "sublink_id", "baseline_method", "waa_method"),
            (1, S, B, W),
        )
        self.v_waa.units = "dB"

        self.v_Aobs = f4var(
            "A_obs_mean",
            ("cml_id", "sublink_id", "baseline_method"),
            (1, S, B),
        )
        self.v_Aobs.units = "dB"

        self.v_Arain = f4var(
            "A_rain_mean",
            ("cml_id", "sublink_id", "baseline_method", "waa_method"),
            (1, S, B, W),
        )
        self.v_Arain.units = "dB"

        self.v_len = self.nc.createVariable("length", "f4", ("cml_id",), fill_value=np.float32(np.nan))
        self.v_len.units = "km"
        self.v_freq = self.nc.createVariable("frequency", "f4", ("cml_id", "sublink_id"), fill_value=np.float32(np.nan))
        self.v_freq.units = "GHz"
        self.v_pol = self.nc.createVariable("polarization", str, ("cml_id", "sublink_id"))

        for k, v in attrs.items():
            try:
                self.nc.setncattr(k, v)
            except Exception:
                self.nc.setncattr(k, str(v))

    def write_metadata(self, out_i, length, freq_list, pol_list):
        self.v_len[out_i] = np.float32(length)
        self.v_freq[out_i, :len(freq_list)] = np.asarray(freq_list, dtype="float32")
        self.v_pol[out_i, :len(pol_list)] = np.asarray(pol_list, dtype=object)

    def write_link(
        self,
        out_i,
        R_sub_5,
        R_mean_5,
        radar_5,
        waa_mean,
        A_obs_mean,
        A_rain_mean,
        R_sub=None,
        R_mean=None,
    ):
        self.v_Rs5[out_i, :, :, :, :, :] = R_sub_5.astype("float32")
        self.v_Rm5[out_i, :, :, :, :] = R_mean_5.astype("float32")
        self.v_rad[out_i, :] = radar_5.astype("float32")
        self.v_waa[out_i, :, :, :] = waa_mean.astype("float32")
        self.v_Aobs[out_i, :, :] = A_obs_mean.astype("float32")
        self.v_Arain[out_i, :, :, :] = A_rain_mean.astype("float32")
        if self.write_1min:
            self.v_Rs[out_i, :, :, :, :, :] = R_sub.astype("float32")
            self.v_Rm[out_i, :, :, :, :] = R_mean.astype("float32")

    def close(self):
        self.nc.close()


# ============================================================
# Config and pipeline
# ============================================================
@dataclass
class PipelineConfig:
    cml_nc: str
    radar_nc: str
    output_nc: str
    log_csv: str
    tstart: str
    tend: str
    radar_wet_thr_5min: float
    baseline_pad_hours: float
    min_link_km: float
    min_valid_samples: int
    n_average_last_dry: int
    baseline_methods: Tuple[str, ...]
    kr_models: Tuple[str, ...]
    waa_methods: Tuple[str, ...]
    waa_gamma: float
    waa_l_antenna: float
    waa_delta: float
    pastorek_a_max: float
    pastorek_zeta: float
    pastorek_d: float
    schleiss_waa_max: float
    schleiss_delta_t: float
    schleiss_tau: float
    only_cml_id: Optional[str]
    max_links: Optional[int]
    overwrite: bool
    chunk_time: int
    complevel: int
    write_1min: bool


class RainfallEstimationPipeline:
    def __init__(self, cfg: PipelineConfig):
        self.cfg = cfg
        self.logs: List[Dict] = []

    def _backend(self):
        return PycomlinkBackend(
            baseline_methods=self.cfg.baseline_methods,
            kr_models=self.cfg.kr_models,
            waa_methods=self.cfg.waa_methods,
            n_average_last_dry=self.cfg.n_average_last_dry,
            waa_gamma=self.cfg.waa_gamma,
            waa_l_antenna=self.cfg.waa_l_antenna,
            waa_delta=self.cfg.waa_delta,
            pastorek_a_max=self.cfg.pastorek_a_max,
            pastorek_zeta=self.cfg.pastorek_zeta,
            pastorek_d=self.cfg.pastorek_d,
            schleiss_waa_max=self.cfg.schleiss_waa_max,
            schleiss_delta_t=self.cfg.schleiss_delta_t,
            schleiss_tau=self.cfg.schleiss_tau,
        )

    def run(self):
        t_start_run = time.perf_counter()
        print("Loading datasets lazily...")
        ds_c = xr.open_dataset(self.cfg.cml_nc)
        ds_r = xr.open_dataset(self.cfg.radar_nc)

        try:
            if "cml_id" in ds_c.coords:
                ds_c = ds_c.assign_coords(cml_id=ds_c["cml_id"].astype(str))
            if "cml_id" in ds_r.coords:
                ds_r = ds_r.assign_coords(cml_id=ds_r["cml_id"].astype(str))

            sublink_dim = find_dim(ds_c, ["sublink_id", "sublink"])
            if sublink_dim is None:
                raise ValueError("Cannot find sublink dimension in CML dataset.")

            cml_ids_all = ds_c["cml_id"].values.astype(str)
            radar_ids = ds_r["cml_id"].values.astype(str)
            radar_pos = {rid: i for i, rid in enumerate(radar_ids)}

            full_time = pd.DatetimeIndex(ds_c["time"].values)
            radar_time = pd.DatetimeIndex(ds_r["time"].values)
            t0 = pd.Timestamp(self.cfg.tstart)
            t1 = pd.Timestamp(self.cfg.tend)
            p0 = t0 - pd.Timedelta(hours=float(self.cfg.baseline_pad_hours))

            proc0, proc1 = time_slice_indices(full_time, p0, t1)
            out0, out1 = time_slice_indices(full_time, t0, t1)
            if proc1 <= proc0 or out1 <= out0:
                raise ValueError("Empty CML time window. Check --tstart and --tend.")

            proc_time = full_time[proc0:proc1]
            out_time = full_time[out0:out1]
            out_rel0 = out0 - proc0
            out_rel1 = out1 - proc0

            rproc0, rproc1 = time_slice_indices(radar_time, p0, t1)
            rout0, rout1 = time_slice_indices(radar_time, t0, t1)
            if rproc1 <= rproc0 or rout1 <= rout0:
                raise ValueError("Empty radar time window. Check --tstart/--tend and radar file.")

            radar_proc_time = radar_time[rproc0:rproc1]
            radar_out_time = radar_time[rout0:rout1]

            sublink_ids = ds_c[sublink_dim].values if sublink_dim in ds_c.coords else np.arange(ds_c.sizes[sublink_dim])
            n_sublink = len(sublink_ids)
            has_trsl = variable_exists(ds_c, "trsl")

            selected_indices = self._select_links(ds_c, cml_ids_all, radar_pos)
            selected_ids = [cml_ids_all[i] for i in selected_indices]

            print(f"  Total CMLs             : {len(cml_ids_all)}")
            print(f"  Radar CML ids          : {len(radar_ids)}")
            print(f"  Selected CMLs          : {len(selected_indices)}")
            print(f"  Sublinks per link      : {n_sublink}")
            print(f"  Processing window      : {proc_time[0]} -> {proc_time[-1]} ({len(proc_time)} min)")
            print(f"  Output window          : {out_time[0]} -> {out_time[-1]} ({len(out_time)} min)")
            print("  Radar convention       : exact event script, radar_mmph = radar_raw * 12")
            print(f"  Baseline methods       : {', '.join(self.cfg.baseline_methods)}")
            print(f"  k-R models             : {', '.join(self.cfg.kr_models)}")
            print(f"  WAA methods            : {', '.join(self.cfg.waa_methods)}")
            print(f"  Write 1-min output     : {self.cfg.write_1min}")

            backend = self._backend()

            attrs = {
                "title": "CML rainfall estimates per sublink, baseline method, k-R model, and WAA method",
                "processing": (
                    "per-sublink TRSL + radar wet/dry + baseline_method + WAA method + k-R model; "
                    "A_obs=max(TRSL-baseline,0); link mean computed after sublink rain-rate estimation"
                ),
                "baseline_methods": ",".join(self.cfg.baseline_methods),
                "kr_models": ",".join(self.cfg.kr_models),
                "waa_methods": ",".join(self.cfg.waa_methods),
                "custom_last_dry_baseline": "closest previous finite dry TRSL; no averaging over baseline pad",
                "tstart": str(self.cfg.tstart),
                "tend": str(self.cfg.tend),
                "baseline_pad_hours": float(self.cfg.baseline_pad_hours),
                "radar_unit_convention": "exact event-script convention: radar_mmph = radar_raw * 12",
                "radar_wet_thr_5min": float(self.cfg.radar_wet_thr_5min),
                "n_average_last_dry": int(self.cfg.n_average_last_dry),
                "source_cml_file": str(self.cfg.cml_nc),
                "source_radar_file": str(self.cfg.radar_nc),
            }

            writer = StreamingNetCDFWriter(
                selected_ids,
                sublink_ids,
                self.cfg.baseline_methods,
                self.cfg.kr_models,
                self.cfg.waa_methods,
                out_time,
                radar_out_time,
                self.cfg.output_nc,
                attrs,
                write_1min=self.cfg.write_1min,
                overwrite=self.cfg.overwrite,
                complevel=self.cfg.complevel,
                chunk_time=self.cfg.chunk_time,
            )

            processed = 0
            for out_i, c_index in enumerate(tqdm(selected_indices, desc="CML links")):
                link_log = self._process_one_link(
                    ds_c, ds_r, c_index, out_i, cml_ids_all[c_index], radar_pos,
                    sublink_dim, n_sublink, has_trsl, proc0, proc1, out_rel0, out_rel1,
                    rproc0, rproc1, rout0, rout1, proc_time, radar_proc_time,
                    out_time, radar_out_time, backend, writer,
                )
                self.logs.extend(link_log)
                if any(x.get("status") == "processed" for x in link_log):
                    processed += 1

            writer.close()
            pd.DataFrame(self.logs).to_csv(self.cfg.log_csv, index=False)

            elapsed = time.perf_counter() - t_start_run
            print("\nProcessing summary")
            print(f"  Links with at least one processed sublink : {processed}")
            print(f"  Log CSV                                  : {self.cfg.log_csv}")
            print(f"  Output NetCDF                            : {self.cfg.output_nc}")
            print(f"  Elapsed time                             : {elapsed:.2f} s")
        finally:
            ds_c.close()
            ds_r.close()

    def _select_links(self, ds_c, cml_ids_all, radar_pos):
        selected = []
        for i, cid in enumerate(cml_ids_all):
            if self.cfg.only_cml_id is not None and str(cid) != str(self.cfg.only_cml_id):
                continue
            if cid not in radar_pos:
                continue
            length = safe_float(ds_c["length"].isel({"cml_id": i}).values)
            if not np.isfinite(length) or length < self.cfg.min_link_km:
                continue
            selected.append(i)
            if self.cfg.max_links is not None and len(selected) >= int(self.cfg.max_links):
                break
        if not selected:
            raise RuntimeError("No CML links selected. Check CML/radar IDs, length filter, or --only-cml-id.")
        return selected

    def _process_one_link(
        self, ds_c, ds_r, c_index, out_i, cid, radar_pos,
        sublink_dim, n_sublink, has_trsl,
        proc0, proc1, out_rel0, out_rel1,
        rproc0, rproc1, rout0, rout1,
        proc_time, radar_proc_time, out_time, radar_out_time,
        backend, writer,
    ):
        logs = []
        B = len(self.cfg.baseline_methods)
        K = len(self.cfg.kr_models)
        W = len(self.cfg.waa_methods)

        length_km = safe_float(ds_c["length"].isel({"cml_id": c_index}).values)
        rpos = radar_pos[str(cid)]

        radar_raw_proc = ds_r["rr"].isel({"cml_id": rpos, "time": slice(rproc0, rproc1)}).values
        radar_mm_h_proc = radar_values_to_mm_h(radar_raw_proc)
        wet_5 = np.isfinite(radar_mm_h_proc) & (radar_mm_h_proc > self.cfg.radar_wet_thr_5min)
        wet_1min_proc = forward_fill_radar_wet_to_cml_time(radar_proc_time, proc_time, wet_5)

        radar_raw_out = ds_r["rr"].isel({"cml_id": rpos, "time": slice(rout0, rout1)}).values
        radar_mm_h_out = radar_values_to_mm_h(radar_raw_out).astype("float32")

        if has_trsl:
            trsl_da2 = ds_c["trsl"].isel({"cml_id": c_index, "time": slice(proc0, proc1)})
        else:
            tx_da2 = ds_c["tx"].isel({"cml_id": c_index, "time": slice(proc0, proc1)})
            rx_da2 = ds_c["rx"].isel({"cml_id": c_index, "time": slice(proc0, proc1)})
            trsl_da2 = tx_da2 - rx_da2

        trsl_both = normalize_2d_time_sublink(trsl_da2, sublink_dim)

        if self.cfg.write_1min:
            R_sub = np.full((n_sublink, B, K, W, len(out_time)), np.nan, dtype="float32")
        else:
            R_sub = None

        R_sub_5 = np.full((n_sublink, B, K, W, len(radar_out_time)), np.nan, dtype="float32")
        waa_mean = np.full((n_sublink, B, W), np.nan, dtype="float32")
        A_obs_mean = np.full((n_sublink, B), np.nan, dtype="float32")
        A_rain_mean = np.full((n_sublink, B, W), np.nan, dtype="float32")

        freq_list, pol_list = [], []

        for sl in range(n_sublink):
            base_log = {
                "cml_id": str(cid),
                "cml_index": int(c_index),
                "sublink_id": int(sl),
                "length_km": length_km,
            }
            try:
                trsl = trsl_both[:, sl].astype("float64")
                valid = int(np.isfinite(trsl).sum())
                base_log["valid_samples"] = valid
                if valid < self.cfg.min_valid_samples:
                    logs.append({**base_log, "status": "skip_no_signal"})
                    continue

                f_raw = safe_float(ds_c["frequency"].isel({"cml_id": c_index, sublink_dim: sl}).values)
                f_GHz = f_raw / 1000.0 if f_raw > 1000 else f_raw
                pol = pol_to_HV(ds_c["polarization"].isel({"cml_id": c_index, sublink_dim: sl}).values)
                freq_list.append(f_GHz)
                pol_list.append(pol)

                base_log["frequency_GHz"] = f_GHz
                base_log["polarization"] = pol
                base_log["wet_fraction_1min"] = float(np.nanmean(wet_1min_proc))
                base_log["dry_samples_1min"] = int(np.sum(~wet_1min_proc))

                if np.sum(~wet_1min_proc) < self.cfg.n_average_last_dry:
                    logs.append({**base_log, "status": "skip_no_dry_baseline"})
                    continue

                est = backend.estimate(trsl, wet_1min_proc, proc_time, f_GHz, pol, length_km)

                for bi, b_method in enumerate(self.cfg.baseline_methods):
                    A_obs_mean[sl, bi] = np.nanmean(est.A_obs[b_method][out_rel0:out_rel1])

                    for wi, w_method in enumerate(self.cfg.waa_methods):
                        waa_mean[sl, bi, wi] = np.nanmean(est.waa[(b_method, w_method)][out_rel0:out_rel1])
                        A_rain_mean[sl, bi, wi] = np.nanmean(est.A_rain[(b_method, w_method)][out_rel0:out_rel1])

                        for ki, kr_model in enumerate(self.cfg.kr_models):
                            key = (b_method, kr_model, w_method)
                            R_out = est.R[key][out_rel0:out_rel1]
                            R_out = np.where(np.isfinite(R_out) & (R_out >= 0), R_out, np.nan).astype("float32")

                            if self.cfg.write_1min:
                                R_sub[sl, bi, ki, wi, :] = R_out
                            R_sub_5[sl, bi, ki, wi, :] = resample_1min_rate_to_5min_rate(
                                R_out, out_time, radar_out_time
                            )

                            q = finite_quantiles(R_out)
                            logs.append({
                                **base_log,
                                "baseline_method": b_method,
                                "kr_model": kr_model,
                                "waa_method": w_method,
                                "status": "processed",
                                "rr_finite": int(np.isfinite(R_out).sum()),
                                "rr_positive": int(np.sum(np.isfinite(R_out) & (R_out > 0))),
                                "rr_min": q["min"],
                                "rr_p50": q["p50"],
                                "rr_p95": q["p95"],
                                "rr_p99": q["p99"],
                                "rr_max": q["max"],
                                "waa_mean_db": float(waa_mean[sl, bi, wi]),
                                "A_obs_mean_db": float(A_obs_mean[sl, bi]),
                                "A_rain_mean_db": float(A_rain_mean[sl, bi, wi]),
                            })

            except Exception as exc:
                logs.append({**base_log, "status": "error", "error_message": f"{type(exc).__name__}: {exc}"})
                continue

        while len(freq_list) < n_sublink:
            sl = len(freq_list)
            f_raw = safe_float(ds_c["frequency"].isel({"cml_id": c_index, sublink_dim: sl}).values)
            freq_list.append(f_raw / 1000.0 if f_raw > 1000 else f_raw)
            pol_list.append(pol_to_HV(ds_c["polarization"].isel({"cml_id": c_index, sublink_dim: sl}).values))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            R_mean_5 = np.nanmean(R_sub_5, axis=0).astype("float32")
            R_mean = np.nanmean(R_sub, axis=0).astype("float32") if self.cfg.write_1min else None

        writer.write_metadata(out_i, length_km, freq_list, pol_list)
        writer.write_link(
            out_i,
            R_sub_5=R_sub_5,
            R_mean_5=R_mean_5,
            radar_5=radar_mm_h_out,
            waa_mean=waa_mean,
            A_obs_mean=A_obs_mean,
            A_rain_mean=A_rain_mean,
            R_sub=R_sub,
            R_mean=R_mean,
        )
        return logs




"""
CML Rainfall -> Gridded Map (Pastorek 2021 / pycml_constant / ITU_2003), STREAMED 128-GB SAFE
================================================================================================

Same scientific workflow as the original fast script:
  * select R_mean_5min for Pastorek 2021 / pycml_constant / ITU_2003
  * attach CML endpoint coordinates from a separate coordinates NetCDF when needed
  * interpolate CML rain-rate fields to a regular grid using IDW
  * optionally interpolate R_radar on the same CML-path/grid geometry
  * write gridded NetCDF and plot total rainfall maps

Main optimisation compared with the original version
----------------------------------------------------
The original script created full arrays in RAM:

    R_grids       = (all_time, y, x)
    R_radar_grids = (all_time, y, x)

That can easily kill the kernel, especially when radar is interpolated after CML.
This version processes TIME CHUNKS and writes each chunk directly to NetCDF.
Only one chunk is in memory at a time. It also chunks the IDW calculation over
GRID CELLS to avoid large temporary arrays of shape:

    batch x grid_cells x nnear

Recommended run on a 128-GB node
--------------------------------
python cml_grid_pastorek_itu03_idw_streamed_128gb.py \
  --input-nc  cml_rainrate_baseline_kr_waa_methods.nc \
  --coords-nc cml_preprocessed_v3.nc \
  --output-nc cml_rainfall_gridded_pastorek_itu03_streamed.nc \
  --time-agg per_timestep \
  --with-radar \
  --time-chunk 288 \
  --batch 32 \
  --cell-chunk 20000

If memory is still high, reduce in this order:
  --time-chunk 144
  --batch 16
  --cell-chunk 10000

Use --no-radar to write only CML.
Use --time-agg total if you only need the total map and not every 5-min grid.
The streamed modes use one clean tqdm progress bar; chunk timings are shown in the bar postfix.
"""


import argparse
import gc
import inspect
import time
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    def tqdm(x, **kwargs):
        return x


# ============================================================
# Small utilities
# ============================================================
def gb(nbytes: float) -> float:
    return float(nbytes) / 1024.0 ** 3


def print_mem_estimate(label: str, n_time: int, ny: int, nx: int) -> None:
    print(f"  [memory] {label}: output chunk ~= {gb(n_time * ny * nx * 4):.2f} GB float32")


def build_grid(lon: np.ndarray, lat: np.ndarray, resolution: float) -> Tuple[np.ndarray, np.ndarray]:
    """Grid construction matching the original arange-based meshgrid."""
    xc = np.arange(lon.min(), lon.max(), resolution)
    yc = np.arange(lat.min(), lat.max(), resolution)
    return np.meshgrid(xc, yc)


# ============================================================
# FAST IDW engine, now also chunked over grid cells
# ============================================================
def idw_fast_series(
    lon: np.ndarray,
    lat: np.ndarray,
    Z: np.ndarray,
    xg2d: np.ndarray,
    yg2d: np.ndarray,
    nnear: int,
    p: int,
    max_distance: Optional[float],
    batch: int = 32,
    cell_chunk: int = 20000,
    workers: int = -1,
) -> np.ndarray:
    """
    Interpolate Z (n_t, n_links) -> (n_t, ny, nx) float32.

    This reproduces pycomlink-like IDW math with exclude_nan=True:
      * weights are 1 / distance**p over nnear nearest valid links;
      * weights are normalized to sum to one;
      * exact grid/link coincidence uses the link value directly;
      * cells beyond max_distance are NaN.

    Memory improvement:
      The original fast version used one einsum over all grid cells, which can
      create a temporary of size batch x grid_cells x nnear. Here we process
      grid cells in chunks controlled by --cell-chunk.
    """
    from scipy.spatial import cKDTree

    Z = np.asarray(Z, dtype="float32")
    n_t, _ = Z.shape
    shape = xg2d.shape
    ny, nx = shape
    ncell = ny * nx
    q = np.c_[xg2d.ravel(), yg2d.ravel()]
    out = np.full((n_t, ny, nx), np.nan, dtype="float32")

    finite = np.isfinite(Z)
    keys = np.packbits(finite, axis=1)
    groups: dict[bytes, list[int]] = {}
    for i in range(n_t):
        groups.setdefault(keys[i].tobytes(), []).append(i)

    cell_chunk = max(1, int(cell_chunk))
    batch = max(1, int(batch))

    for _, idx_list in groups.items():
        idxs = np.asarray(idx_list, dtype=np.int64)
        valid = finite[idxs[0]]
        nv = int(valid.sum())
        if nv == 0:
            continue

        tree = cKDTree(
            np.c_[lon[valid], lat[valid]],
            balanced_tree=False,
            compact_nodes=False,
        )
        k = min(int(nnear), nv)
        d, ix = tree.query(q, k=k, workers=workers)
        if k == 1:
            d = d[:, None]
            ix = ix[:, None]

        with np.errstate(divide="ignore", invalid="ignore"):
            w = 1.0 / d ** p
            wn = w / w.sum(axis=1, keepdims=True)

        near0 = d[:, 0] < 1e-10
        far = (d[:, 0] > max_distance) if max_distance is not None else None
        ix0 = ix[:, 0]

        for b0 in range(0, len(idxs), batch):
            bb = idxs[b0:b0 + batch]
            zb = Z[bb][:, valid]  # (B, nv)

            vals = np.full((len(bb), ncell), np.nan, dtype="float32")
            for c0 in range(0, ncell, cell_chunk):
                c1 = min(c0 + cell_chunk, ncell)
                ixc = ix[c0:c1]
                wnc = wn[c0:c1].astype("float32", copy=False)

                # Temporary controlled by: B x cell_chunk x k
                vals[:, c0:c1] = np.einsum(
                    "bck,ck->bc",
                    zb[:, ixc],
                    wnc,
                    optimize=True,
                ).astype("float32", copy=False)

            if near0.any():
                vals[:, near0] = zb[:, ix0[near0]]
            if far is not None:
                vals[:, far] = np.nan

            out[bb] = vals.reshape((len(bb), ny, nx))

        # Free large per-mask arrays before next mask group
        del tree, d, ix, w, wn, near0, far, ix0
        gc.collect()

    return out


# ============================================================
# Optional pycomlink engine, kept for self-check/fallback
# ============================================================
def make_idw(nnear, p, exclude_nan, max_distance):
    import pycomlink as pycml
    return pycml.spatial.interpolator.IdwKdtreeInterpolator(
        nnear=nnear, p=p, exclude_nan=exclude_nan, max_distance=max_distance
    )


def _supports_grid_kwargs(idw) -> bool:
    try:
        params = inspect.signature(idw.__call__).parameters
        return "xgrid" in params and "ygrid" in params
    except (ValueError, TypeError):
        return False


def idw_pycomlink_series(
    Z,
    lon,
    lat,
    xg2d,
    yg2d,
    resolution,
    nnear,
    p,
    max_distance,
    skip_all_nan=True,
    desc="pycomlink",
) -> np.ndarray:
    idw = make_idw(nnear, p, True, max_distance)
    supports = _supports_grid_kwargs(idw)
    shape = xg2d.shape
    out = np.full((Z.shape[0],) + shape, np.nan, dtype="float32")
    for i in tqdm(range(Z.shape[0]), desc=desc):
        z = np.asarray(Z[i], dtype="float64")
        if skip_all_nan and not np.isfinite(z).any():
            continue
        if supports:
            g = idw(x=lon, y=lat, z=z, xgrid=xg2d, ygrid=yg2d)
        else:
            g = idw(x=lon, y=lat, z=z, resolution=resolution)
        out[i] = np.asarray(g, dtype="float32")
    return out


def interpolate_all(Z, lon, lat, xg2d, yg2d, args, engine, desc="Interpolating"):
    if engine == "fast":
        return idw_fast_series(
            lon,
            lat,
            Z,
            xg2d,
            yg2d,
            args.nnear,
            args.p,
            args.max_dist,
            batch=args.batch,
            cell_chunk=args.cell_chunk,
            workers=args.workers,
        )
    return idw_pycomlink_series(
        Z,
        lon,
        lat,
        xg2d,
        yg2d,
        args.resolution,
        args.nnear,
        args.p,
        args.max_dist,
        desc=desc,
    )


# ============================================================
# Coordinates
# ============================================================
def load_link_coordinates(ds: xr.Dataset, args) -> xr.Dataset:
    needed = [args.site_a_lat, args.site_a_lon, args.site_b_lat, args.site_b_lon]
    if all(v in ds.variables for v in needed):
        print("  Coordinates: found in the estimation file.")
        return ds
    if not args.coords_nc:
        raise ValueError("Estimation file has no site coordinates; pass --coords-nc.")

    print(f"  Coordinates: reading from {args.coords_nc}")
    cds = xr.open_dataset(args.coords_nc)
    if "cml_id" in cds.coords:
        cds = cds.assign_coords(cml_id=cds["cml_id"].astype(str))
    missing = [v for v in needed if v not in cds.variables]
    if missing:
        raise KeyError(f"Coordinate file is missing variables: {missing}")

    sub = cds[needed].reindex(cml_id=ds["cml_id"].values)
    n_missing = int(np.isnan(sub[args.site_a_lat].values).sum())
    if n_missing:
        print(f"  WARNING: {n_missing} links have no coordinates and will be dropped.")
    for v in needed:
        ds[v] = ("cml_id", sub[v].values)
    cds.close()
    return ds


# ============================================================
# Streaming NetCDF writer
# ============================================================
def _datetime64_to_seconds_since_epoch(times: np.ndarray) -> np.ndarray:
    return pd.DatetimeIndex(times).astype("int64") // 10**9


class StreamingGridWriter:
    """Small netCDF4 writer for chunk-by-chunk output."""

    def __init__(self, output_nc: str, times, xgrid, ygrid, args, write_radar: bool):
        try:
            from netCDF4 import Dataset
        except Exception as e:
            raise ImportError(
                "This streamed writer needs netCDF4. Install it in your environment with:\n"
                "  pip install netCDF4\n"
                "or\n"
                "  conda install -c conda-forge netcdf4"
            ) from e

        self.path = str(output_nc)
        self.root = Dataset(self.path, "w", format="NETCDF4")
        self.root.createDimension("time", len(times))
        self.root.createDimension("y", len(ygrid))
        self.root.createDimension("x", len(xgrid))

        tv = self.root.createVariable("time", "i8", ("time",))
        tv[:] = _datetime64_to_seconds_since_epoch(times)
        tv.units = "seconds since 1970-01-01 00:00:00"
        tv.calendar = "proleptic_gregorian"
        tv.long_name = "time"

        xv = self.root.createVariable("x", "f8", ("x",))
        yv = self.root.createVariable("y", "f8", ("y",))
        xv[:] = xgrid
        yv[:] = ygrid
        xv.units = "degrees_east"
        yv.units = "degrees_north"
        xv.long_name = "longitude"
        yv.long_name = "latitude"

        self.root.title = "Gridded CML rainfall (Pastorek 2021, pycml_constant, ITU_2003)"
        self.root.interpolation = f"IDW (nnear={args.nnear}, p={args.p}, max_dist={args.max_dist})"
        self.root.engine = args.engine
        self.root.grid_resolution_deg = float(args.resolution)
        self.root.baseline_method = args.baseline_method
        self.root.kr_model = args.kr_model
        self.root.waa_method = args.waa_method
        self.root.source = args.input_nc
        self.root.note = "Written incrementally in time chunks to avoid holding all grids in RAM."

        ny, nx = len(ygrid), len(xgrid)
        tchunk = max(1, min(int(args.time_chunk), len(times)))
        ychunk = min(ny, int(args.y_chunk)) if args.y_chunk else ny
        xchunk = min(nx, int(args.x_chunk)) if args.x_chunk else nx
        chunksizes = (tchunk, ychunk, xchunk)

        zlib = not args.no_compression
        self.R = self.root.createVariable(
            "R",
            "f4",
            ("time", "y", "x"),
            zlib=zlib,
            complevel=args.complevel,
            shuffle=True,
            chunksizes=chunksizes,
            fill_value=np.float32(np.nan),
        )
        self.R.units = "mm h-1"
        self.R.long_name = "CML rain rate (Pastorek/ITU_2003, IDW)"

        self.R_radar = None
        if write_radar:
            self.R_radar = self.root.createVariable(
                "R_radar",
                "f4",
                ("time", "y", "x"),
                zlib=zlib,
                complevel=args.complevel,
                shuffle=True,
                chunksizes=chunksizes,
                fill_value=np.float32(np.nan),
            )
            self.R_radar.units = "mm h-1"
            self.R_radar.long_name = "Radar rain rate interpolated with the same IDW geometry"

    def write(self, name: str, start: int, data: np.ndarray) -> None:
        end = start + data.shape[0]
        if name == "R":
            self.R[start:end, :, :] = data.astype("float32", copy=False)
        elif name == "R_radar":
            if self.R_radar is None:
                raise ValueError("R_radar variable was not created.")
            self.R_radar[start:end, :, :] = data.astype("float32", copy=False)
        else:
            raise KeyError(name)

    def close(self):
        self.root.close()


# ============================================================
# Plotting
# ============================================================
def plot_cml_lines(ds, ax, args, color="k", linewidth=0.5, alpha=0.4):
    ax.plot(
        [ds[args.site_a_lon].values, ds[args.site_b_lon].values],
        [ds[args.site_a_lat].values, ds[args.site_b_lat].values],
        color=color,
        linewidth=linewidth,
        alpha=alpha,
    )


def plot_total_map(ds, xgrid, ygrid, R_total, R_radar_total, args, out_png):
    bounds = np.array([0.1, 1, 2, 5, 10, 15, 20, 30, 40, 50])
    norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend="both")
    two = R_radar_total is not None
    fig, axes = plt.subplots(1, 2 if two else 1, figsize=(14 if two else 8, 6))
    ax1, ax2 = (axes[0], axes[1]) if two else (axes, None)

    pc1 = ax1.pcolormesh(xgrid, ygrid, R_total, shading="nearest", cmap="turbo", norm=norm)
    plot_cml_lines(ds, ax1, args)
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    ax1.set_aspect("equal")
    ax1.set_title("CML rainfall total (Pastorek 2021, pycml_constant, ITU_2003)")
    fig.colorbar(pc1, ax=ax1, label="Rainfall (mm)", shrink=0.8)

    if two:
        pc2 = ax2.pcolormesh(xgrid, ygrid, R_radar_total, shading="nearest", cmap="turbo", norm=norm)
        plot_cml_lines(ds, ax2, args)
        ax2.set_xlabel("Longitude")
        ax2.set_ylabel("Latitude")
        ax2.set_aspect("equal")
        ax2.set_title("Radar rainfall total (via CML paths)")
        fig.colorbar(pc2, ax=ax2, label="Rainfall (mm)", shrink=0.8)

    plt.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved -> {out_png}")


def _grid_attrs(args, engine) -> dict:
    return {
        "title": "Gridded CML rainfall (Pastorek 2021, pycml_constant, ITU_2003)",
        "interpolation": f"IDW (nnear={args.nnear}, p={args.p}, max_dist={args.max_dist})",
        "engine": engine,
        "grid_resolution_deg": args.resolution,
        "baseline_method": args.baseline_method,
        "kr_model": args.kr_model,
        "waa_method": args.waa_method,
        "source": args.input_nc,
    }


# ============================================================
# Total-only mode
# ============================================================
def run_total_mode(R, ds, good_idx, lon, lat, xg2d, yg2d, xgrid, ygrid, args):
    """Original total-only behaviour, but avoids all-time gridded arrays."""
    print("\nInterpolating total rainfall sum ...")
    Rsum = (R.isel(cml_id=good_idx).sum(dim="time_5min").values.astype("float32")) * (5.0 / 60.0)
    Rg = interpolate_all(Rsum[None, :], lon, lat, xg2d, yg2d, args, args.engine, "total")[0]

    data_vars = {
        "R_total": (
            ["y", "x"],
            Rg.astype("float32"),
            {"units": "mm", "long_name": "Total rainfall (CML, Pastorek/ITU_2003)"},
        )
    }
    ds_grid = xr.Dataset(data_vars, coords={"x": xgrid, "y": ygrid}, attrs=_grid_attrs(args, args.engine))

    enc = {
        "R_total": {
            "zlib": not args.no_compression,
            "complevel": args.complevel,
            "dtype": "float32",
        }
    }
    print("\nSaving total-only NetCDF ...")
    ds_grid.to_netcdf(args.output_nc, encoding=enc)
    print(f"  Saved -> {args.output_nc}")
    return Rg, None


# ============================================================
# Chunked per-timestep mode
# ============================================================
def get_chunk_array(da: xr.DataArray, good_idx: np.ndarray, start: int, end: int) -> np.ndarray:
    """Load one time chunk as contiguous float32 with shape (time, cml_id)."""
    z = da.isel(time_5min=slice(start, end), cml_id=good_idx).transpose("time_5min", "cml_id").values
    return np.ascontiguousarray(z, dtype="float32")


def run_per_timestep_streamed(R, ds, good_idx, lon, lat, xg2d, yg2d, xgrid, ygrid, args):
    times = ds["time_5min"].values
    ntime = len(times)
    ny, nx = xg2d.shape
    write_radar = bool(args.with_radar and "R_radar" in ds.variables)

    if args.with_radar and "R_radar" not in ds.variables:
        print("  WARNING: --with-radar requested but R_radar not found. Writing only CML.")

    print("\nStreaming per-timestep interpolation to NetCDF ...")
    print(f"  Steps: {ntime}")
    print(f"  Time chunk: {args.time_chunk}")
    print(f"  IDW batch: {args.batch}")
    print(f"  IDW cell chunk: {args.cell_chunk}")
    print_mem_estimate("one output variable per time chunk", min(args.time_chunk, ntime), ny, nx)

    writer = StreamingGridWriter(args.output_nc, times, xgrid, ygrid, args, write_radar=write_radar)

    R_total = np.zeros((ny, nx), dtype="float64")
    R_radar_total = np.zeros((ny, nx), dtype="float64") if write_radar else None
    factor = 5.0 / 60.0

    try:
        # One clean progress bar only. Per-chunk timing is shown in the bar postfix.
        chunk_iter = range(0, ntime, args.time_chunk)
        pbar = tqdm(
            chunk_iter,
            desc="gridding",
            unit="chunk",
            dynamic_ncols=True,
            leave=True,
        )

        for start in pbar:
            end = min(start + args.time_chunk, ntime)

            # CML chunk
            ti = time.perf_counter()
            Z = get_chunk_array(R, good_idx, start, end)
            G = interpolate_all(Z, lon, lat, xg2d, yg2d, args, args.engine, "CML chunk")
            writer.write("R", start, G)
            R_total += np.nansum(G, axis=0) * factor
            cml_s = time.perf_counter() - ti
            del Z, G
            gc.collect()

            radar_s = None
            # Radar chunk, processed AFTER freeing CML grid chunk.
            if write_radar:
                ti = time.perf_counter()
                Zr = get_chunk_array(ds["R_radar"], good_idx, start, end)
                Gr = interpolate_all(Zr, lon, lat, xg2d, yg2d, args, args.engine, "radar chunk")
                writer.write("R_radar", start, Gr)
                R_radar_total += np.nansum(Gr, axis=0) * factor
                radar_s = time.perf_counter() - ti
                del Zr, Gr
                gc.collect()

            if radar_s is None:
                pbar.set_postfix_str(f"{start}:{end}/{ntime} | CML {cml_s:.1f}s")
            else:
                pbar.set_postfix_str(
                    f"{start}:{end}/{ntime} | CML {cml_s:.1f}s | Radar {radar_s:.1f}s"
                )

        writer.close()
        print(f"\n  Saved -> {args.output_nc}")

    except Exception:
        writer.close()
        raise

    return R_total.astype("float32"), None if R_radar_total is None else R_radar_total.astype("float32")


# ============================================================
# Optional hourly mode
# ============================================================
def run_hourly_streamed(R, ds, good_idx, lon, lat, xg2d, yg2d, xgrid, ygrid, args):
    """
    Hourly mode kept available. It first creates hourly link fields with xarray,
    then streams the hourly grids. Radar hourly is not implemented here because
    the original script did not create hourly radar output.
    """
    print("\nAggregating CML to hourly link fields ...")
    R_h = R.isel(cml_id=good_idx).resample(time_5min="1h").sum() / 12.0
    times = R_h["time_5min"].values
    ntime = len(times)
    ny, nx = xg2d.shape

    writer = StreamingGridWriter(args.output_nc, times, xgrid, ygrid, args, write_radar=False)
    R_total = np.zeros((ny, nx), dtype="float64")

    try:
        pbar = tqdm(
            range(0, ntime, args.time_chunk),
            desc="hourly gridding",
            unit="chunk",
            dynamic_ncols=True,
            leave=True,
        )
        for start in pbar:
            end = min(start + args.time_chunk, ntime)
            ti = time.perf_counter()
            Z = np.ascontiguousarray(
                R_h.isel(time_5min=slice(start, end)).transpose("time_5min", "cml_id").values,
                dtype="float32",
            )
            G = interpolate_all(Z, lon, lat, xg2d, yg2d, args, args.engine, "hourly CML chunk")
            writer.write("R", start, G)
            R_total += np.nansum(G, axis=0)
            chunk_s = time.perf_counter() - ti
            del Z, G
            gc.collect()
            pbar.set_postfix_str(f"{start}:{end}/{ntime} | CML {chunk_s:.1f}s")
        writer.close()
        print(f"\n  Saved -> {args.output_nc}")
    except Exception:
        writer.close()
        raise

    return R_total.astype("float32"), None


# ============================================================
# Self-check
# ============================================================
def run_selfcheck(R, ds, good_idx, lon, lat, xg2d, yg2d, args):
    n = int(args.selfcheck)
    print(f"\n[selfcheck] comparing fast vs pycomlink on {n} timesteps ...")
    try:
        Zc = np.ascontiguousarray(
            R.isel(cml_id=good_idx, time_5min=slice(0, n)).transpose("time_5min", "cml_id").values,
            dtype="float32",
        )
    except Exception as e:
        print(f"  [selfcheck] could not build subset: {e}")
        return
    try:
        Rf = idw_fast_series(
            lon,
            lat,
            Zc,
            xg2d,
            yg2d,
            args.nnear,
            args.p,
            args.max_dist,
            batch=args.batch,
            cell_chunk=args.cell_chunk,
            workers=args.workers,
        )
        Rp = idw_pycomlink_series(
            Zc,
            lon,
            lat,
            xg2d,
            yg2d,
            args.resolution,
            args.nnear,
            args.p,
            args.max_dist,
            desc="selfcheck",
        )
    except Exception as e:
        print(f"  [selfcheck] pycomlink unavailable or failed ({e}); skipping.")
        return

    a = np.nan_to_num(Rf, nan=-9e9)
    b = np.nan_to_num(Rp, nan=-9e9)
    nan_ok = np.array_equal(np.isnan(Rf), np.isnan(Rp))
    maxdiff = float(np.max(np.abs(a - b)))
    print(f"  [selfcheck] max|diff|={maxdiff:.3e}  nan-pattern-identical={nan_ok}")
    if nan_ok and maxdiff < 1e-3:
        print("  [selfcheck] PASS: engines agree; the fast engine is safe to use.")
    else:
        print("  [selfcheck] WARNING: engines differ; inspect pycomlink version/settings.")






# ============================================================
# MERGED ESTIMATION + GRIDDED OUTPUT + GAUGE IDW
# ============================================================
# The blocks above are copied from the two tested scripts.  The code below only
# orchestrates them: first the original estimation pipeline writes the 5-min
# CML/Radar estimation file, then CML/Radar/Gauge fields are streamed to a grid.
# Gridding can be parallelized over independent time chunks while a single main
# process writes the NetCDF safely.

from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED


# ----------------------------
# Generic coordinate helpers
# ----------------------------
def _pick_name(ds: xr.Dataset, candidates, required=True, label="variable"):
    for name in candidates:
        if name and name in ds.variables:
            return name
    if required:
        raise KeyError(f"Could not find {label}; tried: {candidates}. Available: {list(ds.variables)}")
    return None


def _find_station_dim_for_var(ds: xr.Dataset, var_name: str, preferred=None):
    da = ds[var_name]
    if preferred is not None:
        if preferred not in da.dims:
            raise ValueError(f"Gauge station dim {preferred!r} is not in {var_name}.dims={da.dims}")
        return preferred
    for d in ["WSI", "station", "station_id", "gauge", "site", "id"]:
        if d in da.dims and d != "time":
            return d
    other = [d for d in da.dims if d != "time"]
    if not other:
        raise ValueError(f"Could not infer station dimension for {var_name}.dims={da.dims}")
    return other[0]


def _to_seconds_since_epoch(times: np.ndarray) -> np.ndarray:
    return (pd.DatetimeIndex(times).view("int64") // 10**9).astype("int64")


def _median_step_hours(times: np.ndarray, fallback_hours: float) -> float:
    idx = pd.DatetimeIndex(times)
    if len(idx) < 2:
        return float(fallback_hours)
    d = np.diff(idx.view("int64")) / 1e9 / 3600.0
    d = d[np.isfinite(d) & (d > 0)]
    if d.size == 0:
        return float(fallback_hours)
    return float(np.nanmedian(d))


# ----------------------------
# Merged gridded writer
# ----------------------------
class MergedGridWriter:
    def __init__(self, output_nc: str, times, xgrid, ygrid, args, write_radar: bool, gauge_times=None):
        try:
            from netCDF4 import Dataset
        except Exception as e:
            raise ImportError(
                "This streamed writer needs netCDF4. Install it with:\n"
                "  pip install netCDF4\n"
                "or:\n"
                "  conda install -c conda-forge netcdf4"
            ) from e

        self.path = str(output_nc)
        out = Path(self.path)
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists():
            if args.overwrite_grid:
                out.unlink()
            else:
                raise FileExistsError(f"Gridded output exists: {out}. Use --overwrite-grid.")

        self.root = Dataset(self.path, "w", format="NETCDF4")
        self.root.createDimension("time", len(times))
        self.root.createDimension("y", len(ygrid))
        self.root.createDimension("x", len(xgrid))
        self.has_gauge = gauge_times is not None and len(gauge_times) > 0
        if self.has_gauge:
            self.root.createDimension("time_gauge", len(gauge_times))

        tv = self.root.createVariable("time", "i8", ("time",))
        tv[:] = _to_seconds_since_epoch(times)
        tv.units = "seconds since 1970-01-01 00:00:00"
        tv.calendar = "proleptic_gregorian"
        tv.long_name = "CML/radar 5-min time"

        if self.has_gauge:
            tg = self.root.createVariable("time_gauge", "i8", ("time_gauge",))
            tg[:] = _to_seconds_since_epoch(gauge_times)
            tg.units = "seconds since 1970-01-01 00:00:00"
            tg.calendar = "proleptic_gregorian"
            tg.long_name = "rain-gauge native/resampled time"

        xv = self.root.createVariable("x", "f8", ("x",))
        yv = self.root.createVariable("y", "f8", ("y",))
        xv[:] = xgrid
        yv[:] = ygrid
        xv.units = "degrees_east"
        yv.units = "degrees_north"
        xv.long_name = "longitude"
        yv.long_name = "latitude"

        self.root.title = "Merged rainfall estimation and gridding: CML, radar and rain gauges"
        self.root.interpolation = f"IDW (nnear={args.nnear}, p={args.p}, max_dist={args.max_dist})"
        self.root.engine = "fast chunked IDW"
        self.root.grid_resolution_deg = float(args.resolution)
        self.root.baseline_method = args.grid_baseline_method
        self.root.kr_model = args.grid_kr_model
        self.root.waa_method = args.grid_waa_method
        self.root.source_estimation_file = str(args.estimation_nc)
        self.root.source_cml_file = str(args.cml_nc)
        self.root.source_radar_file = str(args.radar_nc)
        if args.with_gauge:
            self.root.source_gauge_file = str(args.gauge_nc)
        self.root.note = (
            "CML and radar are on dimension time. Gauge interpolation is on time_gauge "
            "to avoid changing the native gauge temporal sampling."
        )

        ny, nx = len(ygrid), len(xgrid)
        tchunk = max(1, min(int(args.grid_time_chunk), len(times)))
        ychunk = min(ny, int(args.y_chunk)) if args.y_chunk else ny
        xchunk = min(nx, int(args.x_chunk)) if args.x_chunk else nx
        chunksizes = (tchunk, ychunk, xchunk)
        zlib = not args.no_compression

        self.R = self.root.createVariable(
            "R", "f4", ("time", "y", "x"),
            zlib=zlib, complevel=args.complevel, shuffle=True,
            chunksizes=chunksizes, fill_value=np.float32(np.nan),
        )
        self.R.units = "mm h-1"
        self.R.long_name = "CML rain rate gridded from R_mean_5min"

        self.R_radar = None
        if write_radar:
            self.R_radar = self.root.createVariable(
                "R_radar", "f4", ("time", "y", "x"),
                zlib=zlib, complevel=args.complevel, shuffle=True,
                chunksizes=chunksizes, fill_value=np.float32(np.nan),
            )
            self.R_radar.units = "mm h-1"
            self.R_radar.long_name = "Radar rain rate interpolated to the same grid"

        self.R_gauge = None
        if self.has_gauge:
            gtchunk = max(1, min(int(args.gauge_time_chunk), len(gauge_times)))
            self.R_gauge = self.root.createVariable(
                "R_gauge", "f4", ("time_gauge", "y", "x"),
                zlib=zlib, complevel=args.complevel, shuffle=True,
                chunksizes=(gtchunk, ychunk, xchunk), fill_value=np.float32(np.nan),
            )
            self.R_gauge.units = "mm h-1"
            self.R_gauge.long_name = "Rain-gauge rain rate interpolated to the same grid"
            self.R_gauge.gauge_amount_conversion = (
                "raw gauge values multiplied by 6 to convert 10-min amount to mm h-1"
                if args.gauge_is_10min_amount else
                "raw gauge values assumed to already be mm h-1"
            )

    def write_time_chunk(self, start: int, R_grid: np.ndarray, R_radar_grid: Optional[np.ndarray] = None):
        end = start + R_grid.shape[0]
        self.R[start:end, :, :] = R_grid.astype("float32", copy=False)
        if R_radar_grid is not None:
            if self.R_radar is None:
                raise ValueError("R_radar variable was not created.")
            self.R_radar[start:end, :, :] = R_radar_grid.astype("float32", copy=False)

    def write_gauge_chunk(self, start: int, R_gauge_grid: np.ndarray):
        if self.R_gauge is None:
            raise ValueError("R_gauge variable was not created.")
        end = start + R_gauge_grid.shape[0]
        self.R_gauge[start:end, :, :] = R_gauge_grid.astype("float32", copy=False)

    def close(self):
        self.root.close()


def prepare_gauge_metadata(args, xgrid, ygrid):
    if not args.with_gauge:
        return None
    if not args.gauge_nc:
        print("  WARNING: --with-gauge requested but --gauge-nc is empty. Skipping gauges.")
        return None
    if not Path(args.gauge_nc).exists():
        print(f"  WARNING: Gauge file not found: {args.gauge_nc}. Skipping gauges.")
        return None

    gd = xr.open_dataset(args.gauge_nc)
    try:
        if args.gauge_var not in gd.variables:
            raise KeyError(f"Gauge variable {args.gauge_var!r} not found in {args.gauge_nc}")
        time_dim = args.gauge_time_dim or _pick_name(gd, ["time", "t"], label="gauge time coordinate")
        station_dim = _find_station_dim_for_var(gd, args.gauge_var, args.gauge_station_dim)
        lon_name = args.gauge_lon or _pick_name(gd, ["lon", "longitude", "station_lon", "station_longitude", "x"], label="gauge longitude")
        lat_name = args.gauge_lat or _pick_name(gd, ["lat", "latitude", "station_lat", "station_latitude", "y"], label="gauge latitude")

        gauge_time_all = pd.DatetimeIndex(gd[time_dim].values)
        # Use the gridding window for gauges. If not provided, fall back to
        # the estimation window so one --tstart/--tend controls the full run.
        tmin = pd.Timestamp(args.grid_tstart or args.tstart)
        tmax = pd.Timestamp(args.grid_tend or args.tend)
        t0 = int(np.searchsorted(gauge_time_all.values, np.datetime64(tmin), "left"))
        t1 = int(np.searchsorted(gauge_time_all.values, np.datetime64(tmax), "right"))
        gauge_times = gauge_time_all[t0:t1]
        if len(gauge_times) == 0:
            print("  WARNING: No gauge time steps overlap the requested period. Skipping gauges.")
            return None

        glon = np.asarray(gd[lon_name].values, dtype="float64")
        glat = np.asarray(gd[lat_name].values, dtype="float64")
        if glon.ndim != 1 or glat.ndim != 1:
            raise ValueError("Gauge lon/lat must be 1-D station coordinates for this gridding script.")
        if glon.size != gd.sizes[station_dim] or glat.size != gd.sizes[station_dim]:
            raise ValueError(
                f"Gauge lon/lat length mismatch. len(lon)={glon.size}, len(lat)={glat.size}, "
                f"station_dim={station_dim} size={gd.sizes[station_dim]}"
            )

        pad = float(args.gauge_bbox_pad if args.gauge_bbox_pad is not None else args.max_dist)
        keep = (
            np.isfinite(glon) & np.isfinite(glat) &
            (glon >= float(np.nanmin(xgrid)) - pad) &
            (glon <= float(np.nanmax(xgrid)) + pad) &
            (glat >= float(np.nanmin(ygrid)) - pad) &
            (glat <= float(np.nanmax(ygrid)) + pad)
        )
        keep_idx = np.flatnonzero(keep)
        if keep_idx.size == 0:
            print("  WARNING: No gauges inside the CML grid bounding box. Skipping gauges.")
            return None

        print(f"  Gauges kept for IDW: {keep_idx.size}/{glon.size}")
        return {
            "gauge_nc": str(args.gauge_nc),
            "var": args.gauge_var,
            "time_dim": time_dim,
            "station_dim": station_dim,
            "time_start_index": t0,
            "time_end_index": t1,
            "times": gauge_times.values,
            "keep_idx": keep_idx.astype("int64"),
            "lon": glon[keep_idx].astype("float64"),
            "lat": glat[keep_idx].astype("float64"),
            "is_10min_amount": bool(args.gauge_is_10min_amount),
        }
    finally:
        gd.close()


def _load_gauge_chunk(meta: dict, rel_start: int, rel_end: int) -> np.ndarray:
    gd = xr.open_dataset(meta["gauge_nc"])
    try:
        abs0 = int(meta["time_start_index"] + rel_start)
        abs1 = int(meta["time_start_index"] + rel_end)
        da = gd[meta["var"]].isel({meta["time_dim"]: slice(abs0, abs1), meta["station_dim"]: meta["keep_idx"]})
        z = da.transpose(meta["time_dim"], meta["station_dim"]).values
        z = np.ascontiguousarray(z, dtype="float32")
        z[~np.isfinite(z)] = np.nan
        z[z < 0] = np.nan
        if meta["is_10min_amount"]:
            z *= np.float32(6.0)
        return z
    finally:
        gd.close()


def _grid_merged_job(job: dict) -> dict:
    kind = job["kind"]
    t0 = time.perf_counter()
    if kind == "cml_radar":
        ds = xr.open_dataset(job["estimation_nc"])
        try:
            if "cml_id" in ds.coords:
                ds = ds.assign_coords(cml_id=ds["cml_id"].astype(str))

            R = ds[job["rain_var"]]
            sel = {}
            for dim, val in [("baseline_method", job["baseline_method"]), ("kr_model", job["kr_model"]), ("waa_method", job["waa_method"])]:
                if dim in R.dims:
                    sel[dim] = val
            if sel:
                R = R.sel(sel)

            # start/end are output-relative indices used by the writer.
            # file_start/file_end are absolute indices inside the estimation file.
            start, end = int(job["start"]), int(job["end"])
            file_start = int(job.get("file_start", start))
            file_end = int(job.get("file_end", end))
            good_idx = job["good_idx"]
            Z = R.isel(cml_id=good_idx, time_5min=slice(file_start, file_end)).transpose("time_5min", "cml_id").values
            Z = np.ascontiguousarray(Z, dtype="float32")
            G = idw_fast_series(job["lon"], job["lat"], Z, job["xg2d"], job["yg2d"], job["nnear"], job["p"], job["max_dist"], batch=job["batch"], cell_chunk=job["cell_chunk"], workers=job["idw_query_workers"])
            total = np.nansum(G, axis=0).astype("float64") * float(job["factor_hours"])
            del Z

            Gr = None
            radar_total = None
            if job["write_radar"] and "R_radar" in ds.variables:
                Zr = ds["R_radar"].isel(cml_id=good_idx, time_5min=slice(file_start, file_end)).transpose("time_5min", "cml_id").values
                Zr = np.ascontiguousarray(Zr, dtype="float32")
                Gr = idw_fast_series(job["lon"], job["lat"], Zr, job["xg2d"], job["yg2d"], job["nnear"], job["p"], job["max_dist"], batch=job["batch"], cell_chunk=job["cell_chunk"], workers=job["idw_query_workers"])
                radar_total = np.nansum(Gr, axis=0).astype("float64") * float(job["factor_hours"])
                del Zr

            return {"kind": kind, "start": start, "end": end, "file_start": file_start, "file_end": file_end, "R": G, "R_radar": Gr, "R_total": total, "R_radar_total": radar_total, "elapsed": time.perf_counter() - t0}
        finally:
            ds.close()

    if kind == "gauge":
        start, end = int(job["start"]), int(job["end"])
        Zg = _load_gauge_chunk(job["gauge_meta"], start, end)
        Gg = idw_fast_series(job["gauge_meta"]["lon"], job["gauge_meta"]["lat"], Zg, job["xg2d"], job["yg2d"], job["nnear"], job["p"], job["max_dist"], batch=job["batch"], cell_chunk=job["cell_chunk"], workers=job["idw_query_workers"])
        total = np.nansum(Gg, axis=0).astype("float64") * float(job["factor_hours"])
        return {"kind": kind, "start": start, "end": end, "R_gauge": Gg, "R_gauge_total": total, "elapsed": time.perf_counter() - t0}

    raise ValueError(f"Unknown grid job kind: {kind}")


def plot_total_map_merged(ds_links, xgrid, ygrid, R_total, R_radar_total, R_gauge_total, args, out_png):
    available = [("CML total", R_total)]
    if R_radar_total is not None:
        available.append(("Radar total", R_radar_total))
    if R_gauge_total is not None:
        available.append(("Rain-gauge total", R_gauge_total))

    bounds = np.array([0.1, 1, 2, 5, 10, 15, 20, 30, 40, 50], dtype=float)
    norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend="both")
    n = len(available)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 6), squeeze=False)
    axes = axes.ravel()
    for ax, (title, field) in zip(axes, available):
        pc = ax.pcolormesh(xgrid, ygrid, field, shading="nearest", cmap="turbo", norm=norm)
        try:
            plot_cml_lines(ds_links, ax, args)
        except Exception:
            pass
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_aspect("equal")
        ax.set_title(title)
        fig.colorbar(pc, ax=ax, label="Rainfall (mm)", shrink=0.8)
    plt.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved -> {out_png}")


def run_estimation_stage(args):
    if args.skip_estimation:
        print("\n[stage 1/2] Skipping estimation; using existing file:")
        print(f"  {args.estimation_nc}")
        if not Path(args.estimation_nc).exists():
            raise FileNotFoundError(f"--skip-estimation was used but file does not exist: {args.estimation_nc}")
        return

    print("\n[stage 1/2] Running original rainfall-estimation pipeline ...")
    cfg = PipelineConfig(
        cml_nc=args.cml_nc, radar_nc=args.radar_nc, output_nc=args.estimation_nc, log_csv=args.log_csv,
        tstart=args.tstart, tend=args.tend, radar_wet_thr_5min=args.radar_wet_thr_5min,
        baseline_pad_hours=args.baseline_pad_hours, min_link_km=args.min_link_km,
        min_valid_samples=args.min_valid_samples, n_average_last_dry=args.n_average_last_dry,
        baseline_methods=tuple(dict.fromkeys(args.baseline_methods)), kr_models=tuple(dict.fromkeys(args.kr_models)),
        waa_methods=tuple(dict.fromkeys(args.waa_methods)), waa_gamma=args.waa_gamma,
        waa_l_antenna=args.waa_l_antenna, waa_delta=args.waa_delta,
        pastorek_a_max=args.pastorek_a_max, pastorek_zeta=args.pastorek_zeta, pastorek_d=args.pastorek_d,
        schleiss_waa_max=args.schleiss_waa_max, schleiss_delta_t=args.schleiss_delta_t, schleiss_tau=args.schleiss_tau,
        only_cml_id=args.only_cml_id, max_links=args.max_links, overwrite=args.overwrite_estimation,
        chunk_time=args.est_chunk_time, complevel=args.complevel, write_1min=args.write_1min,
    )
    RainfallEstimationPipeline(cfg).run()


def _submit_next(executor, pending: set, job_iter):
    try:
        job = next(job_iter)
    except StopIteration:
        return False
    pending.add(executor.submit(_grid_merged_job, job))
    return True


def _resolve_time_window_indices(times, start_value, end_value, label="time"):
    """Return [i0, i1) indices for an inclusive user time window."""
    idx = pd.DatetimeIndex(times)
    if len(idx) == 0:
        raise ValueError(f"Empty {label} coordinate.")

    t0 = pd.Timestamp(start_value) if start_value is not None else idx[0]
    t1 = pd.Timestamp(end_value) if end_value is not None else idx[-1]
    if t1 < t0:
        raise ValueError(f"Invalid {label} window: end {t1} is before start {t0}.")

    i0 = int(np.searchsorted(idx.values, np.datetime64(t0), "left"))
    i1 = int(np.searchsorted(idx.values, np.datetime64(t1), "right"))
    if i1 <= i0:
        raise ValueError(
            f"No {label} steps overlap requested window {t0} -> {t1}. "
            f"Available range is {idx[0]} -> {idx[-1]}."
        )
    return i0, i1, idx[i0:i1].values, t0, t1


def run_gridding_stage(args):
    if args.skip_gridding:
        print("\n[stage 2/2] Skipping gridding.")
        return

    print("\n[stage 2/2] Streaming gridding for CML / radar / gauges ...")
    t0 = time.perf_counter()

    ds = xr.open_dataset(args.estimation_nc)
    if "cml_id" in ds.coords:
        ds = ds.assign_coords(cml_id=ds["cml_id"].astype(str))

    try:
        if args.coords_nc is None:
            args.coords_nc = args.cml_nc
        ds = load_link_coordinates(ds, args)

        lat_center = ((ds[args.site_a_lat] + ds[args.site_b_lat]) / 2.0).values
        lon_center = ((ds[args.site_a_lon] + ds[args.site_b_lon]) / 2.0).values
        good = np.isfinite(lon_center) & np.isfinite(lat_center)
        good_idx = np.flatnonzero(good).astype("int64")
        if good_idx.size == 0:
            raise RuntimeError("No valid CML coordinates found for gridding.")
        lon = lon_center[good].astype("float64")
        lat = lat_center[good].astype("float64")

        xg2d, yg2d = build_grid(lon, lat, args.resolution)
        xgrid, ygrid = xg2d[0, :], yg2d[:, 0]
        ny, nx = xg2d.shape
        print(f"  Links for IDW: {good_idx.size}/{ds.sizes['cml_id']}")
        print(f"  Grid: {ny} x {nx} = {ny * nx:,} cells")

        R_da = ds[args.grid_rain_var]
        sel = {}
        for dim, val in [("baseline_method", args.grid_baseline_method), ("kr_model", args.grid_kr_model), ("waa_method", args.grid_waa_method)]:
            if dim in R_da.dims:
                sel[dim] = val
        R_sel = R_da.sel(sel) if sel else R_da
        if "time_5min" not in R_sel.dims:
            raise ValueError(f"{args.grid_rain_var} must have dimension time_5min. dims={R_sel.dims}")
        times_all = ds["time_5min"].values
        grid_start_value = args.grid_tstart or args.tstart
        grid_end_value = args.grid_tend or args.tend
        time_i0, time_i1, times, grid_t0, grid_t1 = _resolve_time_window_indices(
            times_all, grid_start_value, grid_end_value, label="CML/radar gridding time"
        )
        ntime = len(times)
        cml_factor = _median_step_hours(times, fallback_hours=5.0 / 60.0)
        write_radar = bool(args.with_radar and "R_radar" in ds.variables)
        if args.with_radar and "R_radar" not in ds.variables:
            print("  WARNING: --with-radar requested but R_radar not found in estimation file.")

        gauge_meta = prepare_gauge_metadata(args, xgrid, ygrid)
        gauge_times = None if gauge_meta is None else gauge_meta["times"]
        gauge_factor = None if gauge_meta is None else _median_step_hours(gauge_times, fallback_hours=10.0 / 60.0)

        print(f"  Gridding window requested: {grid_start_value} -> {grid_end_value}")
        print(f"  Gridding window used     : {pd.Timestamp(times[0])} -> {pd.Timestamp(times[-1])}")
        print(f"  Estimation file indices  : {time_i0}:{time_i1} / {len(times_all)}")
        print(f"  Time steps CML/radar: {ntime}")
        print(f"  Time chunk CML/radar: {args.grid_time_chunk}")
        print(f"  Grid workers: {args.grid_workers}")
        print(f"  IDW batch: {args.batch}")
        print(f"  IDW cell chunk: {args.cell_chunk}")
        print_mem_estimate("one CML output time chunk", min(args.grid_time_chunk, ntime), ny, nx)
        if write_radar:
            print_mem_estimate("one radar output time chunk", min(args.grid_time_chunk, ntime), ny, nx)
        if gauge_meta is not None:
            print(f"  Gauge time steps: {len(gauge_times)}")
            print(f"  Gauge time chunk: {args.gauge_time_chunk}")
            print_mem_estimate("one gauge output time chunk", min(args.gauge_time_chunk, len(gauge_times)), ny, nx)

        writer = MergedGridWriter(args.gridded_nc, times, xgrid, ygrid, args, write_radar, gauge_times=gauge_times)

        R_total = np.zeros((ny, nx), dtype="float64")
        R_radar_total = np.zeros((ny, nx), dtype="float64") if write_radar else None
        R_gauge_total = np.zeros((ny, nx), dtype="float64") if gauge_meta is not None else None

        common = {"estimation_nc": str(args.estimation_nc), "rain_var": args.grid_rain_var,
                  "baseline_method": args.grid_baseline_method, "kr_model": args.grid_kr_model, "waa_method": args.grid_waa_method,
                  "good_idx": good_idx, "lon": lon, "lat": lat, "xg2d": xg2d, "yg2d": yg2d,
                  "nnear": int(args.nnear), "p": int(args.p), "max_dist": None if args.max_dist is None else float(args.max_dist),
                  "batch": int(args.batch), "cell_chunk": int(args.cell_chunk), "idw_query_workers": int(args.idw_query_workers),
                  "write_radar": write_radar, "factor_hours": cml_factor}

        jobs = []
        for start in range(0, ntime, int(args.grid_time_chunk)):
            end = min(start + int(args.grid_time_chunk), ntime)
            job = dict(common)
            job.update({
                "kind": "cml_radar",
                "start": start,
                "end": end,
                "file_start": time_i0 + start,
                "file_end": time_i0 + end,
            })
            jobs.append(job)
        if gauge_meta is not None:
            g_common = {"kind": "gauge", "gauge_meta": gauge_meta, "xg2d": xg2d, "yg2d": yg2d,
                        "nnear": int(args.nnear), "p": int(args.p), "max_dist": None if args.max_dist is None else float(args.max_dist),
                        "batch": int(args.batch), "cell_chunk": int(args.cell_chunk), "idw_query_workers": int(args.idw_query_workers),
                        "factor_hours": gauge_factor}
            for start in range(0, len(gauge_times), int(args.gauge_time_chunk)):
                end = min(start + int(args.gauge_time_chunk), len(gauge_times))
                job = dict(g_common); job.update({"start": start, "end": end}); jobs.append(job)

        def consume(res):
            nonlocal R_total, R_radar_total, R_gauge_total
            if res["kind"] == "cml_radar":
                writer.write_time_chunk(res["start"], res["R"], res["R_radar"])
                R_total += res["R_total"]
                if R_radar_total is not None and res["R_radar_total"] is not None:
                    R_radar_total += res["R_radar_total"]
                msg = f"{res['kind']} out {res['start']}:{res['end']} file {res.get('file_start', res['start'])}:{res.get('file_end', res['end'])} {res['elapsed']:.1f}s"
                del res["R"]
                if res["R_radar"] is not None: del res["R_radar"]
                return msg
            if res["kind"] == "gauge":
                writer.write_gauge_chunk(res["start"], res["R_gauge"])
                R_gauge_total += res["R_gauge_total"]
                msg = f"gauge {res['start']}:{res['end']} {res['elapsed']:.1f}s"
                del res["R_gauge"]
                return msg
            return str(res.get("kind"))

        try:
            if int(args.grid_workers) <= 1:
                pbar = tqdm(jobs, desc="gridding", unit="chunk", dynamic_ncols=True)
                for job in pbar:
                    res = _grid_merged_job(job)
                    pbar.set_postfix_str(consume(res)); gc.collect()
            else:
                max_workers = int(args.grid_workers)
                max_pending = max(max_workers, int(args.max_pending_grid_jobs))
                job_iter = iter(jobs); pending = set()
                pbar = tqdm(total=len(jobs), desc="gridding", unit="chunk", dynamic_ncols=True)
                with ProcessPoolExecutor(max_workers=max_workers) as ex:
                    for _ in range(min(max_pending, len(jobs))):
                        if not _submit_next(ex, pending, job_iter): break
                    while pending:
                        done, pending = wait(pending, return_when=FIRST_COMPLETED)
                        for fut in done:
                            res = fut.result(); msg = consume(res)
                            pbar.update(1); pbar.set_postfix_str(msg); gc.collect()
                            _submit_next(ex, pending, job_iter)
                pbar.close()
            writer.close()
        except Exception:
            writer.close(); raise

        print(f"\n  Saved gridded NetCDF -> {args.gridded_nc}")
        if args.plot_total:
            out_png = str(Path(args.gridded_nc).with_suffix("")) + "_total.png"
            plot_total_map_merged(ds, xgrid, ygrid, R_total.astype("float32"),
                                  None if R_radar_total is None else R_radar_total.astype("float32"),
                                  None if R_gauge_total is None else R_gauge_total.astype("float32"), args, out_png)
        print(f"  Gridding elapsed: {time.perf_counter() - t0:.1f}s")
    finally:
        ds.close()


def parse_args_merged():
    ap = argparse.ArgumentParser(description="Merged CML rainfall estimation + streamed parallel IDW gridding for CML/radar/gauges.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--skip-estimation", action="store_true", help="Use an existing --estimation-nc file.")
    ap.add_argument("--skip-gridding", action="store_true", help="Run only the estimation stage.")
    ap.add_argument("--overwrite-estimation", action="store_true", help="Overwrite the estimation NetCDF if it exists.")
    ap.add_argument("--overwrite-grid", action="store_true", help="Overwrite the gridded NetCDF if it exists.")
    ap.add_argument("--cml-nc", default="cml_preprocessed_v3.nc")
    ap.add_argument("--radar-nc", default="RADAR_S2021E2021_5MIN_CZ.nc")
    ap.add_argument("--estimation-nc", default="cml_rainrate_baseline_kr_waa_methods.nc")
    ap.add_argument("--gridded-nc", default="cml_rainfall_gridded_cml_radar_gauge.nc")
    ap.add_argument("--log-csv", default="cml_rainrate_baseline_kr_waa_methods_log.csv")
    ap.add_argument("--tstart", default="2021-04-10 00:00:00",
                    help="Start time for rainfall estimation. Also used as gridding start unless --grid-tstart is given.")
    ap.add_argument("--tend", default="2021-08-16 23:59:00",
                    help="End time for rainfall estimation. Also used as gridding end unless --grid-tend is given.")
    ap.add_argument("--grid-tstart", default=None,
                    help="Optional gridding-only start time. Useful with --skip-estimation on an existing full estimation file.")
    ap.add_argument("--grid-tend", default=None,
                    help="Optional gridding-only end time. Useful with --skip-estimation on an existing full estimation file.")
    ap.add_argument("--baseline-pad-hours", type=float, default=BASELINE_PAD_HOURS_DEFAULT)
    ap.add_argument("--min-link-km", type=float, default=MIN_LINK_KM_DEFAULT)
    ap.add_argument("--min-valid-samples", type=int, default=MIN_VALID_SAMPLES_DEFAULT)
    ap.add_argument("--radar-wet-thr-5min", type=float, default=RADAR_WET_THR_5MIN_DEFAULT)
    ap.add_argument("--n-average-last-dry", type=int, default=N_AVERAGE_LAST_DRY_DEFAULT)
    ap.add_argument("--baseline-methods", nargs="+", default=list(BASELINE_METHODS_DEFAULT), choices=BASELINE_METHODS_ALLOWED)
    ap.add_argument("--kr-models", nargs="+", default=list(KR_MODELS_DEFAULT), choices=KR_MODELS_ALLOWED)
    ap.add_argument("--waa-methods", nargs="+", default=list(WAA_METHODS_DEFAULT), choices=WAA_METHODS_ALLOWED)
    ap.add_argument("--waa-gamma", type=float, default=WAA_GAMMA_DEFAULT)
    ap.add_argument("--waa-l-antenna", type=float, default=WAA_L_ANT_DEFAULT)
    ap.add_argument("--waa-delta", type=float, default=WAA_DELTA_DEFAULT)
    ap.add_argument("--pastorek-a-max", type=float, default=PASTOREK_A_MAX_DEFAULT)
    ap.add_argument("--pastorek-zeta", type=float, default=PASTOREK_ZETA_DEFAULT)
    ap.add_argument("--pastorek-d", type=float, default=PASTOREK_D_DEFAULT)
    ap.add_argument("--schleiss-waa-max", type=float, default=SCHLEISS_WAA_MAX_DEFAULT)
    ap.add_argument("--schleiss-delta-t", type=float, default=SCHLEISS_DELTA_T_DEFAULT)
    ap.add_argument("--schleiss-tau", type=float, default=SCHLEISS_TAU_DEFAULT)
    ap.add_argument("--only-cml-id", default=None)
    ap.add_argument("--max-links", type=int, default=None)
    ap.add_argument("--est-chunk-time", type=int, default=100_000, help="NetCDF chunk time used by estimation writer.")
    ap.add_argument("--write-1min", action="store_true", help="Also write 1-min estimation variables. Large output.")
    ap.add_argument("--grid-rain-var", default="R_mean_5min")
    ap.add_argument("--grid-baseline-method", default="pycml_constant")
    ap.add_argument("--grid-kr-model", default="ITU_2003")
    ap.add_argument("--grid-waa-method", default="pastorek_2021")
    ap.add_argument("--coords-nc", default=None, help="Coordinate file. Defaults to --cml-nc.")
    ap.add_argument("--with-radar", dest="with_radar", action="store_true", default=True)
    ap.add_argument("--no-radar", dest="with_radar", action="store_false")
    ap.add_argument("--with-gauge", dest="with_gauge", action="store_true", default=True)
    ap.add_argument("--no-gauge", dest="with_gauge", action="store_false")
    ap.add_argument("--gauge-nc", default="CHMI_10MIN_v1.nc")
    ap.add_argument("--gauge-var", default="rr")
    ap.add_argument("--gauge-time-dim", default=None)
    ap.add_argument("--gauge-station-dim", default=None)
    ap.add_argument("--gauge-lon", default=None)
    ap.add_argument("--gauge-lat", default=None)
    ap.add_argument("--gauge-is-10min-amount", dest="gauge_is_10min_amount", action="store_true", default=True)
    ap.add_argument("--gauge-is-rate", dest="gauge_is_10min_amount", action="store_false")
    ap.add_argument("--gauge-bbox-pad", type=float, default=None)
    ap.add_argument("--gauge-time-chunk", type=int, default=144)
    ap.add_argument("--nnear", type=int, default=15)
    ap.add_argument("--p", type=int, default=2)
    ap.add_argument("--max-dist", type=float, default=0.3)
    ap.add_argument("--resolution", type=float, default=0.01)
    ap.add_argument("--grid-time-chunk", type=int, default=288)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--cell-chunk", type=int, default=20000)
    ap.add_argument("--grid-workers", type=int, default=1, help="Parallel gridding worker processes. Use 2-4 on 128 GB; keep --idw-query-workers 1.")
    ap.add_argument("--max-pending-grid-jobs", type=int, default=4, help="Limits queued chunk jobs to control RAM.")
    ap.add_argument("--idw-query-workers", type=int, default=1, help="Threads inside scipy cKDTree.query per process. Avoid -1 when --grid-workers > 1.")
    ap.add_argument("--y-chunk", type=int, default=None)
    ap.add_argument("--x-chunk", type=int, default=None)
    ap.add_argument("--complevel", type=int, default=4)
    ap.add_argument("--no-compression", action="store_true")
    ap.add_argument("--plot-total", dest="plot_total", action="store_true", default=True)
    ap.add_argument("--no-plot", dest="plot_total", action="store_false")
    ap.add_argument("--site-a-lat", default="site_a_latitude")
    ap.add_argument("--site-a-lon", default="site_a_longitude")
    ap.add_argument("--site-b-lat", default="site_b_latitude")
    ap.add_argument("--site-b-lon", default="site_b_longitude")
    args = ap.parse_args()
    if args.grid_workers > 1 and args.idw_query_workers != 1:
        print("WARNING: --grid-workers > 1 with --idw-query-workers != 1 can oversubscribe CPUs/RAM.")
    return args


def main():
    args = parse_args_merged()
    total_t0 = time.perf_counter()
    run_estimation_stage(args)
    run_gridding_stage(args)
    print(f"\nAll done in {time.perf_counter() - total_t0:.1f}s")


if __name__ == "__main__":
    main()
