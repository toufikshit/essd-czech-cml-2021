#!/usr/bin/env python3
"""
Fast batch selection and plotting of convective CML-vs-radar instances
=====================================================================

This script finds strong radar events, resamples the gridded CML rain-rate
product to the native radar temporal resolution, and saves one comparison
figure per selected time instance.

Each output figure contains exactly two vertically stacked panels:
    (a) CML rain rate
    (b) Native radar rain rate

Scientific/event-selection logic
--------------------------------
1. Discover daily native-radar files recursively.
2. Read the native radar time axis and infer its temporal resolution.
3. Resample the CML rain-rate field to the same radar resolution using a mean,
   because both inputs are rain rates in mm h-1.
4. Detect candidate convective instances from radar using:
       radar >= 10 mm h-1
   with configurable minimum convective-pixel count and peak intensity.
5. Group nearby candidate times into separate storm events.
6. Select several temporally separated high-ranking instances from multiple
   events.
7. Save one ESSD-style PNG and PDF per instance plus a CSV summary.

Speed/memory choices
--------------------
* Radar files are opened once and cached.
* Event scanning uses only cropped radar arrays and optionally a scan stride.
* CML is resampled once per radar resolution and reused.
* Only selected timestamps are rendered.
* Cartopy background is intentionally restrained for fast publication output.
* Figures are closed immediately after saving.

Important units
---------------
RADAR_SCALE_TO_MMPH must convert the native radar variable to mm h-1.
For the supplied RADOLAN/YW convention, the existing workflow used 12.0.

Example
-------
python batch_cml_vs_radar_convective_events_essd_vertical.py \
  --gridded-nc cml_rainfall_gridded_selected_period.nc \
  --cml-geom-nc cml_preprocessed_v3.nc \
  --radar-root dwd_radolan \
  --convective-threshold 10 \
  --min-convective-pixels 100 \
  --instances-per-event 3 \
  --max-events 6
"""

from __future__ import annotations

import argparse
import gc
import glob
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import AutoMinorLocator

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    def tqdm(x, **kwargs):
        return x


# ============================================================
# ESSD-like publication style
# ============================================================
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

DEFAULT_BOUNDS = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 50, 80]
RADAR_VAR_CANDIDATES = ["YW", "RR", "rainfall_rate", "R", "rainrate", "precip"]


# ============================================================
# Data classes
# ============================================================
@dataclass
class RadarFileInfo:
    date: pd.Timestamp
    path: str


@dataclass
class Candidate:
    time: pd.Timestamp
    radar_path: str
    extent: Tuple[float, float, float, float]
    convective_pixels: int
    wet_pixels: int
    radar_max: float
    radar_p99: float
    radar_mean_convective: float
    cml_wet_pixels: int
    joint_wet_pixels: int
    cml_available_pixels: int
    cml_available_fraction: float
    links_in_window: int
    cml_max: float
    cml_p99: float
    score: float


@dataclass
class Event:
    event_id: int
    start: pd.Timestamp
    end: pd.Timestamp
    candidates: List[Candidate]


# ============================================================
# Basic helpers
# ============================================================
def pick_var(ds: xr.Dataset, candidates: Sequence[str]) -> str:
    for name in candidates:
        if name in ds.data_vars:
            return name
    raise KeyError(f"None of {list(candidates)} found. Available: {list(ds.data_vars)}")


def find_time_dim(da: xr.DataArray) -> Optional[str]:
    return next((d for d in da.dims if "time" in d.lower()), None)


def infer_step_minutes(times: pd.DatetimeIndex, fallback: float = 5.0) -> float:
    if len(times) < 2:
        return float(fallback)
    delta_s = np.diff(times.view("int64")) / 1e9
    delta_s = delta_s[np.isfinite(delta_s) & (delta_s > 0)]
    if delta_s.size == 0:
        return float(fallback)
    return float(np.median(delta_s) / 60.0)


def discover_radar_files(root: str, pattern: str) -> List[RadarFileInfo]:
    files = sorted(glob.glob(os.path.join(root, "**", pattern), recursive=True))
    out: List[RadarFileInfo] = []
    for path in files:
        m = re.search(r"(\d{8})", os.path.basename(path))
        if m:
            out.append(RadarFileInfo(pd.Timestamp(m.group(1)), path))
    return out


def parse_date_filter(values: Optional[Sequence[str]]) -> Optional[set]:
    if not values:
        return None
    return {pd.Timestamp(v).date() for v in values}


def dry_mask(a: np.ndarray, threshold: float) -> np.ndarray:
    out = np.asarray(a, dtype=float).copy()
    out[out < threshold] = np.nan
    return out


def safe_quantile(a: np.ndarray, q: float) -> float:
    vals = np.asarray(a, dtype=float)
    vals = vals[np.isfinite(vals)]
    return float(np.nanquantile(vals, q)) if vals.size else np.nan


def choose_tick_step(extent: Sequence[float]) -> float:
    span = max(extent[1] - extent[0], extent[3] - extent[2])
    if span <= 1.5:
        return 0.25
    if span <= 3:
        return 0.5
    return 1.0


# ============================================================
# Radar handling
# ============================================================
class RadarVariableNotFoundError(RuntimeError):
    """Raised when a NetCDF file has coordinates/metadata but no radar rain field."""


class RadarCache:
    def __init__(self, var_candidates: Sequence[str], scale_to_mmph: float):
        self.var_candidates = list(var_candidates)
        self.scale = float(scale_to_mmph)
        self._datasets: Dict[str, xr.Dataset] = {}
        self._metadata: Dict[str, Tuple[str, Optional[str], pd.DatetimeIndex]] = {}

    def open(self, path: str) -> Tuple[xr.Dataset, xr.DataArray, Optional[str], pd.DatetimeIndex]:
        if path not in self._datasets:
            # Radar files are scanned sequentially. Keep only the current file open
            # to avoid exhausting file descriptors and RAM during year-long scans.
            self.close()
            ds = xr.open_dataset(path)
            try:
                var = pick_var(ds, self.var_candidates)
            except KeyError as exc:
                available = list(ds.data_vars)
                ds.close()
                raise RadarVariableNotFoundError(
                    f"No radar rainfall variable in {path!r}. "
                    f"Available data variables: {available}"
                ) from exc

            da = ds[var]
            tdim = find_time_dim(da)
            times = pd.DatetimeIndex(da[tdim].values) if tdim else pd.DatetimeIndex([])
            self._datasets[path] = ds
            self._metadata[path] = (var, tdim, times)
        ds = self._datasets[path]
        var, tdim, times = self._metadata[path]
        return ds, ds[var], tdim, times

    def field_at(self, path: str, target_time: pd.Timestamp) -> Tuple[np.ndarray, np.ndarray, np.ndarray, pd.Timestamp]:
        ds, da, tdim, _ = self.open(path)
        if tdim is not None:
            selected = da.sel({tdim: target_time}, method="nearest")
            actual = pd.Timestamp(selected[tdim].values)
        else:
            selected = da
            actual = pd.Timestamp(target_time)

        field = np.asarray(selected.values, dtype="float32") * np.float32(self.scale)
        field = np.where(np.isfinite(field), field, np.nan).astype("float32", copy=False)

        latname = next((c for c in ["lat", "latitude"] if c in ds.coords or c in ds.variables), None)
        lonname = next((c for c in ["lon", "longitude"] if c in ds.coords or c in ds.variables), None)
        if latname and lonname:
            lon = np.asarray(ds[lonname].values)
            lat = np.asarray(ds[latname].values)
            return lon, lat, field, actual

        if {"x", "y"}.issubset(set(ds.coords) | set(ds.variables)):
            X, Y = np.meshgrid(np.asarray(ds["x"].values), np.asarray(ds["y"].values))
            radolan = ccrs.Stereographic(
                central_latitude=90.0,
                central_longitude=10.0,
                true_scale_latitude=60.0,
                globe=ccrs.Globe(
                    ellipse=None,
                    semimajor_axis=6370040,
                    semiminor_axis=6370040,
                ),
            )
            pts = ccrs.PlateCarree().transform_points(radolan, X, Y)
            return pts[..., 0], pts[..., 1], field, actual

        raise RuntimeError(f"No lat/lon or projected x/y coordinates in {path}")

    def close(self):
        for ds in self._datasets.values():
            ds.close()
        self._datasets.clear()
        self._metadata.clear()


# ============================================================
# Spatial helpers
# ============================================================
def extent_mask(lon: np.ndarray, lat: np.ndarray, extent: Sequence[float]) -> np.ndarray:
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    return (
        np.isfinite(lon) & np.isfinite(lat) &
        (lon >= extent[0]) & (lon <= extent[1]) &
        (lat >= extent[2]) & (lat <= extent[3])
    )


def derive_dense_cml_extent(
    lon_a: np.ndarray,
    lat_a: np.ndarray,
    lon_b: np.ndarray,
    lat_b: np.ndarray,
    half_lon: float,
    half_lat: float,
    bin_size: float = 0.1,
) -> List[float]:
    lon_c = 0.5 * (lon_a + lon_b)
    lat_c = 0.5 * (lat_a + lat_b)
    valid = np.isfinite(lon_c) & np.isfinite(lat_c)
    lon_c = lon_c[valid]
    lat_c = lat_c[valid]
    if lon_c.size == 0:
        raise ValueError("No finite CML midpoint coordinates.")

    xedges = np.arange(lon_c.min(), lon_c.max() + bin_size, bin_size)
    yedges = np.arange(lat_c.min(), lat_c.max() + bin_size, bin_size)
    H, xe, ye = np.histogram2d(lon_c, lat_c, bins=[xedges, yedges])
    ix, iy = np.unravel_index(np.argmax(H), H.shape)
    cx = 0.5 * (xe[ix] + xe[ix + 1])
    cy = 0.5 * (ye[iy] + ye[iy + 1])
    return [cx - half_lon, cx + half_lon, cy - half_lat, cy + half_lat]


def rain_centroid(field: np.ndarray, lon: np.ndarray, lat: np.ndarray, threshold: float) -> Optional[Tuple[float, float]]:
    m = np.isfinite(field) & (field >= threshold) & np.isfinite(lon) & np.isfinite(lat)
    if not np.any(m):
        return None
    w = np.where(m, field, 0.0)
    sw = float(np.nansum(w))
    if sw <= 0:
        return None
    return float(np.nansum(lon * w) / sw), float(np.nansum(lat * w) / sw)


# ============================================================
# Joint CML-radar spatial selection
# ============================================================
def network_search_extent(lon_a, lat_a, lon_b, lat_b, pad=0.05):
    lon = np.concatenate([np.asarray(lon_a), np.asarray(lon_b)])
    lat = np.concatenate([np.asarray(lat_a), np.asarray(lat_b)])
    ok = np.isfinite(lon) & np.isfinite(lat)
    if not np.any(ok):
        raise ValueError("No finite CML endpoint coordinates.")
    return [float(lon[ok].min()-pad), float(lon[ok].max()+pad),
            float(lat[ok].min()-pad), float(lat[ok].max()+pad)]


def link_density_on_grid(xg, yg, lon_a, lat_a, lon_b, lat_b):
    """Count CML midpoints in the nearest regular-grid cell."""
    lon_c = 0.5 * (np.asarray(lon_a, dtype=float) + np.asarray(lon_b, dtype=float))
    lat_c = 0.5 * (np.asarray(lat_a, dtype=float) + np.asarray(lat_b, dtype=float))
    out = np.zeros((len(yg), len(xg)), dtype=np.int32)
    ok = np.isfinite(lon_c) & np.isfinite(lat_c)
    if not np.any(ok):
        return out
    ix = np.searchsorted(xg, lon_c[ok])
    iy = np.searchsorted(yg, lat_c[ok])
    ix = np.clip(ix, 1, len(xg)-1)
    iy = np.clip(iy, 1, len(yg)-1)
    ix -= (np.abs(lon_c[ok]-xg[ix-1]) <= np.abs(lon_c[ok]-xg[ix])).astype(int)
    iy -= (np.abs(lat_c[ok]-yg[iy-1]) <= np.abs(lat_c[ok]-yg[iy])).astype(int)
    np.add.at(out, (iy, ix), 1)
    return out


def box_sum_same(a: np.ndarray, half_y: int, half_x: int) -> np.ndarray:
    """Fast centered rectangular-window sum using an integral image."""
    a = np.asarray(a, dtype=np.float64)
    ny, nx = a.shape
    py, px = max(0, int(half_y)), max(0, int(half_x))
    padded = np.pad(a, ((py, py), (px, px)), mode="constant")
    ii = np.pad(padded, ((1, 0), (1, 0)), mode="constant").cumsum(0).cumsum(1)
    h, w = 2*py + 1, 2*px + 1
    out = ii[h:, w:] - ii[:-h, w:] - ii[h:, :-w] + ii[:-h, :-w]
    return out[:ny, :nx]


class RadarToCMLGridMapper:
    """Nearest-neighbour radar-to-CML-grid mapping, cached by radar grid signature."""
    def __init__(self, xg, yg):
        self.X, self.Y = np.meshgrid(np.asarray(xg), np.asarray(yg))
        self.target = np.column_stack([self.X.ravel(), self.Y.ravel()])
        self.cache = {}

    @staticmethod
    def _signature(lon, lat):
        lon = np.asarray(lon); lat = np.asarray(lat)
        return (lon.shape, lat.shape,
                float(np.nanmin(lon)), float(np.nanmax(lon)),
                float(np.nanmin(lat)), float(np.nanmax(lat)))

    def map(self, lon, lat, field):
        lon = np.asarray(lon); lat = np.asarray(lat); field = np.asarray(field)
        if lon.ndim == 1 and lat.ndim == 1:
            lon, lat = np.meshgrid(lon, lat)
        sig = self._signature(lon, lat)
        if sig not in self.cache:
            valid = np.isfinite(lon) & np.isfinite(lat)
            src_flat = np.flatnonzero(valid.ravel())
            tree = cKDTree(np.column_stack([lon.ravel()[src_flat], lat.ravel()[src_flat]]))
            _, jj = tree.query(self.target, k=1, workers=-1)
            self.cache[sig] = src_flat[jj]
        idx = self.cache[sig]
        return field.ravel()[idx].reshape(self.X.shape).astype("float32", copy=False)


def best_joint_window(cml_field, radar_on_cml, xg, yg, link_density, search_extent,
                      half_lon, half_lat, cml_rain_threshold, radar_conv_threshold,
                      dry_threshold, min_links, min_avail_fraction,
                      min_cml_wet_pixels, min_radar_conv_pixels, min_joint_pixels):
    """Find the fixed-size window with strong radar+CML rain and dense CML support."""
    dx = float(np.median(np.diff(xg))); dy = float(np.median(np.diff(yg)))
    hx = max(1, int(round(half_lon / abs(dx))))
    hy = max(1, int(round(half_lat / abs(dy))))
    X, Y = np.meshgrid(xg, yg)
    search = ((X >= search_extent[0]) & (X <= search_extent[1]) &
              (Y >= search_extent[2]) & (Y <= search_extent[3]))

    cml_avail = np.isfinite(cml_field)
    cml_wet = cml_avail & (cml_field >= cml_rain_threshold)
    radar_finite = np.isfinite(radar_on_cml)
    radar_wet = radar_finite & (radar_on_cml >= dry_threshold)
    radar_conv = radar_finite & (radar_on_cml >= radar_conv_threshold)
    joint = cml_wet & radar_conv

    win_cells = box_sum_same(np.ones_like(cml_field, dtype=np.uint8), hy, hx)
    avail_n = box_sum_same(cml_avail, hy, hx)
    cml_wet_n = box_sum_same(cml_wet, hy, hx)
    radar_wet_n = box_sum_same(radar_wet, hy, hx)
    radar_conv_n = box_sum_same(radar_conv, hy, hx)
    joint_n = box_sum_same(joint, hy, hx)
    links_n = box_sum_same(link_density, hy, hx)
    avail_frac = np.divide(avail_n, win_cells, out=np.zeros_like(avail_n), where=win_cells > 0)

    # Intensity terms are evaluated only where both sensors indicate rain.
    joint_radar_sum = box_sum_same(np.where(joint, radar_on_cml, 0.0), hy, hx)
    joint_cml_sum = box_sum_same(np.where(joint, cml_field, 0.0), hy, hx)

    valid = (search & (links_n >= min_links) & (avail_frac >= min_avail_fraction) &
             (cml_wet_n >= min_cml_wet_pixels) &
             (radar_conv_n >= min_radar_conv_pixels) &
             (joint_n >= min_joint_pixels))
    if not np.any(valid):
        return None

    # Joint overlap dominates; density and intensity break ties.
    score = (1000.0 * joint_n + 20.0 * np.minimum(cml_wet_n, radar_conv_n) +
             5.0 * links_n + 2.0 * joint_radar_sum + 2.0 * joint_cml_sum +
             0.2 * radar_wet_n)
    score = np.where(valid, score, -np.inf)
    iy, ix = np.unravel_index(np.nanargmax(score), score.shape)
    extent = (float(xg[ix]-half_lon), float(xg[ix]+half_lon),
              float(yg[iy]-half_lat), float(yg[iy]+half_lat))

    xmask = (xg >= extent[0]) & (xg <= extent[1])
    ymask = (yg >= extent[2]) & (yg <= extent[3])
    sl = np.ix_(ymask, xmask)
    csub = cml_field[sl]; rsub = radar_on_cml[sl]
    jsub = joint[sl]
    convsub = radar_conv[sl]
    wetsub = radar_wet[sl]
    return {
        "extent": extent,
        "score": float(score[iy, ix]),
        "links": int(round(links_n[iy, ix])),
        "available": int(round(avail_n[iy, ix])),
        "available_fraction": float(avail_frac[iy, ix]),
        "cml_wet": int(round(cml_wet_n[iy, ix])),
        "radar_conv": int(round(radar_conv_n[iy, ix])),
        "radar_wet": int(round(radar_wet_n[iy, ix])),
        "joint": int(round(joint_n[iy, ix])),
        "radar_max": float(np.nanmax(rsub)) if np.isfinite(rsub).any() else np.nan,
        "radar_p99": safe_quantile(rsub[np.isfinite(rsub)], 0.99),
        "radar_mean_conv": float(np.nanmean(rsub[convsub])) if np.any(convsub) else np.nan,
        "cml_max": float(np.nanmax(csub)) if np.isfinite(csub).any() else np.nan,
        "cml_p99": safe_quantile(csub[np.isfinite(csub)], 0.99),
    }


# ============================================================
# Temporal alignment
# ============================================================
def resample_cml_to_radar_resolution(
    cml_da: xr.DataArray,
    radar_minutes: float,
    label: str = "left",
    closed: str = "left",
) -> xr.DataArray:
    """Resample CML rain rate to radar resolution using arithmetic mean."""
    time_dim = find_time_dim(cml_da)
    if time_dim is None:
        raise ValueError(f"Cannot find time dimension in CML variable dims {cml_da.dims}")
    if time_dim != "time":
        cml_da = cml_da.rename({time_dim: "time"})

    freq_seconds = int(round(float(radar_minutes) * 60.0))
    if freq_seconds <= 0:
        raise ValueError(f"Invalid radar resolution: {radar_minutes} min")
    freq = f"{freq_seconds}s"

    cml_times = pd.DatetimeIndex(cml_da.time.values)
    cml_minutes = infer_step_minutes(cml_times)
    if np.isclose(cml_minutes, radar_minutes, rtol=0, atol=1e-6):
        print(f"  CML already has radar resolution: {radar_minutes:g} min")
        return cml_da

    print(f"  Resampling CML rate: {cml_minutes:g}-min -> {radar_minutes:g}-min using mean")
    return cml_da.resample(time=freq, label=label, closed=closed).mean(skipna=True)


# ============================================================
# Event scanning and grouping
# ============================================================
def scan_candidates(
    radar_files: Sequence[RadarFileInfo],
    cache: RadarCache,
    cml_aligned: xr.DataArray,
    cml_times: pd.DatetimeIndex,
    xg: np.ndarray,
    yg: np.ndarray,
    link_density: np.ndarray,
    mapper: RadarToCMLGridMapper,
    search_extent: Sequence[float],
    convective_threshold: float,
    cml_rain_threshold: float,
    wet_threshold: float,
    min_convective_pixels: int,
    min_cml_wet_pixels: int,
    min_joint_wet_pixels: int,
    min_links_in_window: int,
    min_cml_availability_fraction: float,
    min_peak_mmph: float,
    zoom_half_lon: float,
    zoom_half_lat: float,
    scan_stride: int,
    time_tolerance_min: float,
    date_filter: Optional[set],
) -> List[Candidate]:
    candidates: List[Candidate] = []
    skipped_no_rain_var = 0
    skipped_unreadable = 0
    skipped_no_cml_time = 0
    radar_prescreen_rejected = 0
    joint_rejected = 0

    pbar = tqdm(radar_files, desc="Scanning joint radar-CML events", unit="day")
    for info in pbar:
        if date_filter is not None and info.date.date() not in date_filter:
            continue
        try:
            ds, da, tdim, times = cache.open(info.path)
        except RadarVariableNotFoundError:
            skipped_no_rain_var += 1
            continue
        except Exception:
            skipped_unreadable += 1
            continue

        sample_times = [info.date] if tdim is None or len(times) == 0 else list(times[::max(1, int(scan_stride))])
        radar_search_mask = None

        for target_time in sample_times:
            lon, lat, field, actual = cache.field_at(info.path, pd.Timestamp(target_time))
            if radar_search_mask is None:
                radar_search_mask = extent_mask(lon, lat, search_extent)
            rdomain = field[radar_search_mask]
            rfinite = rdomain[np.isfinite(rdomain)]
            if rfinite.size == 0:
                continue

            # Cheap radar-only prescreen. CML is loaded only for potentially useful times.
            if float(np.nanmax(rfinite)) < min_peak_mmph or int(np.sum(rfinite >= convective_threshold)) < min_joint_wet_pixels:
                radar_prescreen_rejected += 1
                continue

            pos = cml_times.get_indexer([pd.Timestamp(actual)], method="nearest")[0]
            if pos < 0:
                skipped_no_cml_time += 1
                continue
            cml_time = cml_times[pos]
            mismatch = abs((cml_time - pd.Timestamp(actual)).total_seconds()) / 60.0
            if mismatch > time_tolerance_min:
                skipped_no_cml_time += 1
                continue

            cml_field = np.asarray(cml_aligned.isel(time=int(pos)).values, dtype="float32")
            radar_on_cml = mapper.map(lon, lat, field)
            best = best_joint_window(
                cml_field=cml_field,
                radar_on_cml=radar_on_cml,
                xg=xg, yg=yg,
                link_density=link_density,
                search_extent=search_extent,
                half_lon=zoom_half_lon,
                half_lat=zoom_half_lat,
                cml_rain_threshold=cml_rain_threshold,
                radar_conv_threshold=convective_threshold,
                dry_threshold=wet_threshold,
                min_links=min_links_in_window,
                min_avail_fraction=min_cml_availability_fraction,
                min_cml_wet_pixels=min_cml_wet_pixels,
                min_radar_conv_pixels=min_convective_pixels,
                min_joint_pixels=min_joint_wet_pixels,
            )
            if best is None:
                joint_rejected += 1
                continue

            candidates.append(Candidate(
                time=pd.Timestamp(actual),
                radar_path=info.path,
                extent=best["extent"],
                convective_pixels=best["radar_conv"],
                wet_pixels=best["radar_wet"],
                radar_max=best["radar_max"],
                radar_p99=best["radar_p99"],
                radar_mean_convective=best["radar_mean_conv"],
                cml_wet_pixels=best["cml_wet"],
                joint_wet_pixels=best["joint"],
                cml_available_pixels=best["available"],
                cml_available_fraction=best["available_fraction"],
                links_in_window=best["links"],
                cml_max=best["cml_max"],
                cml_p99=best["cml_p99"],
                score=best["score"],
            ))

        pbar.set_postfix_str(
            f"joint={len(candidates)} | noCML={skipped_no_cml_time} | rejected={joint_rejected}"
        )

    print(
        f"  Scan summary: {len(candidates)} joint candidates; "
        f"{skipped_no_rain_var} files without rain, {skipped_unreadable} bad files, "
        f"{skipped_no_cml_time} radar times without matching CML, "
        f"{radar_prescreen_rejected} radar-prescreen rejects, {joint_rejected} joint-window rejects."
    )

    unique: Dict[pd.Timestamp, Candidate] = {}
    for c in candidates:
        old = unique.get(c.time)
        if old is None or c.score > old.score:
            unique[c.time] = c
    return sorted(unique.values(), key=lambda c: c.time)


def group_candidates_into_events(candidates: Sequence[Candidate], event_gap_minutes: float) -> List[Event]:
    if not candidates:
        return []

    groups: List[List[Candidate]] = [[candidates[0]]]
    for candidate in candidates[1:]:
        gap = (candidate.time - groups[-1][-1].time).total_seconds() / 60.0
        if gap <= event_gap_minutes:
            groups[-1].append(candidate)
        else:
            groups.append([candidate])

    events = [
        Event(i + 1, group[0].time, group[-1].time, list(group))
        for i, group in enumerate(groups)
    ]
    return events


def select_instances_from_events(
    events: Sequence[Event],
    max_events: int,
    instances_per_event: int,
    min_instance_gap_minutes: float,
) -> List[Tuple[int, int, Candidate]]:
    # Rank whole events by their strongest candidate.
    ranked_events = sorted(
        events,
        key=lambda e: max(c.score for c in e.candidates),
        reverse=True,
    )[:max_events]

    selected: List[Tuple[int, int, Candidate]] = []
    for event_rank, event in enumerate(ranked_events, 1):
        ranked_candidates = sorted(event.candidates, key=lambda c: c.score, reverse=True)
        chosen: List[Candidate] = []
        for candidate in ranked_candidates:
            if all(
                abs((candidate.time - old.time).total_seconds()) / 60.0 >= min_instance_gap_minutes
                for old in chosen
            ):
                chosen.append(candidate)
            if len(chosen) >= instances_per_event:
                break

        chosen.sort(key=lambda c: c.time)
        for instance_rank, candidate in enumerate(chosen, 1):
            selected.append((event_rank, instance_rank, candidate))
    return selected


# ============================================================
# Plotting
# ============================================================
def setup_map_axis(ax, extent: Sequence[float], proj, show_xlabels: bool = True):
    ax.set_extent(extent, crs=proj)

    # Restrained and fast map context suitable for journal figures.
    ax.add_feature(cfeature.LAND, facecolor="0.98", zorder=0)
    ax.add_feature(cfeature.BORDERS, linewidth=0.55, edgecolor="0.25", zorder=1)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.45, edgecolor="0.35", zorder=1)
    ax.add_feature(cfeature.LAKES, edgecolor="0.65", facecolor="none", linewidth=0.35, zorder=1)

    ax.gridlines(
        crs=proj,
        draw_labels=False,
        linewidth=0.3,
        color="0.72",
        alpha=0.8,
        linestyle="--",
        zorder=0,
    )

    step = choose_tick_step(extent)
    xticks = np.arange(np.ceil(extent[0] / step) * step, extent[1] + 1e-9, step)
    yticks = np.arange(np.ceil(extent[2] / step) * step, extent[3] + 1e-9, step)
    ax.set_xticks(xticks, crs=proj)
    ax.set_yticks(yticks, crs=proj)
    ax.xaxis.set_major_formatter(LongitudeFormatter(number_format=".1f", degree_symbol="°"))
    ax.yaxis.set_major_formatter(LatitudeFormatter(number_format=".1f", degree_symbol="°"))

    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis="both", which="major", labelsize=8.5, length=4, width=0.8)
    ax.tick_params(axis="both", which="minor", length=2.2, width=0.5)
    ax.tick_params(top=False, right=False)
    if not show_xlabels:
        ax.tick_params(axis="x", labelbottom=False)

    try:
        ax.spines["geo"].set_linewidth(0.8)
    except Exception:
        pass


def plot_links(ax, proj, lon_a, lat_a, lon_b, lat_b, extent):
    # Plot only links intersecting the visual domain.
    keep = (
        np.isfinite(lon_a) & np.isfinite(lat_a) & np.isfinite(lon_b) & np.isfinite(lat_b) &
        (((lon_a >= extent[0]) & (lon_a <= extent[1]) & (lat_a >= extent[2]) & (lat_a <= extent[3])) |
         ((lon_b >= extent[0]) & (lon_b <= extent[1]) & (lat_b >= extent[2]) & (lat_b <= extent[3])))
    )
    for i in np.where(keep)[0]:
        ax.plot(
            [lon_a[i], lon_b[i]],
            [lat_a[i], lat_b[i]],
            transform=proj,
            color="0.25",
            linewidth=0.28,
            alpha=0.35,
            zorder=3,
        )


def make_vertical_instance_figure(
    cml_field: np.ndarray,
    cml_x: np.ndarray,
    cml_y: np.ndarray,
    radar_field: np.ndarray,
    radar_lon: np.ndarray,
    radar_lat: np.ndarray,
    actual_time: pd.Timestamp,
    extent: Sequence[float],
    geom: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    args,
    out_png: Path,
    out_pdf: Path,
):
    proj = ccrs.PlateCarree()
    bounds = np.asarray(args.bounds, dtype=float)
    norm = BoundaryNorm(bounds, ncolors=256, extend="both")
    cmap = mpl.colormaps[args.cmap].copy()
    cmap.set_bad("white")

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(6.3, 8.0),
        subplot_kw={"projection": proj},
        squeeze=False,
    )
    ax_cml, ax_radar = axes[:, 0]
    fig.subplots_adjust(left=0.12, right=0.86, bottom=0.08, top=0.94, hspace=0.16)

    cml_plot = dry_mask(cml_field, args.dry_threshold)
    radar_plot = dry_mask(radar_field, args.dry_threshold)
    lon_a, lat_a, lon_b, lat_b = geom

    setup_map_axis(ax_cml, extent, proj, show_xlabels=False)
    if args.draw_links:
        plot_links(ax_cml, proj, lon_a, lat_a, lon_b, lat_b, extent)
    pc = ax_cml.pcolormesh(
        cml_x,
        cml_y,
        cml_plot,
        transform=proj,
        shading="nearest",
        cmap=cmap,
        norm=norm,
        zorder=2,
        rasterized=True,
    )
    ax_cml.text(
        0.018, 0.975, "(a)",
        transform=ax_cml.transAxes,
        ha="left", va="top",
        fontsize=11, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="0.5", lw=0.55, alpha=0.92),
        zorder=6,
    )

    setup_map_axis(ax_radar, extent, proj, show_xlabels=True)
    if args.draw_links:
        plot_links(ax_radar, proj, lon_a, lat_a, lon_b, lat_b, extent)
    ax_radar.pcolormesh(
        radar_lon,
        radar_lat,
        radar_plot,
        transform=proj,
        shading="nearest",
        cmap=cmap,
        norm=norm,
        zorder=2,
        rasterized=True,
    )
    ax_radar.text(
        0.018, 0.975, "(b)",
        transform=ax_radar.transAxes,
        ha="left", va="top",
        fontsize=11, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="0.5", lw=0.55, alpha=0.92),
        zorder=6,
    )

    # One shared time label, no panel titles inside the maps.
    fig.text(
        0.49, 0.965,
        actual_time.strftime("%d %B %Y, %H:%M UTC"),
        ha="center", va="top", fontsize=10,
    )

    p_top = ax_cml.get_position()
    p_bottom = ax_radar.get_position()
    cax = fig.add_axes([0.885, p_bottom.y0, 0.020, p_top.y1 - p_bottom.y0])
    cb = fig.colorbar(pc, cax=cax, extend="both", ticks=bounds)
    cb.set_label(r"Rain rate (mm h$^{-1}$)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


# ============================================================
# CLI
# ============================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Select joint radar+CML convective events in dense/available CML regions and save vertical ESSD-style comparisons."
    )

    p.add_argument("--gridded-nc", default="cml_rainfall_gridded_selected_period.nc")
    p.add_argument("--cml-var", default="R")
    p.add_argument("--cml-geom-nc", default="cml_preprocessed_v3.nc")
    p.add_argument("--site-a-lon-var", default="site_a_longitude")
    p.add_argument("--site-a-lat-var", default="site_a_latitude")
    p.add_argument("--site-b-lon-var", default="site_b_longitude")
    p.add_argument("--site-b-lat-var", default="site_b_latitude")

    p.add_argument("--radar-root", default="dwd_radolan")
    p.add_argument("--radar-pattern", default="YW_2017.002_*.nc")
    p.add_argument("--radar-scale-to-mmph", type=float, default=12.0)
    p.add_argument("--only-dates", nargs="*", default=None)

    p.add_argument("--extent", nargs=4, type=float, default=None,
                   metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"),
                   help="Fixed plot/event-search extent. Default: densest CML region.")
    p.add_argument("--zoom-half-lon", type=float, default=0.9)
    p.add_argument("--zoom-half-lat", type=float, default=0.75)
    p.add_argument("--recenter-on-event", action="store_true", default=True)
    p.add_argument("--no-recenter-on-event", dest="recenter_on_event", action="store_false")

    p.add_argument("--convective-threshold", type=float, default=10.0,
                   help="Radar threshold defining convective pixels in mm h-1.")
    p.add_argument("--dry-threshold", type=float, default=0.1)
    p.add_argument("--min-convective-pixels", type=int, default=100)
    p.add_argument("--min-peak-mmph", type=float, default=10.0)
    p.add_argument("--cml-rain-threshold", type=float, default=0.1,
                   help="Minimum CML rain rate counted as actual CML rain in the joint window.")
    p.add_argument("--min-cml-wet-pixels", type=int, default=40,
                   help="Minimum CML grid cells with rain in the selected zoom window.")
    p.add_argument("--min-joint-wet-pixels", type=int, default=15,
                   help="Minimum cells where CML is raining and radar is >= convective threshold.")
    p.add_argument("--min-links-in-window", type=int, default=40,
                   help="Minimum number of CML midpoints inside the selected zoom window.")
    p.add_argument("--min-cml-availability-fraction", type=float, default=0.15,
                   help="Minimum finite CML-grid fraction inside the selected zoom window.")
    p.add_argument("--scan-stride", type=int, default=1,
                   help="Scan every Nth native radar time step. Use 1 for full detection, 2 for faster scanning.")

    p.add_argument("--event-gap-min", type=float, default=90.0,
                   help="Candidate times separated by more than this start a new event.")
    p.add_argument("--max-events", type=int, default=6)
    p.add_argument("--instances-per-event", type=int, default=3)
    p.add_argument("--min-instance-gap-min", type=float, default=30.0)

    p.add_argument("--cml-time-label", choices=["left", "right"], default="left")
    p.add_argument("--cml-time-closed", choices=["left", "right"], default="left")
    p.add_argument("--radar-time-tolerance-min", type=float, default=30.0)

    p.add_argument("--bounds", nargs="+", type=float, default=DEFAULT_BOUNDS)
    p.add_argument("--cmap", default="turbo")
    p.add_argument("--draw-links", action="store_true", default=True)
    p.add_argument("--no-draw-links", dest="draw_links", action="store_false")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--out-dir", default="convective_event_instances_essd")

    return p.parse_args()


# ============================================================
# Main
# ============================================================
def main():
    args = parse_args()
    t0 = time.perf_counter()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading CML grid ...")
    ds_cml = xr.open_dataset(args.gridded_nc)
    if args.cml_var not in ds_cml.data_vars:
        raise KeyError(f"CML variable {args.cml_var!r} not found. Available: {list(ds_cml.data_vars)}")
    cml_da = ds_cml[args.cml_var]

    x_name = next((v for v in ["x", "lon", "longitude"] if v in ds_cml.coords), None)
    y_name = next((v for v in ["y", "lat", "latitude"] if v in ds_cml.coords), None)
    if x_name is None or y_name is None:
        raise KeyError("CML gridded file must contain x/y or lon/lat coordinates.")
    xg = np.asarray(ds_cml[x_name].values)
    yg = np.asarray(ds_cml[y_name].values)

    print("Loading CML geometry ...")
    with xr.open_dataset(args.cml_geom_nc) as ds_geom:
        lon_a = np.asarray(ds_geom[args.site_a_lon_var].values, dtype=float)
        lat_a = np.asarray(ds_geom[args.site_a_lat_var].values, dtype=float)
        lon_b = np.asarray(ds_geom[args.site_b_lon_var].values, dtype=float)
        lat_b = np.asarray(ds_geom[args.site_b_lat_var].values, dtype=float)
    geom = (lon_a, lat_a, lon_b, lat_b)

    if args.extent is None:
        base_extent = network_search_extent(lon_a, lat_a, lon_b, lat_b, pad=0.05)
    else:
        base_extent = list(args.extent)
    print(f"Joint search domain: {np.round(base_extent, 3).tolist()}")
    print(f"Output zoom size: {2*args.zoom_half_lon:g} deg lon x {2*args.zoom_half_lat:g} deg lat")
    link_density = link_density_on_grid(xg, yg, lon_a, lat_a, lon_b, lat_b)
    print(f"CML links represented in density grid: {int(link_density.sum())}")

    radar_files = discover_radar_files(args.radar_root, args.radar_pattern)
    date_filter = parse_date_filter(args.only_dates)
    if date_filter is not None:
        radar_files = [r for r in radar_files if r.date.date() in date_filter]
    if not radar_files:
        raise FileNotFoundError(f"No radar files found below {args.radar_root!r} with pattern {args.radar_pattern!r}")
    print(f"Radar files found: {len(radar_files)}")

    cache = RadarCache(RADAR_VAR_CANDIDATES, args.radar_scale_to_mmph)
    try:
        # Infer radar resolution from the first file with a valid time axis.
        radar_minutes = None
        skipped_for_resolution = 0
        for info in radar_files:
            try:
                _, _, _, times = cache.open(info.path)
            except RadarVariableNotFoundError:
                skipped_for_resolution += 1
                continue
            except Exception:
                skipped_for_resolution += 1
                continue
            if len(times) >= 2:
                radar_minutes = infer_step_minutes(times)
                break
        if skipped_for_resolution:
            print(f"  Skipped {skipped_for_resolution} invalid radar file(s) while detecting resolution.")
        if radar_minutes is None:
            radar_minutes = 5.0
        print(f"Native radar temporal resolution: {radar_minutes:g} min")

        cml_aligned = resample_cml_to_radar_resolution(
            cml_da,
            radar_minutes=radar_minutes,
            label=args.cml_time_label,
            closed=args.cml_time_closed,
        )
        cml_times = pd.DatetimeIndex(cml_aligned.time.values)

        # Major speed-up: do not scan radar dates outside the available CML period.
        cml_date_min = cml_times.min().date()
        cml_date_max = cml_times.max().date()
        before_overlap = len(radar_files)
        radar_files = [r for r in radar_files if cml_date_min <= r.date.date() <= cml_date_max]
        print(f"Radar files overlapping CML period: {len(radar_files)}/{before_overlap} "
              f"({cml_date_min} to {cml_date_max})")
        mapper = RadarToCMLGridMapper(xg, yg)

        print("Scanning for joint convective radar-CML instances ...")
        candidates = scan_candidates(
            radar_files=radar_files,
            cache=cache,
            cml_aligned=cml_aligned,
            cml_times=cml_times,
            xg=xg, yg=yg,
            link_density=link_density,
            mapper=mapper,
            search_extent=base_extent,
            convective_threshold=args.convective_threshold,
            cml_rain_threshold=args.cml_rain_threshold,
            wet_threshold=args.dry_threshold,
            min_convective_pixels=args.min_convective_pixels,
            min_cml_wet_pixels=args.min_cml_wet_pixels,
            min_joint_wet_pixels=args.min_joint_wet_pixels,
            min_links_in_window=args.min_links_in_window,
            min_cml_availability_fraction=args.min_cml_availability_fraction,
            min_peak_mmph=args.min_peak_mmph,
            zoom_half_lon=args.zoom_half_lon,
            zoom_half_lat=args.zoom_half_lat,
            scan_stride=args.scan_stride,
            time_tolerance_min=args.radar_time_tolerance_min,
            date_filter=date_filter,
        )
        print(f"Candidate convective times: {len(candidates)}")

        events = group_candidates_into_events(candidates, args.event_gap_min)
        selected = select_instances_from_events(
            events,
            max_events=args.max_events,
            instances_per_event=args.instances_per_event,
            min_instance_gap_minutes=args.min_instance_gap_min,
        )
        print(f"Detected events: {len(events)}")
        print(f"Selected instances: {len(selected)}")

        summary_rows = []
        for event_rank, instance_rank, candidate in tqdm(selected, desc="Rendering instances", unit="figure"):
            # Match CML to the actual radar timestamp after resampling.
            if len(cml_times) == 0:
                continue
            pos = cml_times.get_indexer([candidate.time], method="nearest")[0]
            if pos < 0:
                continue
            cml_time = cml_times[pos]
            mismatch_min = abs((cml_time - candidate.time).total_seconds()) / 60.0
            if mismatch_min > args.radar_time_tolerance_min:
                print(f"  Skip {candidate.time}: nearest CML time differs by {mismatch_min:.1f} min")
                continue

            cml_field = np.asarray(cml_aligned.isel(time=int(pos)).values, dtype="float32")
            rlon, rlat, radar_field, actual_radar_time = cache.field_at(candidate.radar_path, candidate.time)

            # The zoom was selected jointly from radar convection, actual CML rain,
            # CML data availability, and static link density during event scanning.
            plot_extent = list(candidate.extent)

            stem = (
                f"event{event_rank:02d}_instance{instance_rank:02d}_"
                f"{actual_radar_time.strftime('%Y%m%d_%H%M')}_"
                f"joint{candidate.joint_wet_pixels}_links{candidate.links_in_window}_max{candidate.radar_max:.1f}"
            )
            out_png = out_dir / f"{stem}.png"
            out_pdf = out_dir / f"{stem}.pdf"

            make_vertical_instance_figure(
                cml_field=cml_field,
                cml_x=xg,
                cml_y=yg,
                radar_field=radar_field,
                radar_lon=rlon,
                radar_lat=rlat,
                actual_time=actual_radar_time,
                extent=plot_extent,
                geom=geom,
                args=args,
                out_png=out_png,
                out_pdf=out_pdf,
            )

            summary_rows.append({
                "event_rank": event_rank,
                "instance_rank": instance_rank,
                "radar_time": actual_radar_time,
                "cml_time": cml_time,
                "time_difference_min": mismatch_min,
                "convective_threshold_mmph": args.convective_threshold,
                "convective_pixels": candidate.convective_pixels,
                "wet_pixels": candidate.wet_pixels,
                "cml_wet_pixels": candidate.cml_wet_pixels,
                "joint_wet_pixels": candidate.joint_wet_pixels,
                "cml_available_pixels": candidate.cml_available_pixels,
                "cml_available_fraction": candidate.cml_available_fraction,
                "links_in_window": candidate.links_in_window,
                "cml_max_mmph": candidate.cml_max,
                "cml_p99_mmph": candidate.cml_p99,
                "radar_max_mmph": candidate.radar_max,
                "radar_p99_mmph": candidate.radar_p99,
                "radar_mean_convective_mmph": candidate.radar_mean_convective,
                "score": candidate.score,
                "extent_lon_min": plot_extent[0],
                "extent_lon_max": plot_extent[1],
                "extent_lat_min": plot_extent[2],
                "extent_lat_max": plot_extent[3],
                "png": out_png.name,
                "pdf": out_pdf.name,
            })

            del cml_field, radar_field
            gc.collect()

        summary = pd.DataFrame(summary_rows)
        summary_path = out_dir / "selected_convective_instances.csv"
        summary.to_csv(summary_path, index=False)

        if len(summary):
            print("\nSelected convective instances")
            cols = [
                "event_rank", "instance_rank", "radar_time",
                "joint_wet_pixels", "cml_wet_pixels", "convective_pixels",
                "links_in_window", "radar_max_mmph", "cml_max_mmph", "png",
            ]
            print(summary[cols].to_string(index=False))
        else:
            print("\nNo figures were rendered. Lower --min-convective-pixels or check time coverage.")

        print(f"\nSaved outputs in: {out_dir}")
        print(f"Summary: {summary_path}")
        print(f"Total elapsed: {time.perf_counter() - t0:.1f} s")

    finally:
        cache.close()
        ds_cml.close()


if __name__ == "__main__":
    main()
