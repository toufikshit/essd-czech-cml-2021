#!/usr/bin/env python3
"""
02_blettner_aproct_dataset_figures_grid.py
==========================================

Step 2 of the April--October Blettner-paper refresh.

Inputs
------
1. Blettner/GRAF V2 rainfall product created by step 1.
2. Monthly availability CSV + period summary CSV created by
   availability_apr_oct_step1_compute.py.
3. Radar/gauge/temperature/humidity source files for observation counts.

Outputs
-------
figures/
  cml_characteristics.png/.pdf
  availability_apr_oct_legacy_style.png/.pdf
  availability_apr_oct_compact_summary.png/.pdf

summary/
  april_october_sensor_summary.csv
  april_october_sensor_summary.json

grid/
  blettner_cml_hourly_apr_oct_2021.nc

Scientific choices
------------------
* No rainfall retrieval is recomputed here.
* Directional rainfall in the V2 file is combined only for mapping:
      link rain rate = nanmean(valid directional sublinks)
* Hourly depth is computed from native 1-min rate:
      depth = sum(R_1min / 60)
  and is retained for a link only when --hourly-min-valid-fraction of the
  expected one-minute values exist.
* IDW uses CML link midpoints. Defaults match the established event workflow:
      resolution = 0.01 degree
      nnear      = 15
      power      = 2
      max_dist   = 0.3 degree
* The availability calculation itself is NOT redone here; this script plots the
  CSV from the already verified availability script.
"""

from __future__ import annotations

import argparse
import json
import math
import time
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, BoundaryNorm
from matplotlib.ticker import AutoMinorLocator
from scipy.spatial import cKDTree

try:
    from netCDF4 import Dataset, date2num
except Exception as exc:
    raise SystemExit(
        "netCDF4 is required. Install with: conda install -c conda-forge netcdf4\n"
        f"Import error: {exc}"
    )

try:
    from tqdm import tqdm
except Exception:
    def tqdm(x, **kwargs):
        return x


ROW_ORDER = ("CML", "Radar", "Gauge", "Temperature", "Humidity")
ROW_LABELS = (
    "CML",
    "Radar",
    "Rain gauge",
    "Air temperature",
    "Relative humidity",
)
NO_DATA_COLOR = "#eeeeee"


# =============================================================================
# Generic helpers
# =============================================================================

def ensure_file(path: str | Path, label: str) -> Path:
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"{label} not found: {p}")
    return p


def decode_strings(values) -> np.ndarray:
    a = np.asarray(values)
    if a.dtype.kind == "S":
        if a.ndim == 1:
            return np.asarray(
                [x.decode("utf-8", errors="ignore").replace("\x00", "").strip()
                 for x in a],
                dtype=object,
            )
        rows = []
        for row in a:
            rows.append(
                "".join(
                    x.decode("utf-8", errors="ignore")
                    if isinstance(x, (bytes, np.bytes_))
                    else str(x)
                    for x in np.ravel(row)
                ).replace("\x00", "").strip()
            )
        return np.asarray(rows, dtype=object)
    return np.asarray([str(x).replace("\x00", "").strip() for x in a.reshape(-1)], dtype=object)


def find_time_dim(obj) -> str:
    for d in obj.dims:
        if d == "time":
            return d
    for d in obj.dims:
        if "time" in d.lower():
            return d
    raise KeyError(f"No time dimension in {obj.dims}")


def count_finite_chunked(da: xr.DataArray, time_dim: str, t0, t1, chunk=4096) -> int:
    times = pd.DatetimeIndex(da[time_dim].values)
    i0 = int(np.searchsorted(times.values, np.datetime64(pd.Timestamp(t0)), "left"))
    i1 = int(np.searchsorted(times.values, np.datetime64(pd.Timestamp(t1)), "right"))
    total = 0
    for s in range(i0, i1, chunk):
        e = min(i1, s + chunk)
        a = np.asarray(da.isel({time_dim: slice(s, e)}).values)
        if np.ma.isMaskedArray(a):
            a = a.filled(np.nan)
        total += int(np.count_nonzero(np.isfinite(a)))
    return total


def get_entity_dim(da: xr.DataArray, time_dim: str):
    candidates = [d for d in da.dims if d != time_dim]
    if not candidates:
        return None
    # Prefer familiar station/link names.
    for cand in ("cml_id", "station", "WSI", "station_id"):
        if cand in candidates:
            return cand
    return candidates[-1]


def choose_main_numeric_var(ds: xr.Dataset, candidates):
    for name in candidates:
        if name in ds.data_vars:
            return name
    for name, da in ds.data_vars.items():
        if np.issubdtype(da.dtype, np.number) and any("time" in d.lower() for d in da.dims):
            return name
    raise KeyError(f"No suitable variable found. Available: {list(ds.data_vars)}")


# =============================================================================
# Publication rc
# =============================================================================

def setup_rc():
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 10,
        "axes.labelsize": 11,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "axes.linewidth": 0.9,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


# =============================================================================
# Availability figures: same old style + one alternative candidate
# =============================================================================

def availability_matrix(csv_path: Path):
    df = pd.read_csv(csv_path)
    df["month"] = pd.to_datetime(df["month"], format="%Y-%m")
    months = pd.DatetimeIndex(sorted(df["month"].unique()))
    M = np.full((len(ROW_ORDER), len(months)), np.nan, dtype=float)
    for i, sensor in enumerate(ROW_ORDER):
        sub = df[df["sensor"] == sensor].set_index("month")
        for j, m in enumerate(months):
            if m in sub.index:
                v = sub.loc[m, "availability_percent"]
                if isinstance(v, pd.Series):
                    v = v.iloc[0]
                M[i, j] = float(v)
    return df, months, M


def avail_text_color(v, cmap, norm):
    r, g, b, _ = cmap(norm(v))
    lum = 0.2126*r + 0.7152*g + 0.0722*b
    return "black" if lum > 0.52 else "white"


def draw_availability_heatmap(ax, M, months, cmap_name="cividis", decimals=0):
    cmap = mpl.colormaps[cmap_name].copy()
    cmap.set_bad(NO_DATA_COLOR)
    norm = Normalize(0, 100)
    mesh = ax.pcolormesh(
        np.arange(M.shape[1] + 1),
        np.arange(M.shape[0] + 1),
        np.ma.masked_invalid(M),
        shading="flat",
        cmap=cmap,
        norm=norm,
        edgecolors="white",
        linewidth=1.0,
        antialiased=True,
    )
    ax.set_xlim(0, M.shape[1])
    ax.set_ylim(M.shape[0], 0)
    ax.set_xticks(np.arange(M.shape[1]) + 0.5)
    ax.set_xticklabels([m.strftime("%b") for m in months])
    ax.set_yticks(np.arange(M.shape[0]) + 0.5)
    ax.set_yticklabels(ROW_LABELS)
    ax.tick_params(axis="x", length=4, width=0.9, pad=3, direction="out")
    ax.tick_params(axis="y", length=0, pad=7)

    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isfinite(v):
                ax.text(
                    j + 0.5, i + 0.5, f"{v:.{decimals}f}",
                    ha="center", va="center", fontsize=9,
                    color=avail_text_color(v, cmap, norm),
                )
    for sp in ax.spines.values():
        sp.set_color("0.25")
        sp.set_linewidth(0.9)
    return mesh


def plot_availability_variants(csv_path, summary_path, out_dir, dpi=600):
    setup_rc()
    _, months, M = availability_matrix(csv_path)
    summary = pd.read_csv(summary_path)

    # A: same publication heatmap style, Apr-Oct only.
    fig = plt.figure(figsize=(9.0, 4.2))
    ax = fig.add_axes([0.18, 0.30, 0.79, 0.63])
    cax = fig.add_axes([0.32, 0.105, 0.50, 0.052])
    mesh = draw_availability_heatmap(ax, M, months)
    ax.text(
        M.shape[1] / 2, -0.34, str(months[0].year),
        transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=10,
    )
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal",
                      ticks=[0, 20, 40, 60, 80, 100])
    cb.set_label("Network-wide monthly availability (%)")
    cb.ax.tick_params(direction="out", length=4, width=0.9)
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"availability_apr_oct_legacy_style.{ext}",
                    dpi=dpi if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)

    # B: alternative. Preserve heatmap language, use right side for sample context.
    fig = plt.figure(figsize=(10.3, 4.2))
    ax = fig.add_axes([0.16, 0.30, 0.58, 0.63])
    sax = fig.add_axes([0.765, 0.30, 0.21, 0.63])
    cax = fig.add_axes([0.25, 0.105, 0.47, 0.052])
    mesh = draw_availability_heatmap(ax, M, months)

    S = summary.set_index("sensor")
    sax.set_xlim(0, 1)
    sax.set_ylim(len(ROW_ORDER), 0)
    sax.axis("off")
    sax.text(0.03, -0.25, "Entities", ha="left", va="center",
             fontsize=9, fontweight="bold")
    sax.text(0.97, -0.25, "Apr–Oct\navailability", ha="right", va="center",
             fontsize=9, fontweight="bold", linespacing=0.95)
    for i, sensor in enumerate(ROW_ORDER):
        y = i + 0.5
        if sensor in S.index:
            row = S.loc[sensor]
            sax.text(0.03, y, f"{int(row['n_entities']):,}",
                     ha="left", va="center", fontsize=10)
            sax.text(0.97, y, f"{float(row['period_availability_percent']):.1f}%",
                     ha="right", va="center", fontsize=10)
        sax.axhline(i, color="0.85", lw=0.7)
    sax.axhline(len(ROW_ORDER), color="0.85", lw=0.7)

    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal",
                      ticks=[0, 20, 40, 60, 80, 100])
    cb.set_label("Network-wide monthly availability (%)")
    cb.ax.tick_params(direction="out", length=4, width=0.9)
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"availability_apr_oct_compact_summary.{ext}",
                    dpi=dpi if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


# =============================================================================
# CML characteristics: preserve old 3-panel language
# =============================================================================

def flatten_freq_length(ds: xr.Dataset):
    if "frequency_GHz" in ds:
        f = np.asarray(ds["frequency_GHz"].values, dtype=float)
    elif "frequency" in ds:
        f = np.asarray(ds["frequency"].values, dtype=float)
        med = np.nanmedian(np.abs(f))
        if med > 1000:
            f = f / 1000.0
    else:
        raise KeyError("No frequency_GHz or frequency variable")

    if "length_km" in ds:
        L = np.asarray(ds["length_km"].values, dtype=float)
    elif "length" in ds:
        L = np.asarray(ds["length"].values, dtype=float)
        if np.nanmedian(np.abs(L)) > 100:
            L = L / 1000.0
    else:
        raise KeyError("No length_km or length variable")

    # Broadcast link-level length to directional frequency if needed.
    if f.ndim == 2 and L.ndim == 1:
        if f.shape[-1] == L.size:
            L = np.broadcast_to(L[None, :], f.shape)
        elif f.shape[0] == L.size:
            L = np.broadcast_to(L[:, None], f.shape)
    elif L.ndim == 2 and f.ndim == 1:
        if L.shape[-1] == f.size:
            f = np.broadcast_to(f[None, :], L.shape)
        elif L.shape[0] == f.size:
            f = np.broadcast_to(f[:, None], L.shape)

    ff = f.reshape(-1)
    ll = L.reshape(-1)
    m = np.isfinite(ff) & np.isfinite(ll) & (ff > 0) & (ll > 0)
    return ff[m], ll[m]


def plot_characteristics(v2_nc: Path, out_dir: Path, dpi=600):
    setup_rc()
    with xr.open_dataset(v2_nc) as ds:
        f, L = flatten_freq_length(ds)

    fig, axes = plt.subplots(1, 3, figsize=(13.7, 4.0),
                             gridspec_kw={"width_ratios": [1, 1, 1.22]})
    fig.subplots_adjust(left=0.065, right=0.93, bottom=0.18, top=0.93, wspace=0.28)

    bar_color = "#ef5a3a"

    # Discrete frequency bands; half-GHz bins preserve narrow licensed bands.
    fbins = np.arange(0, max(90.0, np.nanmax(f) + 1), 1.0)
    axes[0].hist(f, bins=fbins, color=bar_color, edgecolor="white", linewidth=0.45)
    axes[0].set_xlabel("Frequency (GHz)")
    axes[0].set_ylabel("Number of sub-links")
    axes[0].set_xlim(0, 90)

    lbins = np.arange(0, max(30.0, np.nanmax(L) + 1.0), 1.0)
    axes[1].hist(L, bins=lbins, color=bar_color, edgecolor="white", linewidth=0.45)
    axes[1].set_xlabel("Path length (km)")
    axes[1].set_ylabel("Number of links")
    axes[1].set_xlim(0, 30)

    hb = axes[2].hexbin(
        f, L,
        gridsize=(30, 18),
        extent=(0, 90, 0, 30),
        mincnt=1,
        cmap="YlGnBu",
        linewidths=0.35,
        edgecolors="0.45",
    )
    axes[2].set_xlabel("Frequency (GHz)")
    axes[2].set_ylabel("Path length (km)")
    axes[2].set_xlim(0, 90)
    axes[2].set_ylim(0, 30)

    for label, ax in zip(("(a)", "(b)", "(c)"), axes):
        ax.text(0.01, 1.04, label, transform=ax.transAxes,
                ha="left", va="bottom", fontweight="bold", fontsize=11)
        ax.tick_params(direction="out", length=5, width=0.9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    cax = fig.add_axes([0.945, 0.23, 0.014, 0.62])
    cb = fig.colorbar(hb, cax=cax)
    cb.set_label("Number of sub-links", labelpad=7)
    cb.ax.tick_params(direction="out", length=3, width=0.8, labelsize=8.5)

    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"cml_characteristics.{ext}",
                    dpi=dpi if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


# =============================================================================
# Sensor / period summary
# =============================================================================

def summarize_source(name, path, var_candidates, start, end):
    p = ensure_file(path, name)
    with xr.open_dataset(p, decode_cf=True) as ds:
        var = choose_main_numeric_var(ds, var_candidates)
        da = ds[var]
        tdim = find_time_dim(da)
        times = pd.DatetimeIndex(da[tdim].values)
        mask = (times >= pd.Timestamp(start)) & (times <= pd.Timestamp(end))
        ntime = int(mask.sum())
        edim = get_entity_dim(da, tdim)
        n_entities = int(ds.sizes[edim]) if edim else 1
        finite = count_finite_chunked(da, tdim, start, end)
        possible = int(ntime * n_entities)
        # If the variable has extra dimensions (e.g. sublinks), possible should
        # reflect the full variable shape over the selected time interval.
        mult = 1
        for d in da.dims:
            if d != tdim and d != edim:
                mult *= int(ds.sizes[d])
        possible *= mult
        return {
            "sensor": name,
            "file": str(p),
            "variable": var,
            "n_entities": n_entities,
            "n_timestamps": ntime,
            "finite_observations": finite,
            "possible_observations": possible,
            "availability_percent": 100.0 * finite / possible if possible else np.nan,
        }


def summarize_v2(v2_nc, start, end):
    p = ensure_file(v2_nc, "Blettner V2 rainfall")
    with xr.open_dataset(p) as ds:
        if "R_cml" not in ds:
            raise KeyError(f"R_cml missing from {p}")
        da = ds["R_cml"]
        times = pd.DatetimeIndex(ds["time"].values)
        i0 = int(np.searchsorted(times.values, np.datetime64(pd.Timestamp(start)), "left"))
        i1 = int(np.searchsorted(times.values, np.datetime64(pd.Timestamp(end)), "right"))
        ntime = i1 - i0
        nlinks = int(ds.sizes["cml_id"])
        nsublinks = int(ds.sizes["sublink_id"])

        if "qc_final_sublink_pass" in ds:
            qsp = np.asarray(ds["qc_final_sublink_pass"].values)
            passing_sublinks = int(np.count_nonzero(qsp == 1))
        else:
            passing_sublinks = nlinks * nsublinks

        finite = 0
        for s in range(i0, i1, 10080):
            e = min(i1, s + 10080)
            a = np.asarray(da.isel(time=slice(s, e)).values)
            finite += int(np.count_nonzero(np.isfinite(a)))

        return {
            "sensor": "Blettner CML rainfall",
            "file": str(p),
            "variable": "R_cml",
            "n_entities": nlinks,
            "n_sublink_slots": nlinks * nsublinks,
            "n_passing_sublinks": passing_sublinks,
            "n_timestamps": ntime,
            "finite_observations": finite,
            "possible_observations": ntime * nlinks * nsublinks,
            "availability_percent": (
                100.0 * finite / (ntime * nlinks * nsublinks)
                if ntime and nlinks and nsublinks else np.nan
            ),
        }


def write_summary(args, out_dir):
    rows = [
        summarize_v2(args.v2_nc, args.start, args.end),
        summarize_source("Radar path", args.radar_nc, ("rr", "R_radar", "rainfall_rate"),
                         args.start, args.end),
        summarize_source("Rain gauge", args.gauge_nc, ("rr", "R", "precipitation"),
                         args.start, args.end),
        summarize_source("Air temperature", args.temperature_nc, ("T", "temperature"),
                         args.start, args.end),
        summarize_source("Relative humidity", args.humidity_nc, ("H", "humidity"),
                         args.start, args.end),
    ]
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "april_october_sensor_summary.csv", index=False)
    (out_dir / "april_october_sensor_summary.json").write_text(
        json.dumps(rows, indent=2, default=str), encoding="utf-8"
    )
    print("\nAPRIL-OCTOBER DATA SUMMARY")
    print(df.to_string(index=False))
    return df


# =============================================================================
# Hourly gridding
# =============================================================================

class HourlyGridWriter:
    def __init__(self, path, times, x, y, resolution, nnear, power, max_dist,
                 complevel=2):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if self.path.exists():
            self.path.unlink()
        self.nc = Dataset(self.path, "w", format="NETCDF4")
        self.nc.createDimension("time", len(times))
        self.nc.createDimension("y", len(y))
        self.nc.createDimension("x", len(x))

        tv = self.nc.createVariable("time", "f8", ("time",))
        tv.units = "seconds since 1970-01-01 00:00:00 UTC"
        tv.calendar = "standard"
        tv[:] = np.asarray(pd.DatetimeIndex(times).view("int64") / 1e9, dtype=float)

        xv = self.nc.createVariable("x", "f8", ("x",))
        yv = self.nc.createVariable("y", "f8", ("y",))
        xv[:] = x
        yv[:] = y
        xv.units = "degrees_east"
        yv.units = "degrees_north"

        self.R = self.nc.createVariable(
            "R_1h",
            "f4",
            ("time", "y", "x"),
            zlib=True,
            complevel=int(complevel),
            shuffle=True,
            chunksizes=(1, min(160, len(y)), min(240, len(x))),
            fill_value=np.float32(np.nan),
        )
        self.R.units = "mm"
        self.R.long_name = "1-hour CML precipitation accumulation"
        self.R.coordinates = "time y x"

        self.nvalid = self.nc.createVariable(
            "n_valid_links",
            "i2",
            ("time", "y", "x"),
            zlib=True,
            complevel=1,
            shuffle=True,
            chunksizes=(1, min(160, len(y)), min(240, len(x))),
            fill_value=np.int16(-1),
        )
        self.nvalid.long_name = "Number of finite neighboring CML links contributing to IDW"

        self.nc.setncattr("period", "2021-04-01 through 2021-10-31")
        self.nc.setncattr("mapping_input", "mean of valid directional sublink rainfall")
        self.nc.setncattr(
            "hourly_aggregation",
            "sum of one-minute rain rate divided by 60; strict per-link valid-fraction criterion",
        )
        self.nc.setncattr(
            "interpolation",
            f"IDW midpoint nnear={nnear}, power={power}, max_dist={max_dist} degrees, resolution={resolution} degrees",
        )

    def write(self, i, field, nvalid):
        self.R[i, :, :] = field.astype(np.float32)
        self.nvalid[i, :, :] = nvalid.astype(np.int16)

    def close(self):
        self.nc.close()


def build_static_idw(mid_lon, mid_lat, xgrid, ygrid, nnear, max_dist):
    X, Y = np.meshgrid(xgrid, ygrid)
    pts = np.column_stack([mid_lon, mid_lat])
    tree = cKDTree(pts)
    q = np.column_stack([X.ravel(), Y.ravel()])
    dist, idx = tree.query(q, k=min(int(nnear), len(pts)), workers=-1)

    if idx.ndim == 1:
        idx = idx[:, None]
        dist = dist[:, None]

    within = dist <= float(max_dist)
    return X, Y, idx.astype(np.int32), dist.astype(np.float32), within


def idw_field(values, idx, dist, within, power):
    v = np.asarray(values, dtype=np.float64)
    z = v[idx]
    valid = within & np.isfinite(z)

    # Zero distance: exact link midpoint gets exact value.
    zero = valid & (dist <= 1e-12)
    has_zero = np.any(zero, axis=1)

    w = np.zeros_like(dist, dtype=np.float64)
    nonzero = valid & ~zero
    w[nonzero] = 1.0 / np.power(dist[nonzero], float(power))

    num = np.sum(w * np.where(valid, z, 0.0), axis=1)
    den = np.sum(w, axis=1)
    out = np.full(len(idx), np.nan, dtype=np.float32)
    good = den > 0
    out[good] = (num[good] / den[good]).astype(np.float32)

    if np.any(has_zero):
        # Mean in the very unlikely case multiple midpoints are exactly identical.
        zz = np.where(zero, z, np.nan)
        out[has_zero] = np.nanmean(zz[has_zero], axis=1).astype(np.float32)

    nvalid = np.sum(valid, axis=1).astype(np.int16)
    return out, nvalid


def build_hourly_grid(args, grid_path):
    v2 = ensure_file(args.v2_nc, "Blettner V2 rainfall")
    with xr.open_dataset(v2) as ds:
        if "R_cml" not in ds:
            raise KeyError(f"{v2} does not contain R_cml")

        required = [
            "site_a_longitude", "site_a_latitude",
            "site_b_longitude", "site_b_latitude",
        ]
        missing = [v for v in required if v not in ds]
        if missing:
            raise KeyError(f"V2 file lacks geometry variables: {missing}")

        lon0 = np.asarray(ds["site_a_longitude"].values, dtype=float)
        lat0 = np.asarray(ds["site_a_latitude"].values, dtype=float)
        lon1 = np.asarray(ds["site_b_longitude"].values, dtype=float)
        lat1 = np.asarray(ds["site_b_latitude"].values, dtype=float)
        mid_lon = 0.5 * (lon0 + lon1)
        mid_lat = 0.5 * (lat0 + lat1)

        valid_geom = np.isfinite(mid_lon) & np.isfinite(mid_lat)
        if not np.any(valid_geom):
            raise RuntimeError("No finite CML geometry")

        pad = float(args.grid_pad)
        xgrid = np.arange(
            np.nanmin(mid_lon[valid_geom]) - pad,
            np.nanmax(mid_lon[valid_geom]) + pad + args.resolution * 0.5,
            float(args.resolution),
        )
        ygrid = np.arange(
            np.nanmin(mid_lat[valid_geom]) - pad,
            np.nanmax(mid_lat[valid_geom]) + pad + args.resolution * 0.5,
            float(args.resolution),
        )

        geom_idx = np.flatnonzero(valid_geom)
        X, Y, nbr_idx_local, nbr_dist, within = build_static_idw(
            mid_lon[valid_geom], mid_lat[valid_geom],
            xgrid, ygrid,
            args.nnear, args.max_dist,
        )
        nbr_idx_global = geom_idx[nbr_idx_local]

        times = pd.DatetimeIndex(ds["time"].values)
        start = pd.Timestamp(args.start)
        end = pd.Timestamp(args.end)
        hour_times = pd.date_range(start.floor("h"), end.floor("h"), freq="1h")

        writer = HourlyGridWriter(
            grid_path, hour_times, xgrid, ygrid,
            args.resolution, args.nnear, args.idw_power, args.max_dist,
            complevel=args.grid_complevel,
        )

        da = ds["R_cml"]
        expected = 60
        required_n = int(math.ceil(expected * args.hourly_min_valid_fraction))

        try:
            for ih, h0 in enumerate(tqdm(hour_times, desc="Hourly CML grid", unit="h")):
                h1 = h0 + pd.Timedelta(hours=1)
                i0 = int(np.searchsorted(times.values, np.datetime64(h0), "left"))
                i1 = int(np.searchsorted(times.values, np.datetime64(h1), "left"))

                if i1 <= i0:
                    writer.write(
                        ih,
                        np.full(X.shape, np.nan, dtype=np.float32),
                        np.zeros(X.shape, dtype=np.int16),
                    )
                    continue

                a = np.asarray(
                    da.isel(time=slice(i0, i1)).transpose("time", "sublink_id", "cml_id").values,
                    dtype=np.float32,
                )
                # Link mean after directional retrieval, only for mapping.
                finite_s = np.isfinite(a)
                count_s = np.sum(finite_s, axis=1)
                sum_s = np.sum(np.where(finite_s, a, 0.0), axis=1, dtype=np.float64)
                link_rate = np.full((a.shape[0], a.shape[2]), np.nan, dtype=np.float32)
                ok_s = count_s > 0
                link_rate[ok_s] = (sum_s[ok_s] / count_s[ok_s]).astype(np.float32)

                finite_t = np.isfinite(link_rate)
                count_t = np.sum(finite_t, axis=0)
                # 1-min rate [mm/h] -> amount [mm] for each minute.
                amount = np.sum(
                    np.where(finite_t, link_rate, 0.0),
                    axis=0,
                    dtype=np.float64,
                ) / 60.0
                amount[count_t < required_n] = np.nan

                flat, nvalid = idw_field(
                    amount,
                    nbr_idx_global,
                    nbr_dist,
                    within,
                    args.idw_power,
                )
                writer.write(ih, flat.reshape(X.shape), nvalid.reshape(X.shape))
        finally:
            writer.close()

    print(f"[SAVED] {grid_path}")



# =============================================================================
# Fixed 21--22 June event: native 5-min CML grid for exact paper-style map
# =============================================================================

class EventGridWriter:
    def __init__(self, path, times, x, y, resolution, nnear, power, max_dist,
                 complevel=2):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if self.path.exists():
            self.path.unlink()
        self.nc = Dataset(self.path, "w", format="NETCDF4")
        self.nc.createDimension("time", len(times))
        self.nc.createDimension("y", len(y))
        self.nc.createDimension("x", len(x))

        tv = self.nc.createVariable("time", "f8", ("time",))
        tv.units = "seconds since 1970-01-01 00:00:00 UTC"
        tv.calendar = "standard"
        tv[:] = np.asarray(pd.DatetimeIndex(times).view("int64") / 1e9, dtype=float)

        xv = self.nc.createVariable("x", "f8", ("x",))
        yv = self.nc.createVariable("y", "f8", ("y",))
        xv[:] = x
        yv[:] = y
        xv.units = "degrees_east"
        yv.units = "degrees_north"

        self.R = self.nc.createVariable(
            "R",
            "f4",
            ("time", "y", "x"),
            zlib=True,
            complevel=int(complevel),
            shuffle=True,
            chunksizes=(1, min(160, len(y)), min(240, len(x))),
            fill_value=np.float32(np.nan),
        )
        self.R.units = "mm h-1"
        self.R.long_name = "5-minute CML rain rate gridded from Blettner/GRAF retrieval"
        self.nc.setncattr("timestamp_label", "start")
        self.nc.setncattr(
            "temporal_aggregation",
            "Directional 1-min rates -> 5-min accumulation -> equivalent rate; "
            "R5 = sum(R1min/60)*12, min_count=1. Directional sublinks are then "
            "averaged only for spatial mapping.",
        )
        self.nc.setncattr(
            "interpolation",
            f"IDW midpoint nnear={nnear}, power={power}, max_dist={max_dist} degrees, "
            f"resolution={resolution} degrees",
        )

    def write(self, i, field):
        self.R[i, :, :] = field.astype(np.float32)

    def close(self):
        self.nc.close()


def build_event_5min_grid(args, event_grid_path):
    """Create the 5-min CML grid needed by the exact multitime event plot."""
    v2 = ensure_file(args.v2_nc, "Blettner V2 rainfall")
    with xr.open_dataset(v2) as ds:
        da = ds["R_cml"].transpose("time", "sublink_id", "cml_id")
        times = pd.DatetimeIndex(ds["time"].values)

        lon0 = np.asarray(ds["site_a_longitude"].values, dtype=float)
        lat0 = np.asarray(ds["site_a_latitude"].values, dtype=float)
        lon1 = np.asarray(ds["site_b_longitude"].values, dtype=float)
        lat1 = np.asarray(ds["site_b_latitude"].values, dtype=float)
        mid_lon = 0.5 * (lon0 + lon1)
        mid_lat = 0.5 * (lat0 + lat1)
        valid_geom = np.isfinite(mid_lon) & np.isfinite(mid_lat)
        geom_idx = np.flatnonzero(valid_geom)

        pad = float(args.grid_pad)
        xgrid = np.arange(
            np.nanmin(mid_lon[valid_geom]) - pad,
            np.nanmax(mid_lon[valid_geom]) + pad + args.resolution * 0.5,
            float(args.resolution),
        )
        ygrid = np.arange(
            np.nanmin(mid_lat[valid_geom]) - pad,
            np.nanmax(mid_lat[valid_geom]) + pad + args.resolution * 0.5,
            float(args.resolution),
        )
        X, Y, nbr_idx_local, nbr_dist, within = build_static_idw(
            mid_lon[valid_geom], mid_lat[valid_geom],
            xgrid, ygrid, args.nnear, args.max_dist,
        )
        nbr_idx_global = geom_idx[nbr_idx_local]

        event_times = pd.date_range(
            pd.Timestamp(args.event_grid_start),
            pd.Timestamp(args.event_grid_end),
            freq="5min",
        )
        writer = EventGridWriter(
            event_grid_path, event_times, xgrid, ygrid,
            args.resolution, args.nnear, args.idw_power, args.max_dist,
            complevel=args.grid_complevel,
        )

        try:
            for it, t0 in enumerate(tqdm(event_times, desc="Event 5-min CML grid", unit="5min")):
                t1 = t0 + pd.Timedelta(minutes=5)
                i0 = int(np.searchsorted(times.values, np.datetime64(t0), "left"))
                i1 = int(np.searchsorted(times.values, np.datetime64(t1), "left"))
                if i1 <= i0:
                    writer.write(it, np.full(X.shape, np.nan, dtype=np.float32))
                    continue

                a = np.asarray(
                    da.isel(time=slice(i0, i1)).values,
                    dtype=np.float32,
                )  # time, sublink, link

                finite = np.isfinite(a)
                count = np.sum(finite, axis=0)  # sublink, link
                # Exact project convention for 1-min rate -> 5-min equivalent rate:
                # sum(R/60) * 12 = sum(R)/5. Missing minutes contribute no amount,
                # while an all-missing interval remains missing (min_count=1).
                r5_sub = np.sum(
                    np.where(finite, a, 0.0),
                    axis=0,
                    dtype=np.float64,
                ) / 5.0
                r5_sub[count == 0] = np.nan

                finite_s = np.isfinite(r5_sub)
                n_s = np.sum(finite_s, axis=0)
                link_rate = np.full(r5_sub.shape[1], np.nan, dtype=np.float32)
                good = n_s > 0
                if np.any(good):
                    link_rate[good] = (
                        np.sum(np.where(finite_s, r5_sub, 0.0), axis=0)[good]
                        / n_s[good]
                    ).astype(np.float32)

                flat, _ = idw_field(
                    link_rate,
                    nbr_idx_global,
                    nbr_dist,
                    within,
                    args.idw_power,
                )
                writer.write(it, flat.reshape(X.shape))
        finally:
            writer.close()

    print(f"[SAVED] {event_grid_path}")

# =============================================================================
# CLI / main
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Apr-Oct dataset figures, summary and hourly Blettner CML gridding."
    )
    base = "/home/jovyan/ESSD_publication_dataset"

    p.add_argument(
        "--v2-nc",
        default=f"{base}/blettner_apr_oct/BLETTNER_GRAF_RAIN_APR_OCT_2021.nc",
    )
    p.add_argument(
        "--availability-csv",
        default=f"{base}/availability_apr_oct/availability_monthly_essd_apr_oct.csv",
    )
    p.add_argument(
        "--availability-summary",
        default=f"{base}/availability_apr_oct/availability_monthly_essd_apr_oct_period_summary.csv",
    )
    p.add_argument(
        "--availability-compute-script",
        default=f"{base}/availability_apr_oct_step1_compute.py",
        help="Sibling script already created for the verified availability calculation.",
    )
    p.add_argument(
        "--availability-cml-nc",
        default=f"{base}/CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--radar-nc",
        default=f"{base}/CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--gauge-nc",
        default=f"{base}/CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--temperature-nc",
        default="/home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_TEMPERATURE_2021.nc",
    )
    p.add_argument(
        "--humidity-nc",
        default="/home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_HUMIDITY_2021.nc",
    )
    p.add_argument(
        "--out-dir",
        default=f"{base}/blettner_apr_oct/period_products",
    )

    p.add_argument("--start", default="2021-04-01 00:00:00")
    p.add_argument("--end", default="2021-10-31 23:59:59")

    p.add_argument("--resolution", type=float, default=0.01)
    p.add_argument("--grid-pad", type=float, default=0.05)
    p.add_argument("--nnear", type=int, default=15)
    p.add_argument("--idw-power", type=float, default=2.0)
    p.add_argument("--max-dist", type=float, default=0.3)
    p.add_argument("--hourly-min-valid-fraction", type=float, default=1.0)
    p.add_argument("--grid-complevel", type=int, default=2)

    # Exact transboundary-paper event. The earliest hourly panel is 19:45,
    # whose 1-h rolling window needs 5-min fields back to 18:50.
    p.add_argument("--event-grid-start", default="2021-06-21 18:50:00")
    p.add_argument("--event-grid-end", default="2021-06-22 02:45:00")

    p.add_argument("--skip-grid", action="store_true")
    p.add_argument("--skip-event-grid", action="store_true")
    p.add_argument("--skip-availability-plots", action="store_true")
    p.add_argument("--skip-characteristics", action="store_true")
    p.add_argument("--dpi", type=int, default=600)

    return p.parse_args()



def ensure_availability_inputs(args):
    csv_path = Path(args.availability_csv).expanduser().resolve()
    summary_path = Path(args.availability_summary).expanduser().resolve()
    if csv_path.is_file() and summary_path.is_file():
        return csv_path, summary_path

    script = Path(args.availability_compute_script).expanduser().resolve()
    if not script.is_file():
        raise FileNotFoundError(
            "Availability CSVs are missing and the computation script was not found:\n"
            f"  {script}\n"
            "Place availability_apr_oct_step1_compute.py in "
            "/home/jovyan/ESSD_publication_dataset or run it manually first."
        )

    out_dir = csv_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-u", str(script),
        "--cml-file", str(Path(args.availability_cml_nc).expanduser().resolve()),
        "--radar-file", str(Path(args.radar_nc).expanduser().resolve()),
        "--gauge-file", str(Path(args.gauge_nc).expanduser().resolve()),
        "--temperature-file", str(Path(args.temperature_nc).expanduser().resolve()),
        "--humidity-file", str(Path(args.humidity_nc).expanduser().resolve()),
        "--start", "2021-04-01",
        "--end", "2021-10-31",
        "--denominator", "expected",
        "--out-dir", str(out_dir),
        "--basename", "availability_monthly_essd_apr_oct",
    ]
    print("\nAvailability CSVs are missing; running the verified calculation:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)

    if not csv_path.is_file() or not summary_path.is_file():
        raise RuntimeError(
            "Availability computation finished but expected outputs are still missing."
        )
    return csv_path, summary_path



def main():
    args = parse_args()
    started = time.perf_counter()

    out = Path(args.out_dir).expanduser().resolve()
    fig_dir = out / "figures"
    summary_dir = out / "summary"
    grid_dir = out / "grid"
    for d in (fig_dir, summary_dir, grid_dir):
        d.mkdir(parents=True, exist_ok=True)

    ensure_file(args.v2_nc, "Blettner V2 rainfall")

    availability_csv = Path(args.availability_csv).expanduser().resolve()
    availability_summary = Path(args.availability_summary).expanduser().resolve()
    if not args.skip_availability_plots:
        availability_csv, availability_summary = ensure_availability_inputs(args)
        plot_availability_variants(
            availability_csv,
            availability_summary,
            fig_dir,
            dpi=args.dpi,
        )

    if not args.skip_characteristics:
        plot_characteristics(Path(args.v2_nc), fig_dir, dpi=args.dpi)

    write_summary(args, summary_dir)

    grid_path = grid_dir / "blettner_cml_hourly_apr_oct_2021.nc"
    if not args.skip_grid:
        build_hourly_grid(args, grid_path)

    event_grid_path = grid_dir / "blettner_cml_event_20210621_5min.nc"
    if not args.skip_event_grid:
        build_event_5min_grid(args, event_grid_path)

    print("\nSAVED FIGURES")
    for p in sorted(fig_dir.glob("*")):
        print(" ", p)
    print("\nElapsed: %.1f s" % (time.perf_counter() - started))


if __name__ == "__main__":
    main()
