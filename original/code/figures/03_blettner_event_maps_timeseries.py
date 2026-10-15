#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create publication-style multitime event figures that match the provided example:
- columns = successive timestamps
- rows    = CML and radar
- no in-map text
- only top time stamps, left row labels, and one shared bottom colorbar

Outputs:
1) native 5-min rainfall-rate panel figure
2) rolling 1-hour accumulation panel figure
3) optional single-timestep two-panel figures using the same corrected colorbar logic
"""

from __future__ import annotations

import argparse
import gc
import math
import time
from pathlib import Path
from typing import Sequence
from types import SimpleNamespace

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap

import cartopy.crs as ccrs
import cartopy.feature as cfeature


DEFAULT_NATIVE_BOUNDS = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 30, 50]
DEFAULT_HOURLY_BOUNDS = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 30, 50]


def essd_colormap():
    # Tuned to resemble the user's example: pale green -> cyan -> blue -> magenta.
    colors = [
        "#e2edc5", "#d3e8c1", "#c0e2be", "#a6dcbc", "#84d6c2",
        "#5ecfd0", "#40bfd8", "#319fd2", "#2d7dc0", "#285faf",
        "#3345a4", "#592ea9", "#8a23b7", "#d600d9",
    ]
    return LinearSegmentedColormap.from_list("essd_example_like", colors, N=256)



def discrete_essd_cmap(bounds: np.ndarray) -> ListedColormap:
    """
    Create exactly one color for each rainfall interval.

    BoundaryNorm is intentionally used without extension bins. Values above
    the highest boundary use the colormap's over color; values below the first
    boundary are masked before plotting.
    """
    bounds = np.asarray(bounds, dtype=float)
    n_intervals = bounds.size - 1
    if n_intervals < 1:
        raise ValueError("At least two boundaries are required.")

    base = essd_colormap()
    colors = base(np.linspace(0.0, 1.0, n_intervals))
    cmap = ListedColormap(colors, name=f"essd_discrete_{n_intervals}")
    cmap.set_bad("white")
    cmap.set_over(colors[-1])
    return cmap


def format_tick_labels(values: Sequence[float]) -> list[str]:
    labels = []
    for value in values:
        value = float(value)
        if value < 1:
            labels.append(f"{value:.1f}")
        elif np.isclose(value, round(value)):
            labels.append(str(int(round(value))))
        else:
            labels.append(f"{value:g}")
    return labels


def compact_major_ticks(bounds: np.ndarray) -> np.ndarray:
    """
    Keep compact two-panel colorbars readable while retaining all class
    boundaries as minor ticks.
    """
    preferred = [bounds[0], 1, 2, 5, 10, 15, 20, 30, bounds[-1]]
    ticks = []
    for value in preferred:
        hit = bounds[np.isclose(bounds, value)]
        if hit.size:
            value = float(hit[0])
            if not ticks or not np.isclose(ticks[-1], value):
                ticks.append(value)
    return np.asarray(ticks, dtype=float)


def add_discrete_horizontal_colorbar(
    *,
    fig,
    mappable,
    cax,
    bounds: np.ndarray,
    label: str,
    labelsize: float,
    ticksize: float,
    label_all_bounds: bool,
):
    """
    Add a colorbar that exactly matches the discrete rainfall classes.

    - one color box per interval
    - all boundaries shown
    - upper extension only
    - no lower extension because values below the first boundary are masked
    """
    bounds = np.asarray(bounds, dtype=float)
    major_ticks = bounds if label_all_bounds else compact_major_ticks(bounds)

    cb = fig.colorbar(
        mappable,
        cax=cax,
        orientation="horizontal",
        boundaries=bounds,
        ticks=major_ticks,
        spacing="uniform",
        drawedges=True,
        extend="max",
        extendfrac=0.035,
    )

    cb.set_label(label, fontsize=labelsize, labelpad=7)
    cb.set_ticklabels(format_tick_labels(major_ticks))
    cb.ax.tick_params(
        axis="x",
        which="major",
        labelsize=ticksize,
        width=0.9,
        length=5,
        pad=3,
    )

    # Preserve every rainfall boundary as an unlabeled minor tick when the
    # compact figure cannot fit every label.
    cb.ax.xaxis.set_ticks(bounds, minor=True)
    cb.ax.tick_params(
        axis="x",
        which="minor",
        width=0.55,
        length=2.5,
    )

    cb.outline.set_linewidth(0.8)
    try:
        cb.dividers.set_linewidth(0.45)
        cb.dividers.set_color("0.62")
    except Exception:
        pass

    return cb



RADAR_VAR_CANDIDATES = ["YW", "RR", "rainfall_rate", "R", "rainrate", "precip"]


def _pick_radar_var(ds, candidates):
    for name in candidates:
        if name in ds.data_vars:
            return name
    raise KeyError(
        f"No radar rainfall variable found. Tried {candidates}; "
        f"available={list(ds.data_vars)}"
    )


def _find_time_dim(da):
    return next((d for d in da.dims if "time" in d.lower()), None)


def infer_step_minutes(times: pd.DatetimeIndex, fallback: float = 5.0) -> float:
    if len(times) < 2:
        return float(fallback)
    delta_s = np.diff(times.view("int64")) / 1e9
    delta_s = delta_s[np.isfinite(delta_s) & (delta_s > 0)]
    if delta_s.size == 0:
        return float(fallback)
    return float(np.median(delta_s) / 60.0)


class RadarCache:
    """Minimal self-contained copy of the verified project radar reader."""
    def __init__(self, var_candidates, scale_to_mmph):
        self.var_candidates = list(var_candidates)
        self.scale = float(scale_to_mmph)
        self._ds = None
        self._path = None
        self._var = None
        self._tdim = None
        self._times = None

    def close(self):
        if self._ds is not None:
            try:
                self._ds.close()
            except Exception:
                pass
        self._ds = self._path = self._var = self._tdim = self._times = None

    def open(self, path):
        path = str(path)
        if self._path != path:
            self.close()
            ds = xr.open_dataset(path)
            var = _pick_radar_var(ds, self.var_candidates)
            da = ds[var]
            tdim = _find_time_dim(da)
            times = pd.DatetimeIndex(da[tdim].values) if tdim else pd.DatetimeIndex([])
            self._ds = ds
            self._path = path
            self._var = var
            self._tdim = tdim
            self._times = times
        return self._ds, self._ds[self._var], self._tdim, self._times

    def field_at(self, path, target_time):
        ds, da, tdim, _ = self.open(path)
        if tdim is not None:
            selected = da.sel({tdim: pd.Timestamp(target_time)}, method="nearest")
            actual = pd.Timestamp(selected[tdim].values)
        else:
            selected = da
            actual = pd.Timestamp(target_time)

        field = np.asarray(selected.values, dtype="float32") * np.float32(self.scale)
        field = np.where(np.isfinite(field), field, np.nan).astype("float32", copy=False)

        latname = next((c for c in ["lat", "latitude"] if c in ds.coords or c in ds.variables), None)
        lonname = next((c for c in ["lon", "longitude"] if c in ds.coords or c in ds.variables), None)
        if latname and lonname:
            return (
                np.asarray(ds[lonname].values),
                np.asarray(ds[latname].values),
                field,
                actual,
            )

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

def load_base_module(base_path: str | None):
    # Kept as a compatibility shim so the plotting logic below is unchanged.
    # This script is now self-contained.
    mod = SimpleNamespace(
        RADAR_VAR_CANDIDATES=RADAR_VAR_CANDIDATES,
        RadarCache=RadarCache,
        infer_step_minutes=infer_step_minutes,
    )
    return mod, Path(__file__).resolve()


def parse_bounds(vals: Sequence[float]) -> np.ndarray:
    arr = np.asarray(vals, dtype=float)
    if arr.ndim != 1 or arr.size < 2:
        raise ValueError("Need at least two bounds.")
    if not np.all(np.diff(arr) > 0):
        raise ValueError("Bounds must be strictly increasing.")
    return arr


def dry_mask(a: np.ndarray, threshold: float) -> np.ndarray:
    out = np.asarray(a, dtype=float).copy()
    out[out < threshold] = np.nan
    return out


def setup_axis(ax, extent, proj):
    ax.set_extent(extent, crs=proj)
    ax.add_feature(cfeature.LAND, facecolor="0.98", zorder=0)
    ax.add_feature(cfeature.BORDERS, linewidth=0.6, edgecolor="0.76", zorder=1)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.45, edgecolor="0.78", zorder=1)
    try:
        ax.spines["geo"].set_linewidth(0.6)
        ax.spines["geo"].set_edgecolor("0.68")
    except Exception:
        try:
            ax.patch.set_linewidth(0.6)
            ax.patch.set_edgecolor("0.68")
        except Exception:
            pass
    try:
        ax.set_xticks([], crs=proj)
        ax.set_yticks([], crs=proj)
    except Exception:
        pass
    ax.gridlines(draw_labels=False, linewidth=0.0, alpha=0.0)


def plot_links(ax, proj, lon_a, lat_a, lon_b, lat_b, extent):
    keep = (
        np.isfinite(lon_a) & np.isfinite(lat_a) & np.isfinite(lon_b) & np.isfinite(lat_b) &
        (
            ((lon_a >= extent[0]) & (lon_a <= extent[1]) & (lat_a >= extent[2]) & (lat_a <= extent[3])) |
            ((lon_b >= extent[0]) & (lon_b <= extent[1]) & (lat_b >= extent[2]) & (lat_b <= extent[3]))
        )
    )
    idx = np.flatnonzero(keep)
    for i in idx:
        ax.plot(
            [lon_a[i], lon_b[i]],
            [lat_a[i], lat_b[i]],
            transform=proj,
            color="0.35",
            linewidth=0.22,
            alpha=0.18,
            zorder=3,
        )


def nearest_time_index(times: pd.DatetimeIndex, t: pd.Timestamp, tolerance: pd.Timedelta | None):
    pos = times.get_indexer([t], method="nearest")[0]
    if pos < 0:
        return -1
    if tolerance is not None:
        dt = abs(times[pos] - t)
        if dt > tolerance:
            return -1
    return int(pos)


def native_time_series(center_time: pd.Timestamp, ncols: int, step_minutes: int, explicit_start: str | None):
    if explicit_start is not None:
        start = pd.Timestamp(explicit_start)
    else:
        offset = (ncols // 2) * step_minutes
        start = pd.Timestamp(center_time) - pd.Timedelta(minutes=offset)
    return [start + pd.Timedelta(minutes=step_minutes * i) for i in range(ncols)]


def rolling_window_times(end_time: pd.Timestamp, step_minutes: int, n_expected: int = 12):
    start = pd.Timestamp(end_time) - pd.Timedelta(minutes=step_minutes * (n_expected - 1))
    return [start + pd.Timedelta(minutes=step_minutes * i) for i in range(n_expected)]


def cml_rolling_hour(cml_5min, cml_times: pd.DatetimeIndex, end_time: pd.Timestamp, min_valid_fraction: float):
    idxs = [nearest_time_index(cml_times, t, pd.Timedelta(minutes=0)) for t in rolling_window_times(end_time, 5, 12)]
    if any(i < 0 for i in idxs):
        return None, 0, 12
    block = np.asarray(cml_5min.isel(time=idxs).values, dtype="float32")  # 12, y, x
    valid = np.isfinite(block)
    n_valid = np.sum(valid, axis=0)
    required = int(math.ceil(12 * min_valid_fraction))
    accum = np.nansum(block / 12.0, axis=0)
    accum[n_valid < required] = np.nan
    return accum, int(np.nanmax(n_valid)), required


def native_radar_file_for_time(
    radar_root: Path,
    t: pd.Timestamp,
    filename_template: str,
) -> Path:
    """
    Resolve the native daily radar file using the actual directory layout:

        <radar_root>/<YYYY>/<YYYY>/<MM>/
            YW_2017.002_<YYYYMMDD>.nc

    Example:
        dwd_radolan/2021/2021/04/YW_2017.002_20210426.nc

    This deliberately avoids the recursive discovery helper because the radar
    tree also contains derived/path-averaged files such as RADKLIM_5min/*CML_RR.nc.
    """
    t = pd.Timestamp(t)
    yyyy = f"{t.year:04d}"
    mm = f"{t.month:02d}"
    yyyymmdd = t.strftime("%Y%m%d")

    filename = filename_template.format(
        YYYY=yyyy,
        MM=mm,
        YYYYMMDD=yyyymmdd,
    )

    expected = radar_root / yyyy / yyyy / mm / filename
    if expected.is_file():
        return expected

    # Safe exact-basename fallback. Exclude known derived/path-average products.
    matches = [
        p for p in radar_root.rglob(filename)
        if p.is_file()
        and "RADKLIM_5min" not in p.parts
        and "CML_RR" not in p.name
    ]

    if not matches:
        raise FileNotFoundError(
            "Native radar file not found.\n"
            f"Expected:\n  {expected}\n"
            f"Exact basename searched recursively:\n  {filename}"
        )

    expected_tail = (yyyy, yyyy, mm, filename)

    def score(p: Path):
        parts = p.parts
        exact_tree = len(parts) >= 4 and tuple(parts[-4:]) == expected_tail
        return (0 if exact_tree else 1, len(parts), str(p))

    matches.sort(key=score)
    chosen = matches[0]

    if len(matches) > 1:
        print(
            f"  WARNING: {len(matches)} native radar files match {filename}; "
            f"using {chosen}"
        )

    return chosen


def radar_field(
    cache,
    radar_root: Path,
    filename_template: str,
    t: pd.Timestamp,
):
    path = native_radar_file_for_time(
        radar_root=radar_root,
        t=pd.Timestamp(t),
        filename_template=filename_template,
    )
    return cache.field_at(str(path), pd.Timestamp(t))


def radar_rolling_hour(
    cache,
    radar_root: Path,
    filename_template: str,
    end_time: pd.Timestamp,
    min_valid_fraction: float,
):
    stack = []
    expected = 12
    for t in rolling_window_times(end_time, 5, expected):
        _, _, field, _ = radar_field(
            cache,
            radar_root,
            filename_template,
            pd.Timestamp(t),
        )
        stack.append(np.asarray(field, dtype="float32"))
    block = np.stack(stack, axis=0)
    valid = np.isfinite(block)
    n_valid = np.sum(valid, axis=0)
    required = int(math.ceil(expected * min_valid_fraction))
    accum = np.nansum(block / 12.0, axis=0)
    accum[n_valid < required] = np.nan
    return accum, int(np.nanmax(n_valid)), required


def make_multitime_panel(
    *,
    left_fields: list[np.ndarray],
    right_fields: list[np.ndarray],
    left_x: np.ndarray,
    left_y: np.ndarray,
    right_x: np.ndarray,
    right_y: np.ndarray,
    times: Sequence[pd.Timestamp],
    extent: Sequence[float],
    bounds: Sequence[float],
    colorbar_label: str,
    out_png: Path,
    out_pdf: Path,
    row_labels=("CML", "RADAR"),
    dpi: int = 600,
    draw_links: bool = False,
    geom=None,
    mask_threshold: float = 0.1,
):
    proj = ccrs.PlateCarree()
    bounds = parse_bounds(bounds)

    cmap = discrete_essd_cmap(bounds)
    # There are len(bounds)-1 intervals and exactly the same number of colors.
    # The upper extension is handled by cmap.set_over() and the colorbar.
    norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=False)

    ncols = len(times)
    fig_w = max(12.0, 1.70 * ncols)
    fig_h = 4.9

    fig, axes = plt.subplots(
        2,
        ncols,
        figsize=(fig_w, fig_h),
        subplot_kw={"projection": proj},
        squeeze=False,
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.995,
        top=0.92,
        bottom=0.21,
        wspace=0.01,
        hspace=0.01,
    )

    lon_a, lat_a, lon_b, lat_b = (
        geom if geom is not None else (None, None, None, None)
    )
    mappable = None

    for col in range(ncols):
        ax_cml = axes[0, col]
        ax_radar = axes[1, col]

        for ax in (ax_cml, ax_radar):
            setup_axis(ax, extent, proj)
            if draw_links and lon_a is not None:
                plot_links(
                    ax,
                    proj,
                    lon_a,
                    lat_a,
                    lon_b,
                    lat_b,
                    extent,
                )

        cml_plot = dry_mask(left_fields[col], mask_threshold)
        radar_plot = dry_mask(right_fields[col], mask_threshold)

        mappable = ax_cml.pcolormesh(
            left_x,
            left_y,
            cml_plot,
            transform=proj,
            shading="nearest",
            cmap=cmap,
            norm=norm,
            zorder=2,
            rasterized=True,
        )

        ax_radar.pcolormesh(
            right_x,
            right_y,
            radar_plot,
            transform=proj,
            shading="nearest",
            cmap=cmap,
            norm=norm,
            zorder=2,
            rasterized=True,
        )

        ax_cml.set_title(
            pd.Timestamp(times[col]).strftime("%H:%M"),
            fontsize=13,
            fontweight="bold",
            pad=8,
        )

    fig.text(
        0.04,
        0.70,
        row_labels[0],
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.04,
        0.33,
        row_labels[1],
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
    )

    # Wide bar, as in the reference figure. Every boundary label fits without
    # merging 8, 9, 10, 11, ... into one unreadable string.
    cax = fig.add_axes([0.15, 0.095, 0.72, 0.040])
    add_discrete_horizontal_colorbar(
        fig=fig,
        mappable=mappable,
        cax=cax,
        bounds=bounds,
        label=colorbar_label,
        labelsize=14,
        ticksize=8.5,
        label_all_bounds=True,
    )

    fig.savefig(
        out_png,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        out_pdf,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def make_two_panel_map(
    *,
    left_field: np.ndarray,
    right_field: np.ndarray,
    left_x: np.ndarray,
    left_y: np.ndarray,
    right_x: np.ndarray,
    right_y: np.ndarray,
    extent: Sequence[float],
    bounds: Sequence[float],
    colorbar_label: str,
    out_png: Path,
    out_pdf: Path,
    dpi: int,
    draw_links: bool,
    geom,
    mask_threshold: float,
):
    proj = ccrs.PlateCarree()
    bounds = parse_bounds(bounds)

    cmap = discrete_essd_cmap(bounds)
    norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=False)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(6.8, 3.65),
        subplot_kw={"projection": proj},
        squeeze=False,
    )
    fig.subplots_adjust(
        left=0.03,
        right=0.99,
        bottom=0.19,
        top=0.98,
        wspace=0.015,
    )

    lon_a, lat_a, lon_b, lat_b = geom
    for ax in axes[0]:
        setup_axis(ax, extent, proj)
        if draw_links:
            plot_links(
                ax,
                proj,
                lon_a,
                lat_a,
                lon_b,
                lat_b,
                extent,
            )

    mappable = axes[0, 0].pcolormesh(
        left_x,
        left_y,
        dry_mask(left_field, mask_threshold),
        transform=proj,
        shading="nearest",
        cmap=cmap,
        norm=norm,
        zorder=2,
        rasterized=True,
    )

    axes[0, 1].pcolormesh(
        right_x,
        right_y,
        dry_mask(right_field, mask_threshold),
        transform=proj,
        shading="nearest",
        cmap=cmap,
        norm=norm,
        zorder=2,
        rasterized=True,
    )

    for ax, tag in zip(axes[0], ["(a)", "(b)"]):
        ax.text(
            0.02,
            0.98,
            tag,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.5,
            fontweight="bold",
            bbox=dict(
                boxstyle="round,pad=0.16",
                fc="white",
                ec="0.80",
                lw=0.35,
                alpha=0.92,
            ),
            zorder=6,
        )

    # The compact two-panel figure cannot fit nineteen printed labels.
    # It uses selected major labels and still marks every boundary.
    cax = fig.add_axes([0.17, 0.075, 0.66, 0.045])
    add_discrete_horizontal_colorbar(
        fig=fig,
        mappable=mappable,
        cax=cax,
        bounds=bounds,
        label=colorbar_label,
        labelsize=10.5,
        ticksize=8.5,
        label_all_bounds=False,
    )

    fig.savefig(
        out_png,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        out_pdf,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(
        description="Create example-style multitime CML/radar event figures."
    )
    p.add_argument("--gridded-nc", default="/home/jovyan/ESSD_publication_dataset/blettner_apr_oct/period_products/grid/blettner_cml_event_20210621_5min.nc")
    p.add_argument("--cml-var", default="R")
    p.add_argument("--cml-row-label", default="CML")
    p.add_argument("--cml-geom-nc", default="/home/jovyan/ESSD_publication_dataset/blettner_apr_oct/BLETTNER_GRAF_RAIN_APR_OCT_2021.nc")
    p.add_argument("--site-a-lon-var", default="site_a_longitude")
    p.add_argument("--site-a-lat-var", default="site_a_latitude")
    p.add_argument("--site-b-lon-var", default="site_b_longitude")
    p.add_argument("--site-b-lat-var", default="site_b_latitude")

    p.add_argument("--radar-root", default="/home/jovyan/dwd_radolan")
    p.add_argument(
        "--radar-pattern",
        default="YW_2017.002_*.nc",
        help="Compatibility option; native files are resolved using the exact daily template.",
    )
    p.add_argument(
        "--radar-filename-template",
        default="YW_2017.002_{YYYYMMDD}.nc",
        help=(
            "Exact native daily radar basename. Expected directory structure: "
            "<radar-root>/<YYYY>/<YYYY>/<MM>/<filename>."
        ),
    )
    p.add_argument("--radar-scale-to-mmph", type=float, default=12.0)

    p.add_argument("--selection-cache-csv", default=None,
                   help="Existing selected_instances_cache.csv from the earlier event run. "
                        "If omitted, this script searches /home/jovyan recursively and verifies event 3 / instance 1.")
    p.add_argument("--event-rank", type=int, default=3)
    p.add_argument("--instance-rank", type=int, default=1)

    p.add_argument("--base-helper-script", default=None)
    p.add_argument("--extent", nargs=4, type=float, default=None,
                   metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"))
    p.add_argument("--zoom-half-lon", type=float, default=0.55)
    p.add_argument("--zoom-half-lat", type=float, default=0.45)

    p.add_argument("--ncols", type=int, default=8)
    p.add_argument("--panel-step-minutes", type=int, default=60,
                   help="Time spacing between successive columns.")
    p.add_argument("--start-time", default="2021-06-21 19:45:00",
                   help="Optional explicit first panel time, e.g. '2021-06-21 21:50:00'.")
    p.add_argument("--min-valid-fraction-hourly", type=float, default=1.0)

    p.add_argument("--native-bounds", nargs="+", type=float, default=DEFAULT_NATIVE_BOUNDS)
    p.add_argument("--hourly-bounds", nargs="+", type=float, default=DEFAULT_HOURLY_BOUNDS)
    p.add_argument("--native-plot-min", type=float, default=0.1)
    p.add_argument("--hourly-plot-min", type=float, default=0.1)

    p.add_argument("--draw-links", action="store_true", default=False)
    p.add_argument("--save-single-figures", action="store_true", default=False)
    p.add_argument("--dpi", type=int, default=600)

    # Multi-sensor event time series, same night as the transboundary map.
    p.add_argument(
        "--v2-nc",
        default="/home/jovyan/ESSD_publication_dataset/blettner_apr_oct/BLETTNER_GRAF_RAIN_APR_OCT_2021.nc",
    )
    p.add_argument(
        "--path-radar-nc",
        default="/home/jovyan/ESSD_publication_dataset/CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--gauge-nc",
        default="/home/jovyan/ESSD_publication_dataset/CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--temperature-nc",
        default="/home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_TEMPERATURE_2021.nc",
    )
    p.add_argument(
        "--humidity-nc",
        default="/home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_HUMIDITY_2021.nc",
    )
    p.add_argument("--timeseries-start", default="2021-06-21 21:50:00")
    p.add_argument("--timeseries-end", default="2021-06-22 02:30:00")
    p.add_argument("--max-gauge-distance-km", type=float, default=1.0)
    p.add_argument("--skip-timeseries", action="store_true")
    p.add_argument("--out-dir", default="/home/jovyan/ESSD_publication_dataset/blettner_apr_oct/event_20210621")
    return p.parse_args()




# =============================================================================
# Same-night multi-sensor event time series from the Blettner V2 product
# =============================================================================

def _decode_1d_strings(values):
    a = np.asarray(values)
    if a.dtype.kind == "S":
        if a.ndim == 1:
            return np.asarray(
                [x.decode("utf-8", errors="ignore").replace("\x00", "").strip()
                 for x in a],
                dtype=object,
            )
        return np.asarray([
            "".join(
                x.decode("utf-8", errors="ignore")
                if isinstance(x, (bytes, np.bytes_)) else str(x)
                for x in np.ravel(row)
            ).replace("\x00", "").strip()
            for row in a
        ], dtype=object)
    return np.asarray([str(x).replace("\x00", "").strip() for x in a.reshape(-1)], dtype=object)


def _station_dim(ds, da):
    tdim = next((d for d in da.dims if "time" in d.lower()), None)
    other = [d for d in da.dims if d != tdim]
    for cand in ("station", "WSI", "station_id"):
        if cand in other:
            return cand
    return other[-1] if other else None


def _station_ids(ds, dim):
    if dim in ds.coords:
        return _decode_1d_strings(ds[dim].values)
    for cand in ("station_id", "WSI", "station", "GH_ID"):
        if cand in ds and dim in ds[cand].dims:
            return _decode_1d_strings(ds[cand].values)
    return np.asarray([str(i) for i in range(ds.sizes[dim])], dtype=object)


def _lat_lon_arrays(ds, dim):
    lat_name = next((n for n in ("lat", "latitude") if n in ds), None)
    lon_name = next((n for n in ("lon", "longitude") if n in ds), None)
    if lat_name is None or lon_name is None:
        raise KeyError(f"Could not find lat/lon in {list(ds.variables)}")
    lat = np.asarray(ds[lat_name].values, dtype=float).reshape(-1)
    lon = np.asarray(ds[lon_name].values, dtype=float).reshape(-1)
    if lat.size != ds.sizes[dim] or lon.size != ds.sizes[dim]:
        raise ValueError("Station coordinate size does not match station dimension")
    return lat, lon


def _point_segment_distance_km(plon, plat, lon0, lat0, lon1, lat1):
    lat_ref = np.deg2rad((lat0 + lat1 + plat) / 3.0)
    sx = 111.320 * np.cos(lat_ref)
    sy = 110.574
    ax, ay = lon0 * sx, lat0 * sy
    bx, by = lon1 * sx, lat1 * sy
    px, py = plon * sx, plat * sy
    vx, vy = bx - ax, by - ay
    den = vx * vx + vy * vy
    if den <= 0:
        return float(np.hypot(px - ax, py - ay))
    u = ((px - ax) * vx + (py - ay) * vy) / den
    u = float(np.clip(u, 0.0, 1.0))
    qx, qy = ax + u * vx, ay + u * vy
    return float(np.hypot(px - qx, py - qy))


def _choose_numeric_precip_var(ds):
    for name in ("rr", "R", "RR", "precipitation", "precipitation_amount"):
        if name in ds.data_vars:
            return name
    for name, da in ds.data_vars.items():
        if np.issubdtype(da.dtype, np.number) and any("time" in d.lower() for d in da.dims):
            return name
    raise KeyError("No precipitation variable found")


def _choose_event_link_and_gauge(v2, gauge, start, end, max_dist_km):
    cml_ids = _decode_1d_strings(v2["cml_id"].values)
    freq = np.asarray(v2["frequency_GHz"].values, dtype=float)
    if freq.shape[0] != v2.sizes["sublink_id"]:
        freq = np.moveaxis(freq, -1, 0)
    eband = np.any(
        ((freq >= 71.0) & (freq <= 76.0)) |
        ((freq >= 81.0) & (freq <= 86.0)),
        axis=0,
    )

    rain = v2["R_cml"].sel(time=slice(start, end)).transpose(
        "time", "sublink_id", "cml_id"
    )
    a = np.asarray(rain.values, dtype=np.float32)
    finite = np.isfinite(a)
    count = np.sum(finite, axis=1)
    sums = np.sum(np.where(finite, a, 0.0), axis=1)
    link_mean = np.full((a.shape[0], a.shape[2]), np.nan, dtype=np.float32)
    ok = count > 0
    link_mean[ok] = (sums[ok] / count[ok]).astype(np.float32)
    score = np.nanmax(link_mean, axis=0)

    gvar = _choose_numeric_precip_var(gauge)
    gda = gauge[gvar]
    gdim = _station_dim(gauge, gda)
    gids = _station_ids(gauge, gdim)
    glat, glon = _lat_lon_arrays(gauge, gdim)

    lon0 = np.asarray(v2["site_a_longitude"].values, dtype=float)
    lat0 = np.asarray(v2["site_a_latitude"].values, dtype=float)
    lon1 = np.asarray(v2["site_b_longitude"].values, dtype=float)
    lat1 = np.asarray(v2["site_b_latitude"].values, dtype=float)

    best = None
    for j in np.flatnonzero(eband & np.isfinite(score)):
        if score[j] <= 0:
            continue
        dists = np.asarray([
            _point_segment_distance_km(
                glon[k], glat[k], lon0[j], lat0[j], lon1[j], lat1[j]
            )
            for k in range(len(gids))
        ])
        k = int(np.nanargmin(dists))
        if not np.isfinite(dists[k]) or dists[k] > max_dist_km:
            continue
        candidate = (float(score[j]), int(j), int(k), float(dists[k]))
        if best is None or candidate[0] > best[0]:
            best = candidate

    if best is None:
        raise RuntimeError(
            f"No E-band Blettner CML with a gauge within {max_dist_km:g} km "
            f"has positive rainfall during {start} -> {end}"
        )

    _, j, k, dist = best
    return {
        "cml_index": j,
        "cml_id": str(cml_ids[j]),
        "gauge_index": k,
        "gauge_id": str(gids[k]),
        "gauge_distance_km": dist,
        "gauge_var": gvar,
        "gauge_dim": gdim,
        "gauge_lon": float(glon[k]),
        "gauge_lat": float(glat[k]),
    }


def _series_for_station(ds, preferred_var, station_id, target_lon, target_lat):
    var = preferred_var if preferred_var in ds.data_vars else None
    if var is None:
        for name, da in ds.data_vars.items():
            if np.issubdtype(da.dtype, np.number) and any("time" in d.lower() for d in da.dims):
                var = name
                break
    if var is None:
        raise KeyError("No time-dependent numeric variable found")

    da = ds[var]
    sdim = _station_dim(ds, da)
    ids = _station_ids(ds, sdim)
    matches = np.flatnonzero(ids == str(station_id))
    if len(matches):
        idx = int(matches[0])
    else:
        lat, lon = _lat_lon_arrays(ds, sdim)
        scale_x = 111.320 * np.cos(np.deg2rad(target_lat))
        d = np.hypot((lon - target_lon) * scale_x, (lat - target_lat) * 110.574)
        idx = int(np.nanargmin(d))

    tdim = next(d for d in da.dims if "time" in d.lower())
    out = da.isel({sdim: idx})
    return pd.DatetimeIndex(out[tdim].values), np.asarray(out.values, dtype=float).reshape(-1)


def _resample_cml_5min(times, values):
    s = pd.Series(np.asarray(values, dtype=float), index=pd.DatetimeIndex(times))
    # project convention: 1-min rate -> 5-min equivalent rate
    return s.resample("5min", label="left", closed="left").sum(min_count=1) / 5.0


def make_event_timeseries(args):
    v2_path = Path(args.v2_nc).expanduser().resolve()
    if not v2_path.is_file():
        raise FileNotFoundError(v2_path)
    gauge_path = Path(args.gauge_nc).expanduser().resolve()
    radar_path = Path(args.path_radar_nc).expanduser().resolve()
    temp_path = Path(args.temperature_nc).expanduser().resolve()
    hum_path = Path(args.humidity_nc).expanduser().resolve()

    start = pd.Timestamp(args.timeseries_start)
    end = pd.Timestamp(args.timeseries_end)

    with xr.open_dataset(v2_path) as v2, \
         xr.open_dataset(gauge_path) as gauge, \
         xr.open_dataset(radar_path) as radar, \
         xr.open_dataset(temp_path) as temp, \
         xr.open_dataset(hum_path) as hum:

        if "trsl_qc" not in v2:
            raise KeyError(
                "The V2 product has no trsl_qc. Re-run step 1 with --write-trsl; "
                "this is required for the top TRSL panel."
            )

        chosen = _choose_event_link_and_gauge(
            v2, gauge, start, end, args.max_gauge_distance_km
        )
        j = chosen["cml_index"]
        cid = chosen["cml_id"]

        print("\nEVENT TIMESERIES SELECTION")
        print(f"  CML             : {cid}")
        print(f"  Gauge           : {chosen['gauge_id']}")
        print(f"  Gauge distance  : {chosen['gauge_distance_km']:.3f} km")

        tr = v2["trsl_qc"].sel(time=slice(start, end)).isel(cml_id=j).transpose(
            "time", "sublink_id"
        )
        tr_t = pd.DatetimeIndex(tr["time"].values)
        tr_v = np.asarray(tr.values, dtype=float)

        rr = v2["R_cml"].sel(time=slice(start, end)).isel(cml_id=j).transpose(
            "time", "sublink_id"
        )
        rr_t = pd.DatetimeIndex(rr["time"].values)
        rr_a = np.asarray(rr.values, dtype=float)
        rr_count = np.sum(np.isfinite(rr_a), axis=1)
        rr_sum = np.sum(np.where(np.isfinite(rr_a), rr_a, 0.0), axis=1)
        rr_mean = np.full(len(rr_t), np.nan, dtype=float)
        ok = rr_count > 0
        rr_mean[ok] = rr_sum[ok] / rr_count[ok]
        rr5 = _resample_cml_5min(rr_t, rr_mean)

        # Path radar.
        rid = _decode_1d_strings(radar["cml_id"].values)
        pos = np.flatnonzero(rid == cid)
        if len(pos) == 0:
            raise KeyError(f"CML {cid} is not present in path-radar file")
        rvar = next((n for n in ("rr", "R_radar", "rainfall_rate") if n in radar), None)
        if rvar is None:
            raise KeyError("No radar rain variable found")
        rda = radar[rvar].isel(cml_id=int(pos[0]))
        rtdim = next(d for d in rda.dims if "time" in d.lower())
        rtimes = pd.DatetimeIndex(rda[rtdim].values)
        rmask = (rtimes >= start) & (rtimes <= end)
        radar_rate = np.asarray(rda.values, dtype=float)[rmask] * 12.0
        radar_times = rtimes[rmask]

        # Gauge precipitation: CHMI 10-min amount -> intensity.
        gda = gauge[chosen["gauge_var"]].isel(
            {chosen["gauge_dim"]: chosen["gauge_index"]}
        )
        gtdim = next(d for d in gda.dims if "time" in d.lower())
        gtimes = pd.DatetimeIndex(gda[gtdim].values)
        gmask = (gtimes >= start) & (gtimes <= end + pd.Timedelta(minutes=10))
        gauge_rate = np.asarray(gda.values, dtype=float)[gmask] * 6.0
        # CHMI timestamps label the preceding 10-minute interval; plot at interval centre.
        gauge_times = gtimes[gmask] - pd.Timedelta(minutes=5)

        tt, temp_v = _series_for_station(
            temp, "T", chosen["gauge_id"], chosen["gauge_lon"], chosen["gauge_lat"]
        )
        ht, hum_v = _series_for_station(
            hum, "H", chosen["gauge_id"], chosen["gauge_lon"], chosen["gauge_lat"]
        )
        tm = (tt >= start) & (tt <= end)
        hm = (ht >= start) & (ht <= end)

        mpl.rcParams.update({
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 0.9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "pdf.fonttype": 42,
        })

        fig, axes = plt.subplots(
            3, 1, figsize=(10.2, 7.7), sharex=True,
            gridspec_kw={"height_ratios": [1.0, 1.05, 1.0], "hspace": 0.08}
        )

        # TRSL
        axes[0].plot(tr_t, tr_v[:, 0], lw=1.6, label="Sublink 1")
        if tr_v.shape[1] > 1:
            axes[0].plot(tr_t, tr_v[:, 1], lw=1.6, label="Sublink 2")
        axes[0].set_ylabel("TRSL [dB]", fontsize=13)
        finite_tr = tr_v[np.isfinite(tr_v)]
        if finite_tr.size:
            lo = max(0.0, 5.0 * np.floor((np.nanmin(finite_tr) - 2.0) / 5.0))
            hi = 5.0 * np.ceil((np.nanmax(finite_tr) + 4.0) / 5.0)
            axes[0].set_ylim(lo, max(lo + 10.0, hi))
        axes[0].legend(loc="upper right", frameon=False, ncol=2, fontsize=11)

        # Rainfall
        axes[1].step(rr5.index, rr5.values, where="post", lw=2.4,
                     color="#1f4e79", label="CML")
        axes[1].step(radar_times, radar_rate, where="post", lw=1.8, ls="--",
                     color="#d97706", label="Radar")
        axes[1].bar(
            gauge_times, gauge_rate,
            width=pd.Timedelta(minutes=10),
            align="center", alpha=0.42,
            color="#9ecae1", edgecolor="none", label="Gauge",
        )
        axes[1].set_ylabel(r"Rain rate [mm h$^{-1}$]", fontsize=13)
        ymax_data = np.nanmax(np.r_[
            rr5.values[np.isfinite(rr5.values)],
            radar_rate[np.isfinite(radar_rate)],
            gauge_rate[np.isfinite(gauge_rate)],
            [1.0],
        ])
        ymax = min(150.0, max(25.0, 25.0 * np.ceil((ymax_data * 1.15) / 25.0)))
        axes[1].set_ylim(0, ymax)
        axes[1].legend(loc="upper right", frameon=False, ncol=3, fontsize=11)

        # T / RH
        ax3 = axes[2]
        ax3r = ax3.twinx()
        ax3.plot(tt[tm], temp_v[tm], color="#a84636", lw=1.8)
        ax3r.plot(ht[hm], hum_v[hm], color="#347858", lw=1.8)
        ax3.set_ylabel("Temperature [°C]", color="#a84636", fontsize=12.5)
        ax3r.set_ylabel("Humidity [%]", color="#347858", fontsize=12.5, rotation=270, labelpad=18)
        ax3.tick_params(axis="y", colors="#a84636")
        ax3r.tick_params(axis="y", colors="#347858")
        ax3.set_ylim(10, 20)
        ax3r.set_ylim(60, 100)

        for ax in axes:
            ax.set_xlim(start, end)
            ax.tick_params(axis="both", which="major", length=6, width=0.9, labelsize=10)
            ax.tick_params(axis="both", which="minor", length=3, width=0.7)
            ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))

        locator = mdates.HourLocator(interval=1)
        axes[-1].xaxis.set_major_locator(locator)
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        axes[-1].set_xlabel("Time [UTC]", fontsize=13)

        # Date labels at every midnight crossed by the event window.
        midnights = pd.date_range(
            start.normalize() + pd.Timedelta(days=1),
            end.normalize(),
            freq="D",
        )

        for midnight in midnights:
            if start <= midnight <= end:
                axes[-1].text(
                    midnight,
                    -0.17,
                    midnight.strftime("%d %b"),
                    transform=axes[-1].get_xaxis_transform(),
                    ha="center",
                    va="top",
                    fontsize=10,
                    bbox=dict(
                        facecolor="white",
                        edgecolor="none",
                        pad=0.4,
                    ),
                )

        out_dir = Path(args.out_dir) / "timeseries"
        out_dir.mkdir(parents=True, exist_ok=True)

        tag = (
            f"{start:%Y%m%d_%H%M}_"
            f"{end:%Y%m%d_%H%M}"
        )

        png = out_dir / f"event_{tag}_multisensor_timeseries.png"
        pdf = out_dir / f"event_{tag}_multisensor_timeseries.pdf"
        fig.savefig(png, dpi=args.dpi, bbox_inches="tight", pad_inches=0.03)
        fig.savefig(pdf, bbox_inches="tight", pad_inches=0.03)
        plt.close(fig)

        pd.DataFrame([chosen]).to_csv(
            out_dir / f"event_{tag}_timeseries_selection.csv", index=False
        )
        print(f"  Timeseries      : {png}")
        return png, pdf

def resolve_selection_cache(path, event_rank, instance_rank):
    """Reuse the exact old event extent instead of reselecting a new event."""
    candidates = []
    if path:
        candidates.append(Path(path).expanduser().resolve())
    else:
        root = Path("/home/jovyan")
        # Prefer the exact filename written by the old multiscale script.
        candidates.extend(root.rglob("selected_instances_cache.csv"))

    checked = []
    for p in candidates:
        if not p.is_file():
            continue
        checked.append(p)
        try:
            df = pd.read_csv(p)
        except Exception:
            continue
        required = {
            "event_rank", "instance_rank", "time",
            "extent_lon_min", "extent_lon_max",
            "extent_lat_min", "extent_lat_max",
        }
        if not required.issubset(df.columns):
            continue
        hit = df[
            (df["event_rank"].astype(int) == int(event_rank)) &
            (df["instance_rank"].astype(int) == int(instance_rank))
        ]
        if hit.empty:
            continue

        # Accept any event explicitly present in the supplied selection cache.
        t = pd.Timestamp(hit.iloc[0]["time"])
        print(
            f"Using event selection cache: {p} | "
            f"event={event_rank}/{instance_rank} | time={t}"
        )
        return p

    shown = "\n  ".join(str(p) for p in checked[:20])
    raise FileNotFoundError(
        "Could not locate a selected_instances_cache.csv containing "
        f"event {event_rank}, instance {instance_rank}.\n"
        "Pass it explicitly with --selection-cache-csv.\n"
        f"Checked:\n  {shown}"
    )



def main():
    args = parse_args()
    t0 = time.perf_counter()

    base, helper_path = load_base_module(args.base_helper_script)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    native_dir = out_dir / "native_5min"
    hourly_dir = out_dir / "hourly_1h"
    native_dir.mkdir(parents=True, exist_ok=True)
    hourly_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 110)
    print("ESSD EXAMPLE-STYLE MULTITIME EVENT PLOTS")
    print("=" * 110)
    print(f"Base helper       : {helper_path}")
    print(f"Gridded CML       : {args.gridded_nc}")
    print(f"CML geometry      : {args.cml_geom_nc}")
    print(f"Radar root        : {args.radar_root}")
    print("Radar layout      : <root>/<YYYY>/<YYYY>/<MM>/YW_2017.002_<YYYYMMDD>.nc")
    print(f"Selection cache   : {args.selection_cache_csv}")
    print(f"Event/instance    : {args.event_rank}/{args.instance_rank}")
    print(f"Columns           : {args.ncols}")
    print(f"Column step       : {args.panel_step_minutes} min")
    print(f"Output            : {out_dir}")
    print("=" * 110)

    args.selection_cache_csv = str(
        resolve_selection_cache(
            args.selection_cache_csv,
            args.event_rank,
            args.instance_rank,
        )
    )
    selected = pd.read_csv(args.selection_cache_csv, parse_dates=["time"])
    row = selected[
        (selected["event_rank"] == args.event_rank) &
        (selected["instance_rank"] == args.instance_rank)
    ]
    if row.empty:
        raise KeyError(f"Event {args.event_rank}, instance {args.instance_rank} not found in cache.")
    row = row.iloc[0]

    event_time = pd.Timestamp(row["time"])
    if args.extent is None:
        extent = [
            float(row["extent_lon_min"]),
            float(row["extent_lon_max"]),
            float(row["extent_lat_min"]),
            float(row["extent_lat_max"]),
        ]
    else:
        extent = list(args.extent)

    print(f"Anchor event time : {event_time}")
    print(f"Map extent        : {extent}")

    ds_cml = xr.open_dataset(args.gridded_nc)
    if args.cml_var not in ds_cml.data_vars:
        raise KeyError(f"CML variable {args.cml_var!r} not in file. Available: {list(ds_cml.data_vars)}")
    cml_5min = ds_cml[args.cml_var]
    x_name = next((v for v in ["x", "lon", "longitude"] if v in ds_cml.coords), None)
    y_name = next((v for v in ["y", "lat", "latitude"] if v in ds_cml.coords), None)
    if x_name is None or y_name is None:
        raise KeyError("CML grid needs x/y or lon/lat coordinates.")
    xg = np.asarray(ds_cml[x_name].values)
    yg = np.asarray(ds_cml[y_name].values)
    cml_times = pd.DatetimeIndex(cml_5min.time.values)

    with xr.open_dataset(args.cml_geom_nc) as ds_geom:
        lon_a = np.asarray(ds_geom[args.site_a_lon_var].values, dtype=float)
        lat_a = np.asarray(ds_geom[args.site_a_lat_var].values, dtype=float)
        lon_b = np.asarray(ds_geom[args.site_b_lon_var].values, dtype=float)
        lat_b = np.asarray(ds_geom[args.site_b_lat_var].values, dtype=float)
    geom = (lon_a, lat_a, lon_b, lat_b)

    radar_root = Path(args.radar_root).expanduser().resolve()
    if not radar_root.is_dir():
        raise FileNotFoundError(f"Radar root does not exist: {radar_root}")

    cache = base.RadarCache(base.RADAR_VAR_CANDIDATES, args.radar_scale_to_mmph)

    try:
        # Determine native radar timing from the exact nested daily file for
        # the selected event. Do not use recursive discovery here.
        event_radar_path = native_radar_file_for_time(
            radar_root,
            event_time,
            args.radar_filename_template,
        )
        print(f"Native radar file  : {event_radar_path}")

        try:
            _, _, _, times = cache.open(str(event_radar_path))
        except Exception as exc:
            raise RuntimeError(
                f"Could not open expected native radar file {event_radar_path}: {exc}"
            ) from exc

        radar_minutes = base.infer_step_minutes(times) if len(times) >= 2 else 5.0
        print(f"Native radar step : {radar_minutes:g} min")
        tol = pd.Timedelta(minutes=float(radar_minutes) / 2.0 + 0.1)

        panel_times = native_time_series(event_time, args.ncols, args.panel_step_minutes, args.start_time)
        print("Panel times       :")
        for t in panel_times:
            print(f"  {pd.Timestamp(t)}")

        native_cml = []
        native_rad = []
        rlon = rlat = None
        for t in panel_times:
            idx = nearest_time_index(cml_times, pd.Timestamp(t), tol)
            if idx < 0:
                raise RuntimeError(f"No matching CML 5-min time for {t}")
            native_cml.append(np.asarray(cml_5min.isel(time=idx).values, dtype="float32"))

            rlon_, rlat_, field_, actual_ = radar_field(
                cache,
                radar_root,
                args.radar_filename_template,
                pd.Timestamp(t),
            )
            if abs(pd.Timestamp(actual_) - pd.Timestamp(t)) > tol:
                raise RuntimeError(f"No matching radar 5-min field for {t}; got {actual_}")
            rlon, rlat = np.asarray(rlon_), np.asarray(rlat_)
            native_rad.append(np.asarray(field_, dtype="float32"))

        hourly_cml = []
        hourly_rad = []
        for t in panel_times:
            acc_cml, _, _ = cml_rolling_hour(cml_5min, cml_times, pd.Timestamp(t), args.min_valid_fraction_hourly)
            if acc_cml is None:
                raise RuntimeError(f"Could not build rolling 1-hour CML accumulation ending at {t}")
            acc_rad, _, _ = radar_rolling_hour(
                cache,
                radar_root,
                args.radar_filename_template,
                pd.Timestamp(t),
                args.min_valid_fraction_hourly,
            )
            hourly_cml.append(acc_cml.astype("float32"))
            hourly_rad.append(acc_rad.astype("float32"))

        # Multitime figures.
        native_png = native_dir / f"event{args.event_rank:02d}_instance{args.instance_rank:02d}_multitime_native5min.png"
        native_pdf = native_dir / f"event{args.event_rank:02d}_instance{args.instance_rank:02d}_multitime_native5min.pdf"
        make_multitime_panel(
            left_fields=native_cml,
            right_fields=native_rad,
            left_x=xg, left_y=yg,
            right_x=rlon, right_y=rlat,
            times=panel_times,
            extent=extent,
            bounds=args.native_bounds,
            colorbar_label=r"Rainfall rate [mm h$^{-1}$]",
            out_png=native_png,
            out_pdf=native_pdf,
            row_labels=(args.cml_row_label, "RADAR"),
            dpi=args.dpi,
            draw_links=args.draw_links,
            geom=geom,
            mask_threshold=args.native_plot_min,
        )

        hourly_png = hourly_dir / f"event{args.event_rank:02d}_instance{args.instance_rank:02d}_multitime_hourly1h.png"
        hourly_pdf = hourly_dir / f"event{args.event_rank:02d}_instance{args.instance_rank:02d}_multitime_hourly1h.pdf"
        make_multitime_panel(
            left_fields=hourly_cml,
            right_fields=hourly_rad,
            left_x=xg, left_y=yg,
            right_x=rlon, right_y=rlat,
            times=panel_times,
            extent=extent,
            bounds=args.hourly_bounds,
            colorbar_label="Hourly rain sum [mm]",
            out_png=hourly_png,
            out_pdf=hourly_pdf,
            row_labels=(args.cml_row_label, "RADAR"),
            dpi=args.dpi,
            draw_links=args.draw_links,
            geom=geom,
            mask_threshold=args.hourly_plot_min,
        )

        if args.save_single_figures:
            sing_native_dir = native_dir / "single_maps"
            sing_hour_dir = hourly_dir / "single_maps"
            sing_native_dir.mkdir(parents=True, exist_ok=True)
            sing_hour_dir.mkdir(parents=True, exist_ok=True)

            for t, lc, rc, lh, rh in zip(panel_times, native_cml, native_rad, hourly_cml, hourly_rad):
                tag = pd.Timestamp(t).strftime("%Y%m%d_%H%M")
                make_two_panel_map(
                    left_field=lc, right_field=rc,
                    left_x=xg, left_y=yg, right_x=rlon, right_y=rlat,
                    extent=extent, bounds=args.native_bounds,
                    colorbar_label=r"Rainfall rate [mm h$^{-1}$]",
                    out_png=sing_native_dir / f"native5min_{tag}.png",
                    out_pdf=sing_native_dir / f"native5min_{tag}.pdf",
                    dpi=args.dpi, draw_links=args.draw_links, geom=geom,
                    mask_threshold=args.native_plot_min,
                )
                make_two_panel_map(
                    left_field=lh, right_field=rh,
                    left_x=xg, left_y=yg, right_x=rlon, right_y=rlat,
                    extent=extent, bounds=args.hourly_bounds,
                    colorbar_label="Hourly rain sum [mm]",
                    out_png=sing_hour_dir / f"hourly1h_{tag}.png",
                    out_pdf=sing_hour_dir / f"hourly1h_{tag}.pdf",
                    dpi=args.dpi, draw_links=args.draw_links, geom=geom,
                    mask_threshold=args.hourly_plot_min,
                )

        summary = pd.DataFrame({
            "panel_time": panel_times,
            "event_rank": args.event_rank,
            "instance_rank": args.instance_rank,
            "native_multitime_png": str(native_png),
            "hourly_multitime_png": str(hourly_png),
        })
        summary.to_csv(out_dir / f"event{args.event_rank:02d}_instance{args.instance_rank:02d}_multitime_summary.csv", index=False)

        print("\nSaved:")
        print(f"  {native_png}")
        print(f"  {native_pdf}")
        print(f"  {hourly_png}")
        print(f"  {hourly_pdf}")
        if args.save_single_figures:
            print(f"  {native_dir / 'single_maps'}")
            print(f"  {hourly_dir / 'single_maps'}")

        if not args.skip_timeseries:
            make_event_timeseries(args)

        print(f"Elapsed: {time.perf_counter() - t0:.1f} s")

    finally:
        cache.close()
        ds_cml.close()
        gc.collect()


if __name__ == "__main__":
    main()
