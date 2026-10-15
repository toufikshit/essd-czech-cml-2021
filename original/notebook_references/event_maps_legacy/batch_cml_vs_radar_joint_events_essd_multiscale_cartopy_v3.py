
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ESSD publication plotting: for each selected joint convective event instance,
save TWO clean publication-ready comparison figures:

1) native 5-min rain-rate comparison
   (a) CML, (b) radar

2) 1-hour accumulated rainfall comparison
   (a) CML, (b) radar

Design choices follow the user's requested style:
- Cartopy required
- clean basemap and small zoom window
- no in-map titles or time boxes
- only panel annotations (a) and (b)
- only the colorbar label
- PNG + PDF output
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import sys
import time
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

import cartopy.crs as ccrs
import cartopy.feature as cfeature


DEFAULT_NATIVE_BOUNDS = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 30, 50]
DEFAULT_HOURLY_BOUNDS = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 30, 50]


def _candidate_base_paths() -> list[Path]:
    here = Path(__file__).resolve().parent
    return [
        here / "batch_cml_vs_radar_joint_events_essd_1h.py",
        Path("/home/jovyan/batch_cml_vs_radar_joint_events_essd_1h.py"),
        Path("/mnt/data/batch_cml_vs_radar_joint_events_essd_1h.py"),
    ]


def load_base_module():
    base_path = None
    for cand in _candidate_base_paths():
        if cand.exists():
            base_path = cand
            break
    if base_path is None:
        raise FileNotFoundError(
            "Could not locate batch_cml_vs_radar_joint_events_essd_1h.py. "
            "Place it in the same directory as this script or in /home/jovyan."
        )

    spec = importlib.util.spec_from_file_location("essd_base_events", str(base_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load base script from {base_path}")
    mod = importlib.util.module_from_spec(spec)

    # Python 3.13 dataclasses resolves annotations through sys.modules while
    # the imported module is executing. Register it before exec_module().
    sys.modules[spec.name] = mod
    try:
        spec.loader.exec_module(mod)
    except Exception:
        sys.modules.pop(spec.name, None)
        raise

    return mod, base_path


def dry_mask(a: np.ndarray, threshold: float) -> np.ndarray:
    out = np.asarray(a, dtype=float).copy()
    out[out < threshold] = np.nan
    return out


def setup_clean_map_axis(ax, extent: Sequence[float], proj):
    ax.set_extent(extent, crs=proj)
    ax.add_feature(cfeature.LAND, facecolor="0.98", zorder=0)
    ax.add_feature(cfeature.BORDERS, linewidth=0.55, edgecolor="0.70", zorder=1)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.45, edgecolor="0.75", zorder=1)
    ax.add_feature(cfeature.LAKES, edgecolor="0.80", facecolor="none", linewidth=0.30, zorder=1)
    ax.outline_patch.set_linewidth(0.6)
    # remove all tick labels / grid labels for the minimalist journal style
    try:
        ax.set_xticks([], crs=proj)
        ax.set_yticks([], crs=proj)
    except Exception:
        pass
    ax.gridlines(draw_labels=False, linewidth=0.0, alpha=0.0)


def plot_links(ax, proj, lon_a, lat_a, lon_b, lat_b, extent):
    keep = (
        np.isfinite(lon_a) & np.isfinite(lat_a) & np.isfinite(lon_b) & np.isfinite(lat_b)
        & (
            ((lon_a >= extent[0]) & (lon_a <= extent[1]) & (lat_a >= extent[2]) & (lat_a <= extent[3]))
            | ((lon_b >= extent[0]) & (lon_b <= extent[1]) & (lat_b >= extent[2]) & (lat_b <= extent[3]))
        )
    )
    idx = np.flatnonzero(keep)
    for i in idx:
        ax.plot(
            [lon_a[i], lon_b[i]],
            [lat_a[i], lat_b[i]],
            transform=proj,
            color="0.25",
            linewidth=0.25,
            alpha=0.22,
            zorder=3,
        )


def make_two_panel_map(
    *,
    left_field: np.ndarray,
    right_field: np.ndarray,
    x_left: np.ndarray,
    y_left: np.ndarray,
    x_right: np.ndarray,
    y_right: np.ndarray,
    extent: Sequence[float],
    bounds: Sequence[float],
    cmap_name: str,
    left_label: str,
    right_label: str,
    colorbar_label: str,
    out_png: Path,
    out_pdf: Path,
    dpi: int,
    draw_links: bool,
    geom,
    field_mask_threshold: float,
):
    proj = ccrs.PlateCarree()
    bounds = np.asarray(bounds, dtype=float)
    norm = BoundaryNorm(bounds, ncolors=256, extend="both")
    cmap = mpl.colormaps[cmap_name].copy()
    cmap.set_bad("white")

    fig, axes = plt.subplots(
        1, 2,
        figsize=(6.8, 3.65),
        subplot_kw={"projection": proj},
        squeeze=False,
    )
    ax1, ax2 = axes[0]
    fig.subplots_adjust(left=0.03, right=0.985, bottom=0.18, top=0.98, wspace=0.02)

    setup_clean_map_axis(ax1, extent, proj)
    setup_clean_map_axis(ax2, extent, proj)

    lon_a, lat_a, lon_b, lat_b = geom
    if draw_links:
        plot_links(ax1, proj, lon_a, lat_a, lon_b, lat_b, extent)
        plot_links(ax2, proj, lon_a, lat_a, lon_b, lat_b, extent)

    left_plot = dry_mask(left_field, field_mask_threshold)
    right_plot = dry_mask(right_field, field_mask_threshold)

    m1 = ax1.pcolormesh(
        x_left, y_left, left_plot,
        transform=proj,
        shading="nearest",
        cmap=cmap,
        norm=norm,
        zorder=2,
        rasterized=True,
    )
    ax2.pcolormesh(
        x_right, y_right, right_plot,
        transform=proj,
        shading="nearest",
        cmap=cmap,
        norm=norm,
        zorder=2,
        rasterized=True,
    )

    for ax, label in [(ax1, left_label), (ax2, right_label)]:
        ax.text(
            0.02, 0.98, label,
            transform=ax.transAxes,
            ha="left", va="top",
            fontsize=10.5, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="0.80", lw=0.35, alpha=0.92),
            zorder=6,
        )

    cax = fig.add_axes([0.23, 0.08, 0.54, 0.030])
    cb = fig.colorbar(m1, cax=cax, orientation="horizontal", ticks=bounds)
    cb.set_label(colorbar_label, fontsize=10.5)
    cb.ax.tick_params(labelsize=8.5)
    cb.outline.set_linewidth(0.7)

    fig.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(
        description="Create clean Cartopy publication figures for both native 5-min rainfall-rate maps and 1-hour accumulated rainfall maps."
    )

    p.add_argument("--gridded-nc", required=True)
    p.add_argument("--cml-var", default="R")
    p.add_argument("--cml-geom-nc", required=True)
    p.add_argument("--site-a-lon-var", default="site_a_longitude")
    p.add_argument("--site-a-lat-var", default="site_a_latitude")
    p.add_argument("--site-b-lon-var", default="site_b_longitude")
    p.add_argument("--site-b-lat-var", default="site_b_latitude")

    p.add_argument("--radar-root", required=True)
    p.add_argument("--radar-pattern", default="YW_2017.002_*.nc")
    p.add_argument("--radar-scale-to-mmph", type=float, default=12.0)
    p.add_argument("--only-dates", nargs="*", default=None)

    p.add_argument("--extent", nargs=4, type=float, default=None,
                   metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"))
    p.add_argument("--zoom-half-lon", type=float, default=0.55)
    p.add_argument("--zoom-half-lat", type=float, default=0.45)

    p.add_argument("--convective-threshold", type=float, default=10.0)
    p.add_argument("--dry-threshold", type=float, default=0.1)
    p.add_argument("--min-convective-pixels", type=int, default=100)
    p.add_argument("--min-peak-mmph", type=float, default=10.0)
    p.add_argument("--cml-rain-threshold", type=float, default=0.1)
    p.add_argument("--min-cml-wet-pixels", type=int, default=40)
    p.add_argument("--min-joint-wet-pixels", type=int, default=15)
    p.add_argument("--min-links-in-window", type=int, default=40)
    p.add_argument("--min-cml-availability-fraction", type=float, default=0.15)
    p.add_argument("--scan-stride", type=int, default=1)

    p.add_argument("--event-gap-min", type=float, default=90.0)
    p.add_argument("--max-events", type=int, default=6)
    p.add_argument("--instances-per-event", type=int, default=3)
    p.add_argument("--min-instance-gap-min", type=float, default=30.0)

    p.add_argument("--cml-time-label", choices=["left", "right"], default="left")
    p.add_argument("--cml-time-closed", choices=["left", "right"], default="left")
    p.add_argument("--radar-time-tolerance-min", type=float, default=30.0)

    p.add_argument("--hourly-min-valid-fraction", type=float, default=1.0)

    p.add_argument("--native-bounds", nargs="+", type=float, default=DEFAULT_NATIVE_BOUNDS)
    p.add_argument("--hourly-bounds", nargs="+", type=float, default=DEFAULT_HOURLY_BOUNDS)
    p.add_argument("--native-plot-min-mmph", type=float, default=0.1)
    p.add_argument("--hourly-plot-min-mm", type=float, default=0.1)
    p.add_argument("--cmap", default="turbo")
    p.add_argument("--draw-links", action="store_true", default=False)
    p.add_argument("--no-draw-links", dest="draw_links", action="store_false")
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--out-dir", required=True)

    args = p.parse_args()
    if not (0.0 < args.hourly_min_valid_fraction <= 1.0):
        p.error("--hourly-min-valid-fraction must be in (0, 1]")
    return args


def main():
    args = parse_args()
    base, base_path = load_base_module()

    t0 = time.perf_counter()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    native_dir = out_dir / "native_5min"
    hourly_dir = out_dir / "hourly_1h"
    native_dir.mkdir(parents=True, exist_ok=True)
    hourly_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 110)
    print("ESSD CARTOPY MULTISCALE EVENT MAPS")
    print("=" * 110)
    print(f"Base helper script : {base_path}")
    print(f"Gridded CML        : {args.gridded_nc}")
    print(f"CML geometry       : {args.cml_geom_nc}")
    print(f"Radar root         : {args.radar_root}")
    print(f"Output             : {out_dir}")
    print(f"Zoom half-size     : lon={args.zoom_half_lon}, lat={args.zoom_half_lat}")
    print("=" * 110)

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
        base_extent = base.network_search_extent(lon_a, lat_a, lon_b, lat_b, pad=0.05)
    else:
        base_extent = list(args.extent)
    link_density = base.link_density_on_grid(xg, yg, lon_a, lat_a, lon_b, lat_b)

    radar_files = base.discover_radar_files(args.radar_root, args.radar_pattern)
    date_filter = base.parse_date_filter(args.only_dates)
    if date_filter is not None:
        radar_files = [r for r in radar_files if r.date.date() in date_filter]
    if not radar_files:
        raise FileNotFoundError(f"No radar files found below {args.radar_root!r} with pattern {args.radar_pattern!r}")

    cache = base.RadarCache(base.RADAR_VAR_CANDIDATES, args.radar_scale_to_mmph)

    try:
        radar_minutes = None
        for info in radar_files:
            try:
                _, _, _, times = cache.open(info.path)
            except Exception:
                continue
            if len(times) >= 2:
                radar_minutes = base.infer_step_minutes(times)
                break
        if radar_minutes is None:
            radar_minutes = 5.0
        print(f"Native radar temporal resolution: {radar_minutes:g} min")

        cml_aligned = base.resample_cml_to_radar_resolution(
            cml_da,
            radar_minutes=radar_minutes,
            label=args.cml_time_label,
            closed=args.cml_time_closed,
        )
        cml_times = pd.DatetimeIndex(cml_aligned.time.values)

        cml_hourly, cml_native_minutes, cml_hour_expected, cml_hour_required = (
            base.resample_cml_to_hourly_accumulation(
                cml_da,
                min_valid_fraction=args.hourly_min_valid_fraction,
            )
        )
        cml_hour_times = pd.DatetimeIndex(cml_hourly.time.values)
        print(
            f"CML hourly accumulation: native={cml_native_minutes:g} min; "
            f"expected={cml_hour_expected}; required={cml_hour_required}"
        )

        cml_date_min = cml_times.min().date()
        cml_date_max = cml_times.max().date()
        radar_files = [r for r in radar_files if cml_date_min <= r.date.date() <= cml_date_max]
        mapper = base.RadarToCMLGridMapper(xg, yg)

        print("Scanning for joint convective radar-CML instances ...")
        candidates = base.scan_candidates(
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
        events = base.group_candidates_into_events(candidates, args.event_gap_min)
        selected = base.select_instances_from_events(
            events,
            max_events=args.max_events,
            instances_per_event=args.instances_per_event,
            min_instance_gap_minutes=args.min_instance_gap_min,
        )
        print(f"Candidate convective times: {len(candidates)}")
        print(f"Detected events: {len(events)}")
        print(f"Selected instances: {len(selected)}")

        rows = []
        for event_rank, instance_rank, candidate in selected:
            # native 5-min CML match
            if len(cml_times) == 0:
                continue
            pos = cml_times.get_indexer([candidate.time], method="nearest")[0]
            if pos < 0:
                continue
            cml_time = pd.Timestamp(cml_times[pos])
            mismatch_min = abs((cml_time - candidate.time).total_seconds()) / 60.0
            if mismatch_min > args.radar_time_tolerance_min:
                continue

            cml_native_field = np.asarray(cml_aligned.isel(time=int(pos)).values, dtype="float32")

            try:
                rlon_native, rlat_native, radar_native_field, radar_actual = cache.field_at(
                    candidate.radar_path,
                    candidate.time,
                )
            except Exception as exc:
                print(f"Skip native plot {candidate.time}: {exc}")
                continue

            hour_start = pd.Timestamp(candidate.time).floor("h")
            hour_end = hour_start + pd.Timedelta(hours=1)
            hour_pos = cml_hour_times.get_indexer([hour_start])[0]
            if hour_pos < 0 or cml_hour_times[hour_pos] != hour_start:
                print(f"Skip hourly plot {candidate.time}: no exact CML hourly bin")
                continue
            cml_hour_field = np.asarray(cml_hourly.isel(time=int(hour_pos)).values, dtype="float32")

            try:
                (
                    rlon_hour, rlat_hour, radar_hour_field, radar_hour_start, radar_hour_end,
                    radar_step_minutes, radar_fields_used, radar_fields_expected,
                ) = cache.hourly_accumulation(
                    candidate.radar_path,
                    hour_start,
                    min_valid_fraction=args.hourly_min_valid_fraction,
                )
            except Exception as exc:
                print(f"Skip hourly radar accumulation {candidate.time}: {exc}")
                continue

            plot_extent = list(candidate.extent)

            stem_native = (
                f"event{event_rank:02d}_instance{instance_rank:02d}_"
                f"native5min_{candidate.time.strftime('%Y%m%d_%H%M')}"
            )
            stem_hour = (
                f"event{event_rank:02d}_instance{instance_rank:02d}_"
                f"hour1_{hour_start.strftime('%Y%m%d_%H00')}"
            )

            native_png = native_dir / f"{stem_native}.png"
            native_pdf = native_dir / f"{stem_native}.pdf"
            hour_png = hourly_dir / f"{stem_hour}.png"
            hour_pdf = hourly_dir / f"{stem_hour}.pdf"

            make_two_panel_map(
                left_field=cml_native_field,
                right_field=radar_native_field,
                x_left=xg,
                y_left=yg,
                x_right=rlon_native,
                y_right=rlat_native,
                extent=plot_extent,
                bounds=args.native_bounds,
                cmap_name=args.cmap,
                left_label="(a)",
                right_label="(b)",
                colorbar_label=r"Rainfall rate [mm h$^{-1}$]",
                out_png=native_png,
                out_pdf=native_pdf,
                dpi=args.dpi,
                draw_links=args.draw_links,
                geom=geom,
                field_mask_threshold=args.native_plot_min_mmph,
            )

            make_two_panel_map(
                left_field=cml_hour_field,
                right_field=radar_hour_field,
                x_left=xg,
                y_left=yg,
                x_right=rlon_hour,
                y_right=rlat_hour,
                extent=plot_extent,
                bounds=args.hourly_bounds,
                cmap_name=args.cmap,
                left_label="(a)",
                right_label="(b)",
                colorbar_label="Hourly rain sum [mm]",
                out_png=hour_png,
                out_pdf=hour_pdf,
                dpi=args.dpi,
                draw_links=args.draw_links,
                geom=geom,
                field_mask_threshold=args.hourly_plot_min_mm,
            )

            rows.append({
                "event_rank": event_rank,
                "instance_rank": instance_rank,
                "event_native_radar_time": candidate.time,
                "native_cml_match_time": cml_time,
                "native_time_difference_min": mismatch_min,
                "hour_start_utc": hour_start,
                "hour_end_utc": hour_end,
                "radar_hour_fields_found": radar_fields_used,
                "radar_hour_expected_fields": radar_fields_expected,
                "joint_wet_pixels": candidate.joint_wet_pixels,
                "links_in_window": candidate.links_in_window,
                "radar_max_mmph": candidate.radar_max,
                "cml_max_mmph": candidate.cml_max,
                "extent_lon_min": plot_extent[0],
                "extent_lon_max": plot_extent[1],
                "extent_lat_min": plot_extent[2],
                "extent_lat_max": plot_extent[3],
                "native_png": native_png.name,
                "native_pdf": native_pdf.name,
                "hourly_png": hour_png.name,
                "hourly_pdf": hour_pdf.name,
            })

            del cml_native_field, radar_native_field, cml_hour_field, radar_hour_field
            gc.collect()

        summary = pd.DataFrame(rows)
        summary_path = out_dir / "selected_convective_instances_multiscale.csv"
        summary.to_csv(summary_path, index=False)

        print(f"\nSaved output root : {out_dir}")
        print(f"Native figures    : {native_dir}")
        print(f"Hourly figures    : {hourly_dir}")
        print(f"Summary           : {summary_path}")
        print(f"Total elapsed     : {time.perf_counter() - t0:.1f} s")

    finally:
        cache.close()
        ds_cml.close()


if __name__ == "__main__":
    main()
