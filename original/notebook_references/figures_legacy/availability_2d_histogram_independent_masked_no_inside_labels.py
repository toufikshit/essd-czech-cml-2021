#!/usr/bin/env python3
"""
Independent data-availability 2D histograms
===========================================

Outputs one independent figure per sensor:
    availability_cml.png / .pdf
    availability_radar.png / .pdf
    availability_gauge.png / .pdf

Important plotting rules implemented here:
    1. CML and Radar use the SAME colormap and the SAME colorbar normalization.
       Therefore, the same color means the same number of links/entries.
    2. Rain gauge uses a COMPLETELY INDEPENDENT colorbar normalization.
       Gauge station counts are not forced to use the CML/Radar scale.
    3. Exact 0% availability is MASKED, not counted in the 0-10% bin.
       So true zero-availability cells stay white.
    4. Empty histogram cells are also white.
    5. Each sensor plot is saved as a separate standalone figure.

The script keeps the memory optimization from the original version:
extra dimensions such as sublink_id are reduced before any transpose is applied.
"""

from pathlib import Path
import time
import warnings

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

warnings.filterwarnings("ignore")

# ============================================================
# SETTINGS
# ============================================================
CML_FILE = "cml_daily_nc_v2_merged.nc"
RADAR_FILE = "RADAR_S2021E2021_5MIN_CZ.nc"
GAUGE_FILE = "CHMI_10MIN_v1.nc"

OUT_DIR = Path("availability_plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

FREQ = "W"  # weekly availability blocks

# Availability bins in percent.
# Rows are: 0-10, 10-20, ..., 90-100
BIN_EDGES = np.arange(0, 101, 10)
N_BINS = len(BIN_EDGES) - 1

# Critical option requested:
# exact 0% availability is not counted in the first bin.
MASK_ZERO_AVAILABILITY = True
ZERO_EPS = 0.0

# Color scale
CMAP_NAME = "magma"
POWER_GAMMA = 0.5

# Figure options
FIGSIZE_SINGLE = (13, 3.8)
DPI = 300

# Optional combined figure. Default False because the requested output is
# independent standalone plots for each sensor.
SAVE_COMBINED = False


# ============================================================
# HELPER FUNCTIONS
# ============================================================
def find_dim(ds, candidates):
    """Find the first matching dimension name."""
    for c in candidates:
        if c in ds.dims:
            return c
    return None


def find_var(ds, candidates):
    """Find the first matching variable or coordinate name."""
    for c in candidates:
        if c in ds.data_vars or c in ds.coords:
            return c
    return None


def get_coord(ds, candidates):
    """Return numeric coordinate values if present, else None."""
    name = find_var(ds, candidates)
    if name is None:
        return None

    arr = ds[name].values
    if arr.dtype.kind in ("S", "U", "O"):
        return None

    try:
        return np.asarray(arr, dtype=float)
    except Exception:
        return None


def build_valid_mask(da, entity_dim, time_dim, reduce_mode="any"):
    """
    Build a boolean valid mask with shape (n_entities, n_time).

    Extra dimensions are collapsed first. This avoids copying the full
    high-dimensional CML signal array via moveaxis/transpose.

    Parameters
    ----------
    da : xarray.DataArray
        Input data.
    entity_dim : str
        Link/station/entity dimension.
    time_dim : str
        Time dimension.
    reduce_mode : {"any", "all"}
        How to collapse extra dimensions.

    Returns
    -------
    valid : ndarray, bool
        Shape (n_entities, n_time).
    """
    dims = list(da.dims)

    if entity_dim not in dims or time_dim not in dims:
        raise ValueError(
            f"DataArray dims {dims} do not contain entity_dim={entity_dim!r} "
            f"and time_dim={time_dim!r}"
        )

    extra_axes = tuple(
        i for i, d in enumerate(dims)
        if d not in (entity_dim, time_dim)
    )

    # This loads the target variable once.
    finite = np.isfinite(da.values)

    # Collapse sublink/channel/grid-like extra axes before transposing.
    if extra_axes:
        if reduce_mode == "all":
            finite = np.all(finite, axis=extra_axes)
        elif reduce_mode == "any":
            finite = np.any(finite, axis=extra_axes)
        else:
            raise ValueError("reduce_mode must be 'any' or 'all'")

    # After reducing extra axes, the remaining axes are entity and time,
    # but their order follows the original DataArray order.
    remaining = [d for d in dims if d in (entity_dim, time_dim)]
    if remaining != [entity_dim, time_dim]:
        finite = finite.T

    return np.ascontiguousarray(finite, dtype=bool)


def compute_grouped_availability(valid, time_vals, freq="W"):
    """
    Compute grouped availability for a valid mask of shape (n_entities, n_time).

    Returns
    -------
    avail : ndarray, float32
        Shape (n_entities, n_time_blocks), values in [0, 1].
    time_axis : pandas.DatetimeIndex
        Start time of each group.
    """
    time_vals = pd.DatetimeIndex(time_vals)
    periods = time_vals.to_period(freq)

    groups = pd.Series(np.arange(len(time_vals))).groupby(periods).indices
    periods_sorted = sorted(groups.keys())

    avail = np.zeros((valid.shape[0], len(periods_sorted)), dtype=np.float32)

    for i, p in enumerate(periods_sorted):
        idx = np.asarray(groups[p], dtype=int)
        avail[:, i] = valid[:, idx].mean(axis=1)

    time_axis = pd.DatetimeIndex([p.start_time for p in periods_sorted])
    return avail, time_axis


# ============================================================
# SENSOR AVAILABILITY
# ============================================================
def compute_cml_availability(filepath, freq="W"):
    """Compute weekly availability for each CML link."""
    print("[CML] Loading dataset...")
    t0 = time.perf_counter()

    with xr.open_dataset(filepath) as ds:
        cml_dim = find_dim(ds, ["cml_id", "cml", "link", "link_id"])
        time_dim = find_dim(ds, ["time", "t"])

        if cml_dim is None or time_dim is None:
            raise ValueError(f"[CML] Could not find CML/time dims. Found: {list(ds.dims)}")

        signal_var = find_var(ds, ["tsl", "rsl", "tx", "rx", "txrx", "signal"])
        if signal_var is None:
            raise ValueError(f"[CML] Could not find signal variable. Found: {list(ds.data_vars)}")

        signal = ds[signal_var]
        print(f"[CML] Using variable {signal_var!r}, shape={signal.shape}, dims={signal.dims}")

        lon_a = get_coord(ds, ["site_a_longitude", "lon_a", "longitude_a", "start_longitude"])
        lat_a = get_coord(ds, ["site_a_latitude", "lat_a", "latitude_a", "start_latitude"])
        lon_b = get_coord(ds, ["site_b_longitude", "lon_b", "longitude_b", "end_longitude"])
        lat_b = get_coord(ds, ["site_b_latitude", "lat_b", "latitude_b", "end_latitude"])

        time_vals = pd.DatetimeIndex(ds[time_dim].values)
        n_cml = ds.sizes[cml_dim]

        print(f"[CML] {n_cml} links, {len(time_vals)} time steps")
        print(f"[CML] Time range: {time_vals.min()} -> {time_vals.max()}")

        # CML link is available if ANY sublink/channel is valid.
        valid = build_valid_mask(signal, cml_dim, time_dim, reduce_mode="any")

        expected_shape = (n_cml, len(time_vals))
        if valid.shape != expected_shape:
            raise ValueError(f"[CML] valid shape {valid.shape} != expected {expected_shape}")

        avail, time_axis = compute_grouped_availability(valid, time_vals, freq=freq)

    print(f"[CML] Done in {time.perf_counter() - t0:.1f}s")
    return avail, time_axis, lon_a, lat_a, lon_b, lat_b


def compute_radar_availability(filepath, freq="W"):
    """Compute weekly availability for radar data."""
    print("[RADAR] Loading dataset...")
    t0 = time.perf_counter()

    with xr.open_dataset(filepath) as ds:
        time_dim = find_dim(ds, ["time", "t"])
        if time_dim is None:
            raise ValueError(f"[RADAR] Could not find time dim. Found: {list(ds.dims)}")

        rr_var = find_var(ds, ["rr", "rain_rate", "rainfall", "precip", "precipitation", "R", "R_radolan"])
        if rr_var is None:
            raise ValueError(f"[RADAR] Could not find radar rain variable. Found: {list(ds.data_vars)}")

        rr = ds[rr_var]
        print(f"[RADAR] Using variable {rr_var!r}, shape={rr.shape}, dims={rr.dims}")

        cml_dim = find_dim(ds, ["cml_id", "cml", "link", "link_id", "id"])
        if cml_dim is not None and cml_dim in rr.dims:
            spatial_dim = cml_dim
        else:
            other_dims = [d for d in rr.dims if d != time_dim]
            spatial_dim = other_dims[0] if other_dims else None

        time_vals = pd.DatetimeIndex(ds[time_dim].values)

        if spatial_dim is None:
            n_entries = 1
            valid = np.isfinite(rr.values).reshape(1, -1)
        else:
            n_entries = ds.sizes[spatial_dim]
            valid = build_valid_mask(rr, spatial_dim, time_dim, reduce_mode="any")

        print(f"[RADAR] {n_entries} entries, {len(time_vals)} time steps")
        print(f"[RADAR] Time range: {time_vals.min()} -> {time_vals.max()}")

        expected_shape = (n_entries, len(time_vals))
        if valid.shape != expected_shape:
            raise ValueError(f"[RADAR] valid shape {valid.shape} != expected {expected_shape}")

        lon = get_coord(ds, ["longitude", "lon", "x", "site_a_longitude", "lon_a"])
        lat = get_coord(ds, ["latitude", "lat", "y", "site_a_latitude", "lat_a"])
        lon_b = get_coord(ds, ["site_b_longitude", "lon_b", "longitude_b"])
        lat_b = get_coord(ds, ["site_b_latitude", "lat_b", "latitude_b"])

        valid_coords = None
        if spatial_dim is not None:
            if (
                lon is not None and lat is not None and
                lon_b is not None and lat_b is not None and
                len(lon) == n_entries and len(lat) == n_entries and
                len(lon_b) == n_entries and len(lat_b) == n_entries
            ):
                valid_coords = (
                    np.isfinite(lon) & np.isfinite(lat) &
                    np.isfinite(lon_b) & np.isfinite(lat_b)
                )
            elif (
                lon is not None and lat is not None and
                len(lon) == n_entries and len(lat) == n_entries
            ):
                valid_coords = np.isfinite(lon) & np.isfinite(lat)

        avail, time_axis = compute_grouped_availability(valid, time_vals, freq=freq)

    print(f"[RADAR] Done in {time.perf_counter() - t0:.1f}s")
    return avail, time_axis, valid_coords


def compute_gauge_availability(filepath, freq="W"):
    """Compute weekly availability for rain gauge stations."""
    print("[GAUGE] Loading dataset...")
    t0 = time.perf_counter()

    with xr.open_dataset(filepath) as ds:
        time_dim = find_dim(ds, ["time", "t"])
        if time_dim is None:
            raise ValueError(f"[GAUGE] Could not find time dim. Found: {list(ds.dims)}")

        rr_var = find_var(ds, ["rr", "rain_rate", "rainfall", "precip", "precipitation", "P", "R"])
        if rr_var is None:
            raise ValueError(f"[GAUGE] Could not find gauge precipitation variable. Found: {list(ds.data_vars)}")

        rr = ds[rr_var]
        print(f"[GAUGE] Using variable {rr_var!r}, shape={rr.shape}, dims={rr.dims}")

        station_dim = find_dim(ds, ["station", "station_id", "id", "gauge", "site", "WSI"])
        if station_dim is not None and station_dim not in rr.dims:
            station_dim = None

        if station_dim is None:
            other_dims = [d for d in rr.dims if d != time_dim]
            station_dim = other_dims[0] if other_dims else None

        time_vals = pd.DatetimeIndex(ds[time_dim].values)

        if station_dim is None:
            n_stations = 1
            valid = np.isfinite(rr.values).reshape(1, -1)
        else:
            n_stations = ds.sizes[station_dim]
            valid = build_valid_mask(rr, station_dim, time_dim, reduce_mode="any")

        print(f"[GAUGE] {n_stations} stations, {len(time_vals)} time steps")
        print(f"[GAUGE] Time range: {time_vals.min()} -> {time_vals.max()}")

        expected_shape = (n_stations, len(time_vals))
        if valid.shape != expected_shape:
            raise ValueError(f"[GAUGE] valid shape {valid.shape} != expected {expected_shape}")

        avail, time_axis = compute_grouped_availability(valid, time_vals, freq=freq)

    print(f"[GAUGE] Done in {time.perf_counter() - t0:.1f}s")
    return avail, time_axis


# ============================================================
# HISTOGRAM
# ============================================================
def compute_2d_histogram(
    avail,
    bin_edges=BIN_EDGES,
    mask_zero_availability=MASK_ZERO_AVAILABILITY,
):
    """
    Compute counts per availability bin per time block.

    Parameters
    ----------
    avail : ndarray
        Shape (n_sensors, n_blocks), values in [0, 1].
    bin_edges : ndarray
        Availability bin edges in percent.
    mask_zero_availability : bool
        If True, exact 0% availability is ignored, so it appears white instead
        of being counted in the 0-10% availability bin.

    Returns
    -------
    H : ndarray, int32
        Shape (n_bins, n_blocks).
    n_zero_masked : int
        Number of exact 0% availability entries ignored.
    """
    n_blocks = avail.shape[1]
    avail_pct = np.asarray(avail, dtype=np.float32) * 100.0

    finite = np.isfinite(avail_pct)
    n_zero_masked = 0

    if mask_zero_availability:
        zero_mask = finite & (avail_pct <= ZERO_EPS)
        n_zero_masked = int(zero_mask.sum())
        finite &= avail_pct > ZERO_EPS

    if not finite.any():
        H = np.zeros((len(bin_edges) - 1, n_blocks), dtype=np.int32)
        return H, n_zero_masked

    valid_avail = avail_pct[finite]
    valid_avail = np.clip(valid_avail, bin_edges[0], bin_edges[-1])

    # General binning. 100% belongs to the last bin.
    bins = np.searchsorted(bin_edges, valid_avail, side="right") - 1
    bins = np.clip(bins, 0, len(bin_edges) - 2).astype(np.int64)

    blocks = np.broadcast_to(np.arange(n_blocks), avail_pct.shape)[finite].astype(np.int64)

    flat = bins * n_blocks + blocks
    H = np.bincount(
        flat,
        minlength=(len(bin_edges) - 1) * n_blocks,
    ).reshape(len(bin_edges) - 1, n_blocks)

    return H.astype(np.int32), n_zero_masked


# ============================================================
# PLOTTING
# ============================================================
def make_cmap():
    """Create colormap where masked/empty values are white."""
    cmap = plt.get_cmap(CMAP_NAME).copy()
    cmap.set_bad("white")    # masked histogram cells
    cmap.set_under("white")  # safety for under-range values
    return cmap


def make_count_norm(vmax, gamma=POWER_GAMMA):
    """
    Make count normalization for positive counts.

    Count 0 is masked before plotting, so the scale starts at count 1.
    """
    vmax = int(max(vmax, 1))
    if vmax <= 1:
        return mcolors.Normalize(vmin=1, vmax=1.001)
    return mcolors.PowerNorm(gamma=gamma, vmin=1, vmax=vmax)


def mask_empty_cells(H):
    """Mask zero-count histogram cells so they plot as white."""
    return np.ma.masked_where(H <= 0, H)


def format_time_ticks(ax, time_axis, n_blocks):
    """Format one independent x-axis."""
    n_ticks = min(10, n_blocks)

    pos = list(np.unique(np.linspace(0, n_blocks - 1, n_ticks, dtype=int)))

    # Avoid duplicate crowded tick near the final right-edge tick.
    edge_guard = max(1, int(n_blocks * 0.03))
    pos = [p for p in pos if p < n_blocks - edge_guard]
    pos.append(n_blocks)

    labels = [time_axis[min(p, n_blocks - 1)].strftime("%Y-%m") for p in pos]

    ax.set_xticks(pos)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Week", fontsize=10)
    ax.set_xlim(0, n_blocks)


def draw_single_panel(ax, H, time_axis, cmap, norm, annotation, cbar_label):
    """Draw one standalone sensor panel."""
    n_blocks = H.shape[1]
    H_plot = mask_empty_cells(H)

    im = ax.imshow(
        H_plot,
        aspect="auto",
        origin="lower",
        cmap=cmap,
        norm=norm,
        extent=[0, n_blocks, BIN_EDGES[0], BIN_EDGES[-1]],
        interpolation="nearest",
    )

    ax.set_yticks(BIN_EDGES)
    ax.set_yticklabels([str(v) for v in BIN_EDGES], fontsize=8)
    ax.set_ylabel("Availability (%)", fontsize=10)

    format_time_ticks(ax, time_axis, n_blocks)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_position(("outward", 8))
    ax.spines["bottom"].set_position(("outward", 8))
    ax.tick_params(axis="both", direction="out", length=4)

    # No in-panel sensor-name annotation.
    # This keeps the plot interior clean while preserving all other styling,
    # axes, colorbars, color scales, and saved file names.
    if annotation:
        ax.annotate(
            annotation,
            xy=(0.012, 0.94),
            xycoords="axes fraction",
            fontsize=10,
            fontweight="bold",
            va="top",
            ha="left",
            bbox=dict(
                boxstyle="round,pad=0.3",
                fc="white",
                ec="0.7",
                alpha=0.90,
            ),
        )

    cax = make_axes_locatable(ax).append_axes("right", size="1.6%", pad=0.12)
    cbar = ax.figure.colorbar(im, cax=cax)
    cbar.set_label(cbar_label, fontsize=9)
    cbar.ax.tick_params(labelsize=8, direction="out")
    cbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))

    return im


def save_sensor_plot(key, H, time_axis, cmap, norm, annotation, cbar_label):
    """Save one standalone sensor figure as PNG and PDF."""
    fig, ax = plt.subplots(1, 1, figsize=FIGSIZE_SINGLE)
    draw_single_panel(ax, H, time_axis, cmap, norm, annotation, cbar_label)

    png_path = OUT_DIR / f"availability_{key}.png"
    pdf_path = OUT_DIR / f"availability_{key}.pdf"

    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    return png_path, pdf_path


def save_combined_plot(panels, cmap):
    """Optional combined 3-panel figure."""
    fig, axes = plt.subplots(
        len(panels),
        1,
        figsize=(13, 11),
        gridspec_kw={"hspace": 0.55},
    )

    for ax, panel in zip(axes, panels):
        draw_single_panel(
            ax=ax,
            H=panel["H"],
            time_axis=panel["time"],
            cmap=cmap,
            norm=panel["norm"],
            annotation=panel["annotation"],
            cbar_label=panel["cbar_label"],
        )

    png_path = OUT_DIR / "availability_all_sensors.png"
    pdf_path = OUT_DIR / "availability_all_sensors.pdf"

    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    return png_path, pdf_path


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    total_start = time.perf_counter()

    print("=" * 70)
    print("INDEPENDENT DATA AVAILABILITY 2D HISTOGRAMS")
    print("=" * 70)

    # Load data and compute weekly availability.
    cml_avail, cml_time, lon_a, lat_a, lon_b, lat_b = compute_cml_availability(CML_FILE, FREQ)
    radar_avail, radar_time, radar_valid_coords = compute_radar_availability(RADAR_FILE, FREQ)
    gauge_avail, gauge_time = compute_gauge_availability(GAUGE_FILE, FREQ)

    # Metadata filtering.
    print("\n" + "=" * 70)
    print("FILTERING BY METADATA")
    print("=" * 70)

    if (
        lon_a is not None and lat_a is not None and
        lon_b is not None and lat_b is not None and
        len(lon_a) == cml_avail.shape[0] and
        len(lat_a) == cml_avail.shape[0] and
        len(lon_b) == cml_avail.shape[0] and
        len(lat_b) == cml_avail.shape[0]
    ):
        cml_mask = (
            np.isfinite(lon_a) & np.isfinite(lat_a) &
            np.isfinite(lon_b) & np.isfinite(lat_b)
        )
        cml_avail_f = cml_avail[cml_mask]
        print(f"[CML] Kept {cml_avail_f.shape[0]}/{cml_avail.shape[0]} links with valid coordinates")
    else:
        cml_avail_f = cml_avail
        print("[CML] No coordinate filtering")

    if radar_valid_coords is not None and len(radar_valid_coords) == radar_avail.shape[0]:
        radar_avail_f = radar_avail[radar_valid_coords]
        print(f"[RADAR] Kept {radar_avail_f.shape[0]}/{radar_avail.shape[0]} entries with valid coordinates")
    else:
        radar_avail_f = radar_avail
        print("[RADAR] No coordinate filtering")

    gauge_avail_f = gauge_avail
    print(f"[GAUGE] Kept {gauge_avail_f.shape[0]} stations")

    # Histograms.
    print("\n" + "=" * 70)
    print("COMPUTING 2D HISTOGRAMS")
    print("=" * 70)

    H_cml, n_zero_cml = compute_2d_histogram(cml_avail_f)
    H_radar, n_zero_radar = compute_2d_histogram(radar_avail_f)
    H_gauge, n_zero_gauge = compute_2d_histogram(gauge_avail_f)

    print(f"[CML]   H shape={H_cml.shape}, max count={H_cml.max()}, zero-availability entries masked={n_zero_cml}")
    print(f"[RADAR] H shape={H_radar.shape}, max count={H_radar.max()}, zero-availability entries masked={n_zero_radar}")
    print(f"[GAUGE] H shape={H_gauge.shape}, max count={H_gauge.max()}, zero-availability entries masked={n_zero_gauge}")

    # Plotting.
    print("\n" + "=" * 70)
    print("PLOTTING")
    print("=" * 70)

    cmap = make_cmap()

    # Shared scale for CML and Radar.
    vmax_cml_radar = max(int(H_cml.max()), int(H_radar.max()), 1)
    norm_cml_radar = make_count_norm(vmax_cml_radar)

    # Independent scale for gauges.
    vmax_gauge = max(int(H_gauge.max()), 1)
    norm_gauge = make_count_norm(vmax_gauge)

    print(f"[COLOR SCALE] CML/Radar shared vmax = {vmax_cml_radar}")
    print(f"[COLOR SCALE] Gauge independent vmax = {vmax_gauge}")

    panels = [
        {
            "key": "cml",
            "H": H_cml,
            "time": cml_time,
            "norm": norm_cml_radar,
            "annotation": "",
            "cbar_label": "No. of links",
        },
        {
            "key": "radar",
            "H": H_radar,
            "time": radar_time,
            "norm": norm_cml_radar,
            "annotation": "",
            "cbar_label": "No. of radar entries",
        },
        {
            "key": "gauge",
            "H": H_gauge,
            "time": gauge_time,
            "norm": norm_gauge,
            "annotation": "",
            "cbar_label": "No. of stations",
        },
    ]

    saved = []
    for panel in panels:
        png_path, pdf_path = save_sensor_plot(
            key=panel["key"],
            H=panel["H"],
            time_axis=panel["time"],
            cmap=cmap,
            norm=panel["norm"],
            annotation=panel["annotation"],
            cbar_label=panel["cbar_label"],
        )
        saved.extend([png_path, pdf_path])
        print(f"[SAVED] {png_path}")
        print(f"[SAVED] {pdf_path}")

    if SAVE_COMBINED:
        png_path, pdf_path = save_combined_plot(panels, cmap)
        saved.extend([png_path, pdf_path])
        print(f"[SAVED] {png_path}")
        print(f"[SAVED] {pdf_path}")

    print(f"\n[TOTAL] {time.perf_counter() - total_start:.2f}s")
    print("Saved files:")
    for path in saved:
        print(f"  - {path}")
