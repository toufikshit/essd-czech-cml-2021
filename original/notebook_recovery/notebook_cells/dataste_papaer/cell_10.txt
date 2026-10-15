#!/usr/bin/env python3
"""
Publication ready monthly availability heatmap
==============================================

This script creates one compact monthly availability figure for:
    - CML
    - Radar
    - Rain gauges

Output:
    availability_plots/availability_monthly_publication.png
    availability_plots/availability_monthly_publication.pdf

Plot design:
    - rows: CML, Radar, Gauge
    - columns: months
    - cell color: overall monthly availability in %
    - colorbar: 0 to 100 %
    - NaN cells are white
    - 0 % cells are white
    - no text is written inside NaN or 0 % cells
    - month labels are Jan, Feb, Mar, ...
    - full period shown: 2018 to 2021
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# ============================================================
# SETTINGS
# ============================================================
CML_FILE = "ESSD_publication_dataset/CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc"
RADAR_FILE = "ESSD_publication_dataset/CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc"
GAUGE_FILE = "ESSD_publication_dataset/CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc"

OUT_DIR = Path("ESSD_publication_dataset")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Monthly availability
FREQ = "M"

# Full period to display in the paper
START_DATE = "2021-01-01"
END_DATE = "2021-12-31"

# Plot style
DPI = 300
FIGSIZE = (14.5, 2.9)

ANNOTATE_VALUES = True
FONT_CELL = 8

CMAP_NAME = "YlGn"


# ============================================================
# GENERAL HELPERS
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


def reduce_valid_to_entity_time(finite, dims, entity_dim, time_dim, reduce_mode="any"):
    """
    Reduce a boolean array to shape (n_entity, n_time).

    Parameters
    ----------
    finite : ndarray
        Boolean array with the same dimension order as dims.
    dims : list[str]
        Dimension names of finite.
    entity_dim : str
        Entity dimension, for example cml_id or station.
    time_dim : str
        Time dimension.
    reduce_mode : {"any", "all"}
        How to reduce extra dimensions such as sublink_id.

    Returns
    -------
    valid : ndarray
        Boolean array with shape (n_entity, n_time).
    """
    if entity_dim not in dims or time_dim not in dims:
        raise ValueError(
            f"Dims {dims} do not contain entity_dim={entity_dim!r} "
            f"and time_dim={time_dim!r}"
        )

    extra_axes = tuple(
        i for i, d in enumerate(dims)
        if d not in (entity_dim, time_dim)
    )

    if extra_axes:
        if reduce_mode == "any":
            finite = np.any(finite, axis=extra_axes)
        elif reduce_mode == "all":
            finite = np.all(finite, axis=extra_axes)
        else:
            raise ValueError("reduce_mode must be 'any' or 'all'")

    remaining = [d for d in dims if d in (entity_dim, time_dim)]

    # Make sure output is entity x time
    if remaining == [entity_dim, time_dim]:
        valid = finite
    elif remaining == [time_dim, entity_dim]:
        valid = finite.T
    else:
        raise ValueError(f"Unexpected remaining dimension order: {remaining}")

    return np.ascontiguousarray(valid, dtype=bool)


def compute_monthly_overall_availability(valid, time_vals, freq="M"):
    """
    Compute one monthly availability percentage for a sensor network.

    Availability is computed over all entities and all time steps in each month.

    Parameters
    ----------
    valid : ndarray
        Boolean array with shape (n_entity, n_time).
    time_vals : array-like
        Time coordinate.
    freq : str
        Grouping frequency. Use "M" for month.

    Returns
    -------
    monthly : ndarray
        Availability in percent for each month.
    month_axis : pandas.DatetimeIndex
        Month start dates.
    """
    time_vals = pd.DatetimeIndex(time_vals)
    periods = time_vals.to_period(freq)

    groups = pd.Series(np.arange(len(time_vals))).groupby(periods).indices
    periods_sorted = sorted(groups.keys())

    monthly = np.full(len(periods_sorted), np.nan, dtype=float)

    for i, p in enumerate(periods_sorted):
        idx = np.asarray(groups[p], dtype=int)
        monthly[i] = valid[:, idx].mean() * 100.0

    month_axis = pd.DatetimeIndex([p.start_time for p in periods_sorted])
    return monthly, month_axis


def make_complete_month_axis(start_date, end_date):
    """Create a complete monthly time axis from start_date to end_date."""
    return pd.date_range(
        pd.Timestamp(start_date).replace(day=1),
        pd.Timestamp(end_date).replace(day=1),
        freq="MS",
    )


def align_to_complete_axis(series_dict, complete_time):
    """
    Align each sensor monthly series to the complete month axis.

    Missing months become NaN.
    """
    aligned = {}

    for key, (values, time_axis) in series_dict.items():
        s = pd.Series(values, index=pd.DatetimeIndex(time_axis))
        s = s.reindex(complete_time)
        aligned[key] = s.values.astype(float)

    return aligned


# ============================================================
# CML AVAILABILITY
# ============================================================
def load_cml_monthly(filepath, freq="M"):
    """
    Load CML data and compute monthly availability.

    If both transmitted and received signal level variables are available,
    a valid CML sample requires both to be finite.
    If only one signal variable is found, that variable is used.
    """
    print("[CML] Loading dataset...")

    with xr.open_dataset(filepath) as ds:
        cml_dim = find_dim(ds, ["cml_id", "cml", "link", "link_id"])
        time_dim = find_dim(ds, ["time", "t"])

        if cml_dim is None or time_dim is None:
            raise ValueError(f"[CML] Could not find CML/time dims. Found dims: {list(ds.dims)}")

        tx_var = find_var(ds, ["tsl", "tx", "transmitted_signal_level"])
        rx_var = find_var(ds, ["rsl", "rx", "received_signal_level"])

        if tx_var is not None and rx_var is not None:
            tx, rx = xr.align(ds[tx_var], ds[rx_var], join="inner")
            finite = np.isfinite(tx.values) & np.isfinite(rx.values)
            dims = list(tx.dims)
            print(f"[CML] Using variables {tx_var!r} and {rx_var!r}")
        else:
            signal_var = find_var(ds, ["trsl", "signal", "txrx", "rsl", "rx", "tsl", "tx"])
            if signal_var is None:
                raise ValueError(f"[CML] Could not find signal variable. Found vars: {list(ds.data_vars)}")
            signal = ds[signal_var]
            finite = np.isfinite(signal.values)
            dims = list(signal.dims)
            print(f"[CML] Using variable {signal_var!r}")

        valid = reduce_valid_to_entity_time(
            finite=finite,
            dims=dims,
            entity_dim=cml_dim,
            time_dim=time_dim,
            reduce_mode="any",
        )

        time_vals = pd.DatetimeIndex(ds[time_dim].values)

        print(f"[CML] Valid mask shape: {valid.shape}")
        print(f"[CML] Time range: {time_vals.min()} -> {time_vals.max()}")

        monthly, month_axis = compute_monthly_overall_availability(
            valid=valid,
            time_vals=time_vals,
            freq=freq,
        )

    return monthly, month_axis


# ============================================================
# RADAR AVAILABILITY
# ============================================================
def load_radar_monthly(filepath, freq="M"):
    """Load path averaged radar data and compute monthly availability."""
    print("[RADAR] Loading dataset...")

    with xr.open_dataset(filepath) as ds:
        time_dim = find_dim(ds, ["time", "t"])
        if time_dim is None:
            raise ValueError(f"[RADAR] Could not find time dim. Found dims: {list(ds.dims)}")

        rr_var = find_var(
            ds,
            ["rr", "rain_rate", "rainfall", "precip", "precipitation", "R", "R_radolan"],
        )
        if rr_var is None:
            raise ValueError(f"[RADAR] Could not find radar variable. Found vars: {list(ds.data_vars)}")

        rr = ds[rr_var]
        dims = list(rr.dims)

        cml_dim = find_dim(ds, ["cml_id", "cml", "link", "link_id", "id"])

        if cml_dim is not None and cml_dim in rr.dims:
            entity_dim = cml_dim
        else:
            other_dims = [d for d in rr.dims if d != time_dim]
            entity_dim = other_dims[0] if other_dims else None

        finite = np.isfinite(rr.values)

        if entity_dim is None:
            valid = finite.reshape(1, -1)
        else:
            valid = reduce_valid_to_entity_time(
                finite=finite,
                dims=dims,
                entity_dim=entity_dim,
                time_dim=time_dim,
                reduce_mode="any",
            )

        time_vals = pd.DatetimeIndex(ds[time_dim].values)

        print(f"[RADAR] Using variable {rr_var!r}")
        print(f"[RADAR] Valid mask shape: {valid.shape}")
        print(f"[RADAR] Time range: {time_vals.min()} -> {time_vals.max()}")

        monthly, month_axis = compute_monthly_overall_availability(
            valid=valid,
            time_vals=time_vals,
            freq=freq,
        )

    return monthly, month_axis


# ============================================================
# GAUGE AVAILABILITY
# ============================================================
def load_gauge_monthly(filepath, freq="M"):
    """Load rain gauge data and compute monthly availability."""
    print("[GAUGE] Loading dataset...")

    with xr.open_dataset(filepath) as ds:
        time_dim = find_dim(ds, ["time", "t"])
        if time_dim is None:
            raise ValueError(f"[GAUGE] Could not find time dim. Found dims: {list(ds.dims)}")

        rr_var = find_var(
            ds,
            ["rr", "rain_rate", "rainfall", "precip", "precipitation", "P", "R"],
        )
        if rr_var is None:
            raise ValueError(f"[GAUGE] Could not find gauge variable. Found vars: {list(ds.data_vars)}")

        rr = ds[rr_var]
        dims = list(rr.dims)

        station_dim = find_dim(ds, ["station", "station_id", "id", "gauge", "site", "WSI"])
        if station_dim is not None and station_dim not in rr.dims:
            station_dim = None

        if station_dim is None:
            other_dims = [d for d in rr.dims if d != time_dim]
            station_dim = other_dims[0] if other_dims else None

        finite = np.isfinite(rr.values)

        if station_dim is None:
            valid = finite.reshape(1, -1)
        else:
            valid = reduce_valid_to_entity_time(
                finite=finite,
                dims=dims,
                entity_dim=station_dim,
                time_dim=time_dim,
                reduce_mode="any",
            )

        time_vals = pd.DatetimeIndex(ds[time_dim].values)

        print(f"[GAUGE] Using variable {rr_var!r}")
        print(f"[GAUGE] Valid mask shape: {valid.shape}")
        print(f"[GAUGE] Time range: {time_vals.min()} -> {time_vals.max()}")

        monthly, month_axis = compute_monthly_overall_availability(
            valid=valid,
            time_vals=time_vals,
            freq=freq,
        )

    return monthly, month_axis


# ============================================================
# PLOTTING
# ============================================================
def plot_publication_heatmap(data_matrix, month_axis, row_labels, out_dir):
    """
    Plot publication ready availability heatmap.

    NaN and 0 values are both masked and plotted as white.
    Only positive finite values are annotated.
    """
    fig, ax = plt.subplots(figsize=FIGSIZE)

    # Mask missing and zero availability.
    data_plot = np.ma.masked_where(
        (~np.isfinite(data_matrix)) | (data_matrix <= 0),
        data_matrix,
    )

    cmap = plt.get_cmap(CMAP_NAME).copy()
    cmap.set_bad("white")

    im = ax.imshow(
        data_plot,
        aspect="auto",
        cmap=cmap,
        vmin=0,
        vmax=100,
        interpolation="nearest",
    )

    n_rows, n_cols = data_matrix.shape

    # Cell text: only write values > 0.
    if ANNOTATE_VALUES:
        for i in range(n_rows):
            for j in range(n_cols):
                val = data_matrix[i, j]
                if np.isfinite(val) and val > 0:
                    text_color = "black" if val < 70 else "white"
                    ax.text(
                        j,
                        i,
                        f"{val:.0f}",
                        ha="center",
                        va="center",
                        fontsize=FONT_CELL,
                        color=text_color,
                    )

    # Y axis
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=11)


    # X axis: Jan, Feb, ...
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(
        [d.strftime("%b") for d in month_axis],
        rotation=0,
        fontsize=8,
    )

    # White grid lines between cells
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Year separators and year labels
    years = np.array([d.year for d in month_axis])
    for j in range(n_cols - 1):
        if years[j] != years[j + 1]:
            ax.axvline(j + 0.5, color="0.45", lw=0.9)

    for year in sorted(set(years)):
        idx = np.where(years == year)[0]
        x_mid = (idx[0] + idx[-1]) / 2.0
        ax.text(
            x_mid,
            -0.85,
            str(year),
            ha="center",
            va="bottom",
            fontsize=10,
            transform=ax.transData,
        )

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Availability (%)", fontsize=10)
    cbar.set_ticks([0, 20, 40, 60, 80, 100])
    cbar.ax.tick_params(labelsize=8)

    # Clean frame
    ax.set_title("")
    ax.tick_params(axis="both", direction="out", length=3)

    plt.tight_layout()

    png_path = out_dir / "availability_monthly_publication.png"
    pdf_path = out_dir / "availability_monthly_publication.pdf"

    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"[SAVED] {png_path}")
    print(f"[SAVED] {pdf_path}")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("MONTHLY SENSOR AVAILABILITY SUMMARY")
    print("=" * 70)

    complete_time = make_complete_month_axis(START_DATE, END_DATE)

    cml_values, cml_time = load_cml_monthly(CML_FILE, FREQ)
    radar_values, radar_time = load_radar_monthly(RADAR_FILE, FREQ)
    gauge_values, gauge_time = load_gauge_monthly(GAUGE_FILE, FREQ)

    aligned = align_to_complete_axis(
        {
            "CML": (cml_values, cml_time),
            "Radar": (radar_values, radar_time),
            "Gauge": (gauge_values, gauge_time),
        },
        complete_time=complete_time,
    )

    row_labels = ["CML", "Radar", "Gauge"]

    data_matrix = np.vstack(
        [
            aligned["CML"],
            aligned["Radar"],
            aligned["Gauge"],
        ]
    )

    print("\nMonthly availability matrix (%):")
    print(pd.DataFrame(data_matrix, index=row_labels, columns=complete_time.strftime("%Y-%m")))

    plot_publication_heatmap(
        data_matrix=data_matrix,
        month_axis=complete_time,
        row_labels=row_labels,
        out_dir=OUT_DIR,
    )
