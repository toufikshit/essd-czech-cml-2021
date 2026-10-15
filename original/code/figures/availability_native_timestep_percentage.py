#!/usr/bin/env python3
"""
Native-timestep network availability percentage

For each sensor, compute at its own native temporal resolution the percentage
of entities with valid data at every timestamp, then plot all sensors on one
common April--October 2021 calendar axis.

Definition
----------
availability(t) = 100 * N_valid(t) / N_total

CML:
  Uses the preprocessed CML signal file, not the rainfall product.
  A physical CML link is valid at time t when at least one sublink has BOTH
  transmitted and received signal values finite.

Radar:
  Percentage of all radar path series with finite data at each native timestep.

Gauge / temperature / humidity:
  Percentage of all stations with finite data at each native timestep.

No sensor is temporally resampled for the availability calculation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import xarray as xr


DEFAULT_CML_FILE = (
    "/home/jovyan/ESSD_publication_dataset/"
    "CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc"
)
DEFAULT_RADAR_FILE = (
    "/home/jovyan/ESSD_publication_dataset/"
    "CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc"
)
DEFAULT_GAUGE_FILE = (
    "/home/jovyan/ESSD_publication_dataset/"
    "CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc"
)
DEFAULT_TEMPERATURE_FILE = (
    "/home/jovyan/opendata.chmi.cz/processed/"
    "CHMI_10MIN_TEMPERATURE_2021.nc"
)
DEFAULT_HUMIDITY_FILE = (
    "/home/jovyan/opendata.chmi.cz/processed/"
    "CHMI_10MIN_HUMIDITY_2021.nc"
)

DEFAULT_START = "2021-04-01"
DEFAULT_END = "2021-10-31"

ROW_ORDER = ("CML", "Radar", "Gauge", "Temperature", "Humidity")
ROW_LABELS = (
    "CML",
    "Radar",
    "Rain gauge",
    "Air temperature",
    "Relative humidity",
)


@dataclass(frozen=True)
class SensorPlan:
    sensor: str
    time_dim: str
    entity_dim: str | None
    variables: tuple[str, ...]
    n_entities: int


def find_dim(ds: xr.Dataset, candidates: Sequence[str]) -> str | None:
    for name in candidates:
        if name in ds.dims:
            return name
    return None


def find_var(ds: xr.Dataset, candidates: Sequence[str]) -> str | None:
    for name in candidates:
        if name in ds.data_vars or name in ds.coords:
            return name
    return None


def require_file(path: str | Path, label: str) -> Path:
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"{label} file does not exist: {p}")
    return p


def infer_native_step(times: pd.DatetimeIndex) -> pd.Timedelta:
    diffs = np.diff(times.asi8)
    diffs = diffs[diffs > 0]
    if diffs.size == 0:
        raise ValueError("Cannot infer native time step.")
    return pd.to_timedelta(int(np.median(diffs)), unit="ns")


def select_time_indices(
    times: pd.DatetimeIndex,
    start: pd.Timestamp,
    end_inclusive: pd.Timestamp,
) -> tuple[int, int]:
    start_ts = pd.Timestamp(start)
    end_ts = pd.Timestamp(end_inclusive)
    if end_ts == end_ts.normalize():
        end_exclusive = end_ts + pd.Timedelta(days=1)
    else:
        end_exclusive = end_ts + pd.Timedelta(nanoseconds=1)

    i0 = int(np.searchsorted(times.values, np.datetime64(start_ts), side="left"))
    i1 = int(np.searchsorted(times.values, np.datetime64(end_exclusive), side="left"))
    return i0, i1


def build_cml_plan(ds: xr.Dataset) -> SensorPlan:
    time_dim = find_dim(ds, ("time", "t"))
    entity_dim = find_dim(ds, ("cml_id", "cml", "link", "link_id"))
    if time_dim is None or entity_dim is None:
        raise ValueError(f"CML dimensions not recognized: {list(ds.dims)}")

    tx = find_var(ds, ("tx", "tsl", "transmitted_signal_level"))
    rx = find_var(ds, ("rx", "rsl", "received_signal_level"))

    if tx is not None and rx is not None:
        variables = (tx, rx)
    else:
        signal = find_var(ds, ("trsl_qc", "trsl", "signal", "txrx"))
        if signal is None:
            raise ValueError(
                "No usable CML signal pair or combined signal variable found. "
                f"Available: {list(ds.data_vars)}"
            )
        variables = (signal,)

    return SensorPlan(
        "CML", time_dim, entity_dim, variables, int(ds.sizes[entity_dim])
    )


def build_generic_plan(
    ds: xr.Dataset,
    sensor: str,
    variable_candidates: Sequence[str],
    entity_candidates: Sequence[str],
) -> SensorPlan:
    time_dim = find_dim(ds, ("time", "t"))
    if time_dim is None:
        raise ValueError(f"{sensor}: time dimension not found.")

    var = find_var(ds, variable_candidates)
    if var is None:
        raise ValueError(
            f"{sensor}: no usable variable found. Available: {list(ds.data_vars)}"
        )

    entity_dim = find_dim(ds, entity_candidates)
    if entity_dim is not None and entity_dim not in ds[var].dims:
        entity_dim = None

    if entity_dim is None:
        others = [d for d in ds[var].dims if d != time_dim]
        entity_dim = others[0] if others else None

    n_entities = int(ds.sizes[entity_dim]) if entity_dim else 1

    return SensorPlan(sensor, time_dim, entity_dim, (var,), n_entities)


def build_radar_plan(ds: xr.Dataset) -> SensorPlan:
    return build_generic_plan(
        ds,
        "Radar",
        ("rr", "rain_rate", "rainfall", "precipitation", "precip", "R", "R_radolan"),
        ("cml_id", "cml", "link", "link_id", "id"),
    )


def build_gauge_plan(ds: xr.Dataset) -> SensorPlan:
    return build_generic_plan(
        ds,
        "Gauge",
        ("rr", "rain_rate", "rainfall", "precipitation", "precip", "P", "R"),
        ("station", "station_id", "WSI", "id", "gauge", "site"),
    )


def build_temperature_plan(ds: xr.Dataset) -> SensorPlan:
    return build_generic_plan(
        ds,
        "Temperature",
        ("T", "temperature", "air_temperature", "temp", "ta"),
        ("station", "station_id", "WSI", "id", "site"),
    )


def build_humidity_plan(ds: xr.Dataset) -> SensorPlan:
    return build_generic_plan(
        ds,
        "Humidity",
        ("H", "RH", "relative_humidity", "humidity", "relativeHumidity"),
        ("station", "station_id", "WSI", "id", "site"),
    )


def cml_valid_entity_time(
    ds: xr.Dataset,
    plan: SensorPlan,
    time_slice: slice,
) -> np.ndarray:
    indexer = {plan.time_dim: time_slice}

    if len(plan.variables) == 2:
        valid = (
            ds[plan.variables[0]].isel(indexer).notnull()
            & ds[plan.variables[1]].isel(indexer).notnull()
        )
    else:
        valid = ds[plan.variables[0]].isel(indexer).notnull()

    extra_dims = [
        d for d in valid.dims
        if d not in (plan.entity_dim, plan.time_dim)
    ]
    if extra_dims:
        valid = valid.any(dim=extra_dims)

    valid = valid.transpose(plan.entity_dim, plan.time_dim)
    return np.asarray(valid.load().values, dtype=bool)


def generic_valid_entity_time(
    ds: xr.Dataset,
    plan: SensorPlan,
    time_slice: slice,
) -> np.ndarray:
    indexer = {plan.time_dim: time_slice}
    valid = ds[plan.variables[0]].isel(indexer).notnull()

    if plan.entity_dim is None:
        extra_dims = [d for d in valid.dims if d != plan.time_dim]
        if extra_dims:
            valid = valid.any(dim=extra_dims)
        valid = valid.expand_dims({"__entity__": [0]})
        entity_dim = "__entity__"
    else:
        entity_dim = plan.entity_dim
        extra_dims = [
            d for d in valid.dims
            if d not in (entity_dim, plan.time_dim)
        ]
        if extra_dims:
            valid = valid.any(dim=extra_dims)

    valid = valid.transpose(entity_dim, plan.time_dim)
    return np.asarray(valid.load().values, dtype=bool)


def compute_native_availability(
    filepath: Path,
    plan_builder,
    start: pd.Timestamp,
    end: pd.Timestamp,
    chunk_time: int,
) -> pd.DataFrame:
    with xr.open_dataset(filepath, decode_times=True, cache=False) as ds:
        plan = plan_builder(ds)
        times_all = pd.DatetimeIndex(ds[plan.time_dim].values)

        if not times_all.is_monotonic_increasing:
            raise ValueError(f"{plan.sensor}: non-monotonic time axis.")

        i0, i1 = select_time_indices(times_all, start, end)
        if i1 <= i0:
            raise ValueError(f"{plan.sensor}: no timestamps in requested period.")

        native_times = times_all[i0:i1]
        step = infer_native_step(native_times)

        print(
            f"[{plan.sensor}] entities={plan.n_entities:,} | "
            f"native={step / pd.Timedelta(minutes=1):g} min | "
            f"{native_times[0]} -> {native_times[-1]} | "
            f"Ntime={len(native_times):,}"
        )

        n_valid = np.empty(len(native_times), dtype=np.int32)

        out_pos = 0
        for s in range(i0, i1, chunk_time):
            e = min(s + chunk_time, i1)
            sl = slice(s, e)

            if plan.sensor == "CML":
                valid = cml_valid_entity_time(ds, plan, sl)
            else:
                valid = generic_valid_entity_time(ds, plan, sl)

            counts = valid.sum(axis=0, dtype=np.int32)
            n_valid[out_pos:out_pos + len(counts)] = counts
            out_pos += len(counts)

        pct = 100.0 * n_valid.astype(float) / float(plan.n_entities)

        return pd.DataFrame(
            {
                "time": native_times,
                "availability_percent": pct,
                "n_valid": n_valid,
                "n_total": plan.n_entities,
            }
        )


def time_edges(times: pd.DatetimeIndex) -> np.ndarray:
    x = mdates.date2num(times.to_pydatetime())
    dx = float(np.median(np.diff(x)))
    return np.concatenate(
        ([x[0] - dx / 2.0], x[:-1] + np.diff(x) / 2.0, [x[-1] + dx / 2.0])
    )


def plot_native_availability(
    results: dict[str, pd.DataFrame],
    start: pd.Timestamp,
    end: pd.Timestamp,
    output_png: Path,
    output_pdf: Path,
    cmap_name: str,
    dpi: int,
) -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 9,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "savefig.facecolor": "white",
        }
    )

    fig, ax = plt.subplots(figsize=(12.2, 3.6))

    cmap = mpl.colormaps[cmap_name].copy()
    cmap.set_bad("#eeeeee")
    norm = Normalize(vmin=0.0, vmax=100.0)

    mesh = None

    for row_i, sensor in enumerate(ROW_ORDER):
        frame = results[sensor]
        times = pd.DatetimeIndex(frame["time"])
        vals = frame["availability_percent"].to_numpy(float)

        mesh = ax.pcolormesh(
            time_edges(times),
            np.array([row_i, row_i + 1], dtype=float),
            vals[np.newaxis, :],
            cmap=cmap,
            norm=norm,
            shading="flat",
            edgecolors="none",
            antialiased=False,
            rasterized=True,
        )

    ax.set_ylim(len(ROW_ORDER), 0)
    ax.set_yticks(np.arange(len(ROW_ORDER)) + 0.5)
    ax.set_yticklabels(ROW_LABELS)
    ax.tick_params(axis="y", length=0, pad=6)

    start_plot = pd.Timestamp(start).normalize()
    end_plot = pd.Timestamp(end).normalize() + pd.Timedelta(days=1)

    ax.set_xlim(
        mdates.date2num(start_plot.to_pydatetime()),
        mdates.date2num(end_plot.to_pydatetime()),
    )

    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.tick_params(axis="x", length=3, pad=4)

    for d in pd.date_range(start_plot.replace(day=1), end_plot, freq="MS"):
        ax.axvline(
            mdates.date2num(d.to_pydatetime()),
            color="white",
            linewidth=0.7,
            zorder=4,
        )

    for y in range(len(ROW_ORDER) + 1):
        ax.axhline(y, color="white", linewidth=1.0, zorder=4)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    midpoint = start_plot + (end_plot - start_plot) / 2
    ax.text(
        mdates.date2num(midpoint.to_pydatetime()),
        1.04,
        str(start_plot.year),
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="bottom",
        fontsize=9,
        color="0.35",
        clip_on=False,
    )

    cbar = fig.colorbar(
        mesh,
        ax=ax,
        orientation="horizontal",
        fraction=0.08,
        pad=0.22,
        aspect=45,
        ticks=[0, 20, 40, 60, 80, 100],
    )
    cbar.set_label("Entities with valid data at native timestep (%)")
    cbar.ax.tick_params(length=3, width=0.8, pad=2)
    cbar.outline.set_linewidth(0.8)

    fig.subplots_adjust(left=0.16, right=0.985, top=0.87, bottom=0.28)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=dpi, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(output_pdf, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)

    print(f"[SAVED] {output_png}")
    print(f"[SAVED] {output_pdf}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument("--cml-file", default=DEFAULT_CML_FILE)
    p.add_argument("--radar-file", default=DEFAULT_RADAR_FILE)
    p.add_argument("--gauge-file", default=DEFAULT_GAUGE_FILE)
    p.add_argument("--temperature-file", default=DEFAULT_TEMPERATURE_FILE)
    p.add_argument("--humidity-file", default=DEFAULT_HUMIDITY_FILE)

    p.add_argument("--start", default=DEFAULT_START)
    p.add_argument("--end", default=DEFAULT_END)

    p.add_argument("--chunk-time", type=int, default=10000)
    p.add_argument(
        "--out-dir",
        default=(
            "/home/jovyan/ESSD_publication_dataset/"
            "availability_native_timestep"
        ),
    )
    p.add_argument(
        "--basename",
        default="availability_native_timestep_percentage",
    )
    p.add_argument("--cmap", default="RdYlGn")
    p.add_argument("--dpi", type=int, default=600)

    return p.parse_args()


def main() -> None:
    args = parse_args()

    start = pd.Timestamp(args.start)
    end = pd.Timestamp(args.end)

    if end < start:
        raise ValueError("--end must not be earlier than --start")
    if args.chunk_time <= 0:
        raise ValueError("--chunk-time must be > 0")

    jobs = [
        ("CML", require_file(args.cml_file, "CML"), build_cml_plan),
        ("Radar", require_file(args.radar_file, "Radar"), build_radar_plan),
        ("Gauge", require_file(args.gauge_file, "Gauge"), build_gauge_plan),
        (
            "Temperature",
            require_file(args.temperature_file, "Temperature"),
            build_temperature_plan,
        ),
        (
            "Humidity",
            require_file(args.humidity_file, "Humidity"),
            build_humidity_plan,
        ),
    ]

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 92)
    print("NATIVE-TIMESTEP SENSOR AVAILABILITY")
    print("=" * 92)
    print(f"Period     : {start:%Y-%m-%d} -> {end:%Y-%m-%d}")
    print("Definition : 100 * valid entities / total entities at every native timestep")
    print("CML rule   : link valid if >=1 sublink has both TX and RX finite")
    print("=" * 92)

    results = {}

    for sensor, filepath, builder in jobs:
        frame = compute_native_availability(
            filepath,
            builder,
            start,
            end,
            args.chunk_time,
        )
        results[sensor] = frame

        csv_path = out_dir / f"{args.basename}_{sensor.lower()}.csv"
        frame.to_csv(csv_path, index=False, float_format="%.6f")
        print(f"[SAVED] {csv_path}")

        a = frame["availability_percent"].to_numpy(float)
        print(
            f"[{sensor}] min={np.nanmin(a):.2f}% | "
            f"median={np.nanmedian(a):.2f}% | "
            f"mean={np.nanmean(a):.2f}% | "
            f"max={np.nanmax(a):.2f}%"
        )

    png_path = out_dir / f"{args.basename}.png"
    pdf_path = out_dir / f"{args.basename}.pdf"

    plot_native_availability(
        results,
        start,
        end,
        png_path,
        pdf_path,
        args.cmap,
        args.dpi,
    )


if __name__ == "__main__":
    main()
