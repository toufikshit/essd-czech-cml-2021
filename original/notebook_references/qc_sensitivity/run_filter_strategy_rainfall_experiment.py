#!/usr/bin/env python3
"""
run_filter_strategy_rainfall_experiment.py
==========================================

End-to-end fair benchmark of TWO already-produced CML filtering strategies.

The experiment is intentionally designed so that the FILTER is the only
scientific difference:

    Legacy-filtered CML dataset
        -> identical rainfall estimation
        -> legacy rainfall product

    Strict per-link/per-sublink/per-observation filtered CML dataset
        -> identical rainfall estimation
        -> strict rainfall product

Then both rainfall products are evaluated against the SAME radar reference.

Two evaluations are produced:

1) COMMON-SAMPLE accuracy
   The exact same (cml_id, time_5min) positions are used for both products:
       finite legacy CML
       & finite strict CML
       & finite radar
       & radar > wet threshold
       & all values inside the selected plotting/evaluation range

   This gives identical N for both filters and isolates rainfall accuracy.

2) OWN-VALID-POPULATION accuracy
   Each filter is evaluated using its own finite CML samples against the SAME
   radar wet criterion. N may differ.

   This shows the accuracy-vs-availability tradeoff.

Rainfall-estimation configuration is fixed for BOTH filters:
    baseline : pycml_constant
    k-R      : ITU_2003
    WAA      : pastorek_2021

The estimation script is run with --skip-gridding, because gridding is not
needed for the filter benchmark and would only make the experiment slower.

Expected filtered inputs
------------------------
Legacy strategy:
    e.g. /home/jovyan/cml_preprocessed_v3.nc
    total-loss variable normally: trsl

Strict strategy:
    e.g. /home/jovyan/CML_S2020E2021_1MIN_CZ_QC_STRICT_MINIMAL.nc
    total-loss variable normally: trsl_qc

The script auto-detects trsl_qc -> trsl -> TX-RX for each filtered input.

Outputs
-------
<out-dir>/
    rainfall_legacy_filter.nc
    rainfall_strict_filter.nc
    rainfall_legacy_filter_log.csv
    rainfall_strict_filter_log.csv

    hexbin_common_samples.png
    hexbin_common_samples.pdf

    hexbin_own_valid_population.png
    hexbin_own_valid_population.pdf

    filter_strategy_metrics.csv
    filter_strategy_experiment_summary.csv

The metrics are exact over all selected samples.
Only the points used to draw the hexbins are randomly subsampled.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shlex
import subprocess
import sys
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
from matplotlib.colors import LogNorm

try:
    from tqdm import tqdm
except Exception:
    def tqdm(x, **kwargs):
        return x


# =============================================================================
# Constants: fixed rainfall framework for the filter experiment
# =============================================================================
BASELINE_METHOD = "pycml_constant"
KR_MODEL = "ITU_2003"
WAA_METHOD = "pastorek_2021"

# Keep these identical for both runs.
BASELINE_PAD_HOURS = 48.0
MIN_LINK_KM = 0.5
MIN_VALID_SAMPLES = 100
N_AVERAGE_LAST_DRY = 1

# Same retrieval wet/dry threshold for both estimation runs.
# This reproduces the default of the supplied V2 estimation framework.
RADAR_WET_THR_ESTIMATION = 0.00833

# Pastorek parameters used by the supplied estimation framework.
PASTOREK_A_MAX = 14.0
PASTOREK_ZETA = 0.55
PASTOREK_D = 0.1


# =============================================================================
# Small utilities
# =============================================================================
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def run_command(cmd: Sequence[str], dry_run: bool = False) -> None:
    print("\n" + "=" * 100)
    print("RUNNING")
    print(" ".join(shlex.quote(str(x)) for x in cmd))
    print("=" * 100, flush=True)

    if dry_run:
        return

    subprocess.run(list(map(str, cmd)), check=True)


def find_dim(obj: xr.DataArray | xr.Dataset, candidates: Sequence[str], required: bool = True) -> Optional[str]:
    dims = list(obj.dims)

    for candidate in candidates:
        if candidate in dims:
            return candidate

    for d in dims:
        dl = d.lower()
        if any(candidate.lower() in dl for candidate in candidates):
            return d

    if required:
        raise KeyError(f"Could not find dimension matching {list(candidates)}. Available dims: {dims}")
    return None


def detect_trsl_var(path: str) -> Optional[str]:
    """
    Detect the preferred total-loss variable.

    Returns:
        "trsl_qc", "trsl", or None.

    None tells the estimation pipeline to fall back to TX-RX.
    """
    with xr.open_dataset(path) as ds:
        if "trsl_qc" in ds.variables:
            return "trsl_qc"
        if "trsl" in ds.variables:
            return "trsl"
        if "tx" in ds.variables and "rx" in ds.variables:
            return None

    raise KeyError(
        f"{path} has neither trsl_qc, trsl, nor both tx and rx."
    )


def coord_values_as_str(ds: xr.Dataset, dim: str) -> np.ndarray:
    if dim in ds.coords:
        return np.asarray(ds[dim].values).astype(str)
    return np.arange(ds.sizes[dim]).astype(str)


def resolve_coord_index(ds: xr.Dataset, dim: str, value: str) -> int:
    values = [str(v) for v in ds[dim].values]
    if value not in values:
        raise ValueError(
            f"Requested {value!r} not found on {dim!r}. Available: {values}"
        )
    return values.index(value)


def select_rainfall_combo(
    ds: xr.Dataset,
    variable: str,
    baseline_method: str,
    kr_model: str,
    waa_method: str,
) -> xr.DataArray:
    if variable not in ds.data_vars:
        raise KeyError(
            f"Rain variable {variable!r} not found. Available: {list(ds.data_vars)}"
        )

    da = ds[variable]
    selections = {}

    for candidates, value in [
        (["baseline_method", "baseline"], baseline_method),
        (["kr_model", "kr", "itu"], kr_model),
        (["waa_method", "waa"], waa_method),
    ]:
        dim = find_dim(da, candidates, required=False)
        if dim is not None:
            selections[dim] = resolve_coord_index(ds, dim, value)

    if selections:
        da = da.isel(selections)

    return da


def common_order_preserving(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    second_set = set(second.tolist())
    return np.asarray([x for x in first.tolist() if x in second_set], dtype=str)


# =============================================================================
# Exact streamed metrics
# =============================================================================
def init_stats() -> Dict[str, float]:
    return {
        "N": 0,
        "sum_c": 0.0,
        "sum_r": 0.0,
        "sum_c2": 0.0,
        "sum_r2": 0.0,
        "sum_cr": 0.0,
        "sum_diff": 0.0,
        "sum_abs_diff": 0.0,
        "sum_diff2": 0.0,
    }


def update_stats(stats: Dict[str, float], cml: np.ndarray, radar: np.ndarray) -> None:
    if cml.size == 0:
        return

    c = np.asarray(cml, dtype="float64")
    r = np.asarray(radar, dtype="float64")
    d = c - r

    stats["N"] += int(c.size)
    stats["sum_c"] += float(np.sum(c))
    stats["sum_r"] += float(np.sum(r))
    stats["sum_c2"] += float(np.sum(c * c))
    stats["sum_r2"] += float(np.sum(r * r))
    stats["sum_cr"] += float(np.sum(c * r))
    stats["sum_diff"] += float(np.sum(d))
    stats["sum_abs_diff"] += float(np.sum(np.abs(d)))
    stats["sum_diff2"] += float(np.sum(d * d))


def finalize_stats(stats: Dict[str, float]) -> Dict[str, float]:
    n = int(stats["N"])

    if n == 0:
        return {
            "N": 0,
            "MBE": np.nan,
            "MAE": np.nan,
            "RMSE": np.nan,
            "PCC": np.nan,
            "mean_cml": np.nan,
            "mean_radar": np.nan,
        }

    mean_c = stats["sum_c"] / n
    mean_r = stats["sum_r"] / n

    mbe = stats["sum_diff"] / n
    mae = stats["sum_abs_diff"] / n
    rmse = math.sqrt(max(stats["sum_diff2"] / n, 0.0))

    cov = stats["sum_cr"] / n - mean_c * mean_r
    var_c = stats["sum_c2"] / n - mean_c * mean_c
    var_r = stats["sum_r2"] / n - mean_r * mean_r

    if n >= 2 and var_c > 0.0 and var_r > 0.0:
        pcc = cov / math.sqrt(var_c * var_r)
    else:
        pcc = np.nan

    return {
        "N": n,
        "MBE": float(mbe),
        "MAE": float(mae),
        "RMSE": float(rmse),
        "PCC": float(pcc) if np.isfinite(pcc) else np.nan,
        "mean_cml": float(mean_c),
        "mean_radar": float(mean_r),
    }


# =============================================================================
# Uniform plotting reservoir using random priorities
# =============================================================================
class PriorityReservoir:
    """
    Keep a uniform random sample of at most `capacity` rows.

    Each incoming sample gets U(0,1). We retain the `capacity` smallest keys.
    This is equivalent to a uniform random sample without replacement over all
    samples seen so far and works efficiently block-by-block.
    """

    def __init__(self, capacity: int, n_columns: int, seed: int):
        self.capacity = max(0, int(capacity))
        self.n_columns = int(n_columns)
        self.rng = np.random.default_rng(seed)
        self.keys = np.empty(0, dtype="float64")
        self.values = np.empty((0, n_columns), dtype="float32")

    def add(self, *columns: np.ndarray) -> None:
        if self.capacity == 0:
            return

        if len(columns) != self.n_columns:
            raise ValueError(
                f"Expected {self.n_columns} columns, got {len(columns)}"
            )

        if len(columns) == 0 or np.asarray(columns[0]).size == 0:
            return

        arrays = [np.asarray(c, dtype="float32").reshape(-1) for c in columns]
        n = arrays[0].size

        if any(a.size != n for a in arrays):
            raise ValueError("Reservoir input columns have different lengths.")

        vals = np.column_stack(arrays)
        keys = self.rng.random(n)

        if self.values.size == 0:
            merged_values = vals
            merged_keys = keys
        else:
            merged_values = np.concatenate([self.values, vals], axis=0)
            merged_keys = np.concatenate([self.keys, keys])

        if merged_keys.size > self.capacity:
            keep = np.argpartition(merged_keys, self.capacity - 1)[: self.capacity]
            self.keys = merged_keys[keep]
            self.values = merged_values[keep]
        else:
            self.keys = merged_keys
            self.values = merged_values

    def arrays(self) -> Tuple[np.ndarray, ...]:
        if self.values.size == 0:
            return tuple(
                np.empty(0, dtype="float32")
                for _ in range(self.n_columns)
            )

        order = np.argsort(self.keys)
        vals = self.values[order]
        return tuple(vals[:, i] for i in range(self.n_columns))


# =============================================================================
# Block reader
# =============================================================================
@dataclass
class RainSource:
    name: str
    ds: xr.Dataset
    rain: xr.DataArray
    cml_dim: str
    time_dim: str
    cml_ids: np.ndarray
    times: pd.DatetimeIndex
    cml_pos: Dict[str, int]
    time_pos: Dict[np.datetime64, int]


def build_rain_source(
    name: str,
    path: str,
    rain_variable: str,
) -> RainSource:
    ds = xr.open_dataset(path)

    rain = select_rainfall_combo(
        ds,
        variable=rain_variable,
        baseline_method=BASELINE_METHOD,
        kr_model=KR_MODEL,
        waa_method=WAA_METHOD,
    )

    sublink_dim = find_dim(
        rain,
        ["sublink_id", "sublink"],
        required=False,
    )
    if sublink_dim is not None and sublink_dim in rain.dims:
        raise ValueError(
            f"{name}: {rain_variable} still contains {sublink_dim}. "
            "This experiment uses R_mean_5min for final link-level rainfall comparison."
        )

    cml_dim = find_dim(rain, ["cml_id", "cml"])
    time_dim = find_dim(rain, ["time_5min", "time"])

    cml_ids = coord_values_as_str(ds, cml_dim)
    times = pd.DatetimeIndex(ds[time_dim].values)

    return RainSource(
        name=name,
        ds=ds,
        rain=rain,
        cml_dim=cml_dim,
        time_dim=time_dim,
        cml_ids=cml_ids,
        times=times,
        cml_pos={v: i for i, v in enumerate(cml_ids)},
        time_pos={np.datetime64(v): i for i, v in enumerate(times.values)},
    )


def read_block(
    da: xr.DataArray,
    cml_dim: str,
    time_dim: str,
    cml_positions: np.ndarray,
    time_positions: np.ndarray,
) -> np.ndarray:
    """
    Load one small block in canonical shape (cml, time).
    Sequential indexing avoids accidental vectorized paired indexing.
    """
    block = da.isel({cml_dim: cml_positions})
    block = block.isel({time_dim: time_positions})
    block = block.transpose(cml_dim, time_dim)
    return np.asarray(block.values, dtype="float32")


# =============================================================================
# Plotting
# =============================================================================
def metric_box_text(metrics: Dict[str, float]) -> str:
    return (
        f"N={int(metrics['N']):,}\n"
        f"RMSE={metrics['RMSE']:.2f}\n"
        f"MAE={metrics['MAE']:.2f}\n"
        f"MBE={metrics['MBE']:.2f}\n"
        f"PCC={metrics['PCC']:.2f}"
    )


def make_two_panel_hexbin(
    radar_left: np.ndarray,
    cml_left: np.ndarray,
    radar_right: np.ndarray,
    cml_right: np.ndarray,
    metrics_left: Dict[str, float],
    metrics_right: Dict[str, float],
    title_left: str,
    title_right: str,
    main_title: str,
    out_png: Path,
    out_pdf: Path,
    min_rain: float,
    max_rain: float,
    gridsize: int,
) -> None:
    if radar_left.size == 0 or radar_right.size == 0:
        raise RuntimeError(
            f"Cannot create {main_title!r}: one plotting sample is empty."
        )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(10.2, 4.5),
        sharex=True,
        sharey=True,
    )

    hbs = []

    for ax, radar, cml, metrics, title in [
        (axes[0], radar_left, cml_left, metrics_left, title_left),
        (axes[1], radar_right, cml_right, metrics_right, title_right),
    ]:
        hb = ax.hexbin(
            radar,
            cml,
            gridsize=int(gridsize),
            mincnt=1,
            extent=(min_rain, max_rain, min_rain, max_rain),
            linewidths=0.0,
        )
        hbs.append(hb)

        ax.plot(
            [min_rain, max_rain],
            [min_rain, max_rain],
            linestyle="--",
            linewidth=1.0,
            color="0.25",
        )
        ax.set_xlim(min_rain, max_rain)
        ax.set_ylim(min_rain, max_rain)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"Radar rain rate [mm h$^{-1}$]")
        ax.set_title(title)

        ax.text(
            0.97,
            0.03,
            metric_box_text(metrics),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.5,
            bbox=dict(
                boxstyle="round,pad=0.25",
                fc="white",
                ec="0.75",
                alpha=0.92,
            ),
        )

    axes[0].set_ylabel(r"CML rain rate [mm h$^{-1}$]")

    # Force identical count normalization in both panels.
    max_count = max(
        float(np.nanmax(hb.get_array())) if hb.get_array().size else 1.0
        for hb in hbs
    )
    max_count = max(max_count, 1.0)

    norm = LogNorm(vmin=1.0, vmax=max_count)
    for hb in hbs:
        hb.set_norm(norm)

    fig.suptitle(main_title, fontsize=12)
    fig.subplots_adjust(right=0.88, top=0.88, wspace=0.12)

    cax = fig.add_axes([0.90, 0.18, 0.02, 0.66])
    cbar = fig.colorbar(hbs[-1], cax=cax)
    cbar.set_label("Count")

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Rainfall estimation stage
# =============================================================================
def build_estimation_command(
    estimator_script: str,
    cml_nc: str,
    trsl_var: Optional[str],
    radar_nc: str,
    output_nc: str,
    log_csv: str,
    tstart: str,
    tend: str,
    overwrite: bool,
) -> List[str]:
    cmd = [
        sys.executable,
        "-u",
        estimator_script,
        "--cml-nc",
        cml_nc,
        "--radar-nc",
        radar_nc,
        "--estimation-nc",
        output_nc,
        "--log-csv",
        log_csv,
        "--tstart",
        tstart,
        "--tend",
        tend,

        # Run estimation only. Gridding is irrelevant to this filter benchmark.
        "--skip-gridding",

        # Exactly one identical rainfall configuration for BOTH filters.
        "--baseline-methods",
        BASELINE_METHOD,
        "--kr-models",
        KR_MODEL,
        "--waa-methods",
        WAA_METHOD,

        # Explicit identical framework settings.
        "--baseline-pad-hours",
        str(BASELINE_PAD_HOURS),
        "--min-link-km",
        str(MIN_LINK_KM),
        "--min-valid-samples",
        str(MIN_VALID_SAMPLES),
        "--n-average-last-dry",
        str(N_AVERAGE_LAST_DRY),
        "--radar-wet-thr-5min",
        str(RADAR_WET_THR_ESTIMATION),
        "--pastorek-a-max",
        str(PASTOREK_A_MAX),
        "--pastorek-zeta",
        str(PASTOREK_ZETA),
        "--pastorek-d",
        str(PASTOREK_D),
    ]

    if trsl_var is not None:
        cmd.extend(["--trsl-var", trsl_var])

    if overwrite:
        cmd.append("--overwrite-estimation")

    return cmd


# =============================================================================
# Fair comparison stage
# =============================================================================
def compare_rainfall_products(args, legacy_rain_nc: str, strict_rain_nc: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    print("\n" + "=" * 100)
    print("FAIR FILTER COMPARISON")
    print("=" * 100)

    legacy = build_rain_source(
        "legacy",
        legacy_rain_nc,
        args.rain_variable,
    )
    strict = build_rain_source(
        "strict",
        strict_rain_nc,
        args.rain_variable,
    )

    try:
        if args.radar_var not in strict.ds.data_vars:
            raise KeyError(
                f"{args.radar_var!r} not found in strict rainfall file. "
                f"Available: {list(strict.ds.data_vars)}"
            )

        radar = strict.ds[args.radar_var]
        radar_cml_dim = find_dim(radar, ["cml_id", "cml"])
        radar_time_dim = find_dim(radar, ["time_5min", "time"])

        radar_ids = coord_values_as_str(strict.ds, radar_cml_dim)
        radar_times = pd.DatetimeIndex(strict.ds[radar_time_dim].values)
        radar_cml_pos = {v: i for i, v in enumerate(radar_ids)}
        radar_time_pos = {
            np.datetime64(v): i
            for i, v in enumerate(radar_times.values)
        }

        # ------------------------------------------------------------
        # Common links
        # ------------------------------------------------------------
        common_ids = common_order_preserving(
            strict.cml_ids,
            legacy.cml_ids,
        )
        common_ids = np.asarray(
            [cid for cid in common_ids if cid in radar_cml_pos],
            dtype=str,
        )

        if common_ids.size == 0:
            raise RuntimeError(
                "No common cml_id values across legacy, strict and radar."
            )

        # ------------------------------------------------------------
        # Common times
        # ------------------------------------------------------------
        legacy_time_set = set(legacy.time_pos.keys())
        radar_time_set = set(radar_time_pos.keys())

        common_times = np.asarray(
            [
                np.datetime64(t)
                for t in strict.times.values
                if (
                    np.datetime64(t) in legacy_time_set
                    and np.datetime64(t) in radar_time_set
                )
            ],
            dtype="datetime64[ns]",
        )

        if args.eval_tstart is not None:
            common_times = common_times[
                common_times >= np.datetime64(pd.Timestamp(args.eval_tstart))
            ]

        if args.eval_tend is not None:
            common_times = common_times[
                common_times <= np.datetime64(pd.Timestamp(args.eval_tend))
            ]

        if common_times.size == 0:
            raise RuntimeError(
                "No common evaluation timestamps remain."
            )

        # Position arrays
        legacy_cpos = np.asarray(
            [legacy.cml_pos[cid] for cid in common_ids],
            dtype=np.int64,
        )
        strict_cpos = np.asarray(
            [strict.cml_pos[cid] for cid in common_ids],
            dtype=np.int64,
        )
        radar_cpos = np.asarray(
            [radar_cml_pos[cid] for cid in common_ids],
            dtype=np.int64,
        )

        legacy_tpos = np.asarray(
            [legacy.time_pos[t] for t in common_times],
            dtype=np.int64,
        )
        strict_tpos = np.asarray(
            [strict.time_pos[t] for t in common_times],
            dtype=np.int64,
        )
        radar_tpos = np.asarray(
            [radar_time_pos[t] for t in common_times],
            dtype=np.int64,
        )

        n_links = common_ids.size
        n_times = common_times.size

        print(f"Common physical links : {n_links:,}")
        print(f"Common 5-min times    : {n_times:,}")
        print(
            f"Evaluation time       : "
            f"{pd.Timestamp(common_times[0])} -> {pd.Timestamp(common_times[-1])}"
        )
        print(
            f"Final wet criterion   : "
            f"{args.radar_var} > {args.radar_wet_thr:g} mm h-1"
        )
        print(
            f"Evaluation range      : "
            f"{args.min_rain:g} to {args.max_rain:g} mm h-1"
        )

        # Exact metric accumulators
        stats_common_legacy = init_stats()
        stats_common_strict = init_stats()

        stats_own_legacy = init_stats()
        stats_own_strict = init_stats()

        # Plot samples
        common_reservoir = PriorityReservoir(
            capacity=args.max_plot_points,
            n_columns=3,  # radar, legacy, strict
            seed=args.random_seed,
        )
        own_legacy_reservoir = PriorityReservoir(
            capacity=args.max_plot_points,
            n_columns=2,  # radar, legacy
            seed=args.random_seed + 1,
        )
        own_strict_reservoir = PriorityReservoir(
            capacity=args.max_plot_points,
            n_columns=2,  # radar, strict
            seed=args.random_seed + 2,
        )

        # Diagnostics
        total_possible_positions = n_links * n_times
        radar_wet_positions = 0
        legacy_finite_on_radar_wet = 0
        strict_finite_on_radar_wet = 0
        common_finite_on_radar_wet = 0

        n_link_blocks = math.ceil(n_links / args.link_chunk)
        n_time_blocks = math.ceil(n_times / args.time_chunk)
        total_blocks = n_link_blocks * n_time_blocks

        pbar = tqdm(
            total=total_blocks,
            desc="Comparing identical blocks",
            unit="block",
            dynamic_ncols=True,
        )

        for c0 in range(0, n_links, args.link_chunk):
            c1 = min(c0 + args.link_chunk, n_links)

            lc = legacy_cpos[c0:c1]
            sc = strict_cpos[c0:c1]
            rc = radar_cpos[c0:c1]

            for t0 in range(0, n_times, args.time_chunk):
                t1 = min(t0 + args.time_chunk, n_times)

                lt = legacy_tpos[t0:t1]
                st = strict_tpos[t0:t1]
                rt = radar_tpos[t0:t1]

                legacy_block = read_block(
                    legacy.rain,
                    legacy.cml_dim,
                    legacy.time_dim,
                    lc,
                    lt,
                )
                strict_block = read_block(
                    strict.rain,
                    strict.cml_dim,
                    strict.time_dim,
                    sc,
                    st,
                )
                radar_block = read_block(
                    radar,
                    radar_cml_dim,
                    radar_time_dim,
                    rc,
                    rt,
                )

                if (
                    legacy_block.shape != strict_block.shape
                    or legacy_block.shape != radar_block.shape
                ):
                    raise RuntimeError(
                        "Aligned block shapes differ: "
                        f"legacy={legacy_block.shape}, "
                        f"strict={strict_block.shape}, "
                        f"radar={radar_block.shape}"
                    )

                radar_base = (
                    np.isfinite(radar_block)
                    & (radar_block >= args.min_rain)
                    & (radar_block <= args.max_rain)
                    & (radar_block > args.radar_wet_thr)
                )

                legacy_valid = (
                    np.isfinite(legacy_block)
                    & (legacy_block >= args.min_rain)
                    & (legacy_block <= args.max_rain)
                )

                strict_valid = (
                    np.isfinite(strict_block)
                    & (strict_block >= args.min_rain)
                    & (strict_block <= args.max_rain)
                )

                # Same evaluation positions for BOTH products.
                common_mask = (
                    radar_base
                    & legacy_valid
                    & strict_valid
                )

                # Each filter's own valid evaluation population.
                own_legacy_mask = radar_base & legacy_valid
                own_strict_mask = radar_base & strict_valid

                radar_wet_positions += int(radar_base.sum())
                legacy_finite_on_radar_wet += int(own_legacy_mask.sum())
                strict_finite_on_radar_wet += int(own_strict_mask.sum())
                common_finite_on_radar_wet += int(common_mask.sum())

                if common_mask.any():
                    rv = radar_block[common_mask]
                    lv = legacy_block[common_mask]
                    sv = strict_block[common_mask]

                    update_stats(
                        stats_common_legacy,
                        lv,
                        rv,
                    )
                    update_stats(
                        stats_common_strict,
                        sv,
                        rv,
                    )

                    # Same sampled positions in the two common-sample panels.
                    common_reservoir.add(
                        rv,
                        lv,
                        sv,
                    )

                if own_legacy_mask.any():
                    rv = radar_block[own_legacy_mask]
                    lv = legacy_block[own_legacy_mask]

                    update_stats(
                        stats_own_legacy,
                        lv,
                        rv,
                    )
                    own_legacy_reservoir.add(
                        rv,
                        lv,
                    )

                if own_strict_mask.any():
                    rv = radar_block[own_strict_mask]
                    sv = strict_block[own_strict_mask]

                    update_stats(
                        stats_own_strict,
                        sv,
                        rv,
                    )
                    own_strict_reservoir.add(
                        rv,
                        sv,
                    )

                pbar.update(1)

        pbar.close()

        # ------------------------------------------------------------
        # Final exact metrics
        # ------------------------------------------------------------
        common_legacy = finalize_stats(stats_common_legacy)
        common_strict = finalize_stats(stats_common_strict)

        own_legacy = finalize_stats(stats_own_legacy)
        own_strict = finalize_stats(stats_own_strict)

        if common_legacy["N"] == 0:
            raise RuntimeError(
                "No common fair wet samples remain."
            )

        assert common_legacy["N"] == common_strict["N"], (
            "Internal fairness error: common-sample N differs."
        )

        print("\n" + "-" * 100)
        print("COMMON-SAMPLE RESULTS: EXACTLY THE SAME (cml_id, time) POSITIONS")
        print("-" * 100)

        print("\nLegacy filter")
        print(pd.Series(common_legacy).to_string())

        print("\nStrict filter")
        print(pd.Series(common_strict).to_string())

        print("\n" + "-" * 100)
        print("OWN-VALID-POPULATION RESULTS")
        print("-" * 100)

        print("\nLegacy filter")
        print(pd.Series(own_legacy).to_string())

        print("\nStrict filter")
        print(pd.Series(own_strict).to_string())

        # ------------------------------------------------------------
        # Metrics table
        # ------------------------------------------------------------
        rows = []

        for evaluation, filter_name, file_path, metrics in [
            (
                "common_samples",
                args.legacy_label,
                legacy_rain_nc,
                common_legacy,
            ),
            (
                "common_samples",
                args.strict_label,
                strict_rain_nc,
                common_strict,
            ),
            (
                "own_valid_population",
                args.legacy_label,
                legacy_rain_nc,
                own_legacy,
            ),
            (
                "own_valid_population",
                args.strict_label,
                strict_rain_nc,
                own_strict,
            ),
        ]:
            rows.append(
                {
                    "evaluation": evaluation,
                    "filter": filter_name,
                    "rain_file": file_path,
                    "rain_variable": args.rain_variable,
                    "baseline_method": BASELINE_METHOD,
                    "kr_model": KR_MODEL,
                    "waa_method": WAA_METHOD,
                    "radar_wet_operator": ">",
                    "radar_wet_threshold_mmph": args.radar_wet_thr,
                    "min_rain_mmph": args.min_rain,
                    "max_rain_mmph": args.max_rain,
                    **metrics,
                }
            )

        metrics_df = pd.DataFrame(rows)

        metrics_path = Path(args.out_dir) / "filter_strategy_metrics.csv"
        metrics_df.to_csv(
            metrics_path,
            index=False,
        )

        # ------------------------------------------------------------
        # Experiment summary
        # ------------------------------------------------------------
        summary = {
            "legacy_filtered_cml": str(args.legacy_cml),
            "strict_filtered_cml": str(args.strict_cml),
            "legacy_rainfall_product": str(legacy_rain_nc),
            "strict_rainfall_product": str(strict_rain_nc),
            "radar_source_file": str(args.radar_nc),
            "rain_variable": args.rain_variable,
            "baseline_method": BASELINE_METHOD,
            "kr_model": KR_MODEL,
            "waa_method": WAA_METHOD,
            "common_physical_links": int(n_links),
            "common_5min_timestamps": int(n_times),
            "evaluation_start": str(pd.Timestamp(common_times[0])),
            "evaluation_end": str(pd.Timestamp(common_times[-1])),
            "total_possible_common_link_time_positions": int(total_possible_positions),
            "radar_wet_positions": int(radar_wet_positions),
            "legacy_valid_radar_wet_positions": int(legacy_finite_on_radar_wet),
            "strict_valid_radar_wet_positions": int(strict_finite_on_radar_wet),
            "common_valid_radar_wet_positions": int(common_finite_on_radar_wet),
            "legacy_own_valid_fraction_of_radar_wet": (
                float(legacy_finite_on_radar_wet / radar_wet_positions)
                if radar_wet_positions > 0
                else np.nan
            ),
            "strict_own_valid_fraction_of_radar_wet": (
                float(strict_finite_on_radar_wet / radar_wet_positions)
                if radar_wet_positions > 0
                else np.nan
            ),
            "common_valid_fraction_of_radar_wet": (
                float(common_finite_on_radar_wet / radar_wet_positions)
                if radar_wet_positions > 0
                else np.nan
            ),
            "common_sample_N": int(common_legacy["N"]),
            "legacy_own_N": int(own_legacy["N"]),
            "strict_own_N": int(own_strict["N"]),
            "radar_wet_threshold_mmph": float(args.radar_wet_thr),
            "min_rain_mmph": float(args.min_rain),
            "max_rain_mmph": float(args.max_rain),
            "plot_sample_cap_per_population": int(args.max_plot_points),
        }

        summary_path = (
            Path(args.out_dir)
            / "filter_strategy_experiment_summary.csv"
        )
        pd.DataFrame([summary]).to_csv(
            summary_path,
            index=False,
        )

        json_path = (
            Path(args.out_dir)
            / "filter_strategy_experiment_summary.json"
        )
        json_path.write_text(
            json.dumps(summary, indent=2),
            encoding="utf-8",
        )

        # ------------------------------------------------------------
        # Common-sample hexbin
        # ------------------------------------------------------------
        common_radar, common_legacy_plot, common_strict_plot = (
            common_reservoir.arrays()
        )

        make_two_panel_hexbin(
            radar_left=common_radar,
            cml_left=common_legacy_plot,
            radar_right=common_radar,  # identical radar positions
            cml_right=common_strict_plot,
            metrics_left=common_legacy,
            metrics_right=common_strict,
            title_left=args.legacy_label,
            title_right=args.strict_label,
            main_title=(
                "Common-sample comparison "
                "(identical CML links, times, and radar observations)"
            ),
            out_png=(
                Path(args.out_dir)
                / "hexbin_common_samples.png"
            ),
            out_pdf=(
                Path(args.out_dir)
                / "hexbin_common_samples.pdf"
            ),
            min_rain=args.min_rain,
            max_rain=args.max_rain,
            gridsize=args.gridsize,
        )

        # ------------------------------------------------------------
        # Own-population hexbin
        # ------------------------------------------------------------
        own_legacy_radar, own_legacy_plot = (
            own_legacy_reservoir.arrays()
        )
        own_strict_radar, own_strict_plot = (
            own_strict_reservoir.arrays()
        )

        make_two_panel_hexbin(
            radar_left=own_legacy_radar,
            cml_left=own_legacy_plot,
            radar_right=own_strict_radar,
            cml_right=own_strict_plot,
            metrics_left=own_legacy,
            metrics_right=own_strict,
            title_left=(
                f"{args.legacy_label}\n"
                f"own valid N={own_legacy['N']:,}"
            ),
            title_right=(
                f"{args.strict_label}\n"
                f"own valid N={own_strict['N']:,}"
            ),
            main_title=(
                "Each filter's own valid radar-wet population"
            ),
            out_png=(
                Path(args.out_dir)
                / "hexbin_own_valid_population.png"
            ),
            out_pdf=(
                Path(args.out_dir)
                / "hexbin_own_valid_population.pdf"
            ),
            min_rain=args.min_rain,
            max_rain=args.max_rain,
            gridsize=args.gridsize,
        )

        print("\nSaved comparison outputs:")
        for path in [
            Path(args.out_dir) / "hexbin_common_samples.png",
            Path(args.out_dir) / "hexbin_common_samples.pdf",
            Path(args.out_dir) / "hexbin_own_valid_population.png",
            Path(args.out_dir) / "hexbin_own_valid_population.pdf",
            metrics_path,
            summary_path,
            json_path,
        ]:
            print(f"  {path}")

        return metrics_df, pd.DataFrame([summary])

    finally:
        legacy.ds.close()
        strict.ds.close()


# =============================================================================
# CLI
# =============================================================================
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Run the same CML rainfall-estimation framework on two filtered "
            "CML datasets and compare the filters fairly against radar."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument(
        "--legacy-cml",
        default="/home/jovyan/cml_preprocessed_v3.nc",
        help=(
            "Filtered CML dataset from the legacy filtering strategy."
        ),
    )
    p.add_argument(
        "--strict-cml",
        default=(
            "/home/jovyan/"
            "CML_S2020E2021_1MIN_CZ_QC_STRICT_MINIMAL.nc"
        ),
        help=(
            "Filtered CML dataset from the strict "
            "per-link/per-sublink/per-observation QC strategy."
        ),
    )
    p.add_argument(
        "--radar-nc",
        default="/home/jovyan/RADAR_S2021E2021_5MIN_CZ.nc",
    )
    p.add_argument(
        "--estimator-script",
        default=(
            "/home/jovyan/"
            "cml_estimate_and_grid_cml_radar_gauge_merged_parallel_timewindow_v2.py"
        ),
    )
    p.add_argument(
        "--out-dir",
        default="/home/jovyan/filter_strategy_rainfall_experiment",
    )

    # The radar reference exists for 2021. Starting directly in 2021 avoids
    # reprocessing the entire 2020 record unnecessarily. The 48-hour baseline
    # pad is still read before the requested start by the estimation framework.
    p.add_argument(
        "--tstart",
        default="2021-01-01 00:00:00",
    )
    p.add_argument(
        "--tend",
        default="2021-12-31 23:59:00",
    )

    p.add_argument(
        "--eval-tstart",
        default=None,
        help=(
            "Optional comparison-only start time. "
            "Defaults to all common output timestamps."
        ),
    )
    p.add_argument(
        "--eval-tend",
        default=None,
        help=(
            "Optional comparison-only end time. "
            "Defaults to all common output timestamps."
        ),
    )

    p.add_argument(
        "--rain-variable",
        default="R_mean_5min",
        help=(
            "Final link-level rainfall variable used for the filter benchmark."
        ),
    )
    p.add_argument(
        "--radar-var",
        default="R_radar",
    )

    p.add_argument(
        "--legacy-label",
        default="Legacy filter",
    )
    p.add_argument(
        "--strict-label",
        default="Strict link/sublink/sample QC",
    )

    p.add_argument(
        "--radar-wet-thr",
        type=float,
        default=0.1,
        help=(
            "Final evaluation keeps radar > threshold."
        ),
    )
    p.add_argument(
        "--min-rain",
        type=float,
        default=0.0,
    )
    p.add_argument(
        "--max-rain",
        type=float,
        default=100.0,
    )

    p.add_argument(
        "--link-chunk",
        type=int,
        default=32,
        help="Common links read per comparison block.",
    )
    p.add_argument(
        "--time-chunk",
        type=int,
        default=2016,
        help="5-min timesteps read per block. 2016 = seven days.",
    )
    p.add_argument(
        "--max-plot-points",
        type=int,
        default=500_000,
        help=(
            "Uniform random sample cap for each hexbin population. "
            "Exact metrics still use every valid sample."
        ),
    )
    p.add_argument(
        "--gridsize",
        type=int,
        default=50,
    )
    p.add_argument(
        "--random-seed",
        type=int,
        default=42,
    )

    p.add_argument(
        "--skip-estimation",
        action="store_true",
        help=(
            "Do not rerun rainfall estimation. Reuse the two expected "
            "rainfall products in --out-dir."
        ),
    )
    p.add_argument(
        "--overwrite-estimation",
        action="store_true",
        help="Overwrite existing rainfall products.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print estimation commands without executing them.",
    )

    return p


# =============================================================================
# Main
# =============================================================================
def main() -> None:
    args = build_parser().parse_args()

    started = time.perf_counter()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    legacy_rain_nc = str(
        out_dir / "rainfall_legacy_filter.nc"
    )
    strict_rain_nc = str(
        out_dir / "rainfall_strict_filter.nc"
    )

    legacy_log_csv = str(
        out_dir / "rainfall_legacy_filter_log.csv"
    )
    strict_log_csv = str(
        out_dir / "rainfall_strict_filter_log.csv"
    )

    # Validate inputs.
    for path, label in [
        (args.legacy_cml, "legacy filtered CML"),
        (args.strict_cml, "strict filtered CML"),
        (args.radar_nc, "radar"),
        (args.estimator_script, "estimation script"),
    ]:
        if not Path(path).exists():
            raise FileNotFoundError(
                f"{label} file does not exist: {path}"
            )

    legacy_trsl = detect_trsl_var(
        args.legacy_cml
    )
    strict_trsl = detect_trsl_var(
        args.strict_cml
    )

    print("=" * 100)
    print("FILTER-STRATEGY RAINFALL EXPERIMENT")
    print("=" * 100)
    print(f"Legacy filtered input : {args.legacy_cml}")
    print(f"Legacy total loss     : {legacy_trsl or 'TX-RX fallback'}")
    print(f"Strict filtered input : {args.strict_cml}")
    print(f"Strict total loss     : {strict_trsl or 'TX-RX fallback'}")
    print(f"Radar                 : {args.radar_nc}")
    print(f"Estimator             : {args.estimator_script}")
    print(f"Output directory      : {out_dir}")
    print()
    print("IDENTICAL RAINFALL CONFIGURATION")
    print(f"  baseline             : {BASELINE_METHOD}")
    print(f"  k-R                  : {KR_MODEL}")
    print(f"  WAA                  : {WAA_METHOD}")
    print(f"  Pastorek A_max       : {PASTOREK_A_MAX}")
    print(f"  Pastorek zeta        : {PASTOREK_ZETA}")
    print(f"  Pastorek d           : {PASTOREK_D}")
    print(f"  baseline pad         : {BASELINE_PAD_HOURS} h")
    print(f"  min link length      : {MIN_LINK_KM} km")
    print(f"  min valid samples    : {MIN_VALID_SAMPLES}")
    print(
        f"  estimation wet thr   : "
        f"{RADAR_WET_THR_ESTIMATION} mm h-1"
    )
    print(f"  final hexbin wet thr : radar > {args.radar_wet_thr} mm h-1")
    print(f"  estimation window    : {args.tstart} -> {args.tend}")
    print("=" * 100)

    if not args.skip_estimation:
        legacy_cmd = build_estimation_command(
            estimator_script=args.estimator_script,
            cml_nc=args.legacy_cml,
            trsl_var=legacy_trsl,
            radar_nc=args.radar_nc,
            output_nc=legacy_rain_nc,
            log_csv=legacy_log_csv,
            tstart=args.tstart,
            tend=args.tend,
            overwrite=args.overwrite_estimation,
        )

        strict_cmd = build_estimation_command(
            estimator_script=args.estimator_script,
            cml_nc=args.strict_cml,
            trsl_var=strict_trsl,
            radar_nc=args.radar_nc,
            output_nc=strict_rain_nc,
            log_csv=strict_log_csv,
            tstart=args.tstart,
            tend=args.tend,
            overwrite=args.overwrite_estimation,
        )

        run_command(
            legacy_cmd,
            dry_run=args.dry_run,
        )
        run_command(
            strict_cmd,
            dry_run=args.dry_run,
        )

        if args.dry_run:
            print(
                "\nDry run complete. "
                "No rainfall files or comparison plots were created."
            )
            return

    else:
        print("\nSkipping rainfall estimation and reusing:")
        print(f"  {legacy_rain_nc}")
        print(f"  {strict_rain_nc}")

        for path in [
            legacy_rain_nc,
            strict_rain_nc,
        ]:
            if not Path(path).exists():
                raise FileNotFoundError(
                    f"--skip-estimation requested but file is missing: {path}"
                )

    compare_rainfall_products(
        args,
        legacy_rain_nc=legacy_rain_nc,
        strict_rain_nc=strict_rain_nc,
    )

    elapsed = time.perf_counter() - started

    print("\n" + "=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)
    print(f"Elapsed: {elapsed:.1f} s ({elapsed / 60.0:.1f} min)")
    print(f"Results: {out_dir}")
    print("=" * 100)


if __name__ == "__main__":
    main()
