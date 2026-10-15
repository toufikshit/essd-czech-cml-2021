#!/usr/bin/env python3
"""
30b_v6_quantization_from_diagnostics.py
=======================================

Link-derived quantization characterization for the EXISTING V6 population,
without requiring the large V6 NetCDF file.

Why this version exists
-----------------------
The exact V6 merged NetCDF may no longer be present, while the V6 diagnostics
CSV still contains the authoritative directional-link population and metadata.
Quantization itself only needs:
  * the exact V6 population (diagnostics CSV),
  * raw 1-min TX/RX,
  * the exact V6 conservative hardware-state detector,
  * P.676 physics.

The V6 WVD retrieval is NOT recomputed or changed.
CHMI humidity is NEVER read.

For every directional link:
  q_TX, q_RX, q_TL are inferred from native nonzero level changes occurring
  inside the same stable V6 hardware state and outside transition buffers.

Then:
  q_gamma = q_TL / L  [dB km-1]

and the exact nonlinear P.676 quantization-limited WVD resolution is

  Delta_rho_q(T,rho)
    = gamma^{-1}[gamma(T,rho) + q_TL/L] - rho.

Outputs:
  * per-series CSV
  * per-state CSV
  * compact per-series NetCDF
  * JSON summary

The standardized P.676 characterization defaults to:
  T = 10...30 degC
  rho = 2...25 g m-3

A reference value is also reported at:
  T = 20 degC
  rho = 10 g m-3

IMPORTANT:
  This is instrumental quantization-limited resolution, NOT total WVD error.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr


def load_module(path, name):
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def clean_time(x):
    t = pd.DatetimeIndex(pd.to_datetime(np.asarray(x)))
    if t.tz is not None:
        t = t.tz_convert("UTC").tz_localize(None)
    return t


def decode_strings(x):
    a = np.asarray(x)
    if a.dtype.kind == "S":
        if a.ndim == 1:
            return np.asarray([
                v.decode("utf-8", errors="ignore").replace("\x00", "").strip()
                for v in a
            ], dtype=object)
        return np.asarray([
            "".join(
                v.decode("utf-8", errors="ignore")
                if isinstance(v, (bytes, np.bytes_)) else str(v)
                for v in np.ravel(row)
            ).replace("\x00", "").strip()
            for row in a
        ], dtype=object)
    return np.asarray([
        str(v).replace("\x00", "").strip()
        for v in a.reshape(-1)
    ], dtype=object)


def v6_state_args():
    """Exact conservative state-detector defaults from V6."""
    return SimpleNamespace(
        state_window_days=7,
        state_min_days_each_side=4,
        state_min_valid_minutes_per_day=120,
        state_lower_quantile=0.05,
        state_large_tl_step_db=3.0,
        state_tx_step_db=2.0,
        state_min_tl_with_tx_db=2.0,
        state_rx_step_db=3.0,
        state_rx_tl_step_db=3.0,
        state_gap_hours=12.0,
        state_gap_tl_step_db=1.5,
        state_merge_days=5.0,
        hardware_transition_buffer_hours=12.0,
    )


def native_changes(values, state, transition, max_change_db, decimals):
    x = np.asarray(values, dtype=float)
    s = np.asarray(state, dtype=np.int16)
    tr = np.asarray(transition, dtype=bool)

    if len(x) < 2:
        return np.empty(0, dtype=float)

    ok = (
        np.isfinite(x[:-1])
        & np.isfinite(x[1:])
        & (s[:-1] >= 0)
        & (s[:-1] == s[1:])
        & ~tr[:-1]
        & ~tr[1:]
    )
    if not ok.any():
        return np.empty(0, dtype=float)

    d = np.abs(x[1:] - x[:-1])[ok]
    d = np.round(d, int(decimals))
    return d[
        np.isfinite(d)
        & (d > 0)
        & (d <= float(max_change_db))
    ]


def estimate_q(changes, args, state_level=False):
    """
    Estimate the largest supported base quantization step.

    A candidate q is accepted when at least `support_fraction` of observed
    positive changes are integer multiples of q within tolerance, and the
    one-step level itself occurs sufficiently often.

    Choosing the largest accepted q avoids trivially tiny divisors.
    """
    d = np.asarray(changes, dtype=float)
    d = d[np.isfinite(d) & (d > 0)]

    min_changes = (
        args.min_changes_state if state_level else args.min_changes
    )
    min_fund = (
        args.min_fundamental_state if state_level else args.min_fundamental
    )

    if len(d) < int(min_changes):
        return dict(
            q=np.nan, n=int(len(d)), unique=0,
            support=np.nan, fundamental=0, detected=0
        )

    u, count = np.unique(d, return_counts=True)
    candidates = []

    for step in u:
        for k in range(1, int(args.max_harmonic) + 1):
            q = float(step) / k
            if args.q_min_db <= q <= args.q_max_db:
                candidates.append(q)

    if not candidates:
        return dict(
            q=np.nan, n=int(len(d)), unique=int(len(u)),
            support=np.nan, fundamental=0, detected=0
        )

    candidates = np.sort(
        np.unique(np.round(np.asarray(candidates), 6))
    )[::-1]

    total = float(count.sum())

    for q in candidates:
        ratio = u / q
        nearest = np.rint(ratio)
        residual = np.abs(ratio - nearest)

        supported = residual <= float(args.relative_tolerance)
        support = float(count[supported].sum() / total)

        fundamental = np.abs(ratio - 1.0) <= float(args.relative_tolerance)
        fundamental_count = int(count[fundamental].sum())

        if (
            support >= float(args.support_fraction)
            and fundamental_count >= int(min_fund)
        ):
            return dict(
                q=float(q),
                n=int(len(d)),
                unique=int(len(u)),
                support=support,
                fundamental=fundamental_count,
                detected=1,
            )

    return dict(
        q=np.nan, n=int(len(d)), unique=int(len(u)),
        support=np.nan, fundamental=0, detected=0
    )


def q_for_signal(values, state, transition, args, state_id=None):
    s = np.asarray(state, dtype=np.int16)
    if state_id is not None:
        s = np.where(s == int(state_id), s, -1).astype(np.int16)

    d = native_changes(
        values, s, transition,
        args.max_native_change_db,
        args.rounding_decimals,
    )
    return estimate_q(d, args, state_level=state_id is not None)


def weighted_median(values, weights):
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    m = np.isfinite(v) & np.isfinite(w) & (w > 0)
    if not m.any():
        return np.nan

    v = v[m]
    w = w[m]
    order = np.argsort(v)
    v = v[order]
    w = w[order]

    cs = np.cumsum(w)
    return float(v[np.searchsorted(cs, 0.5 * cs[-1])])


def gamma_curve(temp_c, temp_grid, gamma_table):
    if (
        not np.isfinite(temp_c)
        or temp_c < temp_grid[0]
        or temp_c > temp_grid[-1]
    ):
        return None

    return np.asarray([
        np.interp(
            float(temp_c),
            temp_grid,
            np.asarray(gamma_table[:, j], dtype=float),
        )
        for j in range(gamma_table.shape[1])
    ])


def delta_rho_exact(
    q_db, length_km, temp_c, rho,
    temp_grid, rho_grid, gamma_table,
):
    if (
        not np.isfinite(q_db)
        or q_db <= 0
        or not np.isfinite(length_km)
        or length_km <= 0
    ):
        return np.nan

    g = gamma_curve(temp_c, temp_grid, gamma_table)
    if (
        g is None
        or np.any(~np.isfinite(g))
        or np.any(np.diff(g) < -1e-10)
    ):
        return np.nan

    if rho < rho_grid[0] or rho > rho_grid[-1]:
        return np.nan

    g0 = float(np.interp(float(rho), rho_grid, g))
    target = g0 + float(q_db) / float(length_km)

    if target > g[-1]:
        return np.nan

    rho_plus = float(np.interp(target, g, rho_grid))
    return float(rho_plus - rho)


def standardized_range(
    q_db, length_km, temp_grid, rho_grid, gamma_table, args
):
    vals = []

    temperatures = np.arange(
        args.std_t_min,
        args.std_t_max + 0.5 * args.std_t_step,
        args.std_t_step,
    )
    rhos = np.arange(
        args.std_rho_min,
        args.std_rho_max + 0.5 * args.std_rho_step,
        args.std_rho_step,
    )

    for tc in temperatures:
        for rho in rhos:
            x = delta_rho_exact(
                q_db, length_km, tc, rho,
                temp_grid, rho_grid, gamma_table,
            )
            if np.isfinite(x) and x >= 0:
                vals.append(x)

    if not vals:
        return np.nan, np.nan, np.nan

    a = np.asarray(vals, dtype=float)
    return (
        float(np.min(a)),
        float(np.median(a)),
        float(np.max(a)),
    )


def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Compute link-derived radio quantization and P.676 WVD "
            "quantization resolution for the exact V6 population."
        ),
    )

    root = "/home/jovyan/ESSD_publication_dataset"

    p.add_argument(
        "--v6-diagnostics-csv",
        default=f"{root}/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv",
        help="Authoritative exact V6 directional-link population.",
    )
    p.add_argument(
        "--raw-cml-nc",
        default=f"{root}/CZ_CML_RAW_1MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--v6-script",
        default=f"{root}/19_wvd_eband_progressive_baseline_v6.py",
    )
    p.add_argument(
        "--physics-helper-script",
        default=f"{root}/13_wvd_raw_60_100_timevarying_T.py",
    )

    p.add_argument(
        "--out-nc",
        default=f"{root}/CZ_CML_WVD_V6_LINK_QUANTIZATION.nc",
    )
    p.add_argument(
        "--series-csv",
        default=f"{root}/CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv",
    )
    p.add_argument(
        "--states-csv",
        default=f"{root}/CZ_CML_WVD_V6_LINK_QUANTIZATION_STATES.csv",
    )

    p.add_argument("--pressure-hpa", type=float, default=980.0)

    # Quantization inference.
    p.add_argument("--q-min-db", type=float, default=0.01)
    p.add_argument("--q-max-db", type=float, default=2.0)
    p.add_argument("--max-native-change-db", type=float, default=5.0)
    p.add_argument("--rounding-decimals", type=int, default=4)
    p.add_argument("--max-harmonic", type=int, default=12)
    p.add_argument("--support-fraction", type=float, default=0.95)
    p.add_argument("--relative-tolerance", type=float, default=0.03)
    p.add_argument("--min-changes", type=int, default=500)
    p.add_argument("--min-fundamental", type=int, default=20)
    p.add_argument("--min-changes-state", type=int, default=150)
    p.add_argument("--min-fundamental-state", type=int, default=8)

    # P.676 lookup.
    p.add_argument("--temp-grid-min", type=float, default=-35.0)
    p.add_argument("--temp-grid-max", type=float, default=45.0)
    p.add_argument("--temp-grid-step", type=float, default=0.25)
    p.add_argument("--rho-grid-max", type=float, default=60.0)
    p.add_argument("--rho-grid-step", type=float, default=0.02)

    # Standardized Fencl-style characterization domain.
    p.add_argument("--std-t-min", type=float, default=10.0)
    p.add_argument("--std-t-max", type=float, default=30.0)
    p.add_argument("--std-t-step", type=float, default=1.0)
    p.add_argument("--std-rho-min", type=float, default=2.0)
    p.add_argument("--std-rho-max", type=float, default=25.0)
    p.add_argument("--std-rho-step", type=float, default=0.5)

    p.add_argument(
        "--length-min-km",
        type=float,
        default=None,
        help="Optional diagnostic subset; omit to characterize all V6 links.",
    )
    p.add_argument(
        "--length-max-km",
        type=float,
        default=None,
    )
    p.add_argument("--series-start", type=int, default=0,
                   help="Start row (0-based) in sorted V6 diagnostics population.")
    p.add_argument("--series-stop", type=int, default=None,
                   help="Stop row (exclusive) in sorted V6 diagnostics population.")
    p.add_argument("--max-series", type=int, default=None,
                   help="Optional cap after applying series-start/series-stop.")
    p.add_argument("--overwrite", action="store_true")

    return p.parse_args()


def main():
    args = parse_args()
    t0 = time.time()

    diag_path = Path(args.v6_diagnostics_csv).expanduser().resolve()
    raw_path = Path(args.raw_cml_nc).expanduser().resolve()

    if not diag_path.is_file():
        raise FileNotFoundError(
            f"V6 diagnostics CSV not found: {diag_path}"
        )
    if not raw_path.is_file():
        raise FileNotFoundError(
            f"Raw CML NetCDF not found: {raw_path}"
        )

    v6mod = load_module(args.v6_script, "v6_quant_state_helper")
    phys = load_module(args.physics_helper_script, "v6_quant_physics_helper")
    state_args = v6_state_args()

    diag = pd.read_csv(diag_path)

    required = {
        "series_index", "cml_id", "sublink_id",
        "frequency_ghz", "length_km",
    }
    missing = required.difference(diag.columns)
    if missing:
        raise KeyError(
            f"V6 diagnostics missing columns: {sorted(missing)}"
        )

    diag = diag.sort_values("series_index").reset_index(drop=True)

    if args.length_min_km is not None:
        diag = diag[
            pd.to_numeric(diag["length_km"], errors="coerce")
            >= float(args.length_min_km)
        ].copy()

    if args.length_max_km is not None:
        diag = diag[
            pd.to_numeric(diag["length_km"], errors="coerce")
            <= float(args.length_max_km)
        ].copy()

    start = max(0, int(args.series_start))
    stop = len(diag) if args.series_stop is None else min(len(diag), int(args.series_stop))
    if stop < start:
        raise ValueError("series-stop must be >= series-start")
    diag = diag.iloc[start:stop].copy()

    if args.max_series is not None:
        diag = diag.iloc[: int(args.max_series)].copy()

    diag = diag.reset_index(drop=True)

    raw = xr.open_dataset(raw_path, cache=False)

    try:
        raw_cids = decode_strings(raw["cml_id"].values)
        raw_sids = decode_strings(raw["sublink_id"].values)

        cmap = {str(v): i for i, v in enumerate(raw_cids)}
        smap = {str(v): i for i, v in enumerate(raw_sids)}
        raw_times = clean_time(raw["time"].values)

        print("=" * 112)
        print("V6 LINK-DERIVED QUANTIZATION CHARACTERIZATION")
        print("=" * 112)
        print("Population source :", diag_path)
        print("Raw CML           :", raw_path)
        print("Directional links :", len(diag))
        print("V6 retrieval      : NOT CHANGED")
        print("Humidity           : NEVER READ")
        print("Quantization       : inferred independently from each link/state")
        print("=" * 112)

        series_rows = []
        state_rows = []

        for jj, row in diag.iterrows():
            tic = time.time()

            series_index = int(row["series_index"])
            cid = str(row["cml_id"])
            sid = str(row["sublink_id"])
            f = float(row["frequency_ghz"])
            L = float(row["length_km"])

            if cid not in cmap:
                raise KeyError(f"CML {cid} not found in raw file")
            if sid not in smap:
                raise KeyError(f"Sublink {sid} not found in raw file")

            ci = cmap[cid]
            si = smap[sid]

            tx = np.asarray(
                raw["tx"].isel(cml_id=ci, sublink_id=si).values,
                dtype=float,
            ).reshape(-1)
            rx = np.asarray(
                raw["rx"].isel(cml_id=ci, sublink_id=si).values,
                dtype=float,
            ).reshape(-1)

            valid = (
                np.isfinite(tx)
                & np.isfinite(rx)
                & (tx > phys.TX_MIN_VALID_EXCLUSIVE)
                & (tx != phys.TX_FILL_CODE)
                & (rx > phys.RX_MIN_VALID_EXCLUSIVE)
                & (rx <= phys.RX_MAX_VALID_INCLUSIVE)
            )

            txv = np.where(valid, tx, np.nan)
            rxv = np.where(valid, rx, np.nan)
            tl = np.where(valid, tx - rx, np.nan)

            events, _ = v6mod.detect_hardware_events(
                raw_times, tx, rx, valid, state_args
            )
            state, transition, boundaries = v6mod.assign_state_minutes(
                raw_times,
                events,
                state_args.hardware_transition_buffer_hours,
            )
            intervals = v6mod.state_intervals(boundaries)

            qtx_all = q_for_signal(
                txv, state, transition, args
            )
            qrx_all = q_for_signal(
                rxv, state, transition, args
            )
            qtl_all = q_for_signal(
                tl, state, transition, args
            )

            state_q = []
            state_w = []

            for state_id, a, b in intervals:
                qtx = q_for_signal(
                    txv, state, transition, args, state_id
                )
                qrx = q_for_signal(
                    rxv, state, transition, args, state_id
                )
                qtl = q_for_signal(
                    tl, state, transition, args, state_id
                )

                if np.isfinite(qtl["q"]):
                    state_q.append(qtl["q"])
                    state_w.append(max(1, qtl["n"]))

                state_rows.append({
                    "series_index": series_index,
                    "cml_id": cid,
                    "sublink_id": sid,
                    "frequency_ghz": f,
                    "length_km": L,
                    "state_id": int(state_id),
                    "state_start": pd.Timestamp(a),
                    "state_end_exclusive": pd.Timestamp(b),

                    "tx_quantization_db": qtx["q"],
                    "tx_quantization_n_changes": qtx["n"],
                    "tx_quantization_support_fraction": qtx["support"],

                    "rx_quantization_db": qrx["q"],
                    "rx_quantization_n_changes": qrx["n"],
                    "rx_quantization_support_fraction": qrx["support"],

                    "tl_quantization_db": qtl["q"],
                    "tl_quantization_n_changes": qtl["n"],
                    "tl_quantization_n_unique_changes": qtl["unique"],
                    "tl_quantization_fundamental_count": qtl["fundamental"],
                    "tl_quantization_support_fraction": qtl["support"],

                    "specific_attenuation_quantization_db_per_km": (
                        qtl["q"] / L
                        if np.isfinite(qtl["q"]) and L > 0
                        else np.nan
                    ),
                })

            # Prefer state-resolved estimates; fall back to all stable states.
            qtl_rep = weighted_median(state_q, state_w)
            if not np.isfinite(qtl_rep):
                qtl_rep = qtl_all["q"]

            (
                temp_grid,
                rho_grid,
                gamma_table,
                monotonic,
                _,
            ) = phys.p676_temperature_lookup(
                f,
                args.pressure_hpa,
                args.temp_grid_min,
                args.temp_grid_max,
                args.temp_grid_step,
                args.rho_grid_max,
                args.rho_grid_step,
            )

            if np.isfinite(qtl_rep) and bool(monotonic):
                qgamma = qtl_rep / L
                dref = delta_rho_exact(
                    qtl_rep, L, 20.0, 10.0,
                    temp_grid, rho_grid, gamma_table,
                )
                dmin, dmed, dmax = standardized_range(
                    qtl_rep, L,
                    temp_grid, rho_grid, gamma_table,
                    args,
                )
            else:
                qgamma = np.nan
                dref = dmin = dmed = dmax = np.nan

            series_rows.append({
                "series_index": series_index,
                "cml_id": cid,
                "sublink_id": sid,
                "frequency_ghz": f,
                "length_km": L,

                "hardware_n_states": len(intervals),
                "hardware_n_events": len(events),

                "tx_quantization_db": qtx_all["q"],
                "tx_quantization_n_changes": qtx_all["n"],
                "tx_quantization_support_fraction": qtx_all["support"],

                "rx_quantization_db": qrx_all["q"],
                "rx_quantization_n_changes": qrx_all["n"],
                "rx_quantization_support_fraction": qrx_all["support"],

                "tl_quantization_db": qtl_rep,
                "tl_quantization_n_changes": qtl_all["n"],
                "tl_quantization_n_unique_changes": qtl_all["unique"],
                "tl_quantization_fundamental_count": qtl_all["fundamental"],
                "tl_quantization_support_fraction": qtl_all["support"],
                "tl_quantization_n_states_detected": len(state_q),

                "specific_attenuation_quantization_db_per_km": qgamma,

                "wvd_quantization_resolution_ref20c_rho10_gm3": dref,
                "wvd_quantization_resolution_standard_min_gm3": dmin,
                "wvd_quantization_resolution_standard_median_gm3": dmed,
                "wvd_quantization_resolution_standard_max_gm3": dmax,
            })

            if (
                jj < 10
                or (jj + 1) % 25 == 0
                or jj + 1 == len(diag)
            ):
                qstr = (
                    f"{qtl_rep:.4g}"
                    if np.isfinite(qtl_rep)
                    else "NA"
                )
                dstr = (
                    f"{dref:.3g}"
                    if np.isfinite(dref)
                    else "NA"
                )
                print(
                    f"[{jj+1:4d}/{len(diag):4d}] "
                    f"{cid} sub={sid} "
                    f"f={f:.2f}GHz L={L:.2f}km "
                    f"qTL={qstr}dB "
                    f"dRho20C10={dstr}g/m3 "
                    f"states_q={len(state_q)}/{len(intervals)} "
                    f"({time.time()-tic:.1f}s)"
                )

        sdf = pd.DataFrame(series_rows)
        stf = pd.DataFrame(state_rows)

        series_path = Path(args.series_csv).expanduser().resolve()
        states_path = Path(args.states_csv).expanduser().resolve()
        out_nc = Path(args.out_nc).expanduser().resolve()

        for pp in (series_path, states_path, out_nc):
            pp.parent.mkdir(parents=True, exist_ok=True)

        if out_nc.exists() and not args.overwrite:
            raise FileExistsError(
                f"{out_nc} exists; use --overwrite"
            )
        if out_nc.exists():
            out_nc.unlink()

        sdf.to_csv(series_path, index=False)
        stf.to_csv(states_path, index=False)

        with netCDF4.Dataset(out_nc, "w", format="NETCDF4") as nc:
            nc.createDimension("series", len(sdf))

            nc.createVariable(
                "series_index", "i4", ("series",)
            )[:] = sdf["series_index"].to_numpy(np.int32)

            nc.createVariable(
                "cml_id", str, ("series",)
            )[:] = sdf["cml_id"].to_numpy(object)

            nc.createVariable(
                "sublink_id", str, ("series",)
            )[:] = sdf["sublink_id"].to_numpy(object)

            for name, units in (
                ("frequency_ghz", "GHz"),
                ("length_km", "km"),
                ("tx_quantization_db", "dB"),
                ("rx_quantization_db", "dB"),
                ("tl_quantization_db", "dB"),
                (
                    "specific_attenuation_quantization_db_per_km",
                    "dB km-1",
                ),
                (
                    "wvd_quantization_resolution_ref20c_rho10_gm3",
                    "g m-3",
                ),
                (
                    "wvd_quantization_resolution_standard_min_gm3",
                    "g m-3",
                ),
                (
                    "wvd_quantization_resolution_standard_median_gm3",
                    "g m-3",
                ),
                (
                    "wvd_quantization_resolution_standard_max_gm3",
                    "g m-3",
                ),
            ):
                v = nc.createVariable(
                    name, "f4", ("series",), fill_value=np.nan
                )
                v[:] = sdf[name].to_numpy(np.float32)
                v.units = units

            nc.title = (
                "Link-derived quantization characterization "
                "for the existing V6 WVD population"
            )
            nc.population_source = str(diag_path)
            nc.source_raw_cml = str(raw_path)
            nc.v6_retrieval_modified = "no"
            nc.humidity_used = "no"
            nc.quantization_is_not_radio_noise = "true"
            nc.quantization_method = (
                "Largest supported base step whose native stable-state "
                "TX/RX/TL changes are integer multiples within tolerance."
            )
            nc.wvd_resolution_method = (
                "Exact P.676 propagation: "
                "gamma(rho+Delta_rho)=gamma(rho)+q_TL/L."
            )
            nc.pressure_hpa = float(args.pressure_hpa)
            nc.reference_temperature_c = 20.0
            nc.reference_wvd_gm3 = 10.0
            nc.standard_temperature_range_c = (
                f"{args.std_t_min:g}-{args.std_t_max:g}"
            )
            nc.standard_wvd_range_gm3 = (
                f"{args.std_rho_min:g}-{args.std_rho_max:g}"
            )

        summary = {
            "population_source": str(diag_path),
            "raw_cml": str(raw_path),
            "n_series": int(len(sdf)),
            "n_states": int(len(stf)),
            "n_series_with_tl_quantization": int(
                np.isfinite(sdf["tl_quantization_db"]).sum()
            ),
            "fraction_with_tl_quantization": float(
                np.isfinite(sdf["tl_quantization_db"]).mean()
            ),
            "retrieval_modified": False,
            "humidity_used": False,
            "elapsed_seconds": float(time.time() - t0),
        }

        out_nc.with_suffix(".summary.json").write_text(
            json.dumps(summary, indent=2),
            encoding="utf-8",
        )

        print("\n" + "=" * 112)
        print("DONE — V6 QUANTIZATION")
        print("=" * 112)
        print("Series CSV :", series_path)
        print("States CSV :", states_path)
        print("NetCDF     :", out_nc)
        print(
            "Detected qTL:",
            f"{summary['n_series_with_tl_quantization']}/"
            f"{summary['n_series']}",
        )
        print("V6 retrieval was NOT changed.")

    finally:
        raw.close()


if __name__ == "__main__":
    main()
