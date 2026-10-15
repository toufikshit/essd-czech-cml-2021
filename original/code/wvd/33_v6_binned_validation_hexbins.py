#!/usr/bin/env python3
"""
33_v6_binned_validation_hexbins.py
==================================

Strict V6-vs-CHMI validation hexbins stratified by LINK-INTERNAL physics:

  * quantization-limited WVD sensitivity/resolution bins (Delta rho_q)
  * E-band frequency range (71-76 vs 81-86 GHz)
  * disjoint path-length bins (2-3, 3-4, 4-5, 5-6, 6-10 km)
  * frequency x length
  * frequency x sensitivity
  * full length x sensitivity grid separately for each frequency range

The V6 retrieval itself is NOT modified. CHMI humidity is used only here as an
external validation reference. Quantization/sensitivity was already derived
without humidity in CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv.

Validation semantics are intentionally identical to
21_validate_wvd_v6_vs_original.py:

  * nearest CHMI station with H and T <= 5 km from the CML path
  * exact CHMI 10-min timestamps; no interpolation
  * 5-min CML -> 10-min only when exactly two finite V6 WVD values are present
    and retrieval_status == 1 for both 5-min samples
  * [t, t+5 min] CML pair is labelled t+10 min
  * each directional sublink is independent

Sensitivity variable used for the main bins
-------------------------------------------
    wvd_quantization_resolution_ref20c_rho10_gm3

This is the exact nonlinear P.676 WVD increment corresponding to one detected
TL quantization step q_TL at the link's actual frequency and path length,
evaluated at the standardized reference state T=20 C, rho=10 g m-3.

Default DISJOINT sensitivity bins [g m-3]:
    <0.20
    0.20-0.25
    0.25-0.30
    0.30-0.40
    0.40-0.50
    >=0.50

These are diagnostic presentation bins. They are NOT selected by maximizing
CHMI PCC/RMSE/bias.

Default outputs
---------------
WVD_V6_QUANTIZATION_BINNED_VALIDATION/
  validation_pairs_<branch>.npz
  validation_population_links.csv
  binned_validation_metrics.csv
  population_counts.csv
  figures/
    01_hexbin_by_sensitivity_<branch>.png/pdf
    02_hexbin_by_frequency_<branch>.png/pdf
    03_hexbin_by_length_<branch>.png/pdf
    04_hexbin_frequency_x_length_<branch>.png/pdf
    05_hexbin_frequency_x_sensitivity_<branch>.png/pdf
    06_hexbin_joint_71_76GHz_<branch>.png/pdf
    07_hexbin_joint_81_86GHz_<branch>.png/pdf

The NPZ stores the exact strict 10-min pair population compactly, so plotting
can be repeated without reopening CHMI/V6 NetCDF files.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import xarray as xr


FREQ_BANDS = [
    ("71-76 GHz", 71.0, 76.0),
    ("81-86 GHz", 81.0, 86.0),
]

LENGTH_BINS = [
    ("2-3 km", 2.0, 3.0),
    ("3-4 km", 3.0, 4.0),
    ("4-5 km", 4.0, 5.0),
    ("5-6 km", 5.0, 6.0),
    ("6-10 km", 6.0, 10.000001),
]

# Disjoint, physically interpretable standardized WVD-resolution bins.
SENS_BINS = [
    ("<0.20", -np.inf, 0.20),
    ("0.20-0.25", 0.20, 0.25),
    ("0.25-0.30", 0.25, 0.30),
    ("0.30-0.40", 0.30, 0.40),
    ("0.40-0.50", 0.40, 0.50),
    (">=0.50", 0.50, np.inf),
]

# Explicit palette close to the previous WVD validation figures.
HEX_COLORS = [
    "#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
    "#ffffbf", "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142",
]
HEX_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "wvd_hex", HEX_COLORS, N=256
)


def load_module(path: Path, name: str):
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def normalize_id(x):
    s = str(x).strip()
    if s.endswith(".0"):
        try:
            return str(int(float(s)))
        except Exception:
            pass
    return s


def resolve_v6_path(requested: Path, root: Path) -> Path:
    requested = requested.expanduser().resolve()
    if requested.is_file():
        return requested

    candidates = []
    exact = root / "CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc"
    if exact.is_file():
        candidates.append(exact)

    for p in sorted(root.glob("*V6*2021*.nc")):
        u = p.name.upper()
        if "QUANT" in u or "VALID" in u:
            continue
        if p not in candidates:
            candidates.append(p)

    if len(candidates) == 1:
        print(f"Requested V6 path absent; auto-selected: {candidates[0]}")
        return candidates[0]

    if not candidates:
        raise FileNotFoundError(
            f"V6 retrieval NetCDF not found. Requested: {requested}\n"
            f"Expected production name: {exact}"
        )

    raise RuntimeError(
        "More than one possible V6 NetCDF found; supply --v6-nc explicitly:\n  "
        + "\n  ".join(str(x) for x in candidates)
    )


def parse_args():
    root = Path("/home/jovyan/ESSD_publication_dataset")
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--root", type=Path, default=root)
    p.add_argument(
        "--v6-nc", type=Path,
        default=root / "CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc",
    )
    p.add_argument(
        "--v6-diagnostics-csv", type=Path,
        default=root / "CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv",
    )
    p.add_argument(
        "--quant-csv", type=Path,
        default=root / "CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv",
    )
    p.add_argument(
        "--raw-cml-nc", type=Path,
        default=root / "CZ_CML_RAW_1MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--humidity-nc", type=Path,
        default=Path("/home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_HUMIDITY_2021.nc"),
    )
    p.add_argument(
        "--temperature-nc", type=Path,
        default=Path("/home/jovyan/opendata.chmi.cz/processed/CHMI_10MIN_TEMPERATURE_2021.nc"),
    )
    p.add_argument(
        "--validation-helper-script", type=Path,
        default=root / "21_validate_wvd_v6_vs_original.py",
    )
    p.add_argument(
        "--station-assignments-csv", type=Path,
        default=root / "WVD_V6_VALIDATION/station_assignments.csv",
        help="Existing exact station mapping from V6 validation; recomputed if absent.",
    )
    p.add_argument(
        "--out-dir", type=Path,
        default=root / "WVD_V6_QUANTIZATION_BINNED_VALIDATION",
    )
    p.add_argument("--branch", choices=["GRAF", "GAUGE"], default="GRAF")
    p.add_argument("--max-station-distance-km", type=float, default=5.0)
    p.add_argument("--scatter-max-gm3", type=float, default=35.0)
    p.add_argument("--gridsize", type=int, default=55)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--reuse-pairs", action="store_true", default=False)
    p.add_argument("--min-panel-pairs", type=int, default=100)
    return p.parse_args()


def load_quant(path: Path):
    q = pd.read_csv(path)
    required = [
        "cml_id", "sublink_id", "frequency_ghz", "length_km",
        "tl_quantization_db",
        "specific_attenuation_quantization_db_per_km",
        "wvd_quantization_resolution_ref20c_rho10_gm3",
        "wvd_quantization_resolution_standard_min_gm3",
        "wvd_quantization_resolution_standard_median_gm3",
        "wvd_quantization_resolution_standard_max_gm3",
    ]
    missing = [c for c in required if c not in q.columns]
    if missing:
        raise KeyError(f"Quantization CSV missing: {missing}")

    q["cml_id"] = q.cml_id.map(normalize_id)
    q["sublink_id"] = q.sublink_id.map(normalize_id)
    for c in required[2:]:
        q[c] = pd.to_numeric(q[c], errors="coerce")

    if "series_index" in q.columns:
        q["series_index"] = pd.to_numeric(q.series_index, errors="coerce")

    return q


def band_code(f):
    f = np.asarray(f, dtype=float)
    out = np.full(f.shape, -1, dtype=np.int8)
    out[(f >= 71.0) & (f <= 76.0)] = 0
    out[(f >= 81.0) & (f <= 86.0)] = 1
    return out


def length_code(L):
    L = np.asarray(L, dtype=float)
    out = np.full(L.shape, -1, dtype=np.int8)
    for i, (_, lo, hi) in enumerate(LENGTH_BINS):
        out[(L >= lo) & (L < hi)] = i
    return out


def sensitivity_code(x):
    x = np.asarray(x, dtype=float)
    out = np.full(x.shape, -1, dtype=np.int8)
    for i, (_, lo, hi) in enumerate(SENS_BINS):
        out[np.isfinite(x) & (x >= lo) & (x < hi)] = i
    return out


def metrics(ref, est):
    ref = np.asarray(ref, dtype=float)
    est = np.asarray(est, dtype=float)
    m = np.isfinite(ref) & np.isfinite(est)
    ref = ref[m]
    est = est[m]
    if len(ref) == 0:
        return dict(N=0, PCC=np.nan, RMSE=np.nan, Bias=np.nan, MAE=np.nan)
    e = est - ref
    pcc = np.nan
    if len(ref) >= 2 and np.std(ref) > 0 and np.std(est) > 0:
        pcc = float(np.corrcoef(ref, est)[0, 1])
    return dict(
        N=int(len(ref)),
        PCC=pcc,
        RMSE=float(np.sqrt(np.mean(e * e))),
        Bias=float(np.mean(e)),
        MAE=float(np.mean(np.abs(e))),
    )


def savefig(fig, stem: Path, dpi: int):
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def panel_hex(ax, ref, est, title, max_wvd, gridsize, min_pairs):
    met = metrics(ref, est)
    hb = None
    if met["N"] >= min_pairs:
        hb = ax.hexbin(
            ref, est,
            gridsize=gridsize,
            extent=(0, max_wvd, 0, max_wvd),
            mincnt=1,
            cmap=HEX_CMAP,
        )
    ax.plot([0, max_wvd], [0, max_wvd], ls="--", lw=0.8, color="0.45")
    ax.set_xlim(0, max_wvd)
    ax.set_ylim(0, max_wvd)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=8.5)
    ax.tick_params(labelsize=7)
    ax.grid(False)
    txt = (
        f"N={met['N']:,}\n"
        f"PCC={met['PCC']:.3f}\n"
        f"RMSE={met['RMSE']:.2f}\n"
        f"Bias={met['Bias']:+.2f}"
    )
    ax.text(
        0.04, 0.96, txt,
        transform=ax.transAxes,
        va="top", ha="left",
        fontsize=6.7,
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", alpha=0.88, lw=0.5),
    )
    return hb, met


def harmonize_hexbin_norm(hbs):
    vals = []
    for hb in hbs:
        if hb is not None:
            a = np.asarray(hb.get_array(), dtype=float)
            if a.size:
                vals.append(np.nanmax(a))
    vmax = max(vals) if vals else 1.0
    vmax = max(1.0, float(vmax))
    norm = LogNorm(vmin=1.0, vmax=vmax)
    for hb in hbs:
        if hb is not None:
            hb.set_norm(norm)
    return norm


def add_shared_colorbar(fig, axes, hbs, label="Count"):
    valid = [h for h in hbs if h is not None]
    if not valid:
        return
    norm = harmonize_hexbin_norm(hbs)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=HEX_CMAP)
    sm.set_array([])
    fig.colorbar(sm, ax=np.asarray(axes).ravel().tolist(), shrink=0.84, pad=0.015, label=label)


def build_or_load_pairs(args, helper, q, v6_path, out_dir):
    branch = args.branch.lower()
    npz_path = out_dir / f"validation_pairs_{branch}.npz"
    links_csv = out_dir / "validation_population_links.csv"

    if args.reuse_pairs and npz_path.is_file() and links_csv.is_file():
        print(f"Reusing exact pair cache: {npz_path}")
        with np.load(npz_path) as z:
            ref = z["ref"].astype(np.float32, copy=False)
            est = z["est"].astype(np.float32, copy=False)
            link_idx = z["link_idx"].astype(np.int32, copy=False)
        links = pd.read_csv(links_csv)
        return ref, est, link_idx, links

    v6 = xr.open_dataset(v6_path, cache=False, decode_times=True)
    hum = xr.open_dataset(args.humidity_nc, cache=False)
    temp = xr.open_dataset(args.temperature_nc, cache=False)

    try:
        if "H" not in hum or "T" not in temp:
            raise RuntimeError("Expected H in humidity file and T in temperature file")

        # Exact V6 series mapping.
        if "series_index" in q.columns and q.series_index.notna().all():
            series_map = {
                (r.cml_id, r.sublink_id): int(r.series_index)
                for r in q.itertuples()
            }
        else:
            vd = helper.load_diag(args.v6_diagnostics_csv)
            series_map = {
                (normalize_id(r.cml_id), normalize_id(r.sublink_id)): int(r.series_index)
                for r in vd.itertuples()
            }

        qkeys = [(r.cml_id, r.sublink_id) for r in q.itertuples()]
        qkeys = [k for k in qkeys if k in series_map]
        if not qkeys:
            raise RuntimeError("No quantization links map to V6 series")

        # Reuse the exact station mapping from the already validated V6 run.
        if args.station_assignments_csv.is_file():
            assignments = pd.read_csv(args.station_assignments_csv)
            assignments["cml_id"] = assignments.cml_id.map(normalize_id)
            print(f"Using existing station assignments: {args.station_assignments_csv}")
        else:
            print("Station assignment CSV absent; recomputing with validation helper.")
            geom = helper.raw_geometry(args.raw_cml_nc)
            stations = helper.common_met_stations(hum, temp)
            cids = list(dict.fromkeys(k[0] for k in qkeys))
            assignments = helper.assign_stations(
                cids, geom, stations, args.max_station_distance_km
            )
            assignments.to_csv(args.station_assignments_csv, index=False)

        amap = assignments.set_index("cml_id")
        t5 = helper.decode_retrieval_time(v6)
        t10, _ = helper.strict_5_to_10(
            t5, np.zeros(len(t5)), np.ones(len(t5), dtype=np.uint8)
        )

        truth_cache = {}

        def truth_for_station(station_id):
            sid = str(station_id)
            if sid in truth_cache:
                return truth_cache[sid]
            hs = pd.Series(
                np.asarray(hum["H"].sel(station=sid).values, dtype=float),
                index=pd.DatetimeIndex(pd.to_datetime(hum["H"].sel(station=sid)["time"].values)),
            )
            ts = pd.Series(
                np.asarray(temp["T"].sel(station=sid).values, dtype=float),
                index=pd.DatetimeIndex(pd.to_datetime(temp["T"].sel(station=sid)["time"].values)),
            )
            hv = hs.reindex(t10).to_numpy(float)
            tv = ts.reindex(t10).to_numpy(float)
            rho = helper.rh_t_to_wvd(hv, tv)
            truth_cache[sid] = rho
            return rho

        method_names = [str(x).lower() for x in helper.METHODS]
        if branch not in method_names:
            raise RuntimeError(f"Branch {branch} absent from helper METHODS={helper.METHODS}")
        mi = method_names.index(branch)

        ref_parts = []
        est_parts = []
        link_parts = []
        link_rows = []

        qidx = q.set_index(["cml_id", "sublink_id"])
        total = len(qkeys)

        for kk, key in enumerate(qkeys, start=1):
            cid, sid = key
            if cid not in amap.index:
                continue
            ar = amap.loc[cid]
            if isinstance(ar, pd.DataFrame):
                ar = ar.iloc[0]
            if not bool(ar.get("within_validation_radius", False)):
                continue

            rowq = qidx.loc[key]
            if isinstance(rowq, pd.DataFrame):
                rowq = rowq.iloc[0]

            vi = series_map[key]
            w5 = np.asarray(
                v6["wvd_internal_anchor"].isel(method=mi, series=vi).values,
                dtype=float,
            )
            s5 = np.asarray(
                v6["retrieval_status"].isel(method=mi, series=vi).values
            )
            tt, w10 = helper.strict_5_to_10(t5, w5, s5)
            if not tt.equals(t10):
                raise RuntimeError("Internal strict 10-min time mismatch")

            truth = truth_for_station(str(ar["station_id"]))
            common = np.isfinite(truth) & np.isfinite(w10)
            n = int(common.sum())
            if n == 0:
                continue

            link_i = len(link_rows)
            ref_parts.append(truth[common].astype(np.float32))
            est_parts.append(w10[common].astype(np.float32))
            link_parts.append(np.full(n, link_i, dtype=np.int32))

            link_rows.append({
                "link_idx": link_i,
                "series_index": vi,
                "cml_id": cid,
                "sublink_id": sid,
                "station_id": str(ar["station_id"]),
                "station_distance_to_path_km": float(ar["station_distance_to_path_km"]),
                "frequency_ghz": float(rowq["frequency_ghz"]),
                "length_km": float(rowq["length_km"]),
                "tl_quantization_db": float(rowq["tl_quantization_db"]),
                "specific_attenuation_quantization_db_per_km": float(rowq["specific_attenuation_quantization_db_per_km"]),
                "wvd_quantization_resolution_ref20c_rho10_gm3": float(rowq["wvd_quantization_resolution_ref20c_rho10_gm3"]),
                "wvd_quantization_resolution_standard_min_gm3": float(rowq["wvd_quantization_resolution_standard_min_gm3"]),
                "wvd_quantization_resolution_standard_median_gm3": float(rowq["wvd_quantization_resolution_standard_median_gm3"]),
                "wvd_quantization_resolution_standard_max_gm3": float(rowq["wvd_quantization_resolution_standard_max_gm3"]),
                "n_pairs": n,
            })

            if kk % 100 == 0 or kk == total:
                print(f"[{kk:4d}/{total:4d}] strict V6-{branch.upper()} validation pairs")

        ref = np.concatenate(ref_parts) if ref_parts else np.empty(0, np.float32)
        est = np.concatenate(est_parts) if est_parts else np.empty(0, np.float32)
        link_idx = np.concatenate(link_parts) if link_parts else np.empty(0, np.int32)
        links = pd.DataFrame(link_rows)

        np.savez_compressed(npz_path, ref=ref, est=est, link_idx=link_idx)
        links.to_csv(links_csv, index=False)

        return ref, est, link_idx, links

    finally:
        v6.close()
        hum.close()
        temp.close()


def add_codes(links):
    out = links.copy()
    out["frequency_code"] = band_code(out.frequency_ghz.to_numpy(float))
    out["length_code"] = length_code(out.length_km.to_numpy(float))
    out["sensitivity_code"] = sensitivity_code(
        out.wvd_quantization_resolution_ref20c_rho10_gm3.to_numpy(float)
    )
    return out


def pair_mask_from_link_mask(link_idx, link_mask):
    lm = np.asarray(link_mask, dtype=bool)
    return lm[link_idx]


def summarize_group(ref, est, link_idx, links, mask_links, group_type, label, extra=None):
    pm = pair_mask_from_link_mask(link_idx, mask_links)
    met = metrics(ref[pm], est[pm])
    x = links[np.asarray(mask_links, dtype=bool)]
    row = {
        "group_type": group_type,
        "group": label,
        "n_directional_links": int(len(x)),
        **met,
        "median_length_km": float(x.length_km.median()) if len(x) else np.nan,
        "median_frequency_ghz": float(x.frequency_ghz.median()) if len(x) else np.nan,
        "median_qtl_db": float(x.tl_quantization_db.median()) if len(x) else np.nan,
        "median_qgamma_db_per_km": float(x.specific_attenuation_quantization_db_per_km.median()) if len(x) else np.nan,
        "median_delta_rho_q_gm3": float(x.wvd_quantization_resolution_ref20c_rho10_gm3.median()) if len(x) else np.nan,
    }
    if extra:
        row.update(extra)
    return row, pm


def plot_single_axis_family(ref, est, link_idx, links, groups, title, stem, args):
    n = len(groups)
    ncols = min(3, n)
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.35 * ncols, 3.25 * nrows),
        sharex=True, sharey=True,
    )
    axes = np.atleast_1d(axes).ravel()
    hbs = []
    rows = []

    for ax, (label, mask_links, extra) in zip(axes, groups):
        row, pm = summarize_group(
            ref, est, link_idx, links, mask_links,
            stem.name, label, extra,
        )
        hb, _ = panel_hex(
            ax, ref[pm], est[pm], label,
            args.scatter_max_gm3, args.gridsize, args.min_panel_pairs,
        )
        hbs.append(hb)
        rows.append(row)

    for ax in axes[len(groups):]:
        ax.axis("off")

    for i, ax in enumerate(axes[:len(groups)]):
        if i % ncols == 0:
            ax.set_ylabel(r"V6 WVD [g m$^{-3}$]")
        if i // ncols == nrows - 1:
            ax.set_xlabel(r"CHMI WVD [g m$^{-3}$]")

    fig.suptitle(title, fontsize=11, y=1.01)
    add_shared_colorbar(fig, axes[:len(groups)], hbs)
    fig.tight_layout(rect=(0, 0, 0.96, 0.98))
    savefig(fig, stem, args.dpi)
    return rows


def plot_grid_family(ref, est, link_idx, links, row_groups, col_groups, title, stem, args, fixed_band=None):
    nr, nc = len(row_groups), len(col_groups)
    fig, axes = plt.subplots(
        nr, nc,
        figsize=(2.95 * nc, 2.72 * nr),
        sharex=True, sharey=True,
        squeeze=False,
    )
    hbs = []
    rows = []

    for ir, (rlabel, rmask, rextra) in enumerate(row_groups):
        for ic, (clabel, cmask, cextra) in enumerate(col_groups):
            mask = np.asarray(rmask, bool) & np.asarray(cmask, bool)
            if fixed_band is not None:
                mask &= links.frequency_code.to_numpy(int) == int(fixed_band)
            label = f"{rlabel} | {clabel}"
            extra = {}
            extra.update(rextra or {})
            extra.update(cextra or {})
            if fixed_band is not None:
                extra["frequency_band"] = FREQ_BANDS[fixed_band][0]

            row, pm = summarize_group(
                ref, est, link_idx, links, mask,
                stem.name, label, extra,
            )
            hb, _ = panel_hex(
                axes[ir, ic], ref[pm], est[pm], clabel if nr > 1 else label,
                args.scatter_max_gm3, args.gridsize, args.min_panel_pairs,
            )
            hbs.append(hb)
            rows.append(row)

            if ir == 0:
                axes[ir, ic].set_title(clabel, fontsize=8.2)
            if ic == 0:
                axes[ir, ic].set_ylabel(
                    f"{rlabel}\nV6 WVD [g m$^{{-3}}$]",
                    fontsize=7.5,
                )
            if ir == nr - 1:
                axes[ir, ic].set_xlabel(
                    r"CHMI WVD [g m$^{-3}$]", fontsize=7.5
                )

    fig.suptitle(title, fontsize=11, y=1.005)
    add_shared_colorbar(fig, axes, hbs)
    fig.tight_layout(rect=(0, 0, 0.97, 0.985))
    savefig(fig, stem, args.dpi)
    return rows


def main():
    args = parse_args()
    t0 = time.time()

    args.root = args.root.expanduser().resolve()
    args.v6_nc = args.v6_nc.expanduser()
    args.v6_diagnostics_csv = args.v6_diagnostics_csv.expanduser().resolve()
    args.quant_csv = args.quant_csv.expanduser().resolve()
    args.raw_cml_nc = args.raw_cml_nc.expanduser().resolve()
    args.humidity_nc = args.humidity_nc.expanduser().resolve()
    args.temperature_nc = args.temperature_nc.expanduser().resolve()
    args.validation_helper_script = args.validation_helper_script.expanduser().resolve()
    args.station_assignments_csv = args.station_assignments_csv.expanduser().resolve()
    args.out_dir = args.out_dir.expanduser().resolve()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    figdir = args.out_dir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    helper = load_module(args.validation_helper_script, "v6_validation_helper")
    v6_path = resolve_v6_path(args.v6_nc, args.root)
    q = load_quant(args.quant_csv)

    print("=" * 112)
    print("V6 BINNED VALIDATION — QUANTIZATION SENSITIVITY x FREQUENCY x LENGTH")
    print("=" * 112)
    print("V6 NetCDF      :", v6_path)
    print("Quantization   :", args.quant_csv)
    print("Branch         :", args.branch)
    print("CHMI humidity  : VALIDATION ONLY")
    print("Retrieval      : UNCHANGED V6")
    print("Sensitivity    : link-derived quantization-limited Delta rho_q")
    print("=" * 112)

    ref, est, link_idx, links = build_or_load_pairs(
        args, helper, q, v6_path, args.out_dir
    )
    links = add_codes(links)
    links.to_csv(args.out_dir / "validation_population_links.csv", index=False)

    if len(ref) == 0:
        raise RuntimeError("No strict V6/CHMI validation pairs assembled")

    print(
        f"Strict pairs: {len(ref):,} | directional links: {len(links):,} | "
        f"q sensitivity available: "
        f"{int((links.sensitivity_code >= 0).sum()):,}/{len(links):,}"
    )

    all_rows = []

    sens_groups = []
    for i, (lab, lo, hi) in enumerate(SENS_BINS):
        m = links.sensitivity_code.to_numpy(int) == i
        sens_groups.append((lab + r" g m$^{-3}$", m, {"sensitivity_bin": lab}))
    all_rows += plot_single_axis_family(
        ref, est, link_idx, links, sens_groups,
        f"V6 {args.branch} validation by quantization-limited WVD sensitivity bin",
        figdir / f"01_hexbin_by_sensitivity_{args.branch.lower()}", args,
    )

    freq_groups = []
    for i, (lab, _, _) in enumerate(FREQ_BANDS):
        m = links.frequency_code.to_numpy(int) == i
        freq_groups.append((lab, m, {"frequency_band": lab}))
    all_rows += plot_single_axis_family(
        ref, est, link_idx, links, freq_groups,
        f"V6 {args.branch} validation by E-band frequency range",
        figdir / f"02_hexbin_by_frequency_{args.branch.lower()}", args,
    )

    length_groups = []
    for i, (lab, _, _) in enumerate(LENGTH_BINS):
        m = links.length_code.to_numpy(int) == i
        length_groups.append((lab, m, {"length_bin": lab}))
    all_rows += plot_single_axis_family(
        ref, est, link_idx, links, length_groups,
        f"V6 {args.branch} validation by disjoint path-length bin",
        figdir / f"03_hexbin_by_length_{args.branch.lower()}", args,
    )

    all_rows += plot_grid_family(
        ref, est, link_idx, links,
        freq_groups, length_groups,
        f"V6 {args.branch}: frequency range x path length",
        figdir / f"04_hexbin_frequency_x_length_{args.branch.lower()}",
        args,
    )

    all_rows += plot_grid_family(
        ref, est, link_idx, links,
        freq_groups, sens_groups,
        f"V6 {args.branch}: frequency range x quantization-limited WVD sensitivity",
        figdir / f"05_hexbin_frequency_x_sensitivity_{args.branch.lower()}",
        args,
    )

    # Full joint diagnostic: within each frequency range, rows=length,
    # columns=sensitivity. This directly answers whether apparent sensitivity
    # effects remain after controlling for path length.
    for bi, (band_lab, _, _) in enumerate(FREQ_BANDS):
        safe = "71_76GHz" if bi == 0 else "81_86GHz"
        all_rows += plot_grid_family(
            ref, est, link_idx, links,
            length_groups, sens_groups,
            f"V6 {args.branch}: {band_lab} — length x quantization sensitivity",
            figdir / f"0{6+bi}_hexbin_joint_{safe}_{args.branch.lower()}",
            args,
            fixed_band=bi,
        )

    metrics_df = pd.DataFrame(all_rows)
    metrics_df.to_csv(args.out_dir / "binned_validation_metrics.csv", index=False)

    # Link-population audit independent of sample counts.
    pop_rows = []
    for bi, (band_lab, _, _) in enumerate(FREQ_BANDS):
        for li, (len_lab, _, _) in enumerate(LENGTH_BINS):
            for si, (sens_lab, _, _) in enumerate(SENS_BINS):
                m = (
                    (links.frequency_code.to_numpy(int) == bi)
                    & (links.length_code.to_numpy(int) == li)
                    & (links.sensitivity_code.to_numpy(int) == si)
                )
                x = links[m]
                pop_rows.append({
                    "frequency_band": band_lab,
                    "length_bin": len_lab,
                    "sensitivity_bin_gm3": sens_lab,
                    "n_directional_links": int(len(x)),
                    "n_pairs_sum": int(x.n_pairs.sum()) if len(x) else 0,
                    "median_qtl_db": float(x.tl_quantization_db.median()) if len(x) else np.nan,
                    "median_qgamma_db_per_km": float(x.specific_attenuation_quantization_db_per_km.median()) if len(x) else np.nan,
                    "median_delta_rho_q_gm3": float(x.wvd_quantization_resolution_ref20c_rho10_gm3.median()) if len(x) else np.nan,
                })
    pd.DataFrame(pop_rows).to_csv(args.out_dir / "population_counts.csv", index=False)

    overall = metrics(ref, est)
    summary = {
        "branch": args.branch,
        "v6_nc": str(v6_path),
        "quantization_csv": str(args.quant_csv),
        "n_directional_links": int(len(links)),
        "n_strict_pairs": int(len(ref)),
        "n_links_with_quantization_sensitivity": int((links.sensitivity_code >= 0).sum()),
        "overall_metrics": overall,
        "sensitivity_bins_gm3": [x[0] for x in SENS_BINS],
        "length_bins_km": [x[0] for x in LENGTH_BINS],
        "frequency_bands": [x[0] for x in FREQ_BANDS],
        "humidity_role": "validation only",
        "retrieval_modified": False,
        "elapsed_seconds": time.time() - t0,
    }
    (args.out_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True), encoding="utf-8"
    )

    print("\n" + "=" * 112)
    print("DONE — V6 BINNED VALIDATION")
    print("=" * 112)
    print(
        f"Overall {args.branch}: N={overall['N']:,} "
        f"PCC={overall['PCC']:.3f} RMSE={overall['RMSE']:.2f} "
        f"Bias={overall['Bias']:+.2f}"
    )
    print("Output:", args.out_dir)
    print("Metrics:", args.out_dir / "binned_validation_metrics.csv")
    print("Figures:", figdir)
    print("Pair cache:", args.out_dir / f"validation_pairs_{args.branch.lower()}.npz")
    print(f"Elapsed: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
