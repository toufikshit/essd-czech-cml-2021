#!/usr/bin/env python3
"""
CML wet-dry detection by rolling standard deviation (RSD) + contingency table
=============================================================================

Rain DETECTION companion to the rainfall-estimation pipeline. It reproduces the
same data conventions (same CML/radar files, same radar unit rule, same wet
threshold and forward-fill) and evaluates the classic rolling-std wet/dry
detector against the radar reference over ALL links and ALL times.

Detector (as requested):
    wet_rsd = tl.rolling(time=WINDOW, center=True).std() > THR      # per sublink
then the SUBLINKS ARE AVERAGED into a single per-link decision (default:
average the two sublink rolling-std series, then threshold; see --sublink-combine).

Here `tl` (total loss) is the CML `trsl` variable (= TSL - RSL). If `trsl` is
absent it is rebuilt as tx - rx, exactly like the estimation pipeline.

Reference (identical to the estimation pipeline):
    radar_mmph = radar_raw * 12
    wet_5min   = radar_mmph > 0.00833
    wet is forward-filled from the 5-min radar grid onto the 1-min CML grid.

Contingency (2x2), aggregated over every valid (link, minute):
    Hit  (TP): CML wet AND radar wet
    Miss (FN): CML dry AND radar wet
    False alarm (FP): CML wet AND radar dry
    Correct negative (TN): CML dry AND radar dry

Outputs:
    contingency_table.csv            2x2 counts (the detection "dataset")
    detection_metrics.csv            POD/FAR/POFD/CSI/ACC/bias/HSS/MCC + fractions
    per_link_detection_metrics.csv   the same per link
    confusion_matrix_rsd.png/.pdf    publishable matrix only, NO title/subtitle/side metrics

Example
-------
python cml_wet_dry_rsd_contingency.py \
  --cml-nc cml_preprocessed_v3.nc \
  --radar-nc RADAR_S2021E2021_5MIN_CZ.nc \
  --out-dir cml_rain_detection \
  --tstart "2021-04-10 00:00:00" --tend "2021-08-16 23:59:00" \
  --rsd-window 60 --rsd-thr 0.8 --sublink-combine rsd_mean
"""

from __future__ import annotations

import argparse
import math
import time
import warnings
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
except Exception:  # pragma: no cover
    def tqdm(x, **kwargs):
        return x

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ============================================================
# Helpers copied from the estimation pipeline (same conventions)
# ============================================================
def find_dim(ds: xr.Dataset, candidates: Sequence[str]) -> Optional[str]:
    for c in candidates:
        if c in ds.dims:
            return c
    return None


def variable_exists(ds: xr.Dataset, name: str) -> bool:
    return name in ds.variables or name in ds.data_vars


def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def time_slice_indices(times: pd.DatetimeIndex, tstart: pd.Timestamp, tend: pd.Timestamp) -> Tuple[int, int]:
    vals = times.values
    i0 = int(np.searchsorted(vals, np.datetime64(tstart), "left"))
    i1 = int(np.searchsorted(vals, np.datetime64(tend), "right"))
    return i0, i1


def radar_values_to_mm_h(raw: np.ndarray, factor: float = 12.0) -> np.ndarray:
    """Keep the exact working-script radar convention: radar_mmph = radar_raw * 12."""
    raw = np.asarray(raw, dtype="float64")
    return raw * float(factor)


def forward_fill_radar_to_cml_time(
    radar_time: pd.DatetimeIndex,
    cml_time: pd.DatetimeIndex,
    values_5min: np.ndarray,
) -> np.ndarray:
    """Forward-fill a 5-min radar array onto the 1-min CML grid (last value at or before).

    Generalised from the pipeline's forward_fill_radar_wet_to_cml_time so it works
    for boolean wet masks AND for the finite-data mask.
    """
    idx = np.searchsorted(radar_time.values, cml_time.values, "right") - 1
    pre_first = idx < 0
    idx = np.clip(idx, 0, None)
    out = values_5min[idx]
    out[pre_first] = False if out.dtype == bool else np.nan
    return out


def normalize_2d_time_sublink(da: xr.DataArray, sublink_dim: str) -> np.ndarray:
    """Return a (time, sublink) float array, matching the pipeline helper."""
    if sublink_dim in da.dims:
        return da.transpose("time", sublink_dim).values
    return da.values[:, None]


# ============================================================
# Detector core (pure numpy/pandas, unit-testable)
# ============================================================
def rolling_std_1d(x: np.ndarray, window: int, center: bool, min_periods: Optional[int], ddof: int) -> np.ndarray:
    """Rolling standard deviation reproducing xarray's
    `da.rolling(time=window, center=center).std()` (ddof=0, full-window min_periods).
    """
    mp = window if min_periods is None else int(min_periods)
    s = pd.Series(np.asarray(x, dtype="float64"))
    return s.rolling(window=int(window), center=bool(center), min_periods=mp).std(ddof=ddof).to_numpy()


def cml_wet_from_trsl(
    trsl_ts: np.ndarray,           # (T, S) total loss over the processing window
    window: int,
    center: bool,
    min_periods: Optional[int],
    ddof: int,
    thr: float,
    combine: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (wet, valid) 1-D boolean arrays of length T from a (T, S) trsl block.

    combine:
      rsd_mean : mean of per-sublink rolling-std, then > thr        (default)
      tl_mean  : rolling-std of the sublink-mean trsl, then > thr
      any      : wet if ANY sublink's rolling-std > thr
      all      : wet if ALL finite sublinks' rolling-std > thr
    """
    trsl_ts = np.asarray(trsl_ts, dtype="float64")
    if trsl_ts.ndim == 1:
        trsl_ts = trsl_ts[:, None]
    T, S = trsl_ts.shape

    if combine == "tl_mean":
        tl_mean = np.nanmean(trsl_ts, axis=1)
        std_c = rolling_std_1d(tl_mean, window, center, min_periods, ddof)
        valid = np.isfinite(std_c)
        wet = np.where(valid, std_c > thr, False)
        return wet, valid

    # per-sublink rolling std
    std_s = np.column_stack([
        rolling_std_1d(trsl_ts[:, s], window, center, min_periods, ddof) for s in range(S)
    ])  # (T, S)

    if combine == "rsd_mean":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            std_c = np.nanmean(std_s, axis=1)
        valid = np.isfinite(std_c)
        wet = np.where(valid, std_c > thr, False)
        return wet, valid

    if combine == "any":
        valid = np.isfinite(std_s).any(axis=1)
        wet_s = std_s > thr
        wet = np.where(np.isfinite(std_s), wet_s, False).any(axis=1)
        return np.where(valid, wet, False), valid

    if combine == "all":
        finite = np.isfinite(std_s)
        valid = finite.any(axis=1)
        # wet only if every finite sublink exceeds thr
        wet_s = (std_s > thr) | (~finite)          # ignore NaN sublinks in the AND
        wet = wet_s.all(axis=1)
        return np.where(valid, wet, False), valid

    raise ValueError(f"Unknown --sublink-combine {combine!r}")


def contingency_counts(cml_wet: np.ndarray, radar_wet: np.ndarray, valid: np.ndarray) -> Tuple[int, int, int, int]:
    """Return (TP, FP, FN, TN) over the valid samples."""
    v = np.asarray(valid, dtype=bool)
    c = np.asarray(cml_wet, dtype=bool) & v
    r = np.asarray(radar_wet, dtype=bool) & v
    tp = int(np.sum(v & c & r))
    fp = int(np.sum(v & c & ~r))
    fn = int(np.sum(v & ~c & r))
    tn = int(np.sum(v & ~c & ~r))
    return tp, fp, fn, tn


def detection_metrics(tp: int, fp: int, fn: int, tn: int) -> Dict[str, float]:
    """Standard categorical verification scores. NaN where undefined."""
    n = tp + fp + fn + tn

    def div(a, b):
        return float(a) / float(b) if b else float("nan")

    pod = div(tp, tp + fn)               # hit rate / recall / probability of detection
    far = div(fp, tp + fp)               # false alarm ratio
    pofd = div(fp, fp + tn)              # probability of false detection
    sr = div(tp, tp + fp)               # success ratio = 1 - FAR
    csi = div(tp, tp + fp + fn)          # critical success index / threat score
    acc = div(tp + tn, n)                # accuracy
    bias = div(tp + fp, tp + fn)         # frequency bias

    hss_den = (tp + fn) * (fn + tn) + (tp + fp) * (fp + tn)
    hss = div(2.0 * (tp * tn - fp * fn), hss_den)

    mcc_den = math.sqrt(float(tp + fp) * float(tp + fn) * float(tn + fp) * float(tn + fn))
    mcc = (tp * tn - fp * fn) / mcc_den if mcc_den > 0 else float("nan")

    return {
        "N": int(n), "TP": int(tp), "FP": int(fp), "FN": int(fn), "TN": int(tn),
        "wet_frac_radar": div(tp + fn, n),
        "wet_frac_cml": div(tp + fp, n),
        "POD": pod, "FAR": far, "POFD": pofd, "SR": sr, "CSI": csi,
        "ACC": acc, "freq_bias": bias, "HSS": hss, "MCC": mcc,
    }


# ============================================================
# Publishable confusion-matrix plot (no title / no subtitle)
# ============================================================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.8,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

CATEGORY = {
    (0, 0): "Hit",             # CML wet, radar wet
    (0, 1): "False alarm",     # CML wet, radar dry
    (1, 0): "Miss",            # CML dry, radar wet
    (1, 1): "Correct negative"  # CML dry, radar dry
}


def make_confusion_plot(
    tp: int, fp: int, fn: int, tn: int,
    metrics: Dict[str, float],
    out_png: Path,
    out_pdf: Path,
    cmap_name: str = "Blues",
    cell_norm: str = "reference",
    show_metrics_box: bool = False,
    show_colorbar: bool = False,
):
    """Clean 2x2 contingency matrix only.

    Rows = CML rolling-std detection: Wet / Dry
    Columns = Radar reference: Wet / Dry

    The plot intentionally does NOT include the side metrics panel.
    Cell annotations show:
        class name
        raw count
        percentage within the radar reference class

    cell_norm controls the cell colour only:
      reference : each radar column normalised to its own total
      total     : fraction of the grand total
      count     : raw counts on a log colour scale
    """
    counts = np.array([[tp, fp],
                       [fn, tn]], dtype="float64")
    total = counts.sum()

    if cell_norm == "reference":
        col_tot = counts.sum(axis=0, keepdims=True)
        color_vals = np.divide(
            counts, col_tot,
            out=np.zeros_like(counts, dtype="float64"),
            where=col_tot > 0
        )
        norm = None
        cbar_label = "Fraction within radar class"
        vmin, vmax = 0.0, 1.0
    elif cell_norm == "total":
        color_vals = counts / total if total > 0 else counts
        norm = None
        cbar_label = "Fraction of all samples"
        vmin, vmax = 0.0, float(np.nanmax(color_vals)) if total > 0 else 1.0
    elif cell_norm == "count":
        color_vals = counts.copy()
        positive = counts[counts > 0]
        norm = LogNorm(
            vmin=max(1.0, positive.min() if positive.size else 1.0),
            vmax=max(1.0, counts.max())
        )
        cbar_label = "Count"
        vmin = vmax = None
    else:
        raise ValueError(f"Unknown cell_norm {cell_norm!r}")

    fig, ax = plt.subplots(figsize=(5.8, 5.2))

    if norm is None:
        im = ax.imshow(color_vals, cmap=cmap_name, vmin=vmin, vmax=vmax, aspect="equal")
    else:
        im = ax.imshow(color_vals, cmap=cmap_name, norm=norm, aspect="equal")

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Wet", "Dry"], fontsize=12)
    ax.set_yticklabels(["Wet", "Dry"], fontsize=12, rotation=90, va="center")

    ax.set_xlabel("Radar reference", fontsize=13, labelpad=8)
    ax.set_ylabel("CML rolling-σ detection", fontsize=13, labelpad=8)

    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xticks([-0.5, 0.5, 1.5], minor=True)
    ax.set_yticks([-0.5, 0.5, 1.5], minor=True)
    ax.grid(which="minor", color="white", linewidth=2.5)
    ax.tick_params(which="minor", length=0)

    col_tot = counts.sum(axis=0)
    for r in range(2):
        for c in range(2):
            cnt = int(counts[r, c])
            pct_ref = (counts[r, c] / col_tot[c] * 100.0) if col_tot[c] > 0 else float("nan")

            cval = color_vals[r, c]
            if norm is None:
                cmax = vmax if (vmax is not None and vmax > 0) else np.nanmax(color_vals)
                frac = (cval / cmax) if (cmax and np.isfinite(cmax) and cmax > 0) else 0.0
            else:
                frac = norm(cval) if cval > 0 else 0.0

            txt_color = "white" if frac > 0.55 else "0.12"

            ax.text(
                c, r - 0.18, CATEGORY[(r, c)],
                ha="center", va="center",
                fontsize=12.0, fontweight="bold",
                color=txt_color
            )
            ax.text(
                c, r + 0.08, f"{cnt:,}",
                ha="center", va="center",
                fontsize=13.5,
                color=txt_color
            )
            ax.text(
                c, r + 0.30,
                f"{pct_ref:.1f}% of radar {'wet' if c == 0 else 'dry'}",
                ha="center", va="center",
                fontsize=9.2,
                color=txt_color
            )

    if show_colorbar:
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(cbar_label, fontsize=10)
        cbar.ax.tick_params(labelsize=9)

    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


# ============================================================
# Driver
# ============================================================
def run(args) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ds_c = xr.open_dataset(args.cml_nc)
    ds_r = xr.open_dataset(args.radar_nc)
    try:
        if "cml_id" in ds_c.coords:
            ds_c = ds_c.assign_coords(cml_id=ds_c["cml_id"].astype(str))
        if "cml_id" in ds_r.coords:
            ds_r = ds_r.assign_coords(cml_id=ds_r["cml_id"].astype(str))

        sublink_dim = find_dim(ds_c, ["sublink_id", "sublink"])
        if sublink_dim is None:
            raise ValueError("Cannot find a sublink dimension in the CML dataset.")

        cml_ids_all = ds_c["cml_id"].values.astype(str)
        radar_ids = ds_r["cml_id"].values.astype(str)
        radar_pos = {rid: i for i, rid in enumerate(radar_ids)}

        full_time = pd.DatetimeIndex(ds_c["time"].values)
        radar_time = pd.DatetimeIndex(ds_r["time"].values)
        t0 = pd.Timestamp(args.tstart)
        t1 = pd.Timestamp(args.tend)
        # Pad the processing window so the centred rolling std is defined at the
        # edges of the output window (need ~window/2 extra minutes on each side).
        pad = pd.Timedelta(minutes=int(args.rsd_window))
        p0 = t0 - pad

        proc0, proc1 = time_slice_indices(full_time, p0, t1 + pad)
        out0, out1 = time_slice_indices(full_time, t0, t1)
        if proc1 <= proc0 or out1 <= out0:
            raise ValueError("Empty CML time window. Check --tstart/--tend.")
        proc_time = full_time[proc0:proc1]
        out_time = full_time[out0:out1]
        out_rel0 = out0 - proc0
        out_rel1 = out1 - proc0

        rout0, rout1 = time_slice_indices(radar_time, t0, t1)
        if rout1 <= rout0:
            raise ValueError("Empty radar time window. Check --tstart/--tend and radar file.")
        radar_out_time = radar_time[rout0:rout1]

        has_trsl = variable_exists(ds_c, "trsl")

        # Link selection: present in radar AND long enough (mirrors the pipeline).
        selected = []
        for i, cid in enumerate(cml_ids_all):
            if args.only_cml_id is not None and str(cid) != str(args.only_cml_id):
                continue
            if cid not in radar_pos:
                continue
            length = safe_float(ds_c["length"].isel({"cml_id": i}).values) if variable_exists(ds_c, "length") else np.inf
            if not np.isfinite(length) or length < args.min_link_km:
                continue
            selected.append(i)
            if args.max_links is not None and len(selected) >= int(args.max_links):
                break
        if not selected:
            raise RuntimeError("No CML links selected. Check IDs, --min-link-km, or --only-cml-id.")

        print("\nRain-detection (rolling-\u03c3) configuration")
        print(f"  CML file            : {args.cml_nc}")
        print(f"  Radar file          : {args.radar_nc}")
        print(f"  tl source           : {'trsl' if has_trsl else 'tx - rx'}")
        print(f"  Selected links      : {len(selected)} / {len(cml_ids_all)}")
        print(f"  Output window        : {out_time[0]} -> {out_time[-1]} ({len(out_time)} min)")
        print(f"  Rolling window/thr  : {args.rsd_window} min, \u03c3 > {args.rsd_thr:g}")
        print(f"  Sublink combine     : {args.sublink_combine}")
        print(f"  Radar convention    : radar_mmph = radar_raw * {args.radar_mmph_factor:g}; "
              f"wet if > {args.radar_wet_thr_5min:g}")

        TP = FP = FN = TN = 0
        per_link_rows: List[Dict] = []

        for c_index in tqdm(selected, desc="CML links"):
            cid = str(cml_ids_all[c_index])
            rpos = radar_pos[cid]

            # --- radar reference on the OUTPUT window, mapped to 1-min ---
            radar_raw_out = ds_r["rr"].isel({"cml_id": rpos, "time": slice(rout0, rout1)}).values
            radar_mmph = radar_values_to_mm_h(radar_raw_out, args.radar_mmph_factor)
            radar_finite5 = np.isfinite(radar_mmph)
            radar_wet5 = radar_finite5 & (radar_mmph > args.radar_wet_thr_5min)
            radar_wet_1m = forward_fill_radar_to_cml_time(radar_out_time, out_time, radar_wet5.copy())
            radar_ok_1m = forward_fill_radar_to_cml_time(radar_out_time, out_time, radar_finite5.copy())

            # --- CML rolling-std detection on the padded window, then slice to output ---
            if has_trsl:
                trsl_da = ds_c["trsl"].isel({"cml_id": c_index, "time": slice(proc0, proc1)})
            else:
                trsl_da = (ds_c["tx"].isel({"cml_id": c_index, "time": slice(proc0, proc1)})
                           - ds_c["rx"].isel({"cml_id": c_index, "time": slice(proc0, proc1)}))
            trsl_ts = normalize_2d_time_sublink(trsl_da, sublink_dim)  # (Tproc, S)

            cml_wet_proc, cml_valid_proc = cml_wet_from_trsl(
                trsl_ts, window=args.rsd_window, center=True,
                min_periods=args.rsd_min_periods, ddof=args.rsd_ddof,
                thr=args.rsd_thr, combine=args.sublink_combine,
            )
            cml_wet = cml_wet_proc[out_rel0:out_rel1]
            cml_valid = cml_valid_proc[out_rel0:out_rel1]

            # length guard so all three arrays align
            m = min(len(cml_wet), len(radar_wet_1m))
            valid = cml_valid[:m] & radar_ok_1m[:m].astype(bool)
            tp, fp, fn, tn = contingency_counts(cml_wet[:m], radar_wet_1m[:m], valid)

            TP += tp; FP += fp; FN += fn; TN += tn
            lm = detection_metrics(tp, fp, fn, tn)
            lm_row = {"cml_id": cid, **lm}
            per_link_rows.append(lm_row)

        # ---- Global contingency + metrics ----
        metrics = detection_metrics(TP, FP, FN, TN)

        cont_df = pd.DataFrame(
            [[TP, FP, TP + FP],
             [FN, TN, FN + TN],
             [TP + FN, FP + TN, TP + FP + FN + TN]],
            index=["CML wet", "CML dry", "total"],
            columns=["radar wet", "radar dry", "total"],
        )
        cont_csv = out_dir / "contingency_table.csv"
        cont_df.to_csv(cont_csv)

        metrics_csv = out_dir / "detection_metrics.csv"
        pd.DataFrame([metrics]).to_csv(metrics_csv, index=False)

        per_link_csv = out_dir / "per_link_detection_metrics.csv"
        pd.DataFrame(per_link_rows).to_csv(per_link_csv, index=False)

        png = out_dir / "confusion_matrix_rsd.png"
        pdf = out_dir / "confusion_matrix_rsd.pdf"
        make_confusion_plot(
            TP, FP, FN, TN, metrics, png, pdf,
            cmap_name=args.cmap, cell_norm=args.cell_norm,
            show_metrics_box=False,
            show_colorbar=args.colorbar,
        )

        print("\nContingency table (all links, all valid minutes)")
        print(cont_df.to_string())
        print("\nDetection skill")
        for k in ["N", "wet_frac_radar", "wet_frac_cml", "POD", "FAR", "POFD",
                  "CSI", "ACC", "freq_bias", "HSS", "MCC"]:
            v = metrics[k]
            print(f"  {k:16s}: {v:,}" if k == "N" else f"  {k:16s}: {v:.4f}")

        print("\nSaved")
        for p in (cont_csv, metrics_csv, per_link_csv, png, pdf):
            print(f"  {p}")

    finally:
        ds_c.close()
        ds_r.close()


# ============================================================
# CLI
# ============================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="CML rolling-std wet/dry detection and radar contingency table."
    )
    p.add_argument("--cml-nc", default="cml_preprocessed_v3.nc")
    p.add_argument("--radar-nc", default="RADAR_S2021E2021_5MIN_CZ.nc")
    p.add_argument("--out-dir", default="cml_rain_detection")
    p.add_argument("--tstart", default="2021-04-10 00:00:00")
    p.add_argument("--tend", default="2021-08-16 23:59:00")

    p.add_argument("--rsd-window", type=int, default=60, help="Rolling window in minutes (samples).")
    p.add_argument("--rsd-thr", type=float, default=0.8, help="Wet if combined rolling std > this.")
    p.add_argument("--rsd-min-periods", type=int, default=None,
                   help="Min samples in the window. Default = full window (xarray default).")
    p.add_argument("--rsd-ddof", type=int, default=0, help="Std ddof; 0 matches xarray .std().")
    p.add_argument("--sublink-combine", choices=["rsd_mean", "tl_mean", "any", "all"],
                   default="rsd_mean",
                   help="How to average/combine sublinks into one per-link decision.")

    p.add_argument("--radar-mmph-factor", type=float, default=12.0,
                   help="radar_mmph = radar_raw * factor (pipeline convention = 12).")
    p.add_argument("--radar-wet-thr-5min", type=float, default=0.00833,
                   help="Radar wet threshold after conversion (pipeline default 0.00833).")

    p.add_argument("--min-link-km", type=float, default=0.5)
    p.add_argument("--only-cml-id", default=None)
    p.add_argument("--max-links", type=int, default=None)

    p.add_argument("--cmap", default="Blues", help="Confusion-matrix colormap.")
    p.add_argument("--cell-norm", choices=["reference", "total", "count"], default="reference",
                   help="Cell COLOUR normalisation (annotations always show raw counts).")
    # Matrix-only plot by default. Metrics are still saved to CSV and printed,
    # but they are not shown on the figure.
    p.add_argument("--colorbar", action="store_true", default=False,
                   help="Add a colorbar to the matrix. Default: no colorbar.")
    return p.parse_args()


def main():
    t0 = time.perf_counter()
    args = parse_args()
    run(args)
    print(f"\nTotal elapsed: {time.perf_counter() - t0:.2f} s")


if __name__ == "__main__":
    main()
