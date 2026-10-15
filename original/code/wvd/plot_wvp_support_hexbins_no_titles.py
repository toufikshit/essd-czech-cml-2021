#!/usr/bin/env python3
"""
Article-style V6 water-vapour validation support hexbins.

Selection applied before plotting:
  * path length >= 4.0 km
  * quantization-limited WVD resolution Delta rho_q <= 0.20 g m-3
  * E-band only (71-76 and/or 81-86 GHz)

Outputs (no pair-level WVD/WVP data are exported):
  * wvp_support_all_E_band.png/.pdf
  * wvp_support_80GHz.png/.pdf
  * wvp_support_70GHz.png/.pdf
  * wvp_support_three_panel.png/.pdf
  * wvp_support_metrics.csv

The script only writes figures and an aggregated metrics table; it never writes
paired WVD values, timestamps, or link-level observation series.

No plot titles are drawn; population identity is carried by the output filename/caption.

The plot geometry follows the rainfall article style: top/right spines hidden,
left/bottom spines shifted outward and bounded from 0 to the plotting maximum,
outward major/minor ticks, 1:1 line, logarithmic hexagon counts, and metrics box.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd


# ------------------------- science/QC settings -------------------------
MIN_LENGTH_KM = 4.0
MAX_DELTA_RHO_Q_GM3 = 0.20
PLOT_MAX_GM3 = 35.0
GRIDSIZE = 55
DPI = 300
CMAP = "jet"                 # same default as the rainfall article plotting script

# The stored quantity is WATER VAPOUR DENSITY, so WVD is scientifically correct.
# If your manuscript deliberately uses another symbol, change only this label.
SYMBOL = "WVD"

BANDS = {
    "all_E_band": (71.0, 86.0, "All E-band"),
    "80GHz": (81.0, 86.0, "81–86 GHz"),
    "70GHz": (71.0, 76.0, "71–76 GHz"),
}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root", type=Path,
        default=Path("/home/jovyan/ESSD_publication_dataset"),
        help="Dataset root containing WVD_V6_QUANTIZATION_BINNED_VALIDATION",
    )
    ap.add_argument(
        "--out-dir", type=Path, default=None,
        help="Output directory (default: <root>/WVP_SUPPORT_HEXBINS)",
    )
    ap.add_argument("--plot-max", type=float, default=PLOT_MAX_GM3)
    ap.add_argument("--gridsize", type=int, default=GRIDSIZE)
    ap.add_argument("--dpi", type=int, default=DPI)
    return ap.parse_args()


def metrics(ref, est):
    ref = np.asarray(ref, dtype=float)
    est = np.asarray(est, dtype=float)
    m = np.isfinite(ref) & np.isfinite(est)
    ref = ref[m]
    est = est[m]
    err = est - ref
    if len(ref) < 2:
        return dict(N=len(ref), PCC=np.nan, RMSE=np.nan, Bias=np.nan)
    return dict(
        N=int(len(ref)),
        PCC=float(np.corrcoef(ref, est)[0, 1]),
        RMSE=float(np.sqrt(np.mean(err * err))),
        Bias=float(np.mean(err)),
    )


def nice_major_step(plot_max):
    # For 35 g m-3 this gives 5 g m-3 major ticks and 2.5 minor ticks.
    if plot_max <= 20:
        return 2.5
    if plot_max <= 40:
        return 5.0
    if plot_max <= 80:
        return 10.0
    return 20.0


def select_population(links, link_idx, ref, est, f_lo, f_hi):
    # Explicit E-band selection. For the all-E panel, exclude the 76-81 GHz gap.
    if f_lo == 71.0 and f_hi == 86.0:
        freq_ok = (
            links["frequency_ghz"].between(71.0, 76.0, inclusive="both")
            | links["frequency_ghz"].between(81.0, 86.0, inclusive="both")
        )
    else:
        freq_ok = links["frequency_ghz"].between(f_lo, f_hi, inclusive="both")

    link_mask = (
        freq_ok
        & (links["length_km"] >= MIN_LENGTH_KM)
        & np.isfinite(links["wvd_quantization_resolution_ref20c_rho10_gm3"])
        & (
            links["wvd_quantization_resolution_ref20c_rho10_gm3"]
            <= MAX_DELTA_RHO_Q_GM3
        )
    ).to_numpy(bool)

    pair_mask = link_mask[link_idx]
    x = np.asarray(ref[pair_mask], dtype=float)   # observed/reference WVD
    y = np.asarray(est[pair_mask], dtype=float)   # CML WVD
    use_links = links.loc[link_mask].copy()
    met = metrics(x, y)
    met["n_sublinks"] = int(len(use_links))
    met["n_links"] = int(use_links["cml_id"].nunique())
    return x, y, use_links, met


def draw_panel(ax, x, y, met, title, plot_max, gridsize, add_xlabel=True, add_ylabel=True):
    hb = ax.hexbin(
        x, y,
        gridsize=gridsize,
        extent=(0.0, plot_max, 0.0, plot_max),
        mincnt=1,
        cmap=CMAP,
        norm=LogNorm(vmin=1),
        linewidths=0.0,
    )

    ax.plot(
        [0, plot_max], [0, plot_max],
        linestyle="--", color="0.25", linewidth=1.5, zorder=5,
    )

    # Article-style detached axes: a little whitespace around the bounded spines.
    pad = 0.02 * plot_max
    ax.set_xlim(-pad, plot_max + pad * 0.20)
    ax.set_ylim(-pad, plot_max + pad * 0.20)
    ax.set_aspect("equal", adjustable="box")

    if add_xlabel:
        ax.set_xlabel(
            rf"{SYMBOL}$_{{\mathrm{{obs}}}}$ [g m$^{{-3}}$]",
            labelpad=10,
        )
    if add_ylabel:
        ax.set_ylabel(
            rf"{SYMBOL}$_{{\mathrm{{CML}}}}$ [g m$^{{-3}}$]",
            labelpad=10,
        )

    major = nice_major_step(plot_max)
    minor = major / 2.0
    ax.xaxis.set_major_locator(MultipleLocator(major))
    ax.yaxis.set_major_locator(MultipleLocator(major))
    ax.xaxis.set_minor_locator(MultipleLocator(minor))
    ax.yaxis.set_minor_locator(MultipleLocator(minor))
    ax.tick_params(which="major", length=9, width=1.5, direction="out", pad=10)
    ax.tick_params(which="minor", length=5, width=1.0, direction="out")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_position(("outward", 8))
    ax.spines["bottom"].set_position(("outward", 8))
    ax.spines["left"].set_bounds(0, plot_max)
    ax.spines["bottom"].set_bounds(0, plot_max)

    txt = (
        f"N={met['N']:,}\n"
        f"Sublinks={met['n_sublinks']:,}\n"
        f"Links={met['n_links']:,}\n"
        f"RMSE={met['RMSE']:.2f} g m$^{{-3}}$\n"
        f"Bias={met['Bias']:+.2f} g m$^{{-3}}$\n"
        f"PCC={met['PCC']:.3f}"
    )
    ax.text(
        0.04, 0.96, txt,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=13,
        family="monospace",
        bbox=dict(
            boxstyle="round,pad=0.22",
            facecolor="white", edgecolor="0.65",
            alpha=0.88, linewidth=1.0,
        ),
        zorder=10,
    )
    return hb


def main():
    args = parse_args()
    root = args.root.expanduser().resolve()
    src = root / "WVD_V6_QUANTIZATION_BINNED_VALIDATION"
    links_csv = src / "validation_population_links.csv"
    pairs_npz = src / "validation_pairs_graf.npz"
    out = (args.out_dir or (root / "WVP_SUPPORT_HEXBINS")).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    links = pd.read_csv(links_csv)
    with np.load(pairs_npz) as z:
        ref = z["ref"].astype(np.float64, copy=False)
        est = z["est"].astype(np.float64, copy=False)
        link_idx = z["link_idx"].astype(np.int64, copy=False)

    selected = {}
    metric_rows = []
    for key, (f_lo, f_hi, title) in BANDS.items():
        x, y, use_links, met = select_population(
            links, link_idx, ref, est, f_lo, f_hi
        )
        selected[key] = (x, y, use_links, met, title)
        metric_rows.append({
            "population": title,
            "min_length_km": MIN_LENGTH_KM,
            "max_delta_rho_q_gm3": MAX_DELTA_RHO_Q_GM3,
            "N_observations": met["N"],
            "n_sublinks": met["n_sublinks"],
            "n_links": met["n_links"],
            "PCC": met["PCC"],
            "RMSE_gm3": met["RMSE"],
            "Bias_gm3": met["Bias"],
        })

    pd.DataFrame(metric_rows).to_csv(out / "wvp_support_metrics.csv", index=False)

    rc = {
        "font.size": 18,
        "axes.labelsize": 22,
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
        "axes.linewidth": 1.5,
        "mathtext.fontset": "dejavusans",
    }

    # Separate figures.
    with plt.rc_context(rc):
        for key, (x, y, use_links, met, title) in selected.items():
            fig, ax = plt.subplots(figsize=(9.2, 9.2))
            hb = draw_panel(
                ax, x, y, met, title,
                args.plot_max, args.gridsize,
                add_xlabel=True, add_ylabel=True,
            )
            cbar = fig.colorbar(hb, ax=ax, pad=0.035, fraction=0.047)
            cbar.set_label("Count per hexagon", rotation=90, labelpad=16)
            cbar.ax.tick_params(which="major", length=8, width=1.4)
            fig.tight_layout()
            fig.savefig(out / f"wvp_support_{key}.png", dpi=args.dpi, bbox_inches="tight")
            fig.savefig(out / f"wvp_support_{key}.pdf", bbox_inches="tight")
            plt.close(fig)

    # Three-panel article figure: all E-band, 81-86 GHz, 71-76 GHz.
    with plt.rc_context(rc):
        fig, axes = plt.subplots(1, 3, figsize=(23.5, 7.7), sharex=True, sharey=True)
        hbs = []
        for i, key in enumerate(["all_E_band", "80GHz", "70GHz"]):
            x, y, use_links, met, title = selected[key]
            hb = draw_panel(
                axes[i], x, y, met, title,
                args.plot_max, args.gridsize,
                add_xlabel=True,
                add_ylabel=(i == 0),
            )
            hbs.append(hb)

        # Use the largest count as a shared logarithmic color scale.
        vmax = max(float(np.nanmax(hb.get_array())) for hb in hbs)
        vmax = 10.0 ** math.ceil(math.log10(max(vmax, 1.0)))
        shared_norm = LogNorm(vmin=1.0, vmax=vmax)
        for hb in hbs:
            hb.set_norm(shared_norm)

        cbar = fig.colorbar(hbs[-1], ax=axes.tolist(), pad=0.018, fraction=0.025)
        cbar.set_label("Count per hexagon", rotation=90, labelpad=16)
        cbar.ax.tick_params(which="major", length=8, width=1.4)

        fig.subplots_adjust(left=0.07, right=0.91, bottom=0.13, top=0.90, wspace=0.19)
        fig.savefig(out / "wvp_support_three_panel.png", dpi=args.dpi, bbox_inches="tight")
        fig.savefig(out / "wvp_support_three_panel.pdf", bbox_inches="tight")
        plt.close(fig)

    print("Selection:")
    print(f"  length >= {MIN_LENGTH_KM:.1f} km")
    print(f"  Delta rho_q <= {MAX_DELTA_RHO_Q_GM3:.2f} g m-3")
    print("\nMetrics:")
    print(pd.DataFrame(metric_rows).to_string(index=False))
    print("\nOutput:", out)
    print("No pair-level WVD/WVP observations were exported.")


if __name__ == "__main__":
    main()
