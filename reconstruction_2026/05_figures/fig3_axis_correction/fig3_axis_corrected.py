"""fig3_axis_corrected.py - Stage 6: Figure 3 with the approved axis correction. Copy of Stage 5B
fig3_cml_panels_abc.py; only nice_ticks() and the y-limit lines of panels (a) and (b) are changed (data, bins,
colours, fonts and layout unchanged).

Original docstring: fig3_cml_panels_abc.py - Stage 5B, Figure 3 (f03_cml_characteristics.png).

Producer: notebook cell 22 of
  C:\\Users\\taouf\\Downloads\\DISDRODB\\extracted\\matched_code\\matched_code\\cml dataset ()runing.ipynb
(identical code in "CML CZECHIA DATASET.ipynb" cell 107 and "cml_dataset_czechia.ipynb" cell 2), read as JSON
and copied verbatim into ..\\extracted_code\\fig3_cell22_cml_dataset_runing.py. The cell wrote cml_panels_abc.png
(dpi 600, bbox tight) + cml_panels_abc.pdf from `ds_c` (cell 21: xr.open_dataset("CML_S2020E2021_1MIN_CZ.nc")).

This script keeps the cell's plotting code byte-for-byte (block between the "BEGIN VERBATIM" / "END VERBATIM"
markers) and changes only:
  1. cell 21 data loading -> `ds_c` is built from the population selected with --population
     (variables used by the cell: frequency(sublink_id, cml_id) in MHz, length(cml_id) in km);
  2. for --population qc: frequency of sublinks with qc_sublink_status != 0 is set to NaN, so the
     sublink-level panels (a) and (c) show retained sublinks only; panel (b) shows the 1987 retained CMLs;
  3. `! pip install matplotlib` (notebook shell magic) dropped; matplotlib backend Agg; plt.show() dropped;
  4. savefig targets -> --out-png / --out-pdf (same dpi, bbox_inches, pad_inches).
It also prints the plotted counts, ticks, hexbin maximum, legend values and short-path statistics.

Populations:
  qc            approved Stage 3 QC REV1 (1987 CMLs, 3731 retained sublinks)           [new figure]
  qc_allsub     QC REV1, all sublinks of the 1987 CMLs (3974 with frequency)           [information]
  raw2574       approved RAW population (full-history RAW minus the 6 R1 CMLs)         [information]
  raw_union_qc  released RAW 2580 + released-QC-only CMLs (2848) - closest available
                proxy of the unreleased original input (verification v2_01)          [reproduction check]
  raw2580       released RAW 2580                                                      [reproduction check]

Usage:
  python -I -B fig3_cml_panels_abc.py --population qc --out-png ..\\FINAL_FIGURES\\f03_cml_characteristics.png
"""
import argparse
import os
import sys

import numpy as np
import netCDF4
import xarray as xr
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
R = r"C:\Users\taouf\ESSD_RELEASE"
p = argparse.ArgumentParser()
p.add_argument("--population", required=True, choices=["qc", "qc_allsub", "raw2574", "raw_union_qc", "raw2580"])
p.add_argument("--qc", default=R + r"\FINAL_RECONSTRUCTION_20261027\stage3\QC_CORRECTED_CANDIDATE_REV1.nc")
p.add_argument("--raw-full", default=R + r"\extracted\full\ESSD_CZ_2021_HANDOFF\data\source_full_history\raw\CZ_CML_RAW_1MIN_20200101_20211231.nc")
p.add_argument("--raw-approved-csv", default=R + r"\FINAL_RECONSTRUCTION_20261027\stage2\RAW_APPROVED_POPULATION.csv")
p.add_argument("--raw-released", default=R + r"\zenodo_staging\primary\CZ_CML_RAW_1MIN_20210401_20211031.nc")
p.add_argument("--qc-released", default=R + r"\zenodo_staging\primary\CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc")
p.add_argument("--out-png", required=True)
p.add_argument("--out-pdf", default=None)
a = p.parse_args()


def meta(path):
    d = netCDF4.Dataset(path, "r")
    ids = np.array([str(x) for x in d.variables["cml_id"][:]])
    fv = d.variables["frequency"]
    assert fv.dimensions == ("sublink_id", "cml_id"), fv.dimensions
    f = np.ma.filled(fv[:].astype(float), np.nan)
    L = np.ma.filled(d.variables["length"][:].astype(float), np.nan)
    st = None
    if "qc_sublink_status" in d.variables:
        st = np.ma.filled(d.variables["qc_sublink_status"][:].astype(int), -1)
    d.close()
    return ids, f, L, st


# ---------------- population (replaces notebook cell 21) ----------------
if a.population in ("qc", "qc_allsub"):
    ids, f, L, st = meta(a.qc)
    if a.population == "qc":
        f = np.where(st == 0, f, np.nan)          # adaptation: retained sublinks only
elif a.population == "raw2574":
    ids, f, L, _ = meta(a.raw_full)
    import csv
    keep_ids = sorted({r["cml_id"] for r in csv.DictReader(open(a.raw_approved_csv, encoding="utf-8"))})
    k = np.isin(ids, keep_ids)
    assert k.sum() == len(keep_ids) == 2574, (k.sum(), len(keep_ids))
    ids, f, L = ids[k], f[:, k], L[k]
elif a.population == "raw2580":
    ids, f, L, _ = meta(a.raw_released)
else:  # raw_union_qc, as in data_verification_20261027/VERIFICATION_SCRIPTS/v2_01_fig3_population.py
    rid, rf, rL, _ = meta(a.raw_released)
    qid, qf, qL, _ = meta(a.qc_released)
    qo = ~np.isin(qid, rid)
    ids = np.concatenate([rid, qid[qo]])
    f = np.concatenate([rf, qf[:, qo]], axis=1)
    L = np.concatenate([rL, qL[qo]])

ds_c = xr.Dataset(
    {"frequency": (("sublink_id", "cml_id"), f.astype(np.float32)),
     "length": (("cml_id",), L.astype(np.float32))},
    coords={"cml_id": ids, "sublink_id": [0, 1]},
)
print(f"[population] {a.population}: n_cml={ids.size}, sublinks with finite frequency={int(np.isfinite(f).sum())}")

# ================= BEGIN VERBATIM (notebook cell 22, minus '! pip install', savefig paths, plt.show) =================
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import warnings

warnings.filterwarnings("ignore")

# ============================================================
# STYLE — crisp, publication-quality
# ============================================================
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "lines.antialiased": True,
    "patch.antialiased": True,
    "text.antialiased": True,
})


# ============================================================
# EXTRACT VARIABLES
# ============================================================
freq_all = ds_c["frequency"].values.flatten() / 1000.0   # MHz → GHz
freq_all = freq_all[~np.isnan(freq_all)]

length = ds_c["length"].values
length = length[~np.isnan(length)]

freq_c = ds_c["frequency"].values / 1000.0
len_c  = ds_c["length"].values

n_sub = freq_c.shape[0]
freq_flat = freq_c.flatten()
len_flat  = np.tile(len_c, n_sub)

mask = ~np.isnan(freq_flat) & ~np.isnan(len_flat)
freq_flat = freq_flat[mask]
len_flat  = len_flat[mask]

# ============================================================
# DATA-DRIVEN RANGES
# ============================================================
freq_min = np.floor(freq_all.min()) - 3
freq_max = np.ceil(freq_all.max()) + 1

len_max = np.ceil(length.max()) + 1

bins_a = np.arange(freq_min, freq_max + 2, 2)
bins_b = np.arange(0, len_max + 1, 1)

# ============================================================
# HELPER — 4 nice ticks over a data range
# ============================================================
def nice_ticks(lo, hi, n=4):
    # Stage 6 axis correction (approved by the author): regular ticks with one step from {1, 2, 2.5, 5} x 10^k,
    # from floor(lo) to ceil(hi), at most n + 1 ticks, so the axis covers the whole data range.
    raw_step = (hi - lo) / (n - 1)
    mag = 10 ** np.floor(np.log10(raw_step))
    for step in np.array([1, 2, 2.5, 5, 10, 20]) * mag:
        start = np.floor(lo / step) * step
        stop = np.ceil(hi / step) * step
        if round((stop - start) / step) + 1 <= n + 1:
            break
    return np.round(np.arange(start, stop + step * 0.5, step), 6)

# ============================================================
# HELPER — disconnected axes (R-style)
# ============================================================
def disconnect_axes(ax, xticks, yticks):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["bottom"].set_bounds(xticks[0], xticks[-1])
    ax.spines["left"].set_bounds(yticks[0], yticks[-1])

    ax.spines["bottom"].set_position(("outward", 3))
    ax.spines["left"].set_position(("outward", 3))

    ax.set_xticks(xticks)
    ax.set_yticks(yticks)

    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5,
        width=0.6,
        labelsize=9,
        top=False,
        right=False,
        left=True,
        bottom=True,
    )

# ============================================================
# FIGURE
# ============================================================
TOMATO = "#E8573A"

fig, axes = plt.subplots(
    1, 3,
    figsize=(14.5, 4.2),
    gridspec_kw={"width_ratios": [1, 1, 1.25]},
)

# ----- (a) Frequency histogram (sub-links) -----------------
ax = axes[0]
ax.hist(freq_all, bins=bins_a, color=TOMATO, edgecolor="white", linewidth=0.4)
ax.set_xlabel("Frequency (GHz)", fontsize=10)
ax.set_ylabel("Number of sub-links", fontsize=10)

xt_a = nice_ticks(freq_min, freq_max, 4)
yt_a = nice_ticks(0, ax.get_ylim()[1], 4)
ax.set_xlim(xt_a[0] - 1, xt_a[-1] + 1)
ax.set_ylim(-yt_a[-1] * 0.02, yt_a[-1])  # slight negative pad; axis ends at the last tick (Stage 6)
disconnect_axes(ax, xt_a, yt_a)
ax.set_title("(a)", fontsize=11, fontweight="bold", loc="left", pad=8)

# ----- (b) Path-length histogram (links) -------------------
ax = axes[1]
ax.hist(length, bins=bins_b, color=TOMATO, edgecolor="white", linewidth=0.4)
ax.set_xlabel("Path length (km)", fontsize=10)
ax.set_ylabel("Number of links", fontsize=10)

xt_b = nice_ticks(0, len_max, 4)
yt_b = nice_ticks(0, ax.get_ylim()[1], 4)
# Start x slightly below 0 so the first bar is fully visible
ax.set_xlim(-0.8, xt_b[-1] + 0.5)
ax.set_ylim(-yt_b[-1] * 0.02, yt_b[-1])  # axis ends at the last tick (Stage 6)
disconnect_axes(ax, xt_b, yt_b)
ax.set_title("(b)", fontsize=11, fontweight="bold", loc="left", pad=8)

# ----- (c) Hexbin: frequency vs path length ----------------
ax = axes[2]
hb = ax.hexbin(
    freq_flat, len_flat,
    gridsize=25,
    cmap="YlGnBu",
    mincnt=1,
    linewidths=0.3,
    edgecolors="0.3",
)
ax.set_xlabel("Frequency (GHz)", fontsize=10)
ax.set_ylabel("Path length (km)", fontsize=10)

xt_c = nice_ticks(freq_min, freq_max, 4)
yt_c = nice_ticks(0, len_max, 4)
ax.set_xlim(xt_c[0] - 1, xt_c[-1] + 1)
# Start y slightly below 0 so hexbins near zero are fully visible
ax.set_ylim(-0.8, yt_c[-1] + 0.5)
disconnect_axes(ax, xt_c, yt_c)
ax.set_title("(c)", fontsize=11, fontweight="bold", loc="left", pad=8)

# --- Discrete legend (matching the reference) ---
counts = hb.get_array()
vmin_c, vmax_c = int(counts.min()), int(counts.max())

n_legend = min(15, vmax_c)
legend_vals = np.linspace(1, vmax_c, n_legend, dtype=int)
legend_vals = np.unique(legend_vals)[::-1]

cmap_c = mpl.colormaps["YlGnBu"]
norm = mpl.colors.Normalize(vmin=1, vmax=vmax_c)

legend_handles = [
    Line2D(
        [0], [0],
        marker="h",
        color="none",
        markerfacecolor=cmap_c(norm(v)),
        markeredgecolor="0.3",
        markeredgewidth=0.3,
        markersize=8,
        label=str(v),
    )
    for v in legend_vals
]

leg = ax.legend(
    handles=legend_handles,
    title="Counts",
    title_fontsize=7.5,
    fontsize=7,
    loc="upper right",
    framealpha=0.9,
    edgecolor="0.6",
    borderpad=0.4,
    labelspacing=0.25,
    handletextpad=0.3,
)
leg.get_frame().set_linewidth(0.5)

ax.text(
    1.0, 1.02,
    "Number of sub-links",
    transform=ax.transAxes,
    fontsize=9,
    ha="right",
    va="bottom",
    fontfamily="serif",
)

plt.tight_layout(w_pad=2.0)
plt.savefig(a.out_png,                                           # changed: was "cml_panels_abc.png"
            dpi=600, bbox_inches="tight", pad_inches=0.05)
plt.savefig(a.out_pdf or os.path.splitext(a.out_png)[0] + ".pdf",  # changed: was "cml_panels_abc.pdf"
            bbox_inches="tight", pad_inches=0.05)
# ================= END VERBATIM =================

# ---------------- report ----------------
ca = np.histogram(freq_all, bins=bins_a)[0]
cb = np.histogram(length, bins=bins_b)[0]
from PIL import Image
print(f"[saved] {a.out_png}  size(px)={Image.open(a.out_png).size}")
print(f"panel (a) n sublinks = {freq_all.size}; fmin {freq_all.min():.3f} fmax {freq_all.max():.3f} GHz; bins from {freq_min} step 2")
print(f"panel (b) n links    = {length.size}; Lmin {length.min():.3f} Lmax {length.max():.3f} km")
print(f"panel (c) n sublinks = {freq_flat.size}; hexbin sum {int(counts.sum())}; max {vmax_c}")
print(f"ticks a.x {xt_a.tolist()} a.y {yt_a.tolist()} | b.x {xt_b.tolist()} b.y {yt_b.tolist()} | c.x {xt_c.tolist()} c.y {yt_c.tolist()}")
print(f"legend values {legend_vals.tolist()}")
print("(a) non-zero bins: " + ", ".join(f"[{lo:g},{lo+2:g})={int(c)}" for lo, c in zip(bins_a[:-1], ca) if c > 0))
print("(b) bins: " + ", ".join(f"[{lo:g},{lo+1:g})={int(c)}" for lo, c in zip(bins_b[:-1], cb)))
Lk = length
print(f"short paths per link (n={Lk.size}): <1 km {100*np.mean(Lk<1):.2f}% ({int((Lk<1).sum())}); <0.1 km {100*np.mean(Lk<0.1):.2f}% ({int((Lk<0.1).sum())}); min {Lk.min():.3f} km")
print(f"short paths per plotted sublink (n={len_flat.size}): <1 km {100*np.mean(len_flat<1):.2f}% ({int((len_flat<1).sum())}); <0.1 km {100*np.mean(len_flat<0.1):.2f}% ({int((len_flat<0.1).sum())})")
eb = ((freq_all >= 71) & (freq_all <= 76)) | ((freq_all >= 81) & (freq_all <= 86))
fl = ds_c["frequency"].values / 1000.0
eb_link = (((fl >= 71) & (fl <= 76)) | ((fl >= 81) & (fl <= 86))).any(0)
print(f"E-band: sublinks {int(eb.sum())} of {freq_all.size} ({100*eb.mean():.1f}%); CMLs with an E-band sublink {int(eb_link.sum())} of {ids.size} ({100*eb_link.mean():.1f}%)")
