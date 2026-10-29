! pip install matplotlib
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
    raw_step = (hi - lo) / (n - 1)
    mag = 10 ** np.floor(np.log10(raw_step))
    candidates = np.array([1, 2, 2.5, 5, 10]) * mag
    step = candidates[np.argmin(np.abs(candidates - raw_step))]
    start = np.ceil(lo / step) * step
    ticks = np.arange(start, hi + step * 0.1, step)
    if len(ticks) > n:
        ticks = ticks[:n]
    elif len(ticks) < n:
        ticks = np.linspace(lo, hi, n)
        ticks = np.round(ticks / mag) * mag
    return ticks

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
ax.set_ylim(-ax.get_ylim()[1] * 0.02, ax.get_ylim()[1])  # slight negative pad
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
ax.set_ylim(-ax.get_ylim()[1] * 0.02, ax.get_ylim()[1])
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
plt.savefig("cml_panels_abc.png",
            dpi=600, bbox_inches="tight", pad_inches=0.05)
plt.savefig("cml_panels_abc.pdf",
            bbox_inches="tight", pad_inches=0.05)
plt.show()
print("Saved → cml_panels_abc.png + .pdf")