r"""fig2_pixel_compare.py - Stage 5B, Figure 2: pixel comparison and the published crop.

1. Compares the reproduced PNG (unmodified producer on the original full-history inputs) with
   (a) the author's Matplotlib output C:\Users\taouf\Downloads\availability_native_timestep_percentage.png and
   (b) the published figure (original_figures\availability_native_timestep_percentage_v1.png), which equals
       the Matplotlib output with its top --crop-top rows (the "2021" year label) removed.
2. --apply-crop IN OUT: writes OUT = IN without its top --crop-top rows (the post-processing step that produced
   the published *_v1.png), after checking that the removed rows contain no plot content below the year label.

Usage: python -I -B fig2_pixel_compare.py --repro PNG [--mpl PNG] [--published PNG] [--crop-top 125]
       python -I -B fig2_pixel_compare.py --apply-crop IN OUT [--crop-top 125]
"""
import argparse, hashlib
import numpy as np
from PIL import Image
p = argparse.ArgumentParser()
p.add_argument("--repro")
p.add_argument("--mpl", default=r"C:\Users\taouf\Downloads\availability_native_timestep_percentage.png")
p.add_argument("--published", default=r"C:\Users\taouf\ESSD_RELEASE\FINAL_RECONSTRUCTION_20261027\stage5b\original_figures\availability_native_timestep_percentage_v1.png")
p.add_argument("--crop-top", type=int, default=125)
p.add_argument("--apply-crop", nargs=2)
a = p.parse_args()
def load(f): return np.asarray(Image.open(f).convert("RGBA")).astype(int)
def sha(f): return hashlib.sha256(open(f, "rb").read()).hexdigest()
if a.apply_crop:
    src, dst = a.apply_crop
    im = Image.open(src); arr = load(src)
    top = arr[:a.crop_top, :, :3]
    nonwhite_rows = np.where((top < 250).any(axis=(1, 2)))[0]
    print(f"[crop] {src} {im.size}; non-white rows inside the removed band: {nonwhite_rows.min() if nonwhite_rows.size else None}-{nonwhite_rows.max() if nonwhite_rows.size else None}")
    cols = np.where((top < 250).any(axis=(0, 2)))[0]
    print(f"[crop] non-white columns in removed band: {cols.min() if cols.size else None}-{cols.max() if cols.size else None} (year label only expected)")
    out = im.crop((0, a.crop_top, im.size[0], im.size[1]))
    out.save(dst, dpi=(600, 600))
    print(f"[crop] saved {dst} {out.size} sha256 {sha(dst)}")
    raise SystemExit
R = load(a.repro); M = load(a.mpl); P = load(a.published)
print(f"[F] repro {a.repro} {R.shape[1]}x{R.shape[0]} sha256 {sha(a.repro)}")
print(f"[F] author mpl output {a.mpl} {M.shape[1]}x{M.shape[0]} sha256 {sha(a.mpl)} info {Image.open(a.mpl).info}")
print(f"[F] published {a.published} {P.shape[1]}x{P.shape[0]} sha256 {sha(a.published)} info {Image.open(a.published).info}")
if R.shape == M.shape:
    d = np.abs(R - M).max(2)
    print(f"[F] repro vs author mpl output: identical pixels {100*(d==0).mean():.4f}%; max abs diff {d.max()}; differing px {(d>0).sum()}")
    if (d > 0).any():
        ys, xs = np.where(d > 0); print(f"     differing bbox x {xs.min()}-{xs.max()} y {ys.min()}-{ys.max()}")
else:
    print("[F] repro and author mpl output differ in size")
c = R[a.crop_top:]
if c.shape == P.shape:
    d = np.abs(c - P).max(2)
    print(f"[F] repro[top {a.crop_top} rows removed] vs published: identical pixels {100*(d==0).mean():.4f}%; max abs diff {d.max()}; differing px {(d>0).sum()}")
    if (d > 0).any():
        ys, xs = np.where(d > 0); print(f"     differing bbox x {xs.min()}-{xs.max()} y {ys.min()}-{ys.max()}")
else:
    print(f"[F] cropped repro {c.shape} vs published {P.shape}: size differs")

# ---- invert the RdYlGn colours of each data row (column-wise, row centre) to availability percent ----
import matplotlib
cm = matplotlib.colormaps["RdYlGn"]
grid = np.linspace(0, 100, 10001)
lut = (np.array([cm(v / 100.0)[:3] for v in grid]) * 255.0)
def invert(rgb):
    d = ((rgb[:, None, :] - lut[None, :, :]) ** 2).sum(2)
    return grid[d.argmin(1)], np.sqrt(d.min(1))
# row bands: find rows of the 5 data rows from white separator lines in the left part of the plot
col = 3000
Rc = R[a.crop_top:]
names = ["CML", "Radar", "Rain gauge", "Air temperature", "Relative humidity"]
whites = np.where((P[:, col, :3] > 250).all(1))[0]
segs = np.split(np.arange(P.shape[0]), np.where(np.diff((P[:, col, :3] > 250).all(1).astype(int)) != 0)[0] + 1)
bands = [s for s in segs if s.size > 100 and not (P[s[0], col, :3] > 250).all()][:5]
xs = np.arange(740, 6740)
for nm, b in zip(names, bands):
    yc = int(b.mean())
    vp, ep = invert(P[yc, xs, :3].astype(float)); vr, er = invert(Rc[yc, xs, :3].astype(float))
    dd = vr - vp
    print(f"[F] row {nm:18s} y={yc}: published vs repro inverted availability: mean diff {dd.mean():+.3f} pp, "
          f"|diff| p50 {np.median(np.abs(dd)):.3f}, p99 {np.percentile(np.abs(dd), 99):.3f}, max {np.abs(dd).max():.3f} pp; "
          f"published mean {vp.mean():.2f}%, repro mean {vr.mean():.2f}% (colour-fit residual max {max(ep.max(), er.max()):.1f})")
