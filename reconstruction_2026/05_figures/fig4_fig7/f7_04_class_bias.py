#!/usr/bin/env python3
"""Stage 5B / Agent 3b - per-panel class bias between a rendered Fig. 7 PNG (upper row) and a CML grid.
Verification only. Usage: python -I f7_04_class_bias.py <png> <grid.nc> [<diag_png_prefix>]"""
import sys; sys.dont_write_bytecode = True
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import f7_02_pixel_check as P
import numpy as np, pandas as pd
from PIL import Image
png, grid = Path(sys.argv[1]), Path(sys.argv[2]); diag = sys.argv[3] if len(sys.argv) > 3 else None
R = P.renderer_functions(); bounds = R['parse_bounds'](R['DEFAULT_HOURLY_BOUNDS'])
pal = np.asarray(R['discrete_essd_cmap'](bounds).colors)[:, :3] * 255
img = np.asarray(Image.open(png).convert('RGB')).astype(float)
rows, cols = P.frames(img); top, bot = rows[0], rows[1]; lefts, rights = cols[0::2], cols[1::2]
row = pd.read_csv(P.CACHE).query('event_rank==4 and instance_rank==1').iloc[0]
ext = [row.extent_lon_min, row.extent_lon_max, row.extent_lat_min, row.extent_lat_max]
times = [pd.Timestamp('2021-07-13 16:05') + pd.Timedelta(hours=i) for i in range(8)]
x, y, F = P.hourly_fields(grid, R, times); res = float(np.median(np.diff(x)))
panels = []
for p in range(8):
    L, Rr = lefts[p], rights[p]; i0, i1 = int(np.ceil(L)) + 4, int(np.floor(Rr)) - 4; j0, j1 = int(np.ceil(top)) + 4, int(np.floor(bot)) - 4
    pc = P.classify_png(img[j0:j1, i0:i1], pal, np.array([[255, 255, 255], [250, 250, 250]], float))
    jj, ii = np.mgrid[j0:j1, i0:i1]
    lon = ext[0] + (ii + .5 - L) / (Rr - L) * (ext[1] - ext[0]); lat = ext[3] - (jj + .5 - top) / (bot - top) * (ext[3] - ext[2])
    ix = np.clip(np.rint((lon - x[0]) / res).astype(int), 0, len(x) - 1); iy = np.clip(np.rint((lat - y[0]) / res).astype(int), 0, len(y) - 1)
    fc = P.classify_field(F[p][iy, ix], bounds); u = pc != -2
    both = (pc >= 0) & (fc >= 0) & u
    print(times[p].strftime('%H:%M'), 'png wet %.3f grid wet %.3f | png-only wet %.3f grid-only wet %.3f | mean(png-grid class | both wet) %+.3f | png>grid %.3f png<grid %.3f' % (
        ((pc >= 0) & u).sum() / u.sum(), ((fc >= 0) & u).sum() / u.sum(), ((pc >= 0) & (fc == -1) & u).sum() / u.sum(),
        ((fc >= 0) & (pc == -1) & u).sum() / u.sum(), (pc[both] - fc[both]).mean() if both.any() else np.nan,
        (pc[both] > fc[both]).mean() if both.any() else np.nan, (pc[both] < fc[both]).mean() if both.any() else np.nan))
    if diag and p in (2, 3, 4):
        out = np.full(pc.shape + (3,), 255, np.uint8)
        for k in range(len(pal)): out[fc == k] = pal[k]
        panels.append(np.concatenate([img[j0:j1, i0:i1].astype(np.uint8), np.full((6, i1 - i0, 3), 0, np.uint8), out], 0))
if diag:
    Image.fromarray(np.concatenate([np.pad(q, ((0, 0), (0, 8), (0, 0)), constant_values=0) for q in panels], 1)).save(diag)
