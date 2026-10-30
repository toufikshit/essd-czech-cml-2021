#!/usr/bin/env python3
"""Stage 5B / Agent 3b - Fig. 7 CML-row pixel check (verification only, no figure).

The radar row of Fig. 7 needs native RADKLIM YW daily files and cartopy,
neither of which exists on this computer.  This script therefore verifies the
CML (upper) row only: it computes the rolling 1-h CML accumulations with the
renderer's UNCHANGED functions (AST-extracted from
original/code/figures/03_blettner_event_maps_timeseries.py: cml_rolling_hour,
rolling_window_times, nearest_time_index, dry_mask, discrete_essd_cmap,
essd_colormap, parse_bounds; DEFAULT_HOURLY_BOUNDS), classifies every
grid value into the renderer's discrete colour classes, and compares the class
with the colour class of each pixel inside the 8 upper panels of a rendered PNG.

Panel geometry: the axes frames (grey spines) are detected in the PNG; the
inside of each frame is the PlateCarree map extent taken from the event cache
(event_rank 4, instance 1).  Pixels whose colour is far from every class
colour and from the dry background (borders, frame, antialiasing) are ignored.

Usage:
  python -I f7_02_pixel_check.py --png <published or evidence png> --grid <grid.nc> [--grid ...]
"""
from __future__ import annotations

import argparse
import ast
import json
import math
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
import xarray as xr
from PIL import Image
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap

REPO = Path(r'C:\Users\taouf\ESSD_RELEASE\essd-czech-cml-2021')
RENDERER = REPO / 'original' / 'code' / 'figures' / '03_blettner_event_maps_timeseries.py'
CACHE = Path(r'C:\Users\taouf\ESSD_RELEASE\extracted\full\ESSD_CZ_2021_HANDOFF\data\event_support'
             r'\historical_event_cache\selected_instances_cache.csv')
NAMES = {'essd_colormap', 'discrete_essd_cmap', 'parse_bounds', 'dry_mask', 'nearest_time_index',
         'rolling_window_times', 'cml_rolling_hour'}


def renderer_functions(path=RENDERER):
    tree = ast.parse(path.read_text(encoding='utf-8'))
    body = [n for n in tree.body if (isinstance(n, ast.FunctionDef) and n.name in NAMES) or
            (isinstance(n, ast.Assign) and any(getattr(t, 'id', '') == 'DEFAULT_HOURLY_BOUNDS' for t in n.targets))]
    ns = {'np': np, 'pd': pd, 'math': math, 'LinearSegmentedColormap': LinearSegmentedColormap,
          'ListedColormap': ListedColormap, 'BoundaryNorm': BoundaryNorm, 'Sequence': list}
    from typing import Sequence  # noqa
    ns['Sequence'] = Sequence
    exec(compile(ast.Module(body=body, type_ignores=[]), str(path), 'exec'), ns)
    missing = NAMES - set(ns)
    if missing:
        raise RuntimeError(f'renderer functions missing: {missing}')
    return ns


def hourly_fields(grid_nc, R, panel_times):
    out = []
    with xr.open_dataset(grid_nc) as ds:
        da = ds['R']
        times = pd.DatetimeIndex(da.time.values)
        for t in panel_times:
            acc, _, _ = R['cml_rolling_hour'](da, times, pd.Timestamp(t), 1.0)
            if acc is None:
                raise RuntimeError(f'no rolling hour at {t}')
            out.append(acc.astype('float32'))
        x = np.asarray(ds['x'].values); y = np.asarray(ds['y'].values)
    return x, y, out


def frames(img):
    g = (np.abs(img[..., 0] - img[..., 1]) < 6) & (np.abs(img[..., 1] - img[..., 2]) < 6) & \
        (img[..., 0] > 140) & (img[..., 0] < 200)
    H, W = g.shape
    def groups(idx):
        out = []; cur = [idx[0]]
        for i in idx[1:]:
            if i == cur[-1] + 1: cur.append(i)
            else: out.append(cur); cur = [i]
        out.append(cur)
        return [(c[0] + c[-1]) / 2 for c in out]
    rows = groups(np.where(g.sum(1) > W * 0.3)[0])
    cols = groups(np.where(g.sum(0) > H * 0.15)[0])
    return rows, cols


def classify_png(rgb, palette, dry_cols, tol=18.0):
    cand = np.vstack([palette, dry_cols])
    d = np.sqrt(((rgb[..., None, :] - cand[None, None, :, :]) ** 2).sum(-1))
    k = d.argmin(-1); dmin = d.min(-1)
    cls = np.where(k < len(palette), k, -1)  # -1 = dry
    cls = np.where(dmin <= tol, cls, -2)     # -2 = ignore
    return cls


def classify_field(v, bounds):
    c = np.digitize(v, bounds) - 1           # 0..len(bounds)-2 ; >=last bound -> len-1
    c = np.where(c >= len(bounds) - 1, len(bounds) - 2, c)   # over colour == last colour
    c = np.where(~np.isfinite(v) | (v < bounds[0]), -1, c)
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--png', type=Path, required=True)
    ap.add_argument('--grid', type=Path, nargs='+', required=True)
    ap.add_argument('--cache', type=Path, default=CACHE)
    ap.add_argument('--first-panel', default='2021-07-13 16:05:00')
    ap.add_argument('--json', type=Path, default=None)
    a = ap.parse_args()
    R = renderer_functions()
    bounds = R['parse_bounds'](R['DEFAULT_HOURLY_BOUNDS'])
    cmap = R['discrete_essd_cmap'](bounds)
    palette = np.asarray(cmap.colors)[:, :3] * 255.0
    dry_cols = np.array([[255, 255, 255], [250, 250, 250]], float)  # bad=white, LAND facecolor 0.98
    row = pd.read_csv(a.cache).query('event_rank==4 and instance_rank==1').iloc[0]
    ext = [row.extent_lon_min, row.extent_lon_max, row.extent_lat_min, row.extent_lat_max]
    panel_times = [pd.Timestamp(a.first_panel) + pd.Timedelta(hours=i) for i in range(8)]
    img = np.asarray(Image.open(a.png).convert('RGB')).astype(float)
    rows, cols = frames(img)
    print('frame rows', rows); print('frame cols', cols)
    top, bot = rows[0], rows[1]
    lefts = cols[0::2]; rights = cols[1::2]
    if len(lefts) != 8 or len(rights) != 8:
        raise RuntimeError(f'expected 8 panels, found cols {cols}')
    results = {'png': str(a.png), 'extent': ext, 'grids': {}}
    for gpath in a.grid:
        x, y, fields = hourly_fields(gpath, R, panel_times)
        res = float(np.median(np.diff(x)))
        tot_n = tot_ok = 0; per = []
        for p, (L, Rr) in enumerate(zip(lefts, rights)):
            i0, i1 = int(math.ceil(L)) + 4, int(math.floor(Rr)) - 4
            j0, j1 = int(math.ceil(top)) + 4, int(math.floor(bot)) - 4
            sub = img[j0:j1, i0:i1]
            pc = classify_png(sub, palette, dry_cols)
            jj, ii = np.mgrid[j0:j1, i0:i1]
            lon = ext[0] + (ii + 0.5 - L) / (Rr - L) * (ext[1] - ext[0])
            lat = ext[3] - (jj + 0.5 - top) / (bot - top) * (ext[3] - ext[2])
            ix = np.clip(np.rint((lon - x[0]) / res).astype(int), 0, len(x) - 1)
            iy = np.clip(np.rint((lat - y[0]) / res).astype(int), 0, len(y) - 1)
            fc = classify_field(fields[p][iy, ix], bounds)
            use = pc != -2
            ok = (pc == fc) & use
            n, k = int(use.sum()), int(ok.sum())
            within1 = int(((np.abs(pc - fc) <= 1) & (pc >= 0) & (fc >= 0) | ok).sum())
            per.append({'panel': str(panel_times[p]), 'pixels_compared': n, 'same_class': k,
                        'agreement': k / n if n else None, 'within_one_class_or_same': within1 / n if n else None})
            tot_n += n; tot_ok += k
        results['grids'][str(gpath)] = {'resolution': res, 'overall_agreement': tot_ok / tot_n, 'panels': per}
        print(f'{gpath.name}: res {res:.3f}  overall same-class agreement {tot_ok/tot_n:.4%}')
        for d in per:
            print('   ', d)
    if a.json:
        a.json.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
