#!/usr/bin/env python3
"""Stage 5B / Agent 3b - Fig. 7 numerical old-vs-new comparison of the CML (upper-row) fields
and of the CML population entering the IDW; plus a DIAGNOSTIC old|new image of the upper row.

Inputs (CLI, local defaults):
  --old-grid / --new-grid : R(time,y,x) 5-min grids on the identical lattice
  --old-rain / --new-rain : 5-min rainfall files that produced them (R_mean_5min)
  --old-coords / --new-coords : geometry files (site_a/b lon/lat)
  --published-png : the published figure (its upper row is cropped into the diagnostic)
Hourly fields are computed with the renderer's unchanged cml_rolling_hour (AST, see f7_02).
Statistics are over grid nodes inside the map extent of event 4 / instance 1.

The diagnostic image is NOT the publication renderer (no cartopy, no radar row);
it uses the renderer's colour classes, BoundaryNorm and the 0.1 mm mask only.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from PIL import Image

import f7_02_pixel_check as P

REC = Path(r'C:\Users\taouf\ESSD_RELEASE\FINAL_RECONSTRUCTION_20261027')
FH = Path(r'C:\Users\taouf\ESSD_RELEASE\extracted\full\ESSD_CZ_2021_HANDOFF\data\source_full_history')
HERE = Path(__file__).resolve().parent.parent


def ids_of(ds):
    return np.asarray(ds['cml_id'].values).astype(str)


def mids(path):
    with xr.open_dataset(path, decode_times=False) as d:
        lon = (np.asarray(d.site_a_longitude.values, float) + np.asarray(d.site_b_longitude.values, float)) / 2
        lat = (np.asarray(d.site_a_latitude.values, float) + np.asarray(d.site_b_latitude.values, float)) / 2
        return dict(zip(ids_of(d), zip(lon, lat)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--old-grid', type=Path, default=HERE / 'work' / 'grid_safe_original_0p01.nc')
    ap.add_argument('--new-grid', type=Path, default=HERE / 'work' / 'grid_safe_stage4_0p01.nc')
    ap.add_argument('--old-rain', type=Path, default=FH / 'rain' / 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_20200101_20211231.nc')
    ap.add_argument('--new-rain', type=Path, default=REC / 'stage4' / 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc')
    ap.add_argument('--old-coords', type=Path, default=FH / 'qc' / 'CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc')
    ap.add_argument('--new-coords', type=Path, default=REC / 'stage3' / 'QC_CORRECTED_CANDIDATE_REV1.nc')
    ap.add_argument('--published-png', type=Path, default=REC / 'stage5b' / 'original_figures' / 'event04_instance01_multitime_hourly1h.png')
    ap.add_argument('--out-dir', type=Path, default=HERE / 'reproduction')
    ap.add_argument('--comparison-png', type=Path, default=HERE / 'comparison' / 'event04_instance01_multitime_hourly1h_old_vs_new.png')
    ap.add_argument('--buffer-deg', type=float, default=0.3, help='IDW max_dist: links that can influence the extent')
    a = ap.parse_args()
    R = P.renderer_functions()
    bounds = R['parse_bounds'](R['DEFAULT_HOURLY_BOUNDS'])
    cmap = R['discrete_essd_cmap'](bounds); norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=False)
    row = pd.read_csv(P.CACHE).query('event_rank==4 and instance_rank==1').iloc[0]
    ext = [row.extent_lon_min, row.extent_lon_max, row.extent_lat_min, row.extent_lat_max]
    times = [pd.Timestamp('2021-07-13 16:05') + pd.Timedelta(hours=i) for i in range(8)]
    xo, yo, Fo = P.hourly_fields(a.old_grid, R, times)
    xn, yn, Fn = P.hourly_fields(a.new_grid, R, times)
    if not (np.array_equal(xo, xn) and np.array_equal(yo, yn)):
        raise RuntimeError('lattices differ')
    inx = (xo >= ext[0]) & (xo <= ext[1]); iny = (yo >= ext[2]) & (yo <= ext[3])
    panel_rows = []
    for t, fo, fn in zip(times, Fo, Fn):
        o = fo[np.ix_(iny, inx)]; n = fn[np.ix_(iny, inx)]
        d = n - o
        def s(f):
            w = np.isfinite(f) & (f >= 0.1)
            return dict(max=float(np.nanmax(f)), mean_all=float(np.nanmean(f)),
                        wet_frac=float(w.sum() / np.isfinite(f).sum()), mean_wet=float(f[w].mean()) if w.any() else np.nan,
                        n_ge_10mm=int((f >= 10).sum()), n_ge_20mm=int((f >= 20).sum()))
            # noqa
        so, sn = s(o), s(n)
        co = P.classify_field(o, bounds); cn = P.classify_field(n, bounds)
        panel_rows.append({'panel_end_utc': str(t), **{f'old_{k}': v for k, v in so.items()},
                           **{f'new_{k}': v for k, v in sn.items()},
                           'mean_diff_mm': float(np.nanmean(d)), 'mean_abs_diff_mm': float(np.nanmean(np.abs(d))),
                           'max_increase_mm': float(np.nanmax(d)), 'max_decrease_mm': float(np.nanmin(d)),
                           'frac_nodes_class_changed': float((co != cn).mean()),
                           'n_nodes': int(o.size)})
    dfp = pd.DataFrame(panel_rows)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    dfp.to_csv(a.out_dir / 'fig7_panel_stats_old_vs_new.csv', index=False)
    print(dfp.to_string())

    # ---- population entering the IDW (links that can influence the extent) ----
    t0, t1 = pd.Timestamp('2021-07-13 15:10'), pd.Timestamp('2021-07-13 23:05')  # 12 steps before 16:05 .. 23:05
    def pop(rain_p, coords_p):
        m = mids(coords_p)
        with xr.open_dataset(rain_p) as r:
            ids = ids_of(r)
            tt = pd.DatetimeIndex(r.time_5min.values); sel = np.flatnonzero((tt >= t0) & (tt <= t1))
            Z = np.asarray(r['R_mean_5min'].isel(time_5min=sel).transpose('time_5min', 'cml_id').values, float)
            extra = {}
            if 'blackout_minutes_5min' in r:
                extra['blackout'] = np.asarray(r['blackout_minutes_5min'].isel(time_5min=sel).values, float).sum(axis=(0, 1))
                extra['interp'] = np.asarray(r['interpolated_minutes_5min'].isel(time_5min=sel).values, float).sum(axis=(0, 1))
        out = {}
        for j, c in enumerate(ids):
            if c not in m:
                continue
            lo, la = m[c]
            near = (ext[0] - a.buffer_deg <= lo <= ext[1] + a.buffer_deg) and (ext[2] - a.buffer_deg <= la <= ext[3] + a.buffer_deg)
            inside = (ext[0] <= lo <= ext[1]) and (ext[2] <= la <= ext[3])
            z = Z[:, j]
            out[c] = dict(lon=lo, lat=la, near=near, inside=inside, n_finite=int(np.isfinite(z).sum()),
                          acc_mm=float(np.nansum(z) / 12.0) if np.isfinite(z).any() else np.nan,
                          peak=float(np.nanmax(z)) if np.isfinite(z).any() else np.nan,
                          blackout_min=float(extra['blackout'][j]) if 'blackout' in extra else 0.0,
                          interp_min=float(extra['interp'][j]) if 'interp' in extra else 0.0)
        return out, len(sel)
    po, nsteps = pop(a.old_rain, a.old_coords)
    pn, _ = pop(a.new_rain, a.new_coords)
    allc = sorted(set(po) | set(pn))
    rows = []
    for c in allc:
        o = po.get(c); n = pn.get(c)
        base = o or n
        rows.append(dict(cml_id=c, lon=base['lon'], lat=base['lat'], inside_extent=base['inside'], within_extent_plus_0p3deg=base['near'],
                         in_old=o is not None, in_new=n is not None,
                         old_n_finite_5min=o['n_finite'] if o else 0, new_n_finite_5min=n['n_finite'] if n else 0,
                         old_event_acc_mm=o['acc_mm'] if o else np.nan, new_event_acc_mm=n['acc_mm'] if n else np.nan,
                         old_peak_mmph=o['peak'] if o else np.nan, new_peak_mmph=n['peak'] if n else np.nan,
                         new_blackout_minutes=n['blackout_min'] if n else 0.0, new_interpolated_minutes=n['interp_min'] if n else 0.0))
    dfc = pd.DataFrame(rows)
    dfc['acc_diff_mm'] = dfc.new_event_acc_mm - dfc.old_event_acc_mm
    dfc.to_csv(a.out_dir / 'fig7_cml_population_old_vs_new.csv', index=False)
    near = dfc[dfc.within_extent_plus_0p3deg]
    summ = {
        'window_5min_steps': int(nsteps), 'window': [str(t0), str(t1)],
        'old_total_cmls_in_idw': int(dfc.in_old.sum()), 'new_total_cmls_in_idw': int(dfc.in_new.sum()),
        'near_extent': {
            'old': int(near.in_old.sum()), 'new': int(near.in_new.sum()),
            'removed': int((near.in_old & ~near.in_new).sum()), 'added': int((~near.in_old & near.in_new).sum()),
            'old_with_any_finite': int((near.old_n_finite_5min > 0).sum()), 'new_with_any_finite': int((near.new_n_finite_5min > 0).sum()),
            'old_all_finite': int((near.old_n_finite_5min == nsteps).sum()), 'new_all_finite': int((near.new_n_finite_5min == nsteps).sum()),
            'common_with_acc_change_gt_0p1mm': int((near.in_old & near.in_new & (near.acc_diff_mm.abs() > 0.1)).sum()),
            'new_cmls_with_blackout_minutes': int((near.new_blackout_minutes > 0).sum()),
            'new_blackout_minutes_total': float(near.new_blackout_minutes.sum()),
            'new_interpolated_minutes_total': float(near.new_interpolated_minutes.sum()),
        },
        'inside_extent': {
            'old': int(dfc[dfc.inside_extent].in_old.sum()), 'new': int(dfc[dfc.inside_extent].in_new.sum()),
            'removed': int((dfc.inside_extent & dfc.in_old & ~dfc.in_new).sum()),
            'added': int((dfc.inside_extent & ~dfc.in_old & dfc.in_new).sum()),
            'old_with_any_finite': int((dfc.inside_extent & (dfc.old_n_finite_5min > 0)).sum()),
            'new_with_any_finite': int((dfc.inside_extent & (dfc.new_n_finite_5min > 0)).sum()),
        },
    }
    big = near[near.in_old & near.in_new].assign(absd=lambda d: d.acc_diff_mm.abs()).sort_values('absd', ascending=False).head(10)
    summ['top10_common_acc_changes'] = big[['cml_id', 'old_event_acc_mm', 'new_event_acc_mm', 'new_blackout_minutes']].to_dict('records')
    rem = near[near.in_old & ~near.in_new].sort_values('old_event_acc_mm', ascending=False).head(10)
    summ['top10_removed_by_old_acc'] = rem[['cml_id', 'old_event_acc_mm', 'old_n_finite_5min']].to_dict('records')
    add = near[~near.in_old & near.in_new].sort_values('new_event_acc_mm', ascending=False).head(10)
    summ['top10_added_by_new_acc'] = add[['cml_id', 'new_event_acc_mm', 'new_n_finite_5min', 'new_blackout_minutes']].to_dict('records')
    (a.out_dir / 'fig7_population_summary.json').write_text(json.dumps(summ, indent=2, default=float) + '\n')
    print(json.dumps(summ, indent=2, default=float))

    # ---- diagnostic image: published upper row | old (reproduced) | new ----
    img = np.asarray(Image.open(a.published_png).convert('RGB'))
    rws, cls = P.frames(img.astype(float))
    crop = img[int(rws[0]) - 40:int(rws[1]) + 6, int(cls[0]) - 6:int(cls[-1]) + 6]
    fig = plt.figure(figsize=(13.6, 6.6))
    gs = fig.add_gridspec(4, 8, height_ratios=[1.15, 1, 1, 0.16], hspace=0.28, wspace=0.03,
                          left=0.08, right=0.995, top=0.95, bottom=0.07)
    axp = fig.add_subplot(gs[0, :]); axp.imshow(crop); axp.axis('off')
    axp.set_title('Published Fig. 7, upper row (crop of the original PNG)', fontsize=10)
    for r_i, (F, lab) in enumerate([(Fo, 'OLD\nreproduced\n(final Graf\n2219 CMLs)'), (Fn, 'NEW\nStage 4\n(1987 CMLs)')]):
        for c_i in range(8):
            ax = fig.add_subplot(gs[1 + r_i, c_i])
            m = R['dry_mask'](F[c_i], 0.1)
            mp = ax.pcolormesh(xo, yo, m, shading='nearest', cmap=cmap, norm=norm, rasterized=True)
            ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect('equal')
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_edgecolor('0.68'); sp.set_linewidth(0.6)
            if r_i == 0:
                ax.set_title(times[c_i].strftime('%H:%M'), fontsize=10, fontweight='bold')
            if c_i == 0:
                ax.text(-0.08, 0.5, lab, transform=ax.transAxes, ha='right', va='center', fontsize=8.5)
    cax = fig.add_subplot(gs[3, 1:7])
    cb = fig.colorbar(mp, cax=cax, orientation='horizontal', boundaries=bounds, ticks=bounds, spacing='uniform',
                      extend='max', extendfrac=0.035, drawedges=True)
    cb.set_ticklabels([f'{b:g}' for b in bounds]); cb.ax.tick_params(labelsize=7)
    cb.set_label('Hourly rain sum [mm]  - DIAGNOSTIC (no cartopy, no radar row; not the publication renderer)', fontsize=9)
    a.comparison_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.comparison_png, dpi=200, facecolor='white')
    print('WROTE', a.comparison_png)


if __name__ == '__main__':
    main()
