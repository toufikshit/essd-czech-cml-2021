#!/usr/bin/env python3
"""Figure 7: rainfall fields on 13 July 2021, 16:05 to 23:05 UTC: CML field (upper row) and RADKLIM YW (lower row).

Step 1, CML grid: the 5-minute CML rainfall (R_mean_5min of scripts/derive_5min.py, 1987 CMLs) at the link
midpoints is interpolated by inverse distance weighting (nnear 15, power 2, maximum distance 0.3 degree) to the
0.01 degree lattice of config/fig7/grid_lattice.json, 15:05 to 23:05 UTC. The functions build_grid, load_midpoints
and idw_fast_series are taken unchanged from original/code/figures/cml_grid_essd_publication_safe_singleopen.py.
Step 2, figure: the unchanged renderer original/code/figures/03_blettner_event_maps_timeseries.py with the
arguments of original/code/figures/06_run_multiple_events_vs_radar.py for event 4 (config/fig7/event_selection.csv);
the radar row is read from the daily RADKLIM YW 2017.002 file of 13 July 2021 in extras.radklim_dir.
Outputs in <output_dir>/figures/fig07: grid_cml_20210713.nc and hourly_1h/event04_instance01_multitime_hourly1h.png.
Usage: python scripts/figures/fig07_event_maps.py [--config config/user_paths.json]
Needs cartopy for the radar row.
"""
from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from essd_cz.config import ORIGINAL, extra, load_config  # noqa: E402
from essd_cz.rain5 import derived_5min_name  # noqa: E402

GRID_SRC = ORIGINAL / 'figures' / 'cml_grid_essd_publication_safe_singleopen.py'
RENDERER = ORIGINAL / 'figures' / '03_blettner_event_maps_timeseries.py'
FIG7 = ROOT / 'config' / 'fig7'


def extract(src: Path):
    """Take build_grid, load_midpoints and idw_fast_series unchanged from the grid builder (AST, no import)."""
    tree = ast.parse(src.read_text(encoding='utf-8'))
    want = {'build_grid', 'load_midpoints', 'idw_fast_series'}
    body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in want]
    if {n.name for n in body} != want:
        raise RuntimeError('functions not found')
    import gc, typing
    ns = {'np': np, 'xr': xr, 'pd': pd, 'gc': gc, 'Path': Path,
          'Tuple': typing.Tuple, 'Optional': typing.Optional, 'annotations': None}
    exec(compile(ast.Module(body=body, type_ignores=[]), str(src), 'exec'), ns)
    return ns


def build_cml_grid(rain5, coords, output, tstart, tend, nnear=15, p=2, max_dist=0.3, resolution=0.01):
    from types import SimpleNamespace
    F = extract(GRID_SRC)
    with xr.open_dataset(rain5) as ds:
        cml_ids = np.asarray(ds['cml_id'].values).astype(str)
        times_all = pd.DatetimeIndex(ds['time_5min'].values)
        i0 = int(np.searchsorted(times_all.values, np.datetime64(pd.Timestamp(tstart)), side='left'))
        i1 = int(np.searchsorted(times_all.values, np.datetime64(pd.Timestamp(tend)), side='right'))
        times = times_all[i0:i1]
        args = SimpleNamespace(site_a_lon='site_a_longitude', site_a_lat='site_a_latitude',
                               site_b_lon='site_b_longitude', site_b_lat='site_b_latitude')
        lon_all, lat_all = F['load_midpoints'](coords, cml_ids, args)
        good_idx = np.flatnonzero(np.isfinite(lon_all) & np.isfinite(lat_all))
        lon, lat = lon_all[good_idx], lat_all[good_idx]
        lattice = json.loads((FIG7 / 'grid_lattice.json').read_text())
        xg2d, yg2d = np.meshgrid(np.asarray(lattice['x']), np.asarray(lattice['y']))
        Z = ds['R_mean_5min'].isel(cml_id=good_idx, time_5min=slice(i0, i1)).transpose('time_5min', 'cml_id').values
        Z = np.ascontiguousarray(Z, dtype='float32')
    print(f'Links for IDW: {good_idx.size}/{len(cml_ids)}; grid {xg2d.shape}; times {times[0]} -> {times[-1]}', flush=True)
    G = F['idw_fast_series'](lon, lat, Z, xg2d, yg2d, nnear=nnear, p=p, max_distance=max_dist,
                             batch=32, cell_chunk=20000, workers=1)
    out = xr.Dataset({'R': (('time', 'y', 'x'), G.astype('float32'))},
                     coords={'time': times.values, 'x': xg2d[0, :], 'y': yg2d[:, 0]})
    out['R'].attrs['units'] = 'mm h-1'
    out.attrs.update(source_rainrate_file=Path(rain5).name, source_coordinate_file=Path(coords).name,
                     interpolation=f'IDW (nnear={nnear}, p={p}, max_dist={max_dist})', grid_resolution_deg=resolution,
                     links_for_idw=int(good_idx.size),
                     producer_functions=f'{GRID_SRC.name}: build_grid, load_midpoints, idw_fast_series (unchanged)')
    out.to_netcdf(output)
    print('WROTE', output, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
    ap.add_argument('--rain5', type=Path, default=None, help='5-minute rainfall (default: output of derive_5min.py)')
    ap.add_argument('--out-dir', type=Path, default=None)
    ap.add_argument('--grid-only', action='store_true', help='build the CML grid only (no RADKLIM needed)')
    a = ap.parse_args()
    cfg = load_config(a.config)
    ds = cfg['datasets']
    rain5 = a.rain5 or cfg['output_dir'] / 'rainfall_5min' / derived_5min_name(ds['rain'].name)
    if not Path(rain5).is_file():
        raise SystemExit(f'{rain5} not found (run scripts/derive_5min.py first)')
    out = (a.out_dir or cfg['output_dir'] / 'figures' / 'fig07').resolve()
    out.mkdir(parents=True, exist_ok=True)
    grid = out / 'grid_cml_20210713.nc'
    if not grid.exists():
        build_cml_grid(rain5, ds['qc'], grid, '2021-07-13 15:05:00', '2021-07-13 23:05:00')
    if a.grid_only:
        return
    cmd = [sys.executable, str(RENDERER), '--gridded-nc', str(grid), '--cml-var', 'R', '--cml-row-label', 'GRAF',
           '--cml-geom-nc', str(ds['qc']), '--radar-root', str(extra(a.config, 'radklim_dir')),
           '--selection-cache-csv', str(FIG7 / 'event_selection.csv'),
           '--event-rank', '4', '--instance-rank', '1', '--start-time', '2021-07-13 16:05:00',
           '--ncols', '8', '--panel-step-minutes', '60', '--dpi', '600', '--skip-timeseries', '--out-dir', str(out)]
    print('RENDERER:', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
