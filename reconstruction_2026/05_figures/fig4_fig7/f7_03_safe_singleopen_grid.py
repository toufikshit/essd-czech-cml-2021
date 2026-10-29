#!/usr/bin/env python3
"""Stage 5B / Agent 3b - Fig. 7 hypothesis grid with the 'safe single-open' IDW.

Candidate producer of the published Fig. 7 CML grid
(cmlonly_graf_nearby_comparison/apr_oct_2021/grids/CZ_CML_GRID_GRAF_5MIN_APR_OCT_2021_IDW.nc,
the GRAF grid that original/code/figures/06_run_multiple_events_vs_radar.py
passes to 03_blettner_event_maps_timeseries.py):

  Downloads/DISDRODB/extracted/matched_code/matched_code/cml_grid_essd_publication_safe_singleopen.py

The three functions build_grid, load_midpoints and idw_fast_series are taken
UNCHANGED from that file by AST extraction (the file is not imported as a
module; argparse/main and the HPC writer are not used).  This script only
reproduces the data flow of its run() for a short time window:
  ids from the rain file -> midpoints from the coordinate file (same ids) ->
  finite-midpoint links -> lattice build_grid(lon, lat, resolution) ->
  idw_fast_series(...) per time step (nnear 15, p 2, max_dist 0.3, res 0.01).
The output schema R(time,y,x) is what 03_blettner_event_maps_timeseries.py reads.

--lattice-from <grid.nc> forces the x/y node vectors of an existing grid (used
to put the Stage 4 field on the IDENTICAL lattice as the original).
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import sys
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
import xarray as xr

SRC = Path(r'C:\Users\taouf\Downloads\DISDRODB\extracted\matched_code\matched_code'
           r'\cml_grid_essd_publication_safe_singleopen.py')


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def extract(src: Path):
    tree = ast.parse(src.read_text(encoding='utf-8'))
    want = {'build_grid', 'load_midpoints', 'idw_fast_series'}
    body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in want]
    if {n.name for n in body} != want:
        raise RuntimeError('functions not found')
    import gc, typing
    ns = {'np': np, 'xr': xr, 'pd': pd, 'gc': gc, 'Path': Path,
          'Tuple': typing.Tuple, 'Optional': typing.Optional, 'annotations': None}
    code = compile(ast.Module(body=body, type_ignores=[]), str(src), 'exec')
    exec(code, ns)
    return ns


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--src', type=Path, default=SRC)
    ap.add_argument('--input-nc', type=Path, required=True)
    ap.add_argument('--coords-nc', type=Path, required=True)
    ap.add_argument('--output-nc', type=Path, required=True)
    ap.add_argument('--rain-var', default='R_mean_5min')
    ap.add_argument('--tstart', default='2021-07-13 15:05:00')
    ap.add_argument('--tend', default='2021-07-13 23:05:00')
    ap.add_argument('--resolution', type=float, default=0.01)
    ap.add_argument('--nnear', type=int, default=15)
    ap.add_argument('--p', type=int, default=2)
    ap.add_argument('--max-dist', type=float, default=0.3)
    ap.add_argument('--lattice-from', type=Path, default=None)
    a = ap.parse_args()
    if a.output_nc.exists():
        raise FileExistsError(a.output_nc)
    print('source', a.src, 'sha256', sha256(a.src))
    F = extract(a.src)
    with xr.open_dataset(a.input_nc) as ds:
        cml_ids = np.asarray(ds['cml_id'].values).astype(str)
        rain = ds[a.rain_var]
        times_all = pd.DatetimeIndex(ds['time_5min'].values)
        i0 = int(np.searchsorted(times_all.values, np.datetime64(pd.Timestamp(a.tstart)), side='left'))
        i1 = int(np.searchsorted(times_all.values, np.datetime64(pd.Timestamp(a.tend)), side='right'))
        times = times_all[i0:i1]
        args = SimpleNamespace(site_a_lon='site_a_longitude', site_a_lat='site_a_latitude',
                               site_b_lon='site_b_longitude', site_b_lat='site_b_latitude')
        lon_all, lat_all = F['load_midpoints'](a.coords_nc, cml_ids, args)
        good = np.isfinite(lon_all) & np.isfinite(lat_all)
        good_idx = np.flatnonzero(good)
        lon, lat = lon_all[good], lat_all[good]
        if a.lattice_from is None:
            xg2d, yg2d = F['build_grid'](lon, lat, a.resolution)
        else:
            with xr.open_dataset(a.lattice_from) as L:
                xg2d, yg2d = np.meshgrid(np.asarray(L.x.values), np.asarray(L.y.values))
        Z = rain.isel(cml_id=good_idx, time_5min=slice(i0, i1)).transpose('time_5min', 'cml_id').values
        Z = np.ascontiguousarray(Z, dtype='float32')
    print(f'Links for IDW : {good_idx.size}/{len(cml_ids)}; grid {xg2d.shape}; times {times[0]} -> {times[-1]} ({len(times)})')
    print('lattice x0 %.10f y0 %.10f' % (xg2d[0, 0], yg2d[0, 0]))
    G = F['idw_fast_series'](lon, lat, Z, xg2d, yg2d, nnear=a.nnear, p=a.p, max_distance=a.max_dist,
                             batch=32, cell_chunk=20000, workers=1)
    out = xr.Dataset({'R': (('time', 'y', 'x'), G.astype('float32'))},
                     coords={'time': times.values, 'x': xg2d[0, :], 'y': yg2d[:, 0]})
    out['R'].attrs['units'] = 'mm h-1'
    out.attrs.update(source_rainrate_file=str(a.input_nc), source_coordinate_file=str(a.coords_nc),
                     interpolation=f'IDW (nnear={a.nnear}, p={a.p}, max_dist={a.max_dist})',
                     grid_resolution_deg=a.resolution, links_for_idw=int(good_idx.size),
                     producer_functions=f'{a.src.name}: build_grid, load_midpoints, idw_fast_series (unchanged, AST)')
    a.output_nc.parent.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(a.output_nc)
    print('WROTE', a.output_nc, 'sha256', sha256(a.output_nc))


if __name__ == '__main__':
    main()
