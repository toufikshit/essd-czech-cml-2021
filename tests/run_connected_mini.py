#!/usr/bin/env python3
"""Run a tiny, self-contained seven-NetCDF walkthrough; no pytest needed.

Synthetic CML/gauge/weather data, not a test of Graf physical accuracy.
Exercises original availability program and the NEW original-IDW-based 5-min adapter.
"""
from __future__ import annotations
import argparse
import shutil
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from mini_fixture import create
from essd_cz.connected import Project,grid_final_rain,manual_event_selection
from essd_cz.workflows import availability
from essd_cz.notebook_figures import plot_network,plot_characteristics,monthly_from_native_csv
import numpy as np
import xarray as xr
from essd_cz.config import EXPECTED_NAMES


def run(root):
    root=Path(root).expanduser().resolve()
    if root.exists():raise FileExistsError(f'Test target exists; will not delete or overwrite: {root}')
    source=root/'artificial_data'
    create(source)
    for role,nn in [('gauge',3),('humidity',2),('temperature',2)]:
        path=source / EXPECTED_NAMES[role]
        tmp=path.with_suffix('.staging.nc')
        with xr.open_dataset(path) as ds:
            dim='WSI' if role=='gauge' else 'station'
            ds.assign_coords(lat=(dim,np.linspace(49.8,50.2,nn)),lon=(dim,np.linspace(14.1,15.,nn))).to_netcdf(tmp,engine='scipy')
        path.unlink();tmp.replace(path)
    p=Project(source/'data_paths.json',root/'results')
    print('[00] network/metadata and native availability')
    f,_=plot_network(p.paths,p.stage_dir('00')/'network',dpi=60)
    f,_=plot_characteristics(p.paths,p.stage_dir('00')/'characteristics',dpi=60)
    availability(p.paths,p.stage_dir('00')/'native',dpi=60)
    pic,csv,monthly=monthly_from_native_csv(p.stage_dir('00')/'native',p.stage_dir('00')/'monthly',dpi=60)
    p.register('exploration_complete',f,'2-day artificial fixture')
    print('[01] register artificial rain (not original Graf physics)')
    p.source_rainfall()
    print('[02] stage a rain time-series chart from artificial fixture')
    p.register('rain_timeline',pic,'artificial figure handoff')
    print('[03] generate adapted 5-min grid with same original IDW functions')
    grid=p.stage_dir('03')/'mini_idw.nc'
    grid_final_rain(p.get('rainfall'),p.paths['qc'],grid,'2021-04-01 08:00','2021-04-01 09:00',resolution=.04)
    p.register('spatial_grid',grid,'synthetic IDW demo')
    with xr.open_dataset(grid) as ds:
        assert ds.sizes['time']==13 and 'R' in ds and np.isfinite(ds.R.values).any()
    print('[04] register EXPLICIT manual event choice (not historical event discovery)')
    cache=manual_event_selection(p.stage_dir('04')/'example_selection.csv',event_time='2021-04-01 08:00',extent=[14,15.3,49.7,50.4])
    p.register('manual_event_cache',cache,'explicit synthetic demo')
    assert all(x.is_file() for x in (p.get('exploration_complete'),p.get('rainfall'),p.get('spatial_grid'),p.get('manual_event_cache')))
    print('PASS: seven artificial NetCDFs + exploration + native availability + rain handoff + original IDW + map handoff')
    print('Manifest:',p.manifest)

if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--out-dir',type=Path,default=Path('results/mini_connected_demo'))
    a=cli.parse_args()
    run(a.out_dir)
