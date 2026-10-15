"""Regression gates: stage inputs/outputs, ID-order neutrality and unchanged sources."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr
import pytest
from mini_fixture import create
from essd_cz.connected import Project, grid_final_rain, manual_event_selection, _original_idw_functions


def make(tmp_path,scenario='baseline'):
    d=tmp_path/scenario
    create(d,scenario=scenario)
    pr=Project(d/'data_paths.json',tmp_path/'results')
    return d,pr


def test_manifest_independent_of_hpc(tmp_path):
    d,p=make(tmp_path)
    assert 'raw' in p.paths and all(x.is_file() for x in p.paths.values())
    assert p.source_rainfall().is_file()
    again=Project(d/'data_paths.json',tmp_path/'results')
    assert again.get('rainfall')==p.paths['rain']
    assert not str(p.paths['rain']).startswith('/home/jovyan')


def test_idw_original_math_and_missing():
    geometry,field=_original_idw_functions()
    X,Y,idx,dist,within=geometry(np.array([14.,15.]),np.array([50.,50.]),np.array([14.,15.]),np.array([50.]),2,.9)
    z,n=field(np.array([2.,np.nan]),idx,dist,within,2)
    assert np.isclose(z[0],2) and np.isnan(z[1])
    assert np.array_equal(n,[1,0])


def test_grid_produces_native_schema(tmp_path):
    d,p=make(tmp_path)
    p.source_rainfall()
    dest=tmp_path/'grid.nc'
    grid_final_rain(p.paths['rain'],p.paths['qc'],dest,'2021-04-01 08:15','2021-04-01 09:15',resolution=.04,max_grid_cells=10000)
    assert dest.is_file()
    with xr.open_dataset(dest) as ds:
        assert {'R','nvalid'}<=set(ds.data_vars)
        assert ds.R.dims==('time','y','x')
        assert ds.R.attrs['units']=='mm h-1'
        assert ds.sizes['time']==13
        assert np.isfinite(ds.R.values).any()
    with pytest.raises(FileExistsError):
        grid_final_rain(p.paths['rain'],p.paths['qc'],dest,'2021-04-01 08:15','2021-04-01 09:15')
    p.register('spatial_grid',dest,'synthetic test')
    assert Project(d/'data_paths.json',tmp_path/'results').get('spatial_grid')==dest


def test_grid_same_when_rain_cml_ids_reordered(tmp_path):
    d,p=make(tmp_path)
    baseline=tmp_path/'normal.nc'
    grid_final_rain(p.paths['rain'],p.paths['qc'],baseline,'2021-04-01 08:15','2021-04-01 09:15',resolution=.05)
    with xr.open_dataset(p.paths['rain']) as rain:
        reordered=rain.isel(cml_id=[2,0,1]).load()
    try:import netCDF4;engine='netcdf4'  # match mini_fixture: uint8 flags need NETCDF4
    except ImportError:engine='scipy'
    swap=tmp_path/'rain_swapped.nc'; reordered.to_netcdf(swap,engine=engine)
    out=tmp_path/'swap_grid.nc'
    grid_final_rain(swap,p.paths['qc'],out,'2021-04-01 08:15','2021-04-01 09:15',resolution=.05)
    with xr.open_dataset(baseline) as a,xr.open_dataset(out) as b:
        np.testing.assert_allclose(a.R.values,b.R.values,equal_nan=True,rtol=0,atol=0)


def test_manual_cache_is_explicit(tmp_path):
    p=manual_event_selection(tmp_path/'event.csv',event_time='2021-07-13 16:05',extent=[12.0,19.0,48.5,51.5])
    tab=pd.read_csv(p)
    assert tab.iloc[0]['event_rank']==4
    assert 'DEMONSTRATION' in tab.iloc[0]['selection_origin']
    with pytest.raises(FileExistsError):manual_event_selection(p,event_time='2021-07-13 16:05',extent=[12,19,48,51])
