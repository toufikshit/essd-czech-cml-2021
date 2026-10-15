"""Reviewers: adapter guards and CSV aggregation using 2-day simulated data."""
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import pytest
from mini_fixture import create
from essd_cz.notebook_figures import frequency_length_pairs,monthly_from_native_csv,plot_characteristics,station_coordinates,plot_network
from essd_cz.station_workflow import converter_command,run_converter,CONVERTER

def test_frequency_metadata_stays_in_ghz_and_km(tmp_path):
    paths=create(tmp_path/'sim')
    with xr.open_dataset(paths['qc']) as ds:
        f,L=frequency_length_pairs(ds)
    assert len(f)==6
    assert np.isclose(f[0],25.088)
    assert np.isclose(f.max(),73.0)
    assert np.isclose(L.max(),7.5)
    output,counts=plot_characteristics(paths,tmp_path/'plots',dpi=60)
    assert output.is_file() and counts['directional_pairs']==6

def test_missing_station_geo_is_not_silently_faked(tmp_path):
    paths=create(tmp_path/'sim')
    with xr.open_dataset(paths['gauge']) as ds:
        with pytest.raises(ValueError,match='No compatible static station coordinates'):
            station_coordinates(ds)

def test_map_uses_actual_static_station_coordinates(tmp_path):
    paths=create(tmp_path/'sim')
    for role,n in [('gauge',3),('humidity',2),('temperature',2)]:
        src=Path(paths[role]);dest=src.with_suffix('.tmp.nc')
        with xr.open_dataset(src) as ds:
            v=ds.assign_coords(lat=(('WSI' if role=='gauge' else 'station',),np.linspace(50,50.2,n)),lon=(('WSI' if role=='gauge' else 'station',),np.linspace(14.3,14.8,n)))
            v.to_netcdf(dest,engine='scipy')
        src.unlink();dest.rename(src)
    fig,counts=plot_network(paths,tmp_path/'map',dpi=60)
    assert fig.is_file() and counts=={'cml_paths':3,'gauge':3,'temperature':2,'humidity':2}

def test_monthly_means_are_per_native_timestep(tmp_path):
    native=tmp_path/'native';native.mkdir()
    for role in ('cml','radar','gauge','temperature','humidity'):
        pd.DataFrame({'time':['2021-04-01 00:00:00','2021-04-02 00:00:00','2021-05-01 00:00:00'],
                      'availability_percent':[50,100,25]}).to_csv(native/f'availability_native_timestep_percentage_{role}.csv',index=False)
    fig,csv,monthly=monthly_from_native_csv(native,tmp_path/'out',dpi=60)
    assert fig.is_file() and csv.is_file()
    assert monthly.loc['2021-04','cml']==75
    assert monthly.loc['2021-05','radar']==25

def test_chmi_command_guard_and_opt_in(tmp_path):
    out=tmp_path/'new.nc'
    cmd=converter_command('humidity',year=2021,stage='build',root_dir=tmp_path,output=out,elements=['H'])
    assert str(CONVERTER) in cmd and '--elements' in cmd and 'H' in cmd and '--overwrite-nc' not in cmd
    assert run_converter(cmd,approved=False) is None
    out.write_bytes(b'unrelated')
    with pytest.raises(FileExistsError):converter_command('humidity',root_dir=tmp_path,output=out)
    with pytest.raises(ValueError):converter_command('unknown',root_dir=tmp_path)


def test_notebook_json_syntax_and_safe_defaults():
    import nbformat
    base=Path(__file__).resolve().parents[1]/'notebooks'
    files=sorted(base.glob('[0-9][0-9]_*.ipynb'))
    assert len(files)==9
    assert [p.name[:2] for p in files]==[f'{i:02d}' for i in range(9)]
    for file in files:
        nb=nbformat.read(file,as_version=4)
        nbformat.validate(nb)
        for i,cell in enumerate(nb.cells):
            if cell.cell_type=='code':compile(cell.source,f'{file.name}/cell_{i}','exec')
            assert '/home/jovyan/' not in cell.source, 'HPC paths belong in one config, not notebook code'
    assert 'RECOMPUTE = False' in '\n'.join(c.source for c in nbformat.read(base/'01_Graf_Rainfall.ipynb',4).cells)
    source='\n'.join(c.source for c in nbformat.read(base/'02_Rainfall_Validation.ipynb',4).cells)
    assert 'ALLOW_UNVERIFIED_V3 = False' in source
    assert 'RUN_GAUGE_VALIDATION = False' in source
    source='\n'.join(c.source for c in nbformat.read(base/'04_Spatial_Radar_Maps.ipynb',4).cells)
    assert 'RUN_MAPS=False' in source
