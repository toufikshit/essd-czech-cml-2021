"""Tests on small artificial files (no real data): file names, axes, 5-minute derivation, scripts and notebooks."""
import hashlib
from pathlib import Path
import pytest
import numpy as np
import pandas as pd
import xarray as xr
from mini_fixture import create,sha,EXPECTED_NAMES,IDS
from essd_cz.config import load_paths


def test_seven_role_names_axes_units_and_masks(tmp_path):
    paths=create(tmp_path/'fixture')
    names={k:p.name for k,p in ((k,Path(v)) for k,v in paths.items())}
    assert names==EXPECTED_NAMES
    assert load_paths(tmp_path/'fixture'/'data_paths.json')=={k:Path(p) for k,p in paths.items()}
    with xr.open_dataset(paths['qc']) as qc,xr.open_dataset(paths['rain'],mask_and_scale=False) as rain:
        assert qc['trsl'].dims==('time','sublink_id','cml_id')
        assert qc['length'].dims==('sublink_id','cml_id')
        assert qc['frequency'].attrs['units']=='MHz'
        assert rain['R_sublink_1min'].dims==('time_1min','sublink_id','cml_id')
        assert rain['R_sublink_1min'].attrs['units']=='mm h-1'
        assert rain.attrs['kr_model']=='ITU_2005'
        assert qc.sizes=={'time':2880,'sublink_id':2,'cml_id':3}
        assert rain.sizes=={'time_1min':2880,'sublink_id':2,'cml_id':3}
        assert [str(x) for x in qc.cml_id.values]==list(IDS)
        assert [str(x) for x in rain.cml_id.values]==list(IDS)
        assert len(np.unique(qc.sublink_id.values))==2
        assert np.isfinite(qc.trsl.values).sum()<qc.trsl.size
        assert np.isnan(qc.trsl.isel(sublink_id=1,cml_id=1).values).all()
        assert (rain.graf_wet_1min.values==255).any()
        assert np.isnan(rain.R_sublink_1min.values).any()
        assert np.all(np.diff(pd.DatetimeIndex(qc.time.values).asi8)==60_000_000_000)
        assert np.all(np.diff(pd.DatetimeIndex(rain.time_1min.values).asi8)==60_000_000_000)
    for role,n in [('gauge',10),('temperature',10),('humidity',10),('radar',5)]:
        with xr.open_dataset(paths[role]) as ds:
            t=pd.DatetimeIndex(ds.time.values)
            assert np.all(np.diff(t.asi8)==n*60_000_000_000)
            assert len(t)==2880//n
    with xr.open_dataset(paths['gauge']) as ds:
        assert ds.rr.attrs['units']=='mm per 10 min'
    with xr.open_dataset(paths['radar']) as ds:
        assert ds.rr.attrs['units']=='mm per 5 min'


def test_derived_5min_follows_complete_block_semantics(tmp_path):
    """Analysis 5-min values are an exact function of the released 1-min product."""
    from essd_cz import rain5
    paths=create(tmp_path/'fixture')
    with xr.open_dataset(paths['rain'],mask_and_scale=False) as rain:
        assert rain5.is_native_1min(rain)
        derived=rain5.derive_5min_dataset(rain)
        t1=rain['time_1min'].values.copy()
        R1=rain['R_sublink_1min'].values
        V1=(rain['graf_valid_1min'].values==1)&np.isfinite(R1)
        W1=(rain['graf_wet_1min'].values==1)&V1
    assert derived.sizes=={'time_5min':576,'sublink_id':2,'cml_id':3}
    assert np.array_equal(derived.time_5min.values,t1[::5])
    vm=V1.reshape(576,5,2,3).sum(1)
    wm=W1.reshape(576,5,2,3).sum(1)
    complete=vm==5
    # a 5-min rate exists exactly where all five minutes are valid
    assert np.array_equal(np.isfinite(derived.R_sublink_5min.values),complete)
    assert np.array_equal(derived.graf_valid_5min.values==1,complete)
    assert np.array_equal(derived.valid_minutes_5min.values,vm.astype(np.uint8))
    assert np.array_equal(derived.wet_minutes_5min.values,wm.astype(np.uint8))
    assert np.array_equal(derived.graf_wet_5min.values==255,~complete)
    assert np.array_equal((derived.graf_wet_5min.values==1),complete&(wm>0))
    # mean of the five 1-min rates, and physical mean AFTER directional 5-min masking
    means=np.where(V1,R1,0.).reshape(576,5,2,3).sum(axis=1)/5.  # exact only where complete
    assert np.allclose(derived.R_sublink_5min.values[complete],means[complete],rtol=0,atol=1e-6)
    cnt=np.isfinite(derived.R_sublink_5min.values).sum(axis=1)
    some=cnt>0
    expected=np.where(some,np.nansum(derived.R_sublink_5min.values,axis=1)/np.maximum(cnt,1),np.nan)
    assert np.allclose(derived.R_mean_5min.values[some],expected[some],rtol=0,atol=1e-6)
    assert np.isnan(derived.R_mean_5min.values[cnt==0]).all()
    # round trip through the file writer keeps every value and mask
    out=tmp_path/'derived.nc'
    rain5.derive_5min_file(paths['rain'],out)
    with xr.open_dataset(out,mask_and_scale=False) as back:
        a=back['R_sublink_5min'].values;b=derived['R_sublink_5min'].values
        assert np.array_equal(np.isnan(a),np.isnan(b)) and np.array_equal(a[~np.isnan(a)],b[~np.isnan(b)])
        assert np.array_equal(np.asarray(back['graf_wet_5min'].values,np.uint8),derived['graf_wet_5min'].values)


def test_inputs_read_only_and_filename_guards(tmp_path):
    paths=create(tmp_path/'data')
    before={key:sha(p) for key,p in paths.items()}
    with xr.open_dataset(paths['qc']) as qc:
        _=qc.trsl.isel(cml_id=0).values
    assert {key:sha(p) for key,p in paths.items()}==before
    bad=dict(paths);bad['rain']=paths['qc']
    from essd_cz.config import load_paths
    import json,pytest
    (tmp_path/'bad_config.json').write_text(json.dumps({'datasets':bad}))
    with pytest.raises(ValueError,match='wrong file name'):
        load_paths(tmp_path/'bad_config.json')


ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=sorted((ROOT/'scripts').rglob('*.py'))


def test_relative_paths_and_output_dir(tmp_path):
    from essd_cz.config import load_config
    import json
    paths=create(tmp_path/'data')
    cfg_dir=tmp_path/'config';cfg_dir.mkdir()
    rel={k:'../data/'+Path(v).name for k,v in paths.items()}
    (cfg_dir/'user_paths.json').write_text(json.dumps({'datasets':rel,'extras':{'radklim_dir':''},'output_dir':'../out'}))
    cfg=load_config(cfg_dir/'user_paths.json')
    assert cfg['datasets']=={k:Path(v).resolve() for k,v in paths.items()}
    assert cfg['output_dir']==(tmp_path/'out').resolve()
    assert cfg['extras']['radklim_dir'] is None


def test_scripts_compile():
    import py_compile
    assert len(SCRIPTS)>=12
    for script in SCRIPTS:
        py_compile.compile(str(script),doraise=True,cfile=str(Path(__import__('tempfile').gettempdir())/'essd_cz_compile.pyc'))


@pytest.mark.parametrize('name',['check_files.py','derive_5min.py','gauge_validation.py','wvd_validation.py',
                                 'figures/fig02_availability.py','figures/fig05_fig06_validation.py'])
def test_script_help(name):
    import subprocess,sys
    r=subprocess.run([sys.executable,'-B',str(ROOT/'scripts'/name),'--help'],capture_output=True,text=True)
    assert r.returncode==0,r.stderr
    assert 'usage' in r.stdout.lower()


def test_notebooks_use_repository_scripts():
    import json,re
    notebooks=sorted((ROOT/'notebooks').glob('*.ipynb'))
    assert notebooks
    for nb in notebooks:
        cells=json.loads(nb.read_text(encoding='utf-8'))['cells']
        code='\n'.join(''.join(c['source']) for c in cells if c['cell_type']=='code')
        for ref in re.findall(r"scripts/[A-Za-z0-9_/]+\.py",code):
            assert (ROOT/ref).is_file(),f'{nb.name}: {ref} not found'


def test_originals_unchanged():
    import subprocess,sys
    r=subprocess.run([sys.executable,str(ROOT/'tools_verify_originals.py')],capture_output=True,text=True)
    assert r.returncode==0,r.stdout+r.stderr
