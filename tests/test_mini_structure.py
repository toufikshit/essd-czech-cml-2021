"""GATE 1: inspect realistic mini NetCDF structure without performing CML physics."""
import hashlib
from pathlib import Path
import pytest
import numpy as np
import pandas as pd
import xarray as xr
from mini_fixture import create,sha,EXPECTED_NAMES,IDS
from essd_cz.config import load_paths
from essd_cz.rainfall import normalize_id,as_list
from essd_cz.metrics import rain_metrics,flag_metrics


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
        assert as_list(qc.cml_id.values)==list(IDS)
        assert as_list(rain.cml_id.values)==list(IDS)
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


def test_numerical_metrics_reject_mask_flag_and_value_differences():
    a=np.array([0,1.25,4.5,np.nan])
    assert rain_metrics(a,a)['numeric_pass']
    assert not rain_metrics(a,np.array([0,1.25,np.nan,np.nan]))['numeric_pass']
    assert rain_metrics(a,np.array([0,1.25,np.nan,np.nan]))['validity_mismatches']==1
    assert rain_metrics(a,np.array([0,1.25,4.6,np.nan]))['out_of_tolerance']==1
    assert not rain_metrics([np.nan],[np.nan])['numeric_pass']
    assert flag_metrics([0,1,255],[0,1,255])['pass']
    assert not flag_metrics([0,1,255],[0,1,0])['pass']
    assert normalize_id(b'85032_85231_3044\0')=='85032_85231_3044'


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
    with pytest.raises(ValueError,match='wrong filename'):
        load_paths(tmp_path/'bad_config.json')


def test_all_missing_direction_does_not_abort_other_direction(tmp_path):
    """Original producer records one failed direction and writes the other; wrapper must too.

    This uses a fake estimator solely to test error HANDLING, not scientific values.
    """
    pytest.importorskip('numba', reason='science extra not installed (pip install -e .[science])')
    from types import SimpleNamespace
    from essd_cz.rainfall import aggregate_exact
    paths=create(tmp_path/'input')
    def fake_graf(a,*_):
        valid=np.isfinite(a); wet=np.zeros(a.shape,bool)
        counts=valid.sum(axis=0).astype(np.int64)
        return wet,valid,np.zeros(a.shape[1]),np.zeros(a.shape[1]),counts,counts,np.zeros(a.shape[1],np.int64)
    def fake_meta2(ds,name,sd,cd,ns,nl,strings=False):
        return ds[name].transpose(sd,cd).values
    def fake_estimate(tr,wet,valid,times,L,f,pol):
        if np.isfinite(tr).sum()<2:
            return np.full(tr.size,np.nan,np.float32),'<2 valid observations'
        return np.where(np.isfinite(tr),0.,np.nan).astype(np.float32),''
    fake=SimpleNamespace(graf=fake_graf,meta2=fake_meta2,estimate=fake_estimate)
    with xr.open_dataset(paths['qc']) as qc:
        out=aggregate_exact(fake,qc,1,numba_threads=1)
    assert np.isfinite(out['R_sublink_5min'][:,0]).any()
    assert np.isnan(out['R_sublink_5min'][:,1]).all()
    assert np.all(out['graf_wet_5min'][:,1]==255)
    assert len(out['estimate_errors'])==1
