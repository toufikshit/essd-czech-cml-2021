"""GATE 2: integration tests for index/time/flag alignment (synthetic rates only).

These tests monkeypatch ONLY the scientific recomputation boundary; they
DO NOT claim to test Graf/Pastorek physics.  The unchanged real source is
checked separately on HPC by run_mini_science.py.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from mini_fixture import create,sha
from essd_cz import compare as cmp
from essd_cz import rain5


def fake_recompute_factory(canonical):
    def compute(_original,qc,index,numba_threads=2):
        rid=str(qc.cml_id.values[index]);ids=list(map(str,canonical.cml_id.values))
        if rid not in ids:raise AssertionError(rid)
        j=ids.index(rid)
        def get(n):
            arr=canonical[n].isel(cml_id=j).values
            if n in ('graf_threshold_db','graf_rstd_q80_db'):return arr
            if n=='R_mean_5min':return arr
            return arr
        fields=('R_sublink_5min','R_mean_5min','graf_wet_5min','graf_valid_5min',
                'wet_minutes_5min','valid_minutes_5min','graf_threshold_db','graf_rstd_q80_db')
        # canonical original NetCDF uses time,sublink,link axes
        return {'time':pd.DatetimeIndex(qc.time.values)[::5],
                'full_source_minutes':len(qc.time),
                'sublink_ids':list(map(str,qc.sublink_id.values)),
                **{name:get(name) for name in fields}}
    return compute


@pytest.fixture
def canonical(tmp_path,monkeypatch):
    paths=create(tmp_path/'baseline')
    # The published fixture is native 1-minute; the canonical "recomputation"
    # is its exact 5-minute derivation, so an untouched release must PASS and
    # any mutation of the released 1-minute values must FAIL.
    with xr.open_dataset(paths['rain'],mask_and_scale=False) as native:
        ds=rain5.derive_5min_dataset(native)
    monkeypatch.setattr(cmp,'producer',lambda: object())
    monkeypatch.setattr(cmp,'aggregate_exact',fake_recompute_factory(ds))
    monkeypatch.setattr(cmp,'code_sha256',lambda:'FAKE_PURE_COMPARISON_TEST')
    yield paths,ds


def test_baseline_and_order_shuffle_pass(canonical,tmp_path):
    paths,ds=canonical
    report=cmp.compare(paths,tmp_path/'comparison',link_indices=[0,2],days=2,make_plot=False)
    assert report['all_pass']
    assert all(row['compared_5min_times']==576 for row in report['links'])
    other=create(tmp_path/'shuffle',scenario='reordered')
    reordered=cmp.compare(other,tmp_path/'out_shuffle',link_indices=[0,2],days=2,make_plot=False)
    assert reordered['all_pass']


@pytest.mark.parametrize('scenario,expected',[
    ('numeric_corruption','out_of_tolerance'),
    ('mask_corruption','validity_mismatches'),
    ('flag_corruption','n_disagree')])
def test_corruptions_are_detected(canonical,tmp_path,scenario,expected):
    _,ds=canonical
    paths=create(tmp_path/scenario,scenario=scenario)
    r=cmp.compare(paths,tmp_path/('out_'+scenario),link_indices=[0],days=2,make_plot=False)
    assert not r['all_pass']
    fields=r['links'][0]['sublinks'][0]
    assert any(v.get(expected,0)>0 for v in fields.values() if isinstance(v,dict))


def test_timestamp_shift_is_rejected(canonical,tmp_path):
    paths=create(tmp_path/'shift',scenario='time_shift')
    with pytest.raises(ValueError,match='regular 1-minute grid'):
        cmp.compare(paths,tmp_path/'out_shift',link_indices=[0],days=2,make_plot=False)
