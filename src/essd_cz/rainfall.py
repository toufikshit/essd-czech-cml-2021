"""Small, educational functions around the ORIGINAL Graf/ITU2005 producer.

Scientific implementations (graf, estimate, meta2) are loaded UNCHANGED from
original/code/cml_rainfall/cml_graf_itu2005_rebuild.py. For exact numerical
comparison the full QC timeline is passed to those functions. Cropping QC to
14 days BEFORE estimating thresholds/baselines changes the scientific result.
"""
from __future__ import annotations
from pathlib import Path
import importlib.util
import hashlib
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT/'original/code/cml_rainfall/cml_graf_itu2005_rebuild.py'


def producer():
    """Import the original functions, without copying or editing their formulae."""
    spec=importlib.util.spec_from_file_location('essd_original_graf_itu2005', ORIGINAL)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def code_sha256():
    return hashlib.sha256(ORIGINAL.read_bytes()).hexdigest()


def normalize_id(value):
    """Keep structural IDs such as 85032_85231_3044 unchanged."""
    if isinstance(value,(bytes,np.bytes_)):
        return bytes(value).decode('utf-8','replace').replace('\0','').strip()
    if isinstance(value,(np.integer,int)):return str(int(value))
    if isinstance(value,(np.floating,float)) and np.isfinite(value) and float(value).is_integer():
        return str(int(value))
    return str(value).strip()


def as_list(values):
    return [normalize_id(x) for x in np.asarray(values).reshape(-1)]


def select_link_index(dataset, index=None, link_id=None):
    """Select by exact link ID or dataset index, not nearest location."""
    ids=as_list(dataset['cml_id'].values)
    if len(set(ids))!=len(ids):raise ValueError('Duplicate CML IDs')
    if (index is None)==(link_id is None):
        raise ValueError('Specify exactly one of index or link_id')
    if link_id is not None:
        wanted=normalize_id(link_id)
        if wanted not in ids:raise KeyError(f'Unknown physical CML ID {wanted}')
        return ids.index(wanted),wanted
    if not 0<=index<len(ids):raise IndexError(index)
    return index,ids[index]


def aggregate_exact(original, qc, index: int, numba_threads: int=2):
    """Recompute a SINGLE link with the original full-history algorithm.

    Steps: (1) decode geometry and 1-min TRSL, (2) rolling Graf RSTD/Q80,
    (3) constant baseline, (4) Pastorek WAA, (5) ITU2005 k--R,
    (6) aggregate FIVE complete 1-min retrievals into each 5-min sample.
    This reproduces the source script's per-sublink operations. We never
    write back to the input NetCDF or change its link order.
    """
    from numba import set_num_threads
    set_num_threads(max(1,numba_threads))
    nt=qc.sizes['time'];ns=qc.sizes['sublink_id'];nl=qc.sizes['cml_id']
    if ns!=2 or nt%5:raise ValueError(f'Expected two sublinks, 5-minute-aligned time; got {ns}, {nt}')
    times=pd.DatetimeIndex(qc['time'].values)
    if len(times)>1 and not np.all(np.diff(times.asi8)==60_000_000_000):
        raise ValueError('Original QC time is not a regular 1-min axis')
    tr=np.asarray(qc['trsl'].isel(cml_id=index).transpose('time','sublink_id').values,np.float32)
    wet,valid,th,q80,_,_,_=original.graf(np.ascontiguousarray(tr),60,.8,1.12,60,45)
    L=original.meta2(qc,'length','sublink_id','cml_id',ns,nl)[:,index]
    f=original.meta2(qc,'frequency','sublink_id','cml_id',ns,nl)[:,index]*1e6 # source is MHz
    pol=original.meta2(qc,'polarization','sublink_id','cml_id',ns,nl,True)[:,index]
    n5=nt//5
    rates=np.full((n5,ns),np.nan,np.float32)
    flags=np.full((n5,ns),255,np.uint8)
    valid5=np.zeros((n5,ns),np.uint8)
    wet_min=np.zeros((n5,ns),np.uint8)
    valid_min=np.zeros((n5,ns),np.uint8)
    estimate_errors=[]  # Original producer records failed sublinks; other direction survives.
    for sl in range(ns):
        r1,err=original.estimate(tr[:,sl],wet[:,sl],valid[:,sl],times.values,L[sl],f[sl],pol[sl])
        if err:
            estimate_errors.append({'sublink_id':as_list(qc['sublink_id'].values)[sl],
                                    'reason':str(err)})
        complete=valid[:,sl]&np.isfinite(r1)
        raw5=np.asarray(r1).reshape(n5,5)
        v5=complete.reshape(n5,5)
        w5=(wet[:,sl]&complete).reshape(n5,5)
        vc=v5.sum(axis=1).astype(np.uint8)
        wc=w5.sum(axis=1).astype(np.uint8)
        good=vc==5
        rates[good,sl]=np.mean(raw5[good],axis=1,dtype=np.float64).astype(np.float32)
        flags[good,sl]=(wc[good]>0).astype(np.uint8)
        valid5[:,sl]=good.astype(np.uint8)
        wet_min[:,sl]=wc
        valid_min[:,sl]=vc
    count=np.isfinite(rates).sum(axis=1)
    mean=np.full(n5,np.nan,np.float32)
    good=count>0
    mean[good]=(np.nansum(rates,axis=1)[good]/count[good]).astype(np.float32)
    return {
        'time':times[::5], 'R_sublink_5min':rates,'R_mean_5min':mean,
        'graf_wet_5min':flags,'graf_valid_5min':valid5,
        'wet_minutes_5min':wet_min,'valid_minutes_5min':valid_min,
        'graf_threshold_db':th,'graf_rstd_q80_db':q80,
        'sublink_ids':as_list(qc['sublink_id'].values),
        'full_source_minutes':nt,
        'estimate_errors':estimate_errors,
    }
