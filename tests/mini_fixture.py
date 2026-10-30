#!/usr/bin/env python3
"""Make small, artificial (NOT observed) seven-file CML data for the tests.

Uses the same variable and dimension names and units as the released files. The values only test the data
handling (file names, axes, the 5-minute derivation); they have no scientific meaning.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import xarray as xr

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'src'))
from essd_cz.config import EXPECTED_NAMES

IDS=np.array(['85032_85231_3044','093401_90001_0301','10002_11004_0718'],dtype=str)
SUBS=np.array(['0','1'],dtype=str)
START='2021-04-01'
MINUTES=2880  # 2 days; fixed seed for stable regression results


def netcdf_write(ds, target):
    """Prefer real NETCDF4 on HPC; fallback NETCDF3 for offline layout tests."""
    try:import netCDF4 # noqa: F401
    except ImportError: engine='scipy'
    else:engine='netcdf4'
    if engine=='scipy':
        # Classic NetCDF has no native unsigned byte. Preserve flag VALUES, not dtype.
        for key in ('graf_wet_1min','graf_valid_1min',
                    'graf_wet_5min','graf_valid_5min','wet_minutes_5min','valid_minutes_5min'):
            if key in ds:ds[key]=ds[key].astype(np.int16)
    ds.to_netcdf(target,engine=engine)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def create(root:Path):
    """Create isolated 7-file example under root, with no real data or writebacks."""
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    tt=pd.date_range(START,periods=MINUTES,freq='min')
    t5=tt[::5];t10=tt[::10];n5=len(t5)
    rng=np.random.default_rng(20211013)
    signal=np.empty((MINUTES,2,3),np.float32)
    for sl in range(2):
        for c in range(3):
            wave=0.14*np.sin(np.arange(MINUTES)/37.+sl)+0.03*rng.standard_normal(MINUTES)
            wet_peak=np.zeros(MINUTES);wet_peak[500+9*c:650+9*c]=1.9+0.2*sl
            wet_peak[1900:2050]=3.0+0.25*c
            signal[:,sl,c]=(77+1.7*c+wave+wet_peak).astype('f4')
    signal[64:69,0,0]=np.nan # incomplete five-minute block
    signal[1150:1220,:,2]=np.nan # both directional signals absent
    signal[:,1,1]=np.nan # completely absent directional series
    # Real QC contains the combined TRSL and/or TX/RX: make both available here.
    tsl=np.where(np.isfinite(signal),0.0,np.nan).astype(np.float32)
    rsl=(-signal).astype(np.float32)
    lengths=np.array([[1.710,4.2,7.5],[1.710,4.2,7.5]],np.float32)
    freqs=np.array([[25088,39000,73000],[25088,39000,73000]],np.float64)
    pols=np.array([['Vertical']*3,['Horizontal']*3],dtype=str)
    qc=xr.Dataset(
      {'trsl':(('time','sublink_id','cml_id'),signal,{'units':'dB'}),
       'tsl':(('time','sublink_id','cml_id'),tsl,{'units':'dBm'}),
       'rsl':(('time','sublink_id','cml_id'),rsl,{'units':'dBm'}),
       'length':(('sublink_id','cml_id'),lengths,{'units':'km'}),
       'frequency':(('sublink_id','cml_id'),freqs,{'units':'MHz'}),
       'polarization':(('sublink_id','cml_id'),pols),
       'site_a_longitude':(('cml_id',),np.array([14.4,14.6,14.9])),
       'site_a_latitude':(('cml_id',),np.array([50.0,50.1,50.2])),
       'site_b_longitude':(('cml_id',),np.array([14.42,14.63,14.95])),
       'site_b_latitude':(('cml_id',),np.array([50.02,50.13,50.24]))},
      coords={'time':tt,'sublink_id':SUBS,'cml_id':IDS},
      attrs={'preprocessing':'SYNTHETIC: six-step QC output SHAPE, not a QC implementation',
             'scientific_status':'ARTIFICIAL data; not a CHMI or T-Mobile observation'})
    raw=xr.Dataset({'trsl':(('time','sublink_id','cml_id'),
          np.concatenate([signal,signal[:,:,:1]],axis=2))},
          coords={'time':tt,'sublink_id':SUBS,'cml_id':np.r_[IDS,['extra_unfiltered_link']]})
    # The published rain product is NATIVE 1-MINUTE; 5-minute analysis values
    # are derived from it (essd_cz.rain5), exactly like the real release.
    coords={'time_1min':tt,'sublink_id':SUBS,'cml_id':IDS}
    rainfall=xr.Dataset({
       'R_sublink_1min':(('time_1min','sublink_id','cml_id'),np.full((MINUTES,2,3),np.nan,np.float32),{'units':'mm h-1'}),
       'R_mean_1min':(('time_1min','cml_id'),np.full((MINUTES,3),np.nan,np.float32),{'units':'mm h-1'}),
       'graf_wet_1min':(('time_1min','sublink_id','cml_id'),np.full((MINUTES,2,3),255,np.uint8)),
       'graf_valid_1min':(('time_1min','sublink_id','cml_id'),np.zeros((MINUTES,2,3),np.uint8)),
       'graf_threshold_db':(('sublink_id','cml_id'),np.full((2,3),np.nan,np.float32)),
       'graf_rstd_q80_db':(('sublink_id','cml_id'),np.full((2,3),np.nan,np.float32))},
       coords=coords,attrs={'kr_model':'ITU_2005', 'wet_dry_method':'GRAF_RSTD_Q80',
                            'aggregation_to_5min':'arithmetic mean of five complete valid 1-minute directional samples',
                            'scientific_status':'SYNTHETIC fixture; no scientific claim'})
    # Placeholder values: they test ONLY the data handling, not the science.
    if True:
        A=np.zeros((MINUTES,2,3),np.float32)
        A[500:650,:,:]=2.3
        A[1112:1115,0,0]=np.nan   # incomplete 5-minute blocks after derivation
        A[:,:,1]=np.nan           # link without any directional retrieval
        valid=np.isfinite(A)
        rainfall['R_sublink_1min'].values[:]=A
        counts=valid.sum(axis=1)
        totals=np.nansum(A,axis=1)
        means=np.full((MINUTES,3),np.nan,np.float32)
        good=counts>0
        means[good]=(totals[good]/counts[good]).astype(np.float32)
        rainfall['R_mean_1min'].values[:]=means
        rainfall['graf_wet_1min'].values[:]=np.where(valid,(A>0).astype('u1'),255)
        rainfall['graf_valid_1min'].values[:]=valid.astype('u1')
        rainfall['graf_threshold_db'].values[:]=np.full((2,3),0.3,'f4')
        rainfall['graf_rstd_q80_db'].values[:]=np.full((2,3),0.26785715,'f4')
    gw=np.full((len(t10),3),0.,np.float32);gw[30:65,0]=0.5;gw[20,2]=np.nan
    radar=np.zeros((len(t5),2),np.float32);radar[100:125,0]=0.3;radar[15,1]=np.nan
    datasets={
     'raw':raw,'qc':qc,'rain':rainfall,
     'radar':xr.Dataset({'rr':(('time','cml_id'),radar,{'units':'mm per 5 min'})},coords={'time':t5,'cml_id':IDS[[0,2]]}),
     'gauge':xr.Dataset({'rr':(('time','WSI'),gw,{'units':'mm per 10 min'})},coords={'time':t10,'WSI':['G001','G002','G003']}),
     'temperature':xr.Dataset({'T':(('time','station'),np.full((len(t10),2),15.,np.float32))},coords={'time':t10,'station':['T01','T02']}),
     'humidity':xr.Dataset({'H':(('time','station'),np.full((len(t10),2),66.,np.float32))},coords={'time':t10,'station':['H01','H02']})}
    filemap={}
    for role,ds in datasets.items():
        path=root/EXPECTED_NAMES[role]
        if path.exists():raise FileExistsError('Refusing to overwrite fixture: '+str(path))
        netcdf_write(ds,path);filemap[role]=str(path)
    config={'description':'Artificial 2-day schema fixture. NEVER published data.',
            'datasets':filemap}
    (root/'data_paths.json').write_text(json.dumps(config,indent=2)+'\n')
    manifest={role:{'sha256':sha(path),'size':Path(path).stat().st_size} for role,path in filemap.items()}
    (root/'fixture_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('WROTE 7 SYNTHETIC NetCDFs, config, checksums in',root)
    return filemap

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir',type=Path,required=True)
    args=parser.parse_args()
    create(args.out_dir)
