#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
import pycomlink as pycml

try:
    from numba import njit, prange, set_num_threads, get_num_threads
except Exception as exc:
    raise RuntimeError('numba is required; run in tf-env') from exc

VER='ITU2005_REBUILD_1.0'

def pol_hv(x):
    if isinstance(x,(bytes,np.bytes_)): x=bytes(x).decode('utf-8','ignore')
    s=str(x).replace('\x00','').strip().lower()
    if 'v' in s and 'h' not in s: return 'V'
    return 'H'

def call_kr(fn, model, **kw):
    errs=[]
    for key in ('a_b_approximation','approx_type'):
        try: return fn(**kw, **{key:model})
        except (TypeError,ValueError) as e: errs.append(f'{key}:{e}')
    raise RuntimeError('pycomlink rejected model: '+' | '.join(errs))

def estimate(trsl, wet, valid, times, L, f_hz, pol):
    empty=np.full(len(trsl),np.nan,np.float32)
    try:
        valid=np.asarray(valid,bool)&np.isfinite(trsl)
        wet=np.asarray(wet,bool)&valid
        if valid.sum()<2: return empty,'<2 valid observations'
        x=np.where(valid,np.asarray(trsl,float),np.nan)
        xda=xr.DataArray(x,dims='time',coords={'time':times})
        wda=xr.DataArray(wet,dims='time',coords={'time':times})
        b=pycml.processing.baseline.baseline_constant(trsl=xda,wet=wda,n_average_last_dry=1)
        Aobs=xr.where((xda-b)>0,xda-b,0.0)
        waa=pycml.processing.wet_antenna.waa_pastorek_2021_from_A_obs(
            A_obs=Aobs,f_Hz=float(f_hz),pol=pol_hv(pol),L_km=float(L),A_max=14.0,zeta=0.55,d=0.1)
        waa=xr.where(waa<0,0.0,waa); waa=xr.where(waa>Aobs,Aobs,waa)
        Ar=xr.where((Aobs-waa)>0,Aobs-waa,0.0)
        R=call_kr(pycml.processing.k_R_relation.calc_R_from_A,'ITU_2005',
                  A=Ar,L_km=float(L),f_GHz=float(f_hz)/1e9,pol=pol_hv(pol),R_min=0)
        out=np.asarray(getattr(R,'values',R),np.float32)
        out[valid & ~wet]=0.0; out[~valid]=np.nan; out[~np.isfinite(out)]=np.nan; out[out<0]=0
        return out,''
    except Exception as e:
        return empty,f'{type(e).__name__}: {e}'

@njit(parallel=True,cache=False)
def graf(values,window=60,q=0.8,factor=1.12,thr_min=60,det_min=45):
    nt,no=values.shape
    wet=np.zeros((nt,no),np.bool_); valid=np.zeros((nt,no),np.bool_)
    thresholds=np.full(no,np.nan,np.float64); qvals=np.full(no,np.nan,np.float64)
    finite_input=np.zeros(no,np.int64); valid_n=np.zeros(no,np.int64); wet_n=np.zeros(no,np.int64)
    left=window//2; right=window-left
    for j in prange(no):
        cs=np.empty(nt+1,np.float64); cs2=np.empty(nt+1,np.float64); cc=np.empty(nt+1,np.int64)
        cs[0]=0.; cs2[0]=0.; cc[0]=0; nf=0
        for i in range(nt):
            v=values[i,j]
            if np.isfinite(v):
                vv=np.float64(v); cs[i+1]=cs[i]+vv; cs2[i+1]=cs2[i]+vv*vv; cc[i+1]=cc[i]+1; nf+=1
            else:
                cs[i+1]=cs[i]; cs2[i+1]=cs2[i]; cc[i+1]=cc[i]
        finite_input[j]=nf
        st=np.empty(nt,np.float32); nth=0
        for i in range(nt):
            a=i-left
            if a<0: a=0
            b=i+right
            if b>nt: b=nt
            n=cc[b]-cc[a]
            if n>=thr_min and n>1:
                s=cs[b]-cs[a]; s2=cs2[b]-cs2[a]; var=(s2-s*s/n)/(n-1.)
                if var<0: var=0.
                st[i]=np.float32(np.sqrt(var)); nth+=1
            else: st[i]=np.nan
        if nth==0: continue
        fs=np.empty(nth,np.float64); k=0
        for i in range(nt):
            if np.isfinite(st[i]): fs[k]=st[i]; k+=1
        qv=np.quantile(fs,q); th=factor*qv; qvals[j]=qv; thresholds[j]=th
        nv=0; nw=0
        for i in range(nt):
            if det_min==thr_min: sv=st[i]
            else:
                a=i-left
                if a<0: a=0
                b=i+right
                if b>nt: b=nt
                n=cc[b]-cc[a]
                if n>=det_min and n>1:
                    s=cs[b]-cs[a]; s2=cs2[b]-cs2[a]; var=(s2-s*s/n)/(n-1.)
                    if var<0: var=0.
                    sv=np.sqrt(var)
                else: sv=np.nan
            if np.isfinite(sv):
                valid[i,j]=True; nv+=1
                if sv>th: wet[i,j]=True; nw+=1
        valid_n[j]=nv; wet_n[j]=nw
    return wet,valid,thresholds,qvals,finite_input,valid_n,wet_n

def decode(a):
    a=np.asarray(a)
    if a.dtype.kind=='S':
        return np.array([bytes(v).decode('utf-8','ignore').replace('\x00','').strip() for v in a.reshape(-1)],object).reshape(a.shape)
    if a.dtype.kind=='U': return a.astype(object)
    return np.array([str(v) for v in a.reshape(-1)],object).reshape(a.shape)

def meta2(ds,name,sd,cd,ns,nc,strings=False):
    da=ds[name]
    if da.dims==(sd,cd): a=da.values
    elif da.dims==(cd,sd): a=da.transpose(sd,cd).values
    elif da.dims==(cd,): a=np.broadcast_to(da.values[None,:],(ns,nc))
    elif da.dims==(sd,): a=np.broadcast_to(da.values[:,None],(ns,nc))
    else: raise ValueError(f'{name}.dims={da.dims}')
    return decode(a) if strings else np.asarray(a,float)

def make_out(path,ds,t5,cids,sids,L,f,pol,args):
    if path.exists() and not args.overwrite: raise FileExistsError(f'{path} exists; use --overwrite')
    path.parent.mkdir(parents=True,exist_ok=True)
    nc=Dataset(path,'w',format='NETCDF4'); n5=len(t5); ns=len(sids); nl=len(cids)
    nc.createDimension('time_5min',n5); nc.createDimension('sublink_id',ns); nc.createDimension('cml_id',nl)
    tv=nc.createVariable('time_5min','i8',('time_5min',)); tv[:]=(t5.asi8//1_000_000_000).astype(np.int64); tv.units='seconds since 1970-01-01 00:00:00 UTC'; tv.calendar='standard'
    cv=nc.createVariable('cml_id',str,('cml_id',)); cv[:]=np.asarray(cids,object)
    sv=nc.createVariable('sublink_id',str,('sublink_id',)); sv[:]=np.asarray(sids,object)
    lv=nc.createVariable('length','f4',('sublink_id','cml_id')); lv[:]=L.astype(np.float32); lv.units='km'
    fv=nc.createVariable('frequency','f8',('sublink_id','cml_id')); fv[:]=f; fv.units='Hz'
    pv=nc.createVariable('polarization',str,('sublink_id','cml_id')); pv[:]=np.asarray(pol,object)
    for name in ('site_a_latitude','site_a_longitude','site_b_latitude','site_b_longitude'):
        if name in ds and ds[name].dims==(args.cml_dim,):
            a=np.asarray(ds[name].values); v=nc.createVariable(name,'f8',('cml_id',)); v[:]=a
    tc=min(2016,n5); lc=min(args.output_link_chunk,nl); c3=(tc,1,lc); c2=(tc,lc); cm=(1,lc); z=dict(zlib=True,complevel=args.compression_level)
    R=nc.createVariable('R_sublink_5min','f4',('time_5min','sublink_id','cml_id'),fill_value=np.float32(np.nan),chunksizes=c3,**z); R.units='mm h-1'
    Rm=nc.createVariable('R_mean_5min','f4',('time_5min','cml_id'),fill_value=np.float32(np.nan),chunksizes=c2,**z); Rm.units='mm h-1'
    w5=nc.createVariable('graf_wet_5min','u1',('time_5min','sublink_id','cml_id'),fill_value=np.uint8(255),chunksizes=c3,**z); w5.flag_values=np.array([0,1],np.uint8); w5.flag_meanings='dry wet'
    v5=nc.createVariable('graf_valid_5min','u1',('time_5min','sublink_id','cml_id'),chunksizes=c3,**z)
    wm=nc.createVariable('wet_minutes_5min','u1',('time_5min','sublink_id','cml_id'),chunksizes=c3,**z)
    vm=nc.createVariable('valid_minutes_5min','u1',('time_5min','sublink_id','cml_id'),chunksizes=c3,**z)
    th=nc.createVariable('graf_threshold_db','f4',('sublink_id','cml_id'),fill_value=np.float32(np.nan),chunksizes=cm,**z)
    q80=nc.createVariable('graf_rstd_q80_db','f4',('sublink_id','cml_id'),fill_value=np.float32(np.nan),chunksizes=cm,**z)
    nc.title='CML-only rainfall using GRAF_RSTD_Q80 + Pastorek + ITU_2005'; nc.script=Path(__file__).name; nc.script_version=VER
    nc.source_cml_file=str(Path(args.input_nc).resolve()); nc.wet_dry_method='GRAF_RSTD_Q80'; nc.wet_dry_reference_data_used='none'
    nc.rstd_window=60; nc.graf_quantile=0.8; nc.graf_factor=1.12; nc.threshold_min_periods=60; nc.detection_min_periods=45
    nc.baseline_method='pycomlink baseline_constant'; nc.n_average_last_dry=1; nc.waa_method='Pastorek_2021'; nc.pastorek_A_max=14.; nc.pastorek_zeta=.55; nc.pastorek_d=.1; nc.kr_model='ITU_2005'; nc.min_output_valid_fraction=1.; nc.wet_aggregation='any wet minute'
    return nc,dict(R=R,Rm=Rm,w5=w5,v5=v5,wm=wm,vm=vm,th=th,q80=q80)

def main():
    p=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--input-nc',default='/home/jovyan/ESSD_publication_dataset/CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc')
    p.add_argument('--output-nc',default='/home/jovyan/ESSD_publication_dataset/cmlonly_graf_itu2005_pastorek/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_20200101_20211231.nc')
    p.add_argument('--summary-json',default='/home/jovyan/ESSD_publication_dataset/cmlonly_graf_itu2005_pastorek/run_summary.json')
    p.add_argument('--time-dim',default='time'); p.add_argument('--sublink-dim',default='sublink_id'); p.add_argument('--cml-dim',default='cml_id'); p.add_argument('--trsl-var',default='trsl')
    p.add_argument('--frequency-unit',choices=('Hz','MHz','GHz'),default='MHz'); p.add_argument('--length-unit',choices=('m','km'),default='km')
    p.add_argument('--link-batch',type=int,default=8); p.add_argument('--workers',type=int,default=8); p.add_argument('--output-link-chunk',type=int,default=8); p.add_argument('--compression-level',type=int,default=1); p.add_argument('--overwrite',action='store_true')
    a=p.parse_args(); t0=time.time(); set_num_threads(max(1,a.workers)); print('numba threads:',get_num_threads())
    graf(np.array([[1.,np.nan],[1.1,2.],[1.2,2.1],[1.3,2.2]],np.float32),2,.8,1.12,2,2)
    ds=xr.open_dataset(a.input_nc,decode_times=True)
    try:
        nt=ds.sizes[a.time_dim]; ns=ds.sizes[a.sublink_dim]; nl=ds.sizes[a.cml_dim]; times=pd.DatetimeIndex(ds[a.time_dim].values)
        if nt%5: raise RuntimeError(f'time length {nt} not divisible by 5')
        if len(times)>1 and not np.all(np.diff(times.asi8)==60*1_000_000_000): raise RuntimeError('time is not a complete regular 1-min grid')
        if nl!=2219 or ns!=2: raise RuntimeError(f'expected 2219 links x 2 sublinks, got {nl} x {ns}')
        cids=decode(ds[a.cml_dim].values).astype(str); sids=decode(ds[a.sublink_dim].values).astype(str)
        L=meta2(ds,'length',a.sublink_dim,a.cml_dim,ns,nl); L=L/1000 if a.length_unit=='m' else L
        f=meta2(ds,'frequency',a.sublink_dim,a.cml_dim,ns,nl); f=f*({'Hz':1,'MHz':1e6,'GHz':1e9}[a.frequency_unit])
        pol=meta2(ds,'polarization',a.sublink_dim,a.cml_dim,ns,nl,True); t5=times[::5]; n5=nt//5
        print('dims:',dict(time=nt,time_5min=n5,sublink_id=ns,cml_id=nl)); print('GHz range:',np.nanmin(f)/1e9,np.nanmax(f)/1e9)
        out=Path(a.output_nc); nc,v=make_out(out,ds,t5,cids,sids,L,f,pol,a); failures=[]; rows=[]
        try:
            for b0 in range(0,nl,a.link_batch):
                b1=min(b0+a.link_batch,nl); nb=b1-b0; print(f'[{b0}:{b1}/{nl}] elapsed={(time.time()-t0)/60:.1f} min',flush=True)
                tr=np.asarray(ds[a.trsl_var].isel({a.cml_dim:slice(b0,b1)}).transpose(a.time_dim,a.sublink_dim,a.cml_dim).values,np.float32)
                flat=np.ascontiguousarray(tr.reshape(nt,ns*nb)); wet,valid,thr,q80,fin,vn,wn=graf(flat,60,.8,1.12,60,45)
                R5=np.full((n5,ns,nb),np.nan,np.float32); W5=np.full((n5,ns,nb),255,np.uint8); V5=np.zeros((n5,ns,nb),np.uint8); WM=np.zeros((n5,ns,nb),np.uint8); VM=np.zeros((n5,ns,nb),np.uint8)
                def one(si):
                    sl=si//nb; lc=si%nb; r1,err=estimate(flat[:,si],wet[:,si],valid[:,si],times.values,L[sl,b0+lc],f[sl,b0+lc],pol[sl,b0+lc])
                    vv=valid[:,si]&np.isfinite(r1); rb=r1.reshape(n5,5); vb=vv.reshape(n5,5); wb=(wet[:,si]&vv).reshape(n5,5); counts=vb.sum(1).astype(np.uint8); wc=wb.sum(1).astype(np.uint8); complete=counts==5
                    r5=np.full(n5,np.nan,np.float32); r5[complete]=np.mean(rb[complete],axis=1,dtype=np.float64).astype(np.float32); w5=np.full(n5,255,np.uint8); w5[complete]=(wc[complete]>0).astype(np.uint8)
                    return si,sl,lc,r5,w5,complete.astype(np.uint8),wc,counts,err
                with ThreadPoolExecutor(max_workers=min(max(1,a.workers),ns*nb)) as ex:
                    futs=[ex.submit(one,si) for si in range(ns*nb)]
                    for fu in as_completed(futs):
                        si,sl,lc,r5,w5,v5,wm,vm,err=fu.result(); R5[:,sl,lc]=r5; W5[:,sl,lc]=w5; V5[:,sl,lc]=v5; WM[:,sl,lc]=wm; VM[:,sl,lc]=vm
                        rows.append(dict(cml_id=cids[b0+lc],sublink_id=sids[sl],finite_5min=int(np.isfinite(r5).sum()),graf_threshold_db=float(thr[si]) if np.isfinite(thr[si]) else np.nan,graf_q80_db=float(q80[si]) if np.isfinite(q80[si]) else np.nan,finite_input_1min=int(fin[si]),valid_detector_1min=int(vn[si]),wet_detector_1min=int(wn[si])))
                        if err: failures.append(dict(cml_id=cids[b0+lc],sublink_id=sids[sl],error=err))
                v['R'][:,:,b0:b1]=R5; cnt=np.sum(np.isfinite(R5),axis=1); sm=np.nansum(R5,axis=1); rm=np.full((n5,nb),np.nan,np.float32); ok=cnt>0; rm[ok]=(sm[ok]/cnt[ok]).astype(np.float32); v['Rm'][:,b0:b1]=rm
                v['w5'][:,:,b0:b1]=W5; v['v5'][:,:,b0:b1]=V5; v['wm'][:,:,b0:b1]=WM; v['vm'][:,:,b0:b1]=VM; v['th'][:,b0:b1]=thr.reshape(ns,nb).astype(np.float32); v['q80'][:,b0:b1]=q80.reshape(ns,nb).astype(np.float32); nc.sync()
        finally: nc.close()
        sp=Path(a.summary_json); sp.parent.mkdir(parents=True,exist_ok=True); pd.DataFrame(rows).to_csv(sp.with_name('rainfall_retrieval_per_sublink.csv'),index=False)
        summary=dict(script_version=VER,input_cml=str(Path(a.input_nc).resolve()),output=str(out.resolve()),n_links=nl,n_sublinks=ns,n_time_1min=nt,n_time_5min=n5,wet_dry_method='GRAF_RSTD_Q80',rstd_window=60,graf_quantile=.8,graf_factor=1.12,threshold_min_periods=60,detection_min_periods=45,baseline_method='pycomlink baseline_constant',n_average_last_dry=1,waa_method='Pastorek_2021',pastorek_A_max=14.,pastorek_zeta=.55,pastorek_d=.1,kr_model='ITU_2005',min_output_valid_fraction=1.,n_series_failures=len(failures),failures=failures[:100],elapsed_seconds=time.time()-t0)
        sp.write_text(json.dumps(summary,indent=2),encoding='utf-8'); print('\nDONE:',out); print('failures:',len(failures)); print('elapsed min:',(time.time()-t0)/60)
    finally: ds.close()
if __name__=='__main__': main()
