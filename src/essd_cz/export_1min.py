#!/usr/bin/env python3
"""Maintained entry point for the native 1-minute rainfall publication export.

This is the verified MetaCentrum exporter (handoff 2026-10-23, sha256
88013e76920d3da1e86c5dd46925ff95a73b5a497652257ff79f706e87f4fe10, preserved
byte-identical in docs/evidence/hpc_run_2026-10/cml_graf_export_1min_verified.py)
with ONE change: the three data paths are required arguments instead of
hard-coded HPC defaults. No scientific operation differs. It needs the
FULL-HISTORY QC input (Graf thresholds and baseline state use the entire
recorded history) and the validated full-history 5-minute reference, against
which every batch is compared before anything is kept.
"""
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


# Additional exporter for the manuscript release: the original algorithm above is
# unchanged. Work with FULL-HISTORY QC for Graf thresholds and baseline state,
# and write only the Apr--Oct 2021 native 1-minute rainfall results.

def exact_floats(name, computed, stored):
    a = np.asarray(computed)
    b = np.asarray(stored)
    same = (a.shape == b.shape) and np.all((a == b) | (np.isnan(a) & np.isnan(b)))
    if not same:
        if a.shape != b.shape:
            raise AssertionError(f'{name} shape mismatch: {a.shape} != {b.shape}')
        mismatch = ~((a == b) | (np.isnan(a) & np.isnan(b)))
        idx = tuple(int(q) for q in np.argwhere(mismatch)[0])
        raise AssertionError(f'{name} differs at {idx}: regenerated={a[idx]}, validated={b[idx]}, mismatches={np.count_nonzero(mismatch)}')


def make_native_file(path, src, release_times, cids, sids, L, f, pol, args):
    n = len(release_times)
    ns, nl = len(sids), len(cids)
    nc = Dataset(path, 'w', format='NETCDF4')
    try:
        nc.createDimension('time_1min', n)
        nc.createDimension('sublink_id', ns)
        nc.createDimension('cml_id', nl)
        tv = nc.createVariable('time_1min', 'i8', ('time_1min',))
        tv[:] = (release_times.asi8 // 1_000_000_000).astype(np.int64)
        tv.units = 'seconds since 1970-01-01 00:00:00 UTC'
        tv.calendar = 'standard'
        cv = nc.createVariable('cml_id', str, ('cml_id',))
        cv[:] = np.asarray(cids, object)
        sv = nc.createVariable('sublink_id', str, ('sublink_id',))
        sv[:] = np.asarray(sids, object)
        lv = nc.createVariable('length', 'f4', ('sublink_id', 'cml_id'))
        lv[:] = L.astype(np.float32)
        lv.units = 'km'
        fv = nc.createVariable('frequency', 'f8', ('sublink_id', 'cml_id'))
        fv[:] = f
        fv.units = 'Hz'
        pv = nc.createVariable('polarization', str, ('sublink_id', 'cml_id'))
        pv[:] = np.asarray(pol, object)
        for name in ('site_a_latitude', 'site_a_longitude', 'site_b_latitude', 'site_b_longitude'):
            if name in src and src[name].dims == (args.cml_dim,):
                v = nc.createVariable(name, 'f8', ('cml_id',))
                v[:] = np.asarray(src[name].values)

        lc = min(args.output_link_chunk, nl)
        tc = min(2880, n)
        z = dict(zlib=True, shuffle=True, complevel=args.compression_level)
        c3, c2, cm = (tc, 1, lc), (tc, lc), (1, lc)
        R = nc.createVariable('R_sublink_1min', 'f4', ('time_1min', 'sublink_id', 'cml_id'),
                              fill_value=np.float32(np.nan), chunksizes=c3, **z)
        R.units = 'mm h-1'
        R.long_name = 'Rainfall estimate at 1-minute time step for each directional sublink'
        Rm = nc.createVariable('R_mean_1min', 'f4', ('time_1min', 'cml_id'),
                               fill_value=np.float32(np.nan), chunksizes=c2, **z)
        Rm.units = 'mm h-1'
        Rm.long_name = 'Mean of available valid directional sublink rain rates at each minute'
        w = nc.createVariable('graf_wet_1min', 'u1', ('time_1min', 'sublink_id', 'cml_id'),
                              fill_value=np.uint8(255), chunksizes=c3, **z)
        w.flag_values = np.array([0, 1], np.uint8)
        w.flag_meanings = 'dry wet'
        w.comment = '255 signifies invalid or unavailable 1-minute rain-rate estimate'
        v = nc.createVariable('graf_valid_1min', 'u1', ('time_1min', 'sublink_id', 'cml_id'),
                              chunksizes=c3, **z)
        v.flag_values = np.array([0, 1], np.uint8)
        v.flag_meanings = 'invalid valid'
        v.comment = 'Valid detector classification AND finite retrieved rain rate'
        th = nc.createVariable('graf_threshold_db', 'f4', ('sublink_id', 'cml_id'),
                               fill_value=np.float32(np.nan), chunksizes=cm, **z)
        q80 = nc.createVariable('graf_rstd_q80_db', 'f4', ('sublink_id', 'cml_id'),
                                fill_value=np.float32(np.nan), chunksizes=cm, **z)
        nc.title = '1-minute CML-only rainfall using Graf RSTD + Pastorek + ITU 2005'
        nc.script = Path(__file__).name
        nc.script_version = VER + '_1MIN_EXPORT'
        nc.source_cml_file = Path(args.input_nc).name
        nc.reference_5min_file = Path(args.reference_5min).name
        nc.license = 'CC BY-NC 4.0'
        nc.data_provider = 'T-Mobile Czech Republic a.s.'
        nc.temporal_subset = '2021-04-01T00:00:00 to 2021-10-31T23:59:00 UTC'
        nc.wet_dry_method = 'GRAF_RSTD_Q80'
        nc.wet_dry_reference_data_used = 'none'
        nc.rstd_window = 60
        nc.graf_quantile = 0.8
        nc.graf_factor = 1.12
        nc.threshold_min_periods = 60
        nc.detection_min_periods = 45
        nc.baseline_method = 'pycomlink baseline_constant'
        nc.n_average_last_dry = 1
        nc.waa_method = 'Pastorek_2021'
        nc.pastorek_A_max = 14.0
        nc.pastorek_zeta = 0.55
        nc.pastorek_d = 0.1
        nc.kr_model = 'ITU_2005'
        nc.aggregation_to_5min = 'arithmetic mean of five complete valid 1-minute directional samples'
        nc.comment = 'Graf threshold and baseline are calculated on FULL-HISTORY QC, before restricting output to 2021-04-01 through 2021-10-31.'
        return nc, dict(R=R, Rm=Rm, wet=w, valid=v, th=th, q80=q80)
    except Exception:
        nc.close()
        raise


def main():
    p = argparse.ArgumentParser(description='Export native 1-minute rain-rate from unchanged validated Graf/Pastorek/ITU2005 retrieval. Compare aggregated values bit-for-bit to the existing 5-minute file.')
    p.add_argument('--input-nc', required=True, help='FULL-HISTORY QC NetCDF (CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc)')
    p.add_argument('--reference-5min', required=True, help='Validated full-history 5-minute reference NetCDF')
    p.add_argument('--output-nc', required=True, help='Target native 1-minute NetCDF (must not exist)')
    p.add_argument('--time-dim', default='time')
    p.add_argument('--sublink-dim', default='sublink_id')
    p.add_argument('--cml-dim', default='cml_id')
    p.add_argument('--trsl-var', default='trsl')
    p.add_argument('--frequency-unit', choices=('Hz','MHz','GHz'), default='MHz')
    p.add_argument('--length-unit', choices=('m','km'), default='km')
    p.add_argument('--link-batch', type=int, default=8)
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--output-link-chunk', type=int, default=8)
    p.add_argument('--compression-level', type=int, choices=range(1, 10), default=4)
    p.add_argument('--preflight-only', action='store_true', help='Calculate the first batch of links and compare to validated 5-minute data without writing any output file.')
    a = p.parse_args()
    if a.link_batch < 1 or a.workers < 1 or a.output_link_chunk < 1:
        p.error('link-batch, workers and output-link-chunk must be positive')
    start_wall = time.time()
    set_num_threads(a.workers)
    print('Numba threads:', get_num_threads(), flush=True)
    graf(np.array([[1., np.nan], [1.1, 2.], [1.2, 2.1], [1.3, 2.2]], np.float32), 2, .8, 1.12, 2, 2)
    ds = xr.open_dataset(a.input_nc, decode_times=True)
    ref = Dataset(a.reference_5min, 'r')
    ref.set_auto_mask(False)
    nc = None
    tmp = None
    try:
        nt = ds.sizes[a.time_dim]
        ns = ds.sizes[a.sublink_dim]
        nl = ds.sizes[a.cml_dim]
        times = pd.DatetimeIndex(ds[a.time_dim].values)
        if nt % 5 != 0 or len(times) > 1 and not np.all(np.diff(times.asi8) == 60 * 1_000_000_000):
            raise RuntimeError('Input time must be complete, regular 1-minute grid, and divisible by 5.')
        if nl != 2219 or ns != 2:
            raise RuntimeError(f'Expected 2219 physical links x 2 sublinks, got {nl} x {ns}.')
        start_ts = pd.Timestamp('2021-04-01T00:00:00')
        stop_ts = pd.Timestamp('2021-11-01T00:00:00')
        ia = int(np.searchsorted(times.asi8, start_ts.value))
        ib = int(np.searchsorted(times.asi8, stop_ts.value))
        if ia % 5 or ib % 5 or ib-ia != 308160 or times[ia] != start_ts or times[ib-1] != pd.Timestamp('2021-10-31T23:59:00'):
            raise RuntimeError(f'Release cropping not aligned or wrong length: i={ia}:{ib}, n={ib-ia}.')
        i5a, i5b = ia // 5, ib // 5
        sec5 = (times[i5a*5:i5b*5:5].asi8 // 1_000_000_000).astype(np.int64)
        oldtime = np.asarray(ref.variables['time_5min'][:], dtype=np.int64)
        j5a = int(np.searchsorted(oldtime, sec5[0]))
        j5b = j5a + len(sec5)
        if j5b > len(oldtime) or not np.array_equal(oldtime[j5a:j5b], sec5):
            raise RuntimeError('Validated 5-minute reference time axis does not align with selected 1-minute release.')
        cids = decode(ds[a.cml_dim].values).astype(str)
        sids = decode(ds[a.sublink_dim].values).astype(str)
        if not np.array_equal(decode(ref.variables['cml_id'][:]).astype(str), cids):
            raise RuntimeError('Reference cml_id order differs from QC input.')
        if not np.array_equal(decode(ref.variables['sublink_id'][:]).astype(str), sids):
            raise RuntimeError('Reference sublink_id order differs from QC input.')
        L = meta2(ds, 'length', a.sublink_dim, a.cml_dim, ns, nl)
        L = L / 1000 if a.length_unit == 'm' else L
        f = meta2(ds, 'frequency', a.sublink_dim, a.cml_dim, ns, nl)
        f = f * {'Hz': 1, 'MHz': 1e6, 'GHz': 1e9}[a.frequency_unit]
        pol = meta2(ds, 'polarization', a.sublink_dim, a.cml_dim, ns, nl, True)
        print(f'Input: {nt} minutes, release: {ib-ia} minutes, links: {nl}, frequencies: {np.nanmin(f)/1e9:.2f}--{np.nanmax(f)/1e9:.2f} GHz', flush=True)
        out = Path(a.output_nc)
        if not a.preflight_only:
            if out.exists():
                raise FileExistsError(f'Final file exists, refusing to overwrite: {out}')
            out.parent.mkdir(parents=True, exist_ok=True)
            tmp = out.with_name(out.name + '.partial')
            if tmp.exists():
                raise FileExistsError(f'Partial file already exists, inspect before retry: {tmp}')
            nc, v = make_native_file(tmp, ds, times[ia:ib], cids, sids, L, f, pol, a)
        n5 = nt // 5
        fail = []
        total_finite = 0
        batches = 0
        for b0 in range(0, 8 if a.preflight_only else nl, a.link_batch):
            b1 = min(b0 + a.link_batch, nl)
            nb = b1-b0
            print(f'Processing links {b0}:{b1}/{nl}', flush=True)
            tr = np.asarray(ds[a.trsl_var].isel({a.cml_dim: slice(b0,b1)}).transpose(a.time_dim,a.sublink_dim,a.cml_dim).values, np.float32)
            flat = np.ascontiguousarray(tr.reshape(nt, ns*nb))
            wet, valid, thr, q80, fin, vn, wn = graf(flat, 60, .8, 1.12, 60, 45)
            R5 = np.full((i5b-i5a,ns,nb), np.nan, np.float32)
            W5 = np.full((i5b-i5a,ns,nb),255,np.uint8)
            V5 = np.zeros((i5b-i5a,ns,nb),np.uint8)
            WM = np.zeros((i5b-i5a,ns,nb),np.uint8)
            VM = np.zeros((i5b-i5a,ns,nb),np.uint8)
            if not a.preflight_only:
                R1 = np.full((ib-ia,ns,nb), np.nan, np.float32)
                W1 = np.full((ib-ia,ns,nb),255,np.uint8)
                V1 = np.zeros((ib-ia,ns,nb),np.uint8)

            def one(si):
                sl, lc = si//nb, si%nb
                r1, err = estimate(flat[:,si], wet[:,si], valid[:,si], times.values,
                                   L[sl,b0+lc], f[sl,b0+lc], pol[sl,b0+lc])
                vv = valid[:,si] & np.isfinite(r1)
                rb = r1.reshape(n5,5)
                vb = vv.reshape(n5,5)
                wb = (wet[:,si] & vv).reshape(n5,5)
                counts = vb.sum(1).astype(np.uint8)
                wc = wb.sum(1).astype(np.uint8)
                complete = counts == 5
                r5 = np.full(n5,np.nan,np.float32)
                r5[complete] = np.mean(rb[complete],axis=1,dtype=np.float64).astype(np.float32)
                w5 = np.full(n5,255,np.uint8)
                w5[complete] = (wc[complete] > 0).astype(np.uint8)
                return si,sl,lc,r1[ia:ib],(wet[ia:ib,si] & vv[ia:ib]),vv[ia:ib],r5[i5a:i5b],w5[i5a:i5b],complete[i5a:i5b].astype(np.uint8),wc[i5a:i5b],counts[i5a:i5b],err

            with ThreadPoolExecutor(max_workers=min(a.workers,ns*nb)) as pool:
                futures = [pool.submit(one,si) for si in range(ns*nb)]
                for fu in as_completed(futures):
                    si,sl,lc,r1,w1,v1,r5,w5,v5,wm,vm,err = fu.result()
                    R5[:,sl,lc] = r5
                    W5[:,sl,lc] = w5
                    V5[:,sl,lc] = v5
                    WM[:,sl,lc] = wm
                    VM[:,sl,lc] = vm
                    if not a.preflight_only:
                        R1[:,sl,lc] = r1
                        W1[:,sl,lc] = np.where(v1,w1,255).astype(np.uint8)
                        V1[:,sl,lc] = v1.astype(np.uint8)
                    total_finite += int(np.isfinite(r1).sum())
                    if err:
                        fail.append((cids[b0+lc],sids[sl],err))
            sl5 = slice(j5a,j5b)
            exact_floats('R_sublink_5min', R5, ref.variables['R_sublink_5min'][sl5,:,b0:b1])
            for name, arr in [('graf_wet_5min',W5),('graf_valid_5min',V5),('wet_minutes_5min',WM),('valid_minutes_5min',VM)]:
                if not np.array_equal(arr, np.asarray(ref.variables[name][sl5,:,b0:b1])):
                    raise AssertionError(f'{name} does not match validated reference for links {b0}:{b1}')
            for name, arr in [('graf_threshold_db',thr.reshape(ns,nb).astype(np.float32)),('graf_rstd_q80_db',q80.reshape(ns,nb).astype(np.float32))]:
                exact_floats(name, arr, ref.variables[name][:,b0:b1])
            cnt = np.sum(np.isfinite(R5),axis=1)
            sm = np.nansum(R5,axis=1)
            rm5 = np.full((len(sec5),nb),np.nan,np.float32)
            ok = cnt>0
            rm5[ok] = (sm[ok]/cnt[ok]).astype(np.float32)
            exact_floats('R_mean_5min', rm5, ref.variables['R_mean_5min'][sl5,b0:b1])
            if not a.preflight_only:
                v['R'][:,:,b0:b1] = R1
                cnt1 = np.sum(np.isfinite(R1),axis=1)
                sm1 = np.nansum(R1,axis=1)
                mean1 = np.full((ib-ia,nb),np.nan,np.float32)
                ok1 = cnt1>0
                mean1[ok1] = (sm1[ok1]/cnt1[ok1]).astype(np.float32)
                v['Rm'][:,b0:b1] = mean1
                v['wet'][:,:,b0:b1] = W1
                v['valid'][:,:,b0:b1] = V1
                v['th'][:,b0:b1] = thr.reshape(ns,nb).astype(np.float32)
                v['q80'][:,b0:b1] = q80.reshape(ns,nb).astype(np.float32)
                nc.sync()
            batches += 1
            print(f'  VERIFIED against original 5-minute values / flags / thresholds; total finite 1-min directional samples: {total_finite:,}', flush=True)
        if nc is not None:
            nc.close()
            nc = None
            tmp.rename(out)
            summary = {'output': str(out), 'reference_5min': str(Path(a.reference_5min).resolve()),
                       'input_full_history': str(Path(a.input_nc).resolve()), 'n_physical_links':nl,
                       'n_directional_sublinks':ns, 'n_release_minutes':ib-ia,
                       'n_finite_directional_values':total_finite, 'batches':batches,
                       'identical_5min_comparison':'PASSED', 'retrieval_errors': fail[:100],
                       'n_retrieval_errors':len(fail),'elapsed_seconds':time.time()-start_wall}
            out.with_suffix(out.suffix+'.summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
            print('DONE:',out, flush=True)
        else:
            print('PREFLIGHT PASSED: first batch reproduces original 5-minute reference; no data written.',flush=True)
        print(f'Elapsed minutes: {(time.time()-start_wall)/60:.1f}',flush=True)
    finally:
        if nc is not None:
            nc.close()
        ref.close()
        ds.close()

if __name__ == '__main__':
    main()
