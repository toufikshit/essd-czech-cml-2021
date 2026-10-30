#!/usr/bin/env python3
"""Stage 4: native 1-minute CML rainfall from the approved QC candidate QC_CORRECTED_CANDIDATE_REV1.nc.

Copy of the final maintained exporter `essd-czech-cml-2021/src/essd_cz/export_1min.py` (sha256 abf56a31...; the
verified MetaCentrum exporter sha256 88013e76... with only the data paths made arguments). The scientific part
(imports, pol_hv, call_kr, estimate, graf, decode, meta2: Graf RSTD 60 min / q 0.8 / factor 1.12 / 60 and 45
minimum periods; pycomlink baseline_constant with n_average_last_dry = 1; Pastorek 2021 wet antenna, A_max 14,
zeta 0.55, d 0.1; ITU-R P.838-3 2005 k-R; R_min 0) is copied VERBATIM; s4_00_check_verbatim.py proves it.

Changes, all authorised by the author (Stage 4, 2026-10-28):
  * the 2219-CML requirement and the bit-for-bit comparison with the old 5-minute product are removed;
  * replaced by checks of the new input: dimensions, 1-minute time axis exactly 2021-04-01 00:00 .. 2021-10-31
    23:59, unique ids equal to the QC attribute count, finite frequency and length and Vertical polarization for
    retained sublinks, excluded sublinks entirely missing, no NaN-inconsistency of the output;
  * calibration: Graf thresholds and the dry baseline are computed on the QC period April-October 2021 (the QC
    candidate covers only this period); previously they used the full history (2021-01-13..19 and 04-01..11-05);
  * the input trsl is used as approved, including blackout-reconstructed (qc_flag bit 8) and interpolated (bit 16)
    minutes; the output variable `trsl_source_flag_1min` carries these two bits for every minute, so measured and
    reconstructed signal loss stay distinguishable in the rainfall product.

Usage:
  python -I s4_01_rain_1min.py --input-nc ../../stage3/QC_CORRECTED_CANDIDATE_REV1.nc --output-nc ../CZ_CML_RAINRATE_..._1MIN_...nc
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


# ---------------------------------------------------------------------------------------------------------------
# Stage 4 export (replaces the release-specific part of export_1min.py; the algorithm above is unchanged)
# ---------------------------------------------------------------------------------------------------------------
import hashlib
import sys

WIN_START = pd.Timestamp('2021-04-01T00:00:00')
WIN_LAST = pd.Timestamp('2021-10-31T23:59:00')
N_WIN = 308160


def sha256(path, block=1 << 24):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for ch in iter(lambda: fh.read(block), b''):
            h.update(ch)
    return h.hexdigest()


def make_native_file(path, src, release_times, cids, sids, L, f, pol, args, input_sha):
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
        sf = nc.createVariable('trsl_source_flag_1min', 'u1', ('time_1min', 'sublink_id', 'cml_id'),
                               chunksizes=c3, **z)
        sf.flag_masks = np.array([8, 16], np.uint8)
        sf.flag_meanings = 'blackout_total_loss_reconstructed short_gap_interpolated'
        sf.comment = ('Bits 8 and 16 of qc_flag of the QC input at the same minute and sublink. 0 = the total loss used '
                      'at this minute is measured (TX - RX) or missing. Bit 8: total loss reconstructed for a blackout '
                      '(max TX - min RX of the sublink, adaptation of Blettner et al., 2023; RX missing). Bit 16: total '
                      'loss linearly interpolated over a raw gap of at most 5 minutes. The rain rate at a flagged minute '
                      'is computed from that reconstructed or interpolated total loss; the Graf wet/dry decision uses a '
                      'centred 60-minute window and the dry baseline holds the last dry value, so flagged minutes can '
                      'also influence neighbouring minutes.')
        th = nc.createVariable('graf_threshold_db', 'f4', ('sublink_id', 'cml_id'),
                               fill_value=np.float32(np.nan), chunksizes=cm, **z)
        q80 = nc.createVariable('graf_rstd_q80_db', 'f4', ('sublink_id', 'cml_id'),
                                fill_value=np.float32(np.nan), chunksizes=cm, **z)
        nc.title = '1-minute CML-only rainfall using Graf RSTD + Pastorek + ITU 2005 (Stage 4 candidate from the corrected QC)'
        nc.script = Path(__file__).name
        nc.script_sha256 = sha256(__file__)
        nc.script_version = VER + '_1MIN_EXPORT_STAGE4'
        nc.source_cml_file = Path(args.input_nc).name
        nc.source_cml_sha256 = input_sha
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
        nc.comment = ('Graf thresholds and the dry baseline are calculated on the QC period 2021-04-01 to 2021-10-31 '
                      '(the corrected QC covers only this period; the earlier product used the full recorded history '
                      '2021-01-13..19 and 2021-04-01..11-05). Input total loss includes flagged blackout-reconstructed '
                      'and interpolated minutes; see trsl_source_flag_1min. No TX or RX values were generated.')
        return nc, dict(R=R, Rm=Rm, wet=w, valid=v, sf=sf, th=th, q80=q80)
    except Exception:
        nc.close()
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input-nc', required=True, help='Approved QC candidate (QC_CORRECTED_CANDIDATE_REV1.nc)')
    p.add_argument('--output-nc', required=True, help='Target native 1-minute NetCDF (must not exist)')
    p.add_argument('--expected-input-sha256', default='ca8afa95aacb50dacbd5f7d96f776c56a78c63d1a4888e9f721a29cdcfafde65')
    p.add_argument('--time-dim', default='time')
    p.add_argument('--sublink-dim', default='sublink_id')
    p.add_argument('--cml-dim', default='cml_id')
    p.add_argument('--trsl-var', default='trsl')
    p.add_argument('--link-batch', type=int, default=16)
    p.add_argument('--workers', type=int, default=12)
    p.add_argument('--output-link-chunk', type=int, default=8)
    p.add_argument('--compression-level', type=int, choices=range(1, 10), default=4)
    p.add_argument('--only-first', type=int, default=0, help='test: process only the first N CMLs')
    a = p.parse_args()
    start_wall = time.time()
    set_num_threads(a.workers)
    print('Numba threads:', get_num_threads(), flush=True)
    input_sha = sha256(a.input_nc)
    if a.expected_input_sha256 and input_sha != a.expected_input_sha256:
        sys.exit(f'input sha256 {input_sha} differs from the approved candidate {a.expected_input_sha256}')
    graf(np.array([[1., np.nan], [1.1, 2.], [1.2, 2.1], [1.3, 2.2]], np.float32), 2, .8, 1.12, 2, 2)
    ds = xr.open_dataset(a.input_nc, decode_times=True)
    if a.only_first:
        ds = ds.isel({a.cml_dim: slice(0, a.only_first)})
    nc = None
    tmp = None
    try:
        # ---- input checks (replace the release-specific 2219 and reference checks) ----
        nt = ds.sizes[a.time_dim]; ns = ds.sizes[a.sublink_dim]; nl = ds.sizes[a.cml_dim]
        times = pd.DatetimeIndex(ds[a.time_dim].values)
        if nt != N_WIN or ns != 2:
            raise RuntimeError(f'Expected {N_WIN} minutes x 2 sublinks, got {nt} x {ns}')
        if times[0] != WIN_START or times[-1] != WIN_LAST or not np.all(np.diff(times.asi8) == 60 * 1_000_000_000):
            raise RuntimeError('Input time axis is not the regular 1-minute grid 2021-04-01 00:00 .. 2021-10-31 23:59')
        cids = decode(ds[a.cml_dim].values).astype(str)
        sids = decode(ds[a.sublink_dim].values).astype(str)
        if len(set(cids)) != nl:
            raise RuntimeError('Duplicate cml_id in input')
        if not a.only_first and nl != int(ds.attrs['n_cml_retained']):
            raise RuntimeError(f'cml count {nl} differs from n_cml_retained {ds.attrs["n_cml_retained"]}')
        status = np.asarray(ds['qc_sublink_status'].transpose(a.sublink_dim, a.cml_dim).values)
        L = meta2(ds, 'length', a.sublink_dim, a.cml_dim, ns, nl)
        f = meta2(ds, 'frequency', a.sublink_dim, a.cml_dim, ns, nl) * 1e6   # QC stores MHz
        pol = meta2(ds, 'polarization', a.sublink_dim, a.cml_dim, ns, nl, True)
        ret = status == 0
        if not (ret.any(axis=0)).all():
            raise RuntimeError('A CML without a retained sublink is present')
        if not (np.isfinite(f[ret]).all() and np.isfinite(L[ret]).all() and (L[ret] > 0).all()):
            raise RuntimeError('Retained sublink with missing frequency or length')
        if set(np.unique(pol)) != {'Vertical'}:
            raise RuntimeError(f'Unexpected polarization values {set(np.unique(pol))}')
        print(f'Input: {nt} minutes {times[0]}..{times[-1]}, links: {nl}, retained sublinks: {int(ret.sum())}, '
              f'frequencies: {np.nanmin(f)/1e9:.2f}--{np.nanmax(f)/1e9:.2f} GHz, sha256 {input_sha}', flush=True)
        out = Path(a.output_nc)
        if out.exists():
            raise FileExistsError(f'Final file exists, refusing to overwrite: {out}')
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_name(out.name + '.partial')
        if tmp.exists():
            raise FileExistsError(f'Partial file already exists, inspect before retry: {tmp}')
        nc, v = make_native_file(tmp, ds, times, cids, sids, L, f, pol, a, input_sha)
        fail = []; total_finite = 0; total_flagged = {8: 0, 16: 0}; batches = 0
        for b0 in range(0, nl, a.link_batch):
            b1 = min(b0 + a.link_batch, nl); nb = b1 - b0
            print(f'Processing links {b0}:{b1}/{nl}  {(time.time()-start_wall)/60:.1f} min', flush=True)
            tr = np.asarray(ds[a.trsl_var].isel({a.cml_dim: slice(b0, b1)}).transpose(a.time_dim, a.sublink_dim, a.cml_dim).values, np.float32)
            qf = np.asarray(ds['qc_flag'].isel({a.cml_dim: slice(b0, b1)}).transpose(a.time_dim, a.sublink_dim, a.cml_dim).values, np.uint8)
            for lc in range(nb):
                for sl in range(ns):
                    if status[sl, b0 + lc] != 0 and np.isfinite(tr[:, sl, lc]).any():
                        raise RuntimeError(f'Excluded sublink {cids[b0+lc]} s{sl} has finite trsl')
            flat = np.ascontiguousarray(tr.reshape(nt, ns * nb))
            wet, valid, thr, q80, fin, vn, wn = graf(flat, 60, .8, 1.12, 60, 45)
            R1 = np.full((nt, ns, nb), np.nan, np.float32)
            W1 = np.full((nt, ns, nb), 255, np.uint8)
            V1 = np.zeros((nt, ns, nb), np.uint8)

            def one(si):
                sl, lc = si // nb, si % nb
                r1, err = estimate(flat[:, si], wet[:, si], valid[:, si], times.values,
                                   L[sl, b0 + lc], f[sl, b0 + lc], pol[sl, b0 + lc])
                vv = valid[:, si] & np.isfinite(r1)
                return si, sl, lc, r1, (wet[:, si] & vv), vv, err

            with ThreadPoolExecutor(max_workers=min(a.workers, ns * nb)) as pool:
                futures = [pool.submit(one, si) for si in range(ns * nb)]
                for fu in as_completed(futures):
                    si, sl, lc, r1, w1, v1, err = fu.result()
                    R1[:, sl, lc] = r1
                    W1[:, sl, lc] = np.where(v1, w1, 255).astype(np.uint8)
                    V1[:, sl, lc] = v1.astype(np.uint8)
                    total_finite += int(np.isfinite(r1).sum())
                    if err and status[sl, b0 + lc] == 0:
                        fail.append((cids[b0 + lc], sids[sl], err))
            if (np.isfinite(R1) != (V1 == 1)).any():
                raise RuntimeError('finite rain rate and graf_valid_1min disagree')
            src_flag = (qf & 24).astype(np.uint8)
            for bit in (8, 16):
                total_flagged[bit] += int(((src_flag & bit) > 0)[np.isfinite(R1)].sum())
            v['R'][:, :, b0:b1] = R1
            cnt1 = np.sum(np.isfinite(R1), axis=1); sm1 = np.nansum(R1, axis=1)
            mean1 = np.full((nt, nb), np.nan, np.float32); ok1 = cnt1 > 0
            mean1[ok1] = (sm1[ok1] / cnt1[ok1]).astype(np.float32)
            v['Rm'][:, b0:b1] = mean1
            v['wet'][:, :, b0:b1] = W1
            v['valid'][:, :, b0:b1] = V1
            v['sf'][:, :, b0:b1] = src_flag
            v['th'][:, b0:b1] = thr.reshape(ns, nb).astype(np.float32)
            v['q80'][:, b0:b1] = q80.reshape(ns, nb).astype(np.float32)
            nc.sync(); batches += 1
        nc.close(); nc = None
        tmp.rename(out)
        summary = {'output': str(out), 'output_sha256': sha256(out), 'input_qc': str(Path(a.input_nc).resolve()),
                   'input_sha256': input_sha, 'n_physical_links': nl, 'n_directional_sublinks': ns,
                   'n_retained_sublinks': int(ret.sum()), 'n_minutes': nt,
                   'n_finite_directional_values': total_finite,
                   'n_finite_values_at_blackout_reconstructed_minutes': total_flagged[8],
                   'n_finite_values_at_interpolated_minutes': total_flagged[16],
                   'batches': batches, 'retrieval_errors_retained_sublinks': fail[:100],
                   'n_retrieval_errors_retained_sublinks': len(fail), 'elapsed_seconds': time.time() - start_wall}
        out.with_suffix(out.suffix + '.summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
        print(json.dumps(summary, indent=2), flush=True)
        print('DONE:', out, flush=True)
    finally:
        if nc is not None:
            nc.close()
        ds.close()


if __name__ == '__main__':
    main()
