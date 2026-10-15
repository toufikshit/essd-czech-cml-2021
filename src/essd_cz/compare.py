"""Numeric identity check: original full-record recomputation vs published NetCDF.

Do NOT compare a 14-day-only recalculation against a full-history baseline.
Compute on the entire recorded QC time axis first, THEN subset output.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr
from .rainfall import producer,aggregate_exact,normalize_id,as_list,select_link_index,code_sha256
from .metrics import rain_metrics,flag_metrics
from .plots import comparison_plot

FIELDS=('R_sublink_5min','R_mean_5min','graf_wet_5min','graf_valid_5min',
        'wet_minutes_5min','valid_minutes_5min','graf_threshold_db','graf_rstd_q80_db')


def compare(paths, output_dir, *, link_indices=(3,), start='2021-04-01', days=14,
            atol=1e-4, rtol=1e-4, numba_threads=2, make_plot=True):
    """Return a machine-readable outcome; raise on alignment or missing variables.

    This is an acceptance test, not an optimization or replacement algorithm.
    2020--2021 input permits the same full-history Graf threshold as the
    published product. A shorter input may NOT reproduce it identically.
    """
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    begin=pd.Timestamp(start);end=begin+pd.Timedelta(days=days)
    if days<1:raise ValueError('days must be positive')
    model=producer()
    report={'comparison_period':[str(begin),str(end)+' (exclusive)'],
            'original_script_sha256':code_sha256(),
            'method':'original full-history Graf/Pastorek/ITU2005 per-link functions',
            'threshold_context':'full QC source time coordinate, not the comparison subset',
            'links':[]}
    with xr.open_dataset(paths['qc'],decode_times=True) as qc, \
         xr.open_dataset(paths['rain'],decode_times=True,mask_and_scale=False) as rain:
        if 'R_sublink_5min' not in rain:raise KeyError('Published file lacks R_sublink_5min')
        kr=str(rain.attrs.get('kr_model',''))
        if '2005' not in kr:raise ValueError(f'Expected ITU2005 rain product, got kr_model={kr!r}')
        ref_ids=as_list(rain['cml_id'].values)
        if len(set(ref_ids))!=len(ref_ids):raise ValueError('Duplicate published CML IDs')
        ref_times=pd.DatetimeIndex(rain['time_5min'].values)
        if not ref_times.is_monotonic_increasing or ref_times.has_duplicates:
            raise ValueError('Published time axis not unique and increasing')
        ref_sl=as_list(rain['sublink_id'].values)
        if len(set(ref_sl))!=len(ref_sl):raise ValueError('Duplicate published sublink IDs')
        for idx in link_indices:
            index,cml_id=select_link_index(qc,index=int(idx))
            if cml_id not in ref_ids:raise KeyError(f'Published rainfall lacks QC CML ID {cml_id}')
            rid=ref_ids.index(cml_id)
            print(f'COMPUTING full-history physical link index={index}, cml_id={cml_id}',flush=True)
            computed=aggregate_exact(model,qc,index,numba_threads=numba_threads)
            source_t=computed['time']
            selected=np.flatnonzero((source_t>=begin)&(source_t<end))
            if not selected.size:raise ValueError(f'No source 5-min time slots in [{begin},{end})')
            target_t=source_t[selected]
            ri=ref_times.get_indexer(target_t)
            if np.any(ri<0):
                missing=target_t[np.flatnonzero(ri<0)[0]]
                raise ValueError(f'Published rainfall missing timestamp {missing}; do not compare adjacent times')
            r0,r1=int(ri[0]),int(ri[-1])+1
            if not np.array_equal(ri,np.arange(r0,r1)):
                raise ValueError('Reference timestamps not consecutive or in a different order')
            result={'index':index,'cml_id':cml_id,'minutes_used_for_baseline':computed['full_source_minutes'],
                    'compared_5min_times':len(selected),'results':{},'sublinks':[],
                    'estimate_errors':computed.get('estimate_errors',[])}
            all_pass=True;combined_model=[];combined_ref=[]
            for sl_name in computed['sublink_ids']:
                if sl_name not in ref_sl:raise ValueError(f'Published rain missing sublink {sl_name}')
                si=computed['sublink_ids'].index(sl_name)
                rsi=ref_sl.index(sl_name)
                per={'sublink_id':sl_name}
                for field in FIELDS:
                    if field not in rain:
                        per[field]={'skipped':'not present in released NetCDF'}
                        continue
                    src_arr=computed[field]
                    if field=='R_mean_5min':continue
                    if field in ('graf_threshold_db','graf_rstd_q80_db'):
                        recon=np.asarray([src_arr[si]])
                        published=np.asarray([rain[field].isel(cml_id=rid,sublink_id=rsi).values])
                    else:
                        recon=np.asarray(src_arr[selected,si])
                        published=np.asarray(rain[field].isel(cml_id=rid,sublink_id=rsi,time_5min=slice(r0,r1)).values)
                    if field.startswith('R_') or field.endswith('_db'):
                        score=rain_metrics(recon,published,atol,rtol)
                        matched=score['numeric_pass']
                    else:
                        score=flag_metrics(recon,published)
                        matched=score['pass']
                    score['pass']=matched;per[field]=score
                    all_pass=all_pass and matched
                    if field=='R_sublink_5min':
                        combined_model.extend(np.asarray(recon,float).tolist())
                        combined_ref.extend(np.asarray(published,float).tolist())
                        pd.DataFrame({
                            'timestamp_utc':target_t.astype(str),'cml_id':cml_id,
                            'sublink_id':sl_name,'recomputed_mm_h':recon,'published_mm_h':published,
                            'difference_mm_h':recon-published,
                        }).to_csv(output_dir/f'comparison_link_{index}_sublink_{si}.csv',index=False)
                result['sublinks'].append(per)
            if 'R_mean_5min' in rain:
                recomputed=computed['R_mean_5min'][selected]
                published=rain['R_mean_5min'].isel(cml_id=rid,time_5min=slice(r0,r1)).values
                score=rain_metrics(recomputed,published,atol,rtol)
                score['pass']=score['numeric_pass'];result['results']['R_mean_5min']=score
                all_pass=all_pass and score['pass']
            result['pass']=bool(all_pass)
            combined=rain_metrics(combined_model,combined_ref,atol,rtol)
            result['combined_R_sublink_5min']=combined
            if make_plot and combined['n_shared_finite']:
                comparison_plot(combined_model,combined_ref,
                                output_dir/f'comparison_link_{index}.png',
                                f'CML {cml_id} (two directions)',combined['rmse_mm_h'])
            print(f"LINK {index}: {'PASS' if all_pass else 'DIFFERENCES'} "
                  f"shared={combined['n_shared_finite']}, mask_mismatch={combined['validity_mismatches']}, "
                  f"RMSE={combined['rmse_mm_h']}, out_of_tol={combined['out_of_tolerance']}",flush=True)
            report['links'].append(result)
    report['all_pass']=bool(report['links'] and all(x['pass'] for x in report['links']))
    (output_dir/'comparison_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return report
