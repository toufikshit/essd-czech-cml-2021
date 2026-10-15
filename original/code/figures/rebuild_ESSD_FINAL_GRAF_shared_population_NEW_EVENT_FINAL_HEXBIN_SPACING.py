#!/usr/bin/env python3
"""Recreate ESSD final-Graf figures (April--October 2021).

Input rainfall MUST be an actual Graf-based product; the old
CZ_CML_RAINRATE_5MIN_20200101_20211231.nc and *_RAINRATE_NEARBY_* are rejected.
No rainfall retrieval is rerun and no missing data are replaced by zeros.

Plots:
  figures/FINAL_GRAF_hexbin_gauge_202104_202110.png/.pdf
  figures/FINAL_GRAF_hexbin_gauge_WETONLY_202104_202110.png/.pdf
  figures/FINAL_GRAF_confusion_GAUGE_same_population_202104_202110.png/.pdf
  figures/FINAL_GRAF_event_<different selected event dates>.png/.pdf
  figures/FINAL_GRAF_validation_summary.json

Gauge validation uses the nearest CHMI gauge to the CML midpoint within 5 km,
with end-labelled CHMI 10-min accumulation matched to two start-labelled CML
5-min rain rates. The hexbin uses ALL valid gauge intervals (wet and dry), and
also produces a separate gauge-wet-only hexbin at >= 0.1 mm/h.

Rain detection is calculated independently from PREPROCESSED 1-min signal
using the ORIGINAL project's cml_graf_rstd_hess_radar_gauge.py routines.
BOTH the gauge confusion matrix and the gauge hexbin use the IDENTICAL
CML, sublink, gauge, time, rainfall, and complete Graf-validity mask.
Reference is CHMI gauge (not radar); confusion N MUST equal hexbin N.
A second event is selected from April--October 2021, excluding 23--24 June.

Radar RR is stored as 5-min accumulation in mm: multiply by 12 to get mm/h.
No native RADKLIM field is required for this workflow.

Runs where the final Graf product or the original Graf detector is absent MUST
fail rather than silently falling back to historical or NEARBY products.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import warnings
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator, FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd
import xarray as xr

DEFAULT_ROOT = Path('/home/jovyan/ESSD_publication_dataset')
START, END = pd.Timestamp('2021-04-01'), pd.Timestamp('2021-10-31 23:59:59')
MAX_R = 150.0
# Selected publication product: Graf RSTD + ITU-R P.838-3 (2005) + Pastorek WAA.
PREFERRED_GRAF_RAIN = (
    Path('cmlonly_graf_itu2005_pastorek') /
    'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_20200101_20211231.nc'
)


def key(value):
    """Preserve structural CML IDs (e.g. 85032_85231_3044).

    Python float() accepts underscores inside numbers and would silently turn
    85032_85231_3044 into 85032852313044. Never float-convert string IDs.
    Only normalize decimal .0 decorations added when an integer is represented
    as a floating-point numeric coordinate or exported as a numeric-looking ID.
    """
    if isinstance(value, bytes):
        value=value.decode('utf-8','replace')
    if isinstance(value, str):
        result=value.strip()
        if re.fullmatch(r'[+-]?[0-9]+\.0+',result):
            return result.split('.')[0]
        return result
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)) and np.isfinite(value):
        if float(value).is_integer():
            return str(int(value))
    return str(value).strip()


def ids(ds, name='cml_id'):
    return [key(v) for v in ds[name].values]


def positions(ds, name='cml_id'):
    out = dict((v, i) for i, v in enumerate(ids(ds, name)))
    if len(out) != ds.sizes[name]:
        raise ValueError(f'Duplicate {name} values in dataset')
    return out


def original_graf(root):
    p = root/'cml_graf_rstd_hess_radar_gauge.py'
    if not p.is_file():
        raise FileNotFoundError(f'Original Graf RSTD detection implementation is needed: {p}')
    spec = importlib.util.spec_from_file_location('essd_exact_graf', p)
    g = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g)
    required = ['detect_trsl_source', 'graf_rstd_matrix', 'aggregate_to_reference',
                'build_reference_intervals', 'prepare_gauge_mapping']
    absent = [v for v in required if not hasattr(g, v)]
    if absent:
        raise RuntimeError(f'Original Graf implementation lacks expected functions: {absent}')
    return g


def choose_graf_rain(root: Path, provided: Path|None):
    if provided:
        p = provided.expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(p)
        choices = [p]
    else:
        # Select the exact publication product, not the other Graf comparison
        # product. No directory-wide guessing or silent fallback.
        p = (root / PREFERRED_GRAF_RAIN).resolve()
        if not p.is_file():
            raise FileNotFoundError(
                f'Expected final ITU2005/Pastorek Graf rainfall product: {p}\n'
                'Pass --graf-rain-nc with the correct final product to override.'
            )
        choices = [p]
    valid = []
    for p in choices:
        n = p.name.upper()
        if ('GRAF' not in n or 'RAIN' not in n or 'NEARBY' in n
                or n == 'CZ_CML_RAINRATE_5MIN_20200101_20211231.NC'):
            raise ValueError(f'REJECTED non-final-Graf rainfall input: {p}')
        with xr.open_dataset(p, decode_times=False) as ds:
            candidate = [v for v in ('R_sublink_5min', 'rain_rate', 'rainrate', 'R_5min')
                         if v in ds and 'cml_id' in ds[v].dims and 'sublink_id' in ds[v].dims]
            if candidate:
                valid.append((p, ds.sizes.get('cml_id', 0), candidate))
    if not valid:
        raise ValueError(f'No directional rain-rate variable in Graf candidates: {choices}')
    # Exactly one selected candidate by construction.
    p,n,variables = valid[0]
    print(f'FINAL GRAF RAIN SOURCE: {p}\n  cml_id count: {n}; rainfall vars: {variables}', flush=True)
    if n < 1800:
        warnings.warn(f'Graf file has only {n} physical CMLs: manuscript expects 2219; confirm selected release.')
    return p


def select_rain_rate(rain):
    candidates = [v for v in ('R_sublink_5min', 'rain_rate', 'rainrate', 'R_5min')
                  if v in rain and 'sublink_id' in rain[v].dims]
    if not candidates:
        raise ValueError(f'Graf rain file lacks a directional rain-rate variable: {list(rain.data_vars)}')
    da = rain[candidates[0]]
    for dim, preferred in [('baseline_method', 'pycml_constant'),
                           ('kr_model', 'ITU_2005'), ('waa_method', 'pastorek_2021'),
                           ('wet_dry_method', 'graf'), ('wetdry_method', 'graf'),
                           ('wetdry', 'graf')]:
        if dim not in da.dims:
            continue
        val = ids(da, dim)
        match = next((i for i,v in enumerate(val) if v.lower() == preferred.lower()), None)
        if match is None:
            if len(val) == 1:
                match = 0
            else:
                raise ValueError(f'Multiple {dim} options in final Graf product: {val}; cannot select automatically')
        da = da.isel({dim: match})
    expected = {'time_5min','sublink_id','cml_id'}
    if set(da.dims) != expected:
        raise ValueError(f'Expected directional (cml_id,sublink_id,time_5min) rain rate, got {da.dims}')
    units = str(da.attrs.get('units', ''))
    if 'mm' in units and 'min' in units and 'h' not in units:
        raise ValueError(f'Rain-rate variable appears to be an amount, not mm/h: units={units}')
    return da


def save(fig, stem):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix('.png'), dpi=450, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(stem.with_suffix('.pdf'), bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print('SAVED',stem.with_suffix('.png'), 'and PDF', flush=True)


def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.labelsize':11,'axes.spines.right':False,
                         'axes.spines.top':False,'xtick.direction':'out',
                         'ytick.direction':'out','pdf.fonttype':42,'ps.fonttype':42,
                         'axes.linewidth':.8})


def compute_stats(x,y):
    if len(x)==0:
        return {'N':0, 'RMSE':None, 'MBE':None, 'PCC':None}
    x=np.asarray(x,dtype=float); y=np.asarray(y,dtype=float)
    diff=y-x
    return {'N':int(len(x)), 'RMSE':float(np.sqrt(np.mean(diff*diff))),
            'MBE':float(np.mean(diff)),
            'PCC':float(np.corrcoef(x,y)[0,1]) if len(x)>=2 and x.std()>0 and y.std()>0 else None}


def hexbin(x,y,stem):
    """Gauge-validation hexbin with the approved spacing style.

    Keep the original 0--150 hexbin grid and colormap, and create the desired
    visual spacing by moving the left and bottom axes outward.
    """
    style(); fig,ax=plt.subplots(figsize=(6.9,6.5))
    h=ax.hexbin(x,y,gridsize=40,extent=(0,MAX_R,0,MAX_R),mincnt=1,
                norm=LogNorm(vmin=1),cmap='jet',linewidths=0)
    ax.plot([0,MAX_R],[0,MAX_R],'--',lw=1,color='0.26',zorder=3)
    ax.set_xlim(0,MAX_R);ax.set_ylim(0,MAX_R)
    ax.set_aspect('equal',adjustable='box')
    ax.set_xticks([0,30,60,90,120,150]);ax.set_yticks([0,30,60,90,120,150])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='both',direction='out',right=False,top=False)
    ax.spines['right'].set_visible(False);ax.spines['top'].set_visible(False)
    ax.spines['left'].set_position(('outward',12))
    ax.spines['bottom'].set_position(('outward',12))
    ax.set_xlabel(r'RR$_{gauge}$ [mm h$^{-1}$]')
    ax.set_ylabel(r'RR$_{CML}$ [mm h$^{-1}$]')
    stats=compute_stats(x,y)
    label = (f'N={stats["N"]:,}\n'
             f'RMSE={stats["RMSE"]:.2f} mm h$^{{-1}}$\n'
             f'MBE={stats["MBE"]:.2f} mm h$^{{-1}}$\n'
             f'PCC={stats["PCC"]:.3f}') if stats['PCC'] is not None else f'N={stats["N"]:,}'
    ax.text(.025,.972,label,
            transform=ax.transAxes,ha='left',va='top',fontsize=10,
            bbox={'facecolor':'white','edgecolor':'.55','alpha':.88,'boxstyle':'round,pad=.35'})
    cb=fig.colorbar(h,ax=ax,shrink=.98,pad=.035,fraction=.05)
    cb.set_label('Count per hexagon');cb.ax.tick_params(direction='out')
    save(fig,stem)
    return stats

def plot_confusion(tp,fp,fn,tn,stem):
    style()
    fig,ax=plt.subplots(figsize=(8.4,8))
    colors=[[(151/255,185/255,206/255),(227/255,233/255,237/255)],
            [(75/255,147/255,195/255),(11/255,63/255,133/255)]]
    cells=[[('Hit',tp,tp/(tp+fn) if tp+fn else np.nan,'gauge wet'),
            ('False alarm',fp,fp/(fp+tn) if fp+tn else np.nan,'gauge dry')],
           [('Miss',fn,fn/(tp+fn) if tp+fn else np.nan,'gauge wet'),
            ('Correct negative',tn,tn/(fp+tn) if fp+tn else np.nan,'gauge dry')]]
    for r in range(2):
        for c in range(2):
            title,n,ratio,lab=cells[r][c]
            ax.add_patch(Rectangle((c,1-r),1,1,facecolor=colors[r][c],
                                   edgecolor='white',linewidth=6))
            col='white' if r==1 else 'black'
            ax.text(c+.5,1-r+.76,title,ha='center',va='center',fontsize=16 if len(title)>15 else 18,weight='bold',color=col)
            ax.text(c+.5,1-r+.44,f'{n:,}',ha='center',va='center',fontsize=23,color=col)
            ax.text(c+.5,1-r+.17,f'{100*ratio:.1f}% of {lab}',ha='center',va='center',fontsize=13,color=col)
    ax.set(xlim=(0,2),ylim=(0,2));ax.set_aspect('equal')
    ax.set_xticks([.5,1.5],['Wet','Dry'])
    ax.set_yticks([1.5,.5],['Wet','Dry'])
    ax.tick_params(length=0,labelsize=16,pad=10)
    for sp in ax.spines.values():sp.set_visible(False)
    ax.set_xlabel('CHMI gauge reference',fontsize=21,labelpad=17)
    ax.set_ylabel('CML Graf/RSTD detection',fontsize=21,labelpad=17)
    save(fig,stem)


def numeric(da,*dims):
    arr=np.asarray(da.values,dtype=float)
    if set(da.dims)!=set(dims):
        raise ValueError(f'Expected {dims}, found {da.dims}')
    return np.transpose(arr,[da.dims.index(s) for s in dims])


def station_series(ds,var,midlat,midlon,t0,t1):
    field=ds[var]
    tdim='time'; st_dims=[s for s in field.dims if s!='time']
    if len(st_dims)!=1:
        raise ValueError(f'{var} must have exactly one station dim: {field.dims}')
    sd=st_dims[0];n=ds.sizes[sd]
    lat=next((np.asarray(ds[v].values,float).reshape(-1) for v in ('lat','latitude','station_lat')
              if v in ds and ds[v].size==n),None)
    lon=next((np.asarray(ds[v].values,float).reshape(-1) for v in ('lon','longitude','station_lon')
              if v in ds and ds[v].size==n),None)
    if lat is None or lon is None:raise KeyError(f'Missing station coordinates for {var}')
    d=hav(midlat,midlon,lat,lon)
    d[~np.isfinite(d)]=np.inf
    for i in np.argsort(d)[:min(30,n)]:
        s=field.isel({sd:int(i)}).sel(time=slice(t0,t1))
        if np.isfinite(np.asarray(s.values,float)).sum()>=2:
            return pd.DatetimeIndex(s['time'].values),np.asarray(s.values,float),float(d[i])
    raise ValueError(f'No finite {var} stations for event')


def hav(lat,lon,lat2,lon2):
    a=(np.asarray(lat2)-lat)*np.pi/360
    b=(np.asarray(lon2)-lon)*np.pi/360
    h=np.sin(a)**2+np.cos(lat*np.pi/180)*np.cos(np.asarray(lat2)*np.pi/180)*np.sin(b)**2
    return 12742*np.arcsin(np.sqrt(np.clip(h,0,1)))


def resolve_paths(root, a):
    pp={
        'cml':root/'CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc',
        'gauge':root/'CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc',
        'radar':root/'CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc',
        'temp':root/'CHMI_10MIN_TEMPERATURE_2021.nc',
        'hum':root/'CHMI_10MIN_HUMIDITY_2021.nc',
        'rain':choose_graf_rain(root,a.graf_rain_nc)
    }
    met=Path('/home/jovyan/opendata.chmi.cz/processed')
    for k in ('temp','hum'):
        if not pp[k].is_file():pp[k]=met/pp[k].name
    for k,p in pp.items():
        if not p.is_file():raise FileNotFoundError(f'{k}: {p}')
    return pp


def gauge_validation(cml, rain, gauge, g, out, start, end, threshold=0.1):
    """ONE joint mask for all-period hexbin, wet-only hexbin, confusion.

    As in the project's original same-population analysis, Graf detection uses
    the full 1-minute QC time series to establish each sublink's Q80 threshold.
    Gauge truth and 10min retrieved rainfall then share the same exact native
    interval. No detector interpolation, no radar, and no silent fallbacks.
    """
    rate=select_rain_rate(rain).sel(time_5min=slice(start-pd.Timedelta(minutes=10),end))
    rt=pd.DatetimeIndex(rate.time_5min.values)
    times=pd.DatetimeIndex(gauge.time.values)
    keep=np.flatnonzero((times>=start+pd.Timedelta(minutes=10))&(times<=end))
    gt=times[keep]
    i1=rt.get_indexer(gt-pd.Timedelta(minutes=10))
    i2=rt.get_indexer(gt-pd.Timedelta(minutes=5))
    exact=(i1>=0)&(i2>=0)
    if not exact.any():raise ValueError('No exact 10min gauge / 5min final Graf rainfall matches')
    i1=np.clip(i1,0,len(rt)-1);i2=np.clip(i2,0,len(rt)-1)
    st=pd.DatetimeIndex(cml.time.values)
    imap=g.build_reference_intervals(st,gt,'end',10.)
    if int(imap['expected'])!=10:raise ValueError('Expected ten 1min Graf samples per gauge interval')
    rainpos=positions(rain);cpos=positions(cml)
    # Preserve literal QC IDs for the original helper (no float conversion).
    cids=np.asarray(cml['cml_id'].values).astype(str)
    arg=SimpleNamespace(site_a_lon='site_a_longitude',site_a_lat='site_a_latitude',
                        site_b_lon='site_b_longitude',site_b_lat='site_b_latitude',
                        gauge_lon=None,gauge_lat=None,gauge_id=None,max_gauge_dist_km=5)
    info=g.prepare_gauge_mapping(cml=cml,cml_ids=cids,gauge=gauge,gauge_var='rr',args=arg)
    mapping=info['mapping'];mapping=mapping[(mapping.gauge_position>=0)&(mapping.distance_km<=5)].copy()
    print(f'Gauge midpoint matched CMLs <=5km: {len(mapping):,}',flush=True)
    gv=info['rr'];gdim=info['station_dim'];tdim=info['time_dim']
    sid_c=ids(cml,'sublink_id');sid_r=ids(rate,'sublink_id')
    if len(set(sid_c))!=len(sid_c) or set(sid_c)!=set(sid_r):
        raise ValueError(f'Directional sublinks differ: QC {sid_c}; Graf rainfall {sid_r}')
    ri={sid:i for i,sid in enumerate(sid_r)}
    source=g.detect_trsl_source(cml,None)
    cache={};X=[];Y=[];audit=[]
    TP=FP=FN=TN=0;pre_rain=pre_det=0;rain_membership=0
    for k,(raw_cid,station,gid,dist) in enumerate(mapping[['cml_id','gauge_position','gauge_id','distance_km']].itertuples(index=False,name=None),1):
        cid=key(raw_cid)
        if cid not in rainpos or cid not in cpos:continue
        rain_membership+=1
        station=int(station)
        if station not in cache:
            cache[station]=np.asarray(gv.isel({gdim:station,tdim:keep}).values,float).reshape(-1)*6
        ref=cache[station]
        cp=cpos[cid]
        sig=(cml['tx'].isel(cml_id=cp)-cml['rx'].isel(cml_id=cp)) if source is None else cml[source].isel(cml_id=cp)
        sig=numeric(sig,'time','sublink_id')
        wet,valid,*_=g.graf_rstd_matrix(sig,window=60,quantile=.8,factor=1.12,
                                        threshold_min_periods=60,detection_min_periods=45)
        pred,predvalid=g.aggregate_to_reference(wet,valid,imap,1.,'any')
        if pred.shape!=(len(gt),len(sid_c)):raise ValueError(f'Graf detector shape mismatch: {pred.shape}')
        rainfall=numeric(rate.isel(cml_id=rainpos[cid]),'sublink_id','time_5min')
        base=exact & imap['complete'] & np.isfinite(ref) & (ref>=0) & (ref<=MAX_R)
        for i,sid in enumerate(sid_c):
            j=ri[sid]
            a=rainfall[j,i1];b=rainfall[j,i2]
            y=(a+b)*.5
            rainvalid=base & np.isfinite(a)&np.isfinite(b)&np.isfinite(y)&(y>=0)&(y<=MAX_R)
            detvalid=base & predvalid[:,i]
            pre_rain+=int(rainvalid.sum());pre_det+=int(detvalid.sum())
            ok=rainvalid&detvalid
            if not ok.any():continue
            x=ref[ok].astype(np.float32); yy=y[ok].astype(np.float32)
            truth=(x>=threshold)
            predicted=np.asarray(pred[ok,i],bool)
            tp=int(np.count_nonzero(predicted&truth))
            fp=int(np.count_nonzero(predicted&~truth))
            fn=int(np.count_nonzero(~predicted&truth))
            tn=int(np.count_nonzero(~predicted&~truth))
            TP+=tp;FP+=fp;FN+=fn;TN+=tn
            X.append(x);Y.append(yy)
            audit.append({'cml_id':cid,'sublink_id':sid,'gauge_id':key(gid),
                          'distance_km':float(dist),'N':len(x),
                          'TP':tp,'FP':fp,'FN':fn,'TN':tn})
        if k%150==0:print(f'Gauge validated {k}/{len(mapping)} midpoint-matched physical CMLs',flush=True)
    if not X:raise RuntimeError('No joint valid final Graf / Graf RSTD / gauge samples')
    x=np.concatenate(X);y=np.concatenate(Y)
    counts={'TP':TP,'FP':FP,'FN':FN,'TN':TN}
    n=TP+FP+FN+TN
    if n!=len(x) or n!=sum(d['N'] for d in audit):
        raise AssertionError(f'Population mismatch: confusion={n}, hexbin={len(x)}')
    wet=(x>=threshold)
    if int(wet.sum())!=(TP+FN):raise AssertionError('Wet hexbin count != gauge-wet confusion total')
    stats=hexbin(x,y,out/'FINAL_GRAF_hexbin_gauge_202104_202110')
    if not wet.any():raise RuntimeError('No gauge wet intervals in validation')
    wet_stats=hexbin(x[wet],y[wet],out/'FINAL_GRAF_hexbin_gauge_WETONLY_202104_202110')
    plot_confusion(TP,FP,FN,TN,out/'FINAL_GRAF_confusion_GAUGE_same_population_202104_202110')
    den=math.sqrt(float(TP+FP)*float(TP+FN)*float(TN+FP)*float(TN+FN))
    skill={'POD':TP/(TP+FN) if TP+FN else None,
           'FAR':FP/(TP+FP) if TP+FP else None,
           'CSI':TP/(TP+FP+FN) if TP+FP+FN else None,
           'MCC':(TP*TN-FP*FN)/den if den else None}
    result={'all_period':stats,'gauge_wet_ge_0_1':wet_stats,
            'confusion':{**counts,'N':n,**skill},
            'population':{'assigned_midpoint_CMLs':len(mapping),
                          'matched_CMLs_in_final_rainfall':rain_membership,
                          'contributing_physical_CMLs':len({d['cml_id'] for d in audit}),
                          'contributing_directional_sublinks':len(audit),
                          'distinct_gauges':len({d['gauge_id'] for d in audit}),
                          'rainvalid_pre_intersection':pre_rain,
                          'detvalid_pre_intersection':pre_det},
            'reference':'nearest CHMI gauge to CML midpoint <=5km; end-labelled 10min *6; wet>=0.1',
            'rainfall':'two complete final Graf ITU2005/Pastorek 5min bins at t-10 and t-5',
            'detection':'original Graf RSTD, full 1min Q80 * 1.12; all 10min valid',
            'validation_mask':'identical (link,gauge,sublink,10min timestamp) for confusion and all-period hexbin',
            'in_plot_range':'0..150 mm/h inclusive',
            'equality_check':'PASS'}
    pd.DataFrame(audit).to_csv(out/'FINAL_GRAF_gauge_shared_population_audit.csv',index=False)
    print('='*70)
    print(f'FINAL GRAF HEXBIN N={len(x):,}, GAUGE CONFUSION N={n:,}; EQUALITY: PASS')
    print(f'GAUGE WET N={int(wet.sum()):,}; CHECK TP+FN={TP+FN:,}')
    print('All-period:',stats,'Wet-only:',wet_stats,flush=True)
    print('Gauge confusion:',{**counts,**skill},flush=True)
    print('='*70)
    return result


def event_plot(cml, rain, radar, gauge, temp, hum, out, t0, t1,
               preferred=None, require_preferred=False):
    """Four-panel FINAL Graf event for a NEW date, never June 23--24.

    An explicit --event-cml-id requires that ID to pass all checks.
    For automatic dates, select the strongest valid event link instead.
    """
    rate=select_rain_rate(rain).sel(time_5min=slice(t0,t1))
    rainpos=positions(rain);radpos=positions(radar);cpos=positions(cml)
    if 'frequency' not in cml:raise KeyError('CML metadata missing frequency')
    f=numeric(cml['frequency'],'cml_id','sublink_id')
    if 'MHz' in str(cml['frequency'].attrs.get('units','')) or np.nanmedian(f)>1000:f=f/1000
    elinks=np.any(((f>=71)&(f<=76))|((f>=81)&(f<=86)),axis=1)
    t_r=pd.DatetimeIndex(rate.time_5min.values)
    t_ra=pd.DatetimeIndex(radar.time.values)
    rad_sel=np.flatnonzero((t_ra>=t0)&(t_ra<=t1))
    if not len(rad_sel):raise ValueError('No radar timestamps in event')
    t_g=pd.DatetimeIndex(gauge.time.values)
    gauge_sel=np.flatnonzero((t_g>=t0)&(t_g<=t1))
    if not len(gauge_sel):raise ValueError('No gauge timestamps in event')
    st_dim=next(d for d in gauge['rr'].dims if d!='time')
    st_n=gauge.sizes[st_dim]
    lat=next(np.asarray(gauge[v].values,float).reshape(-1) for v in ('lat','latitude','station_lat') if v in gauge and gauge[v].size==st_n)
    lon=next(np.asarray(gauge[v].values,float).reshape(-1) for v in ('lon','longitude','station_lon') if v in gauge and gauge[v].size==st_n)
    a,b=np.asarray(cml['site_a_latitude'].values,float),np.asarray(cml['site_b_latitude'].values,float)
    c,d=np.asarray(cml['site_a_longitude'].values,float),np.asarray(cml['site_b_longitude'].values,float)
    mids_lat=(a+b)/2;mids_lon=(c+d)/2

    preferred=key(preferred) if preferred is not None else None
    if t0 <= pd.Timestamp('2021-06-24 03:00:00') and t1 >= pd.Timestamp('2021-06-23 09:00:00'):
        raise ValueError('Different event requested: exclude the original 23--24 June event')
    if require_preferred and preferred not in cpos:
        raise ValueError(f'Explicit event CML ID {preferred!r} not found in preprocessed CML file')
    rejection={};counts={}
    candidates=[]
    st=pd.DatetimeIndex(cml.time.values)
    sigsel=np.flatnonzero((st>=t0)&(st<=t1))
    for k,cid in enumerate(ids(cml)):
        if require_preferred and cid!=preferred:continue
        if not elinks[k]:
            if cid==preferred:rejection[cid]='not E-band according to QC frequency metadata'
            continue
        reason=None
        if cid not in rainpos:reason='CML ID absent from FINAL Graf rainfall product'
        elif cid not in radpos:reason='CML ID absent from path-radar product'
        if reason:
            counts[reason]=counts.get(reason,0)+1
            if cid==preferred:rejection[cid]=reason
            continue
        distances=hav(float(mids_lat[k]),float(mids_lon[k]),lat,lon)
        distances[~np.isfinite(distances)]=np.inf
        gi=int(distances.argmin())
        if distances[gi]>5:
            reason=f'nearest midpoint gauge {float(distances[gi]):.2f} km away (>5 km)'
        else:
            rainvals=numeric(rate.isel(cml_id=rainpos[cid]),'sublink_id','time_5min')
            valid=np.isfinite(rainvals)
            # Sum only valid sublinks and leave samples with zero valid directions as NaN.
            cnt=valid.sum(axis=0)
            cr=np.full(rainvals.shape[1],np.nan)
            np.divide(np.where(valid,rainvals,0).sum(axis=0),cnt,
                      out=cr,where=cnt>0)
            rv=np.asarray(radar['rr'].isel(cml_id=radpos[cid],time=rad_sel).values,float).reshape(-1)*12
            gv=np.asarray(gauge['rr'].isel({st_dim:gi,'time':gauge_sel}).values,float).reshape(-1)*6
            if np.isfinite(cr).sum()<8:reason='fewer than 8 finite final Graf event rain-rate samples'
            elif np.isfinite(rv).sum()<8:reason='fewer than 8 finite radar event samples'
            elif np.isfinite(gv).sum()<5:reason='fewer than 5 finite gauge event samples'
            elif np.nanmax(cr)<=1:reason='final Graf event rain rate does not exceed 1 mm/h'
            elif np.nanmax(rv)<=1:reason='radar event rain rate does not exceed 1 mm/h'
            elif np.nanmax(gv)<=1:reason='gauge event rain rate does not exceed 1 mm/h'
            else:
                signal=cml['trsl'].isel(cml_id=k,time=sigsel) if 'trsl' in cml else (cml['tx']-cml['rx']).isel(cml_id=k,time=sigsel)
                if np.isfinite(numeric(signal,'time','sublink_id')).sum()<30:
                    reason='fewer than 30 finite raw TRSL samples during event'
        if reason:
            counts[reason]=counts.get(reason,0)+1
            if cid==preferred:rejection[cid]=reason
            continue
        score=(min(float(np.nanmax(cr)),float(np.nanmax(rv)),float(np.nanmax(gv))),-float(distances[gi]))
        candidates.append((score,cid,k,gi,float(distances[gi]),cr,rv,gv))

    if preferred and preferred in rejection:
        print(f'EARLIER EVENT CML {preferred!r} REJECTED: {rejection[preferred]}',flush=True)
    if not candidates:
        print('Event rejection summary:',json.dumps(counts,indent=2),flush=True)
        raise RuntimeError('No eligible FINAL Graf E-band event link in this window. '
                           'The old event CML is not silently reused. '
                           'See the rejection reason above.')
    candidates.sort(key=lambda z:(int(preferred is not None and z[1]==preferred),z[0]),reverse=True)
    for score,cid,k,gi,gkm,cr,rv,gv in candidates:
        try:
            midlat=float(mids_lat[k]);midlon=float(mids_lon[k])
            tt,tv,tdist=station_series(temp,'T',midlat,midlon,t0,t1)
            ht,hv,hdist=station_series(hum,'H',midlat,midlon,t0,t1)
            break
        except ValueError as exc:
            if cid==preferred:print(f'EARLIER EVENT CML {cid!r} has no suitable meteorological station: {exc}',flush=True)
    else:
        raise RuntimeError('No valid event CML has both temperature and humidity station observations')
    if preferred and cid!=preferred:
        raise AssertionError(f'Explicit CML ID {preferred} was not selected; got {cid}')
    print(f'NEW EVENT CML SELECTED: {cid}',flush=True)

    sigda=cml['trsl'].isel(cml_id=k,time=sigsel) if 'trsl' in cml else (cml['tx']-cml['rx']).isel(cml_id=k,time=sigsel)
    signals=numeric(sigda,'time','sublink_id')
    style();fig,axes=plt.subplots(4,1,sharex=True,figsize=(8.6,7.8),
                                  gridspec_kw={'height_ratios':[1.27,1.42,1.,1.]},
                                  constrained_layout=True)
    for m,ax in enumerate(axes):
        ax.text(0,1.02,f'({"abcd"[m]})',transform=ax.transAxes,va='bottom',fontsize=10.5)
        ax.tick_params(which='both',direction='out',right=False,top=False)
        ax.spines['bottom'].set_position(('outward',4))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    for sl in range(signals.shape[1]):
        axes[0].plot(st[sigsel],signals[:,sl],label=f'Sublink {sl+1}',lw=1.1,
                     color=['#233B78','#B65A3A'][sl%2])
    axes[0].set_ylabel('TRSL [dB]');axes[0].legend(frameon=False,loc='upper right',ncol=2,fontsize=8.5)
    axes[1].step(t_r,cr,where='post',label='CML (Graf)',color='#1F4E79',lw=1.6)
    axes[1].step(t_ra[rad_sel],rv,where='post',label='Radar',color='#C05E25',lw=1.45,ls='--')
    axes[1].bar(t_g[gauge_sel]-pd.Timedelta(minutes=10),gv,width=10/(24*60),align='edge',
                label='Gauge',color='#5A9C83',alpha=.36)
    axes[1].set_ylabel('Rain rate [mm h$^{-1}$]')
    axes[1].legend(frameon=False,ncol=3,loc='upper right',fontsize=8.5)
    axes[1].set_ylim(bottom=0)
    axes[2].plot(tt,tv,color='#A64940',lw=1.3)
    axes[2].set_ylabel('Temperature [°C]')
    axes[3].plot(ht,hv,color='#40727D',lw=1.3)
    axes[3].set_ylabel('Humidity [%]')
    if np.isfinite(hv).any():
        # Actual humidity range; never enlarge a 65--75% series to 104%.
        lower=max(0.,float(np.nanmin(hv))-3.)
        upper=min(100.,float(np.nanmax(hv))+3.)
        if lower>=upper:lower=max(0.,upper-1.)
        axes[3].set_ylim(lower,upper)
    # Native x limits end on outward MAJOR ticks, including exact first/last time.
    major=list(pd.date_range(t0.ceil('3h'),t1.floor('3h'),freq='3h'))
    # Remove interior labels nearer than two hours to an endpoint, so the
    # last major tick remains legible (e.g. 00:00 vs 01:00).
    major=sorted(set([t0,*[x for x in major if x-t0>=pd.Timedelta(hours=2)
                                               and t1-x>=pd.Timedelta(hours=2)],t1]))
    axes[-1].set_xticks(mdates.date2num(major))
    axes[-1].xaxis.set_major_formatter(FuncFormatter(
        lambda x,_:mdates.num2date(x).strftime('%d %b %H:%M')
        if mdates.num2date(x).hour==0 or abs(x-mdates.date2num(t0))<1e-8
        else mdates.num2date(x).strftime('%H:%M')))
    axes[-1].set_xlabel('Time [UTC]')
    for ax in axes:
        ax.set_xlim(t0,t1);ax.margins(x=0)
        # Use outward major ticks at BOTH actual y-spine endpoints.
        lo,hi=ax.get_ylim()
        locator=MaxNLocator(nbins=5)
        ticks=locator.tick_values(lo,hi)
        if ax is axes[1]:ticks=ticks[ticks>=0]
        if ax is axes[3]:
            ticks=ticks[(ticks>=0)&(ticks<=100)]
            # A fully saturated observation must remain on the plot.
            if np.isfinite(hv).any() and float(np.nanmax(hv))>ticks[-1]:
                ticks=np.r_[ticks,100.]
        if len(ticks)>1:
            ax.set_ylim(float(ticks[0]),float(ticks[-1]));ax.set_yticks(ticks)
        ax.tick_params(which='both',direction='out',top=False,right=False)
    save(fig,out/f'FINAL_GRAF_event_{t0:%Y%m%d}_{t1:%Y%m%d}')
    result={'cml_id':cid,'midpoint_to_gauge_km':gkm,'temp_station_km':tdist,
            'humidity_station_km':hdist,'event_start':str(t0),'event_end':str(t1),
            'rain_rate_source':'Final Graf ITU2005 Pastorek product',
            'explicit_event_id_requested':preferred,
            'explicit_event_id_rejection_reason':rejection.get(preferred),
            'selection':'different event, strongest eligible E-band Graf/radar/gauge CML'}
    (out/'FINAL_GRAF_event_metadata.json').write_text(json.dumps(result,indent=2)+'\n')
    print('FINAL GRAF event:',result,flush=True)
    return result

def event_windows(gauge,start,end):
    """Choose another real 2021 storm, with 13 July preferred (paper Fig. 7).

    If no E-band link has all five reference streams for that event, try
    strongest alternative CHMI gauge rain days (without using old 23--24 June).
    Ranking identifies candidates only; event_plot verifies all streams.
    """
    a=pd.Timestamp('2021-07-13 12:00');b=pd.Timestamp('2021-07-14 01:00')
    windows=[]
    if a>=start and b<=end:windows.append((a,b,'2021-07-13 event (Figure 7)'))
    da=gauge['rr']
    time=pd.DatetimeIndex(gauge.time.values)
    mask=np.flatnonzero((time>=start)&(time<=end))
    # Load the ~7 month gauge field once and find strong rainy intervals.
    values=np.asarray(da.isel(time=mask).values,dtype=np.float32)
    timeaxis=da.dims.index('time')
    if timeaxis!=0:values=np.moveaxis(values,timeaxis,0)
    with np.errstate(all='ignore'):
        peak=np.where(np.isfinite(values),values,-np.inf).max(axis=1)
    order=np.argsort(peak)[::-1]
    used={pd.Timestamp('2021-07-13').date()}
    for j in order:
        if len(windows)>=17:break
        moment=time[mask[j]]
        if moment.date() in used or moment.date() in (pd.Timestamp('2021-06-23').date(),pd.Timestamp('2021-06-24').date()):continue
        used.add(moment.date())
        t0=moment-pd.Timedelta(hours=6);t1=moment+pd.Timedelta(hours=7)
        if t0>=start and t1<=end:windows.append((t0,t1,f'CHMI strong-rain candidate {moment}'))
    return windows


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    p.add_argument('--graf-rain-nc',type=Path,help='Override the default final ITU2005/Pastorek Graf rainfall NetCDF')
    p.add_argument('--start',default=str(START))
    p.add_argument('--end',default=str(END))
    p.add_argument('--event-start',default=None,help='Explicit NEW event start UTC; defaults to 13 July candidate, then other real storms')
    p.add_argument('--event-end',default=None,help='Explicit NEW event end UTC (requires --event-start)')
    p.add_argument('--event-cml-id',default=None)
    p.add_argument('--skip-confusion',action='store_true',help='Ignored: shared-population confusion is required with hexbin')
    p.add_argument('--skip-event',action='store_true')
    p.add_argument('--only-event',action='store_true',help='Create only the event plot; preserve completed hexbin and confusion files')
    p.add_argument('--only-check-sources',action='store_true')
    a=p.parse_args()
    root=a.root.expanduser().resolve()
    paths=resolve_paths(root,a)
    if a.only_check_sources:
        for k,v in paths.items():print(f'{k}: {v}')
        return
    if a.only_event and a.skip_event:
        raise ValueError('--only-event and --skip-event cannot both be set')
    if a.skip_confusion and not a.only_event:
        raise ValueError('The requested equal-population gauge confusion matrix is mandatory; do not use --skip-confusion')
    g=None if a.only_event else original_graf(root)
    out=root/'figures';out.mkdir(exist_ok=True)
    for k,v in paths.items():print(f'{k}: {v}',flush=True)
    result={'inputs':{k:str(v) for k,v in paths.items()},'analysis_period':[str(a.start),str(a.end)]}
    with xr.open_dataset(paths['cml']) as cml, xr.open_dataset(paths['rain']) as rain,\
         xr.open_dataset(paths['radar']) as radar, xr.open_dataset(paths['gauge']) as gauge,\
         xr.open_dataset(paths['temp']) as temp, xr.open_dataset(paths['hum']) as hum:
        if not a.only_event:
            result['shared_gauge_validation']=gauge_validation(cml,rain,gauge,g,out,pd.Timestamp(a.start),pd.Timestamp(a.end))
        if not a.skip_event:
            if (a.event_start is None)!=(a.event_end is None):
                raise ValueError('--event-start and --event-end must be supplied together')
            windows=[(pd.Timestamp(a.event_start),pd.Timestamp(a.event_end),'explicit user event')] if a.event_start else event_windows(gauge,pd.Timestamp(a.start),pd.Timestamp(a.end))
            last=None
            for t0,t1,reason in windows:
                print(f'TRY DIFFERENT EVENT: {t0} -> {t1}; source: {reason}',flush=True)
                try:
                    result['event']=event_plot(cml,rain,radar,gauge,temp,hum,out,t0,t1,
                                               preferred=a.event_cml_id,
                                               require_preferred=(a.event_cml_id is not None))
                    result['event']['candidate_origin']=reason
                    print(f'SELECTED NEW EVENT: {t0} to {t1}; physical CML {result["event"]["cml_id"]}',flush=True)
                    break
                except RuntimeError as exc:
                    last=exc
                    print(f'NO ELIGIBLE LINK for {reason}: {exc}',flush=True)
            else:
                raise RuntimeError(f'No NEW event was found with suitable E-band/Graf/gauge/radar and meteorological coverage. Last error: {last}')

    if not a.only_event:
        (out/'FINAL_GRAF_validation_summary.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    else:
        print('Event-only run: existing hexbins and confusion matrix were not recalculated or modified.',flush=True)
    print('SUCCESS: final Graf script used only:',paths['rain'])


if __name__=='__main__':main()
