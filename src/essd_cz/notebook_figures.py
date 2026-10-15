"""Small non-scientific figure adapters. Never edit original retrieval algorithms.

Plots map coordinates/metadata and summarizes CSVs produced by the ORIGINAL
availability implementation. Approximate reconstructions are named reviewer_*.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STATIONS = {'gauge': ('Rain gauges', '#e67e22', '^'),
            'temperature': ('Temperature stations', '#2a9d8f', 'o'),
            'humidity': ('Humidity stations', '#5e60ce', 'o')}

def station_coordinates(ds):
    """Return static latitude/longitude only. No station values or dates changed."""
    def find(choices):
        for name in choices:
            if name in ds:
                vals = np.asarray(ds[name].values, dtype=float).ravel()
                if vals.ndim == 1: return vals
        return None
    lat=find(('lat','latitude','station_lat'))
    lon=find(('lon','longitude','station_lon'))
    if lat is None or lon is None or lat.size != lon.size:
        raise ValueError('No compatible static station coordinates; inspect this NetCDF, and do not invent station locations')
    good=np.isfinite(lat)&np.isfinite(lon)&(lat>=-90)&(lat<=90)&(lon>=-180)&(lon<=180)
    return lon[good],lat[good]

def plot_network(paths, out_dir, dpi=200):
    """Map geometry from QC CML endpoints and CHMI station coordinates.

    This does NOT reproduce original national outline or native RADKLIM footprint.
    Original map generator was not recovered; title and filename disclose status.
    """
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    with xr.open_dataset(paths['qc']) as ds:
        names=('site_a_longitude','site_b_longitude','site_a_latitude','site_b_latitude')
        if not all(v in ds for v in names):raise KeyError('CML geometry variables missing: '+str([v for v in names if v not in ds]))
        lon0,lon1,lat0,lat1=[np.asarray(ds[v].values,float).ravel() for v in names]
    good=np.isfinite(lon0)&np.isfinite(lon1)&np.isfinite(lat0)&np.isfinite(lat1)
    fig,ax=plt.subplots(figsize=(9,6))
    ax.plot(np.stack([lon0[good],lon1[good]]),np.stack([lat0[good],lat1[good]]),color='0.24',lw=.42,alpha=.4)
    counts={'cml_paths':int(good.sum())}
    for role,(label,color,marker) in STATIONS.items():
        with xr.open_dataset(paths[role]) as ds:
            lon,lat=station_coordinates(ds)
        counts[role]=len(lon)
        ax.scatter(lon,lat,s=16 if role=='gauge' else 12,marker=marker,color=color,alpha=.83,label=f'{label} ({len(lon)})',zorder=3)
    ax.set(xlabel='Longitude (°E)',ylabel='Latitude (°N)',title='Review reconstruction: CML paths and CHMI observation sites')
    ax.grid(alpha=.15);ax.legend(fontsize=8,frameon=False,loc='best')
    fig.tight_layout()
    target=out_dir/'reviewer_f01_network_geometry.png'
    fig.savefig(target,dpi=dpi);plt.close(fig)
    return target,counts

def frequency_length_pairs(ds):
    """Flatten *directional* frequency–length pairs with explicit dimensional alignment."""
    fkey='frequency_GHz' if 'frequency_GHz' in ds else 'frequency'
    lkey='length_km' if 'length_km' in ds else 'length'
    if fkey not in ds or lkey not in ds:raise KeyError('Need frequency and length metadata')
    f=ds[fkey];l=ds[lkey]
    f,l=xr.broadcast(f,l)
    freq=np.asarray(f.values,float).reshape(-1)
    km=np.asarray(l.values,float).reshape(-1)
    unit=str(ds[fkey].attrs.get('units','')).lower()
    if fkey=='frequency_GHz' or unit=='ghz': pass
    elif unit in ('mhz','megahertz',''):freq=freq/1000.
    else:raise ValueError('Unknown carrier-frequency units: '+unit)
    lu=str(ds[lkey].attrs.get('units','')).lower()
    if lkey=='length_km' or lu in ('km','kilometer','kilometre',''):pass
    elif lu in ('m','meter','metre'):km=km/1000.
    else:raise ValueError('Unknown path-length units: '+lu)
    good=np.isfinite(freq)&np.isfinite(km)&(freq>0)&(km>0)
    return freq[good],km[good]

def plot_characteristics(paths,out_dir,dpi=200):
    """Reconstruct Figure 3 metadata distributions (not exact original styling)."""
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    with xr.open_dataset(paths['qc']) as ds:f,km=frequency_length_pairs(ds)
    fig,axes=plt.subplots(1,3,figsize=(13,3.8))
    axes[0].hist(f,bins=np.arange(0,max(91,int(np.nanmax(f))+2)),color='#e76f51')
    axes[1].hist(km,bins=np.arange(0,max(31,int(np.nanmax(km))+2)),color='#e76f51')
    hb=axes[2].hexbin(f,km,gridsize=30,mincnt=1,cmap='YlGnBu')
    axes[0].set(xlabel='Carrier frequency (GHz)',ylabel='Directional sublinks')
    axes[1].set(xlabel='Path length (km)',ylabel='Directional sublinks')
    axes[2].set(xlabel='Carrier frequency (GHz)',ylabel='Path length (km)')
    fig.colorbar(hb,ax=axes[2],label='Directional sublinks')
    for ax in axes:ax.grid(alpha=.13)
    fig.tight_layout()
    target=out_dir/'reviewer_f03_cml_characteristics.png'
    fig.savefig(target,dpi=dpi);plt.close(fig)
    return target,{'directional_pairs':int(len(f)),'eband_pairs':int(( ((f>=71)&(f<=76)) | ((f>=81)&(f<=86))).sum())}

def monthly_from_native_csv(native_dir,out_dir,dpi=200,basename='availability_native_timestep_percentage'):
    """Monthly arithmetic mean of native-timestep network-wide availability percentages.

    NOT the same as pooling station-time across sensors. Frequencies stay separate.
    """
    native_dir=Path(native_dir);out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    frames=[]
    for role in ('cml','radar','gauge','temperature','humidity'):
        p=native_dir/f'{basename}_{role}.csv'
        if not p.is_file():raise FileNotFoundError('Run the original availability script first: '+str(p))
        x=pd.read_csv(p,usecols=['time','availability_percent'],parse_dates=['time'])
        if x.empty:raise ValueError('Empty native availability CSV: '+str(p))
        if x['time'].duplicated().any():raise ValueError('Duplicate timestamps: '+str(p))
        x['month']=x['time'].dt.to_period('M').astype(str)
        y=x.groupby('month',sort=True)['availability_percent'].mean().rename(role)
        frames.append(y)
    monthly=pd.concat(frames,axis=1).sort_index()
    csv_path=out_dir/'reviewer_f02_monthly_availability.csv';monthly.to_csv(csv_path,float_format='%.6f')
    fig,ax=plt.subplots(figsize=(9,4.5))
    for name in monthly.columns:ax.plot(monthly.index,monthly[name],marker='o',ms=4,label=name.capitalize())
    ax.set(xlabel='Month',ylabel='Mean network availability at native timestep (%)',ylim=(0,105),title='Reconstructed monthly network availability')
    ax.legend(ncol=3,loc='lower left',fontsize=8);ax.grid(alpha=.2)
    fig.tight_layout();target=out_dir/'reviewer_f02_monthly_availability.png'
    fig.savefig(target,dpi=dpi);plt.close(fig)
    return target,csv_path,monthly
