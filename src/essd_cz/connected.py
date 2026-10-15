"""Connected, reviewer-facing ESSD workflow. Original science code stays byte-identical.

One configuration for seven released NetCDF paths; one manifest for derived artifacts.
The new IDW step adapts the original project's unchanged midpoint IDW functions
from older 1-min products to the *final* 5-min Graf link rain rates. Its output is
an adapted review grid, not an asserted byte-identical historic IDW grid.
"""
from __future__ import annotations
import ast
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
from .config import load_paths

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT / 'original' / 'code'


def extras_from_config(config):
    p = Path(config).expanduser().resolve()
    obj = json.loads(p.read_text(encoding='utf-8'))
    out = {}
    for k, v in obj.get('extras', {}).items():
        if not v:
            out[k] = None
        else:
            x = Path(v).expanduser()
            out[k] = (x if x.is_absolute() else (p.parent / x)).resolve()
    return out


class Project:
    """Load source data from any location, pass each stage's files to the next."""
    def __init__(self, config, output=None):
        self.config = Path(config).expanduser().resolve()
        self.paths = load_paths(self.config)
        self.extras = extras_from_config(self.config)
        self.output = Path(output or ROOT / 'results' / 'connected_workflow').expanduser().resolve()
        self.output.mkdir(parents=True, exist_ok=True)
        self.manifest = self.output / 'pipeline_manifest.json'
        self.state = json.loads(self.manifest.read_text()) if self.manifest.exists() else {'stages': {}, 'source_config': str(self.config)}
        if self.state.get('source_config') != str(self.config):
            raise ValueError('This output manifest was created with another configuration. Set a different workspace in config/active_profile.json before switching sources.')

    def _resolve_artifact(self, name):
        value = self.state['stages'].get(name)
        if value is None:
            raise RuntimeError(f'No output registered for {name!r}. Run its previous notebook first.')
        p = self.output / value['relative_path']
        if not p.is_file():
            raise FileNotFoundError(f'Staged artifact {name} is absent: {p}')
        return p.resolve()

    def register(self, name, path, provenance, **params):
        path = Path(path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        try:
            rel = str(path.relative_to(self.output))
        except ValueError:
            # External input not copied; its absolute path is a property of the user's config.
            rel = str(path)
        self.state['stages'][name] = {'relative_path': rel, 'provenance': provenance,
                                      'parameters': params, 'date_utc': datetime.now(timezone.utc).isoformat()}
        self.manifest.write_text(json.dumps(self.state, indent=2) + '\n')
        print(f'[REGISTERED] {name}: {path}\n  source: {provenance}', flush=True)
        return path

    def get(self, name):
        return self._resolve_artifact(name)

    def stage_dir(self, name):
        p = self.output / name
        p.mkdir(parents=True, exist_ok=True)
        return p

    def source_rainfall(self):
        """Stage the published rain product without suggesting it was recomputed."""
        return self.register('rainfall', self.paths['rain'], 'published-final-ITU2005-Pastorek-Graf')

    def recompute_rainfall(self, workers=4, python_executable=None):
        """Opt-in ORIGINAL full-history producer; optional separate scientific Python."""
        from .workflows import rebuild
        out = self.stage_dir('01_rainfall_rebuild')
        target = out / 'recomputed_CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN.nc'
        if target.exists():
            raise FileExistsError(target)
        rebuild(self.paths, out, workers=workers, python_executable=python_executable)
        return self.register('rainfall', target, 'recomputed-original-full-history',
                             workers=workers, science_python=str(python_executable or sys.executable))


def _original_idw_functions():
    """Load the exact two unchanged mathematical functions by AST, not a rewrite.

    The original module also imports netCDF4 and plotting code. Loading only the
    two original function nodes makes this geometry testable without those dependencies.
    Source bytes are verified independently by tools_verify_originals.py.
    """
    source = ORIGINAL / 'figures' / '02_blettner_aproct_dataset_figures_grid.py'
    tree = ast.parse(source.read_text(encoding='utf-8'), filename=str(source))
    chosen = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in {'build_static_idw','idw_field'}]
    if len(chosen) != 2:
        raise RuntimeError('Original IDW functions are absent or renamed; refusing to substitute')
    namespace = {'np': np, 'cKDTree': cKDTree}
    code = compile(ast.Module(body=chosen, type_ignores=[]), str(source), 'exec')
    exec(code, namespace)
    return namespace['build_static_idw'], namespace['idw_field']


def grid_final_rain(rain_nc, qc_nc, output_nc, start, end, *,
                    resolution=.01, pad=.05, nnear=15, power=2., max_dist=.3,
                    extent=None, max_grid_cells=100000):
    """Stream 5-min final CML R_mean rates to lon/lat IDW grid using ORIGINAL weights.

    The source historical gridding used an earlier 1-min Blettner product, so
    this adapter is a documented new implementation of the published 5-min rates.
    Its output (time,y,x,R) is compatible with original spatial map renderer.
    End time is inclusive. No input files are modified. No silent time interpolation.
    """
    output_nc = Path(output_nc).expanduser().resolve()
    if output_nc.exists(): raise FileExistsError(output_nc)
    if resolution <= 0 or nnear < 1 or power <= 0 or max_dist <= 0: raise ValueError('Invalid IDW parameters')
    build_static_idw, idw_field = _original_idw_functions()
    with xr.open_dataset(qc_nc, decode_times=False) as qc, xr.open_dataset(rain_nc) as rain:
        ids_q = [str(v.decode('utf8','replace') if isinstance(v,bytes) else v).rstrip('\x00').strip() for v in qc['cml_id'].values]
        ids_r = [str(v.decode('utf8','replace') if isinstance(v,bytes) else v).rstrip('\x00').strip() for v in rain['cml_id'].values]
        if len(set(ids_q)) != len(ids_q) or len(set(ids_r)) != len(ids_r): raise ValueError('Duplicate structural CML IDs')
        ridx = {c:i for i,c in enumerate(ids_r)}
        source_vars = ('site_a_longitude','site_b_longitude','site_a_latitude','site_b_latitude')
        for name in source_vars:
            if name not in qc: raise KeyError('Missing CML geometry '+name)
        lon0,lon1,lat0,lat1 = [np.asarray(qc[v].values,dtype=float).reshape(-1) for v in source_vars]
        mid_lon,mid_lat = (lon0+lon1)/2,(lat0+lat1)/2
        keep = np.array([i for i,c in enumerate(ids_q) if c in ridx and np.isfinite(mid_lon[i]) and np.isfinite(mid_lat[i])], dtype=int)
        if keep.size == 0: raise ValueError('No matching CML IDs with valid geometry')
        indices = np.array([ridx[ids_q[i]] for i in keep], dtype=int)
        if extent is None:
            x0,x1 = np.min(mid_lon[keep])-pad, np.max(mid_lon[keep])+pad
            y0,y1 = np.min(mid_lat[keep])-pad, np.max(mid_lat[keep])+pad
        else:
            if len(extent)!=4: raise ValueError('extent must be [lon_min, lon_max, lat_min, lat_max]')
            x0,x1,y0,y1 = map(float,extent)
        if not (x0 < x1 and y0 < y1): raise ValueError('Invalid geographic extent')
        xx=np.arange(x0,x1+resolution*.1,resolution)
        yy=np.arange(y0,y1+resolution*.1,resolution)
        if not xx.size or not yy.size or xx.size*yy.size > max_grid_cells:
            raise ValueError(f'Grid too large: {xx.size} × {yy.size}. Increase resolution or set an extent.')
        X,Y,idx,dist,within=build_static_idw(mid_lon[keep],mid_lat[keep],xx,yy,nnear,max_dist)
        if 'R_mean_5min' not in rain: raise KeyError('Final rainfall needs R_mean_5min')
        r=rain['R_mean_5min'].transpose('time_5min','cml_id')
        all_times=pd.DatetimeIndex(rain['time_5min'].values)
        start=pd.Timestamp(start);end=pd.Timestamp(end)
        selected=np.flatnonzero((all_times>=start)&(all_times<=end))
        if not len(selected): raise ValueError('No rainfall times inside selected date interval')
        ts=all_times[selected]
        if len(ts)>1 and not np.all(np.diff(ts.asi8)==300_000_000_000):raise ValueError('Non-contiguous 5-min samples; refusing implicit interpolation')
        print(f'GRID: {len(keep)} CMLs, {len(ts)} time steps, {len(yy)} × {len(xx)} cells; no changes to input',flush=True)
        output_nc.parent.mkdir(parents=True, exist_ok=True)
        try:
            import netCDF4
            fh=netCDF4.Dataset(output_nc,'w',format='NETCDF4')
            netcdf4=True
        except ImportError:
            from scipy.io import netcdf_file
            fh=netcdf_file(str(output_nc),'w')
            netcdf4=False
        try:
            for key,count in [('time',len(ts)),('y',len(yy)),('x',len(xx))]: fh.createDimension(key,count)
            def v(name, dtype, dims, **opts):
                if netcdf4:return fh.createVariable(name,dtype,dims,zlib=True,complevel=2,**opts)
                return fh.createVariable(name,{'f8':'d','f4':'f','i2':'h'}[dtype],dims)
            v('time','f8',('time',))[:] = ts.asi8.astype('float64')/1e9
            fh.variables['time'].units='seconds since 1970-01-01 00:00:00'
            fh.variables['time'].calendar='standard'
            v('x','f8',('x',))[:] = xx
            v('y','f8',('y',))[:] = yy
            fh.variables['x'].units='degrees_east';fh.variables['y'].units='degrees_north'
            R=v('R','f4',('time','y','x'))
            R.units='mm h-1'
            Nv=v('nvalid','i2',('time','y','x'))
            fh.title='Reviewer adaptation: midpoint IDW of final Graf 5-min CML rainfall'
            fh.processing_note='Unchanged original IDW functions, final R_mean_5min; NOT byte-identical historic 1-min-grid workflow'
            fh.spatial_idw_parameters=json.dumps(dict(resolution=resolution,nnear=nnear,power=power,max_dist=max_dist,pad=pad))
            for oi,ti in enumerate(selected):
                rr=np.asarray(r.isel(time_5min=int(ti)).values,dtype='f4')[indices]
                field,nvalid=idw_field(rr,idx,dist,within,power)
                R[oi,:,:]=field.reshape(X.shape)
                Nv[oi,:,:]=nvalid.reshape(X.shape)
        except BaseException:
            fh.close()
            output_nc.unlink(missing_ok=True)
            raise
        else:
            fh.close()
    return output_nc


def manual_event_selection(csv, *, event_time, extent, event_rank=4, instance_rank=1):
    """An EXPLICIT hand-selected demonstration event, never claimed original event discovery."""
    if len(extent)!=4: raise ValueError('extent = lonmin lonmax latmin latmax')
    dest=Path(csv).expanduser().resolve()
    if dest.exists():raise FileExistsError(dest)
    dest.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame([dict(event_rank=int(event_rank),instance_rank=int(instance_rank),time=str(pd.Timestamp(event_time)),start_time=str(pd.Timestamp(event_time)),extent_lon_min=float(extent[0]),extent_lon_max=float(extent[1]),extent_lat_min=float(extent[2]),extent_lat_max=float(extent[3]),selection_origin='EXPLICIT_REVIEWER_DEMONSTRATION_NOT_HISTORIC_EVENT_DISCOVERY')]).to_csv(dest,index=False)
    return dest


def run_native_radar_map(project, *, first_panel='2021-07-13 16:05:00', ncols=8,
                         step_minutes=60, event_rank=4, instance_rank=1,
                         extent=None, cache=None, dpi=200):
    """Use original map renderer, passing grid generated in previous notebook."""
    grid=project.get('spatial_grid')
    root=project.extras.get('native_radklim_root')
    if root is None or not root.is_dir():raise FileNotFoundError('Set extras.native_radklim_root in chosen config')
    script=ORIGINAL/'figures'/'03_blettner_event_maps_timeseries.py'
    if cache is None:
        cache=project.extras.get('historical_event_cache')
    if not cache or not Path(cache).is_file():
        if extent is None:raise FileNotFoundError('Provide historical event cache OR explicit demonstration extent (lon_min,lon_max,lat_min,lat_max)')
        cache=manual_event_selection(project.stage_dir('04_maps')/'explicit_demo_event_cache.csv',event_time=first_panel,extent=extent,event_rank=event_rank,instance_rank=instance_rank)
        print('NOTE: Manual selection does not reproduce historical event discovery!',flush=True)
    else:
        cache=Path(cache)
    # Don't allow unbounded recursive fallback in the original by always passing cache explicitly.
    last_panel=pd.Timestamp(first_panel)+pd.Timedelta(minutes=(ncols-1)*step_minutes)
    with xr.open_dataset(grid) as ds:
        tt=pd.DatetimeIndex(ds['time'].values)
        if tt.min()>pd.Timestamp(first_panel)-pd.Timedelta(hours=1) or tt.max()<last_panel:
            raise ValueError('Previous notebook grid does not cover entire 1-h rolling map window. Extend grid time first.')
    try:import cartopy  # noqa:F401
    except ImportError as ex:raise RuntimeError('Install cartopy in HPC environment for original map script') from ex
    out=project.stage_dir('04_maps')
    cmd=[sys.executable,str(script),'--gridded-nc',str(grid),'--cml-var','R',
        '--cml-row-label','GRAF','--cml-geom-nc',str(project.paths['qc']),
        '--radar-root',str(root),'--selection-cache-csv',str(cache),
        '--event-rank',str(event_rank),'--instance-rank',str(instance_rank),
        '--start-time',str(first_panel),'--ncols',str(ncols),
        '--panel-step-minutes',str(step_minutes),'--dpi',str(dpi),
        '--skip-timeseries','--out-dir',str(out)]
    print('ORIGINAL MAP:', ' '.join(cmd),flush=True)
    subprocess.run(cmd,check=True)
    project.register('map_hourly',out/'hourly_1h'/f'event{event_rank:02d}_instance{instance_rank:02d}_multitime_hourly1h.png','original renderer with adapted 5-min IDW input',event_cache=str(cache))
    return out


def stage_wvd(project, *, recompute=False, max_series=3, pressure_hpa=980.):
    """Use existing V6 output by default; optional original small E-band retrieval."""
    out=project.stage_dir('05_wvd')
    if not recompute:
        p=project.extras.get('wvd_v6_nc')
        if p is None or not p.is_file():raise FileNotFoundError('Provide extras.wvd_v6_nc or set recompute=True for original V6')
        return project.register('wvd_v6',p,'precomputed V6 supporting demonstration')
    script=ORIGINAL/'wvd'/'19_wvd_eband_progressive_baseline_v6.py'
    helper=ORIGINAL/'wvd'/'13_wvd_raw_60_100_independent.py'
    dest=out/'reviewer_wvd_v6_small.nc'
    if dest.exists():raise FileExistsError(dest)
    prefix=out/'reviewer_wvd_v6_small'
    cmd=[sys.executable,str(script),'--helper-script',str(helper),'--cml-nc',str(project.paths['raw']),
         '--gauge-nc',str(project.paths['gauge']),'--temperature-nc',str(project.paths['temperature']),
         '--out-nc',str(dest),'--diagnostics-csv',str(prefix)+'_diagnostics.csv',
         '--states-csv',str(prefix)+'_states.csv','--events-csv',str(prefix)+'_events.csv',
         '--reinit-csv',str(prefix)+'_reinit.csv','--max-series',str(max_series),
         '--pressure-hpa',str(pressure_hpa)]
    print('ORIGINAL V6:', ' '.join(cmd),flush=True)
    subprocess.run(cmd,check=True)
    return project.register('wvd_v6',dest,'original V6 small-series run',max_series=max_series,pressure_hpa=pressure_hpa)


def wvd_sensitivity(project, *, branch='GRAF', max_station_distance_km=5., dpi=200):
    """Original WVD binned sensitivity; next notebook consumes saved pair cache."""
    quant=project.extras.get('wvd_quantization_csv')
    diag=project.extras.get('wvd_v6_diagnostics_csv')
    if quant is None or not quant.is_file():raise FileNotFoundError('Set extras.wvd_quantization_csv; quantization diagnostic is NOT in seven NetCDFs')
    if diag is None or not diag.is_file():raise FileNotFoundError('Set extras.wvd_v6_diagnostics_csv')
    v6=project.get('wvd_v6')
    out=project.stage_dir('WVD_V6_QUANTIZATION_BINNED_VALIDATION')
    script=ORIGINAL/'wvd'/'33_v6_binned_validation_hexbins.py'
    helper=ORIGINAL/'wvd'/'21_validate_wvd_v6_vs_original.py'
    cmd=[sys.executable,str(script),'--root',str(project.output),'--v6-nc',str(v6),
        '--v6-diagnostics-csv',str(diag),'--quant-csv',str(quant),
        '--raw-cml-nc',str(project.paths['raw']),'--humidity-nc',str(project.paths['humidity']),
        '--temperature-nc',str(project.paths['temperature']),'--validation-helper-script',str(helper),
        '--station-assignments-csv',str(out/'station_assignments.csv'),'--out-dir',str(out),
        '--branch',branch,'--max-station-distance-km',str(max_station_distance_km),'--dpi',str(dpi)]
    print('ORIGINAL WVD SENSITIVITY:', ' '.join(cmd),flush=True)
    subprocess.run(cmd,check=True)
    project.register('wvd_binned_metrics',out/'binned_validation_metrics.csv','original 33 V6 sensitivity',branch=branch)
    project.register('wvd_link_population',out/'validation_population_links.csv','original 33 station matched links')
    return out


def plot_paper_wvd(project, *, dpi=300, gridsize=55, plot_max=35.):
    """Original paper-population WVD hexbin with exact L and quantization filters."""
    project.get('wvd_binned_metrics')
    out=project.stage_dir('07_wvd_paper_population')
    script=ORIGINAL/'wvd'/'plot_wvp_support_hexbins_no_titles.py'
    cmd=[sys.executable,str(script),'--root',str(project.output),'--out-dir',str(out),
         '--dpi',str(dpi),'--gridsize',str(gridsize),'--plot-max',str(plot_max)]
    print('ORIGINAL WVD PAPER HEXBIN:', ' '.join(cmd),flush=True)
    subprocess.run(cmd,check=True)
    project.register('wvd_paper_metrics',out/'wvp_support_metrics.csv',
                     'original paper population L>=4 km, delta-rho-q<=0.20 g/m3, E-band, GRAF')
    return out
