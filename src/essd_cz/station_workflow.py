"""Safe command builder and NetCDF metadata inspection for CHMI converter v5.

Does not reimplement CHMI parsing, value cleaning, or metadata joins.
"""
from __future__ import annotations
from pathlib import Path
import subprocess,sys
from .rainfall import ROOT

CONVERTER=ROOT/'original/code/chmi/chmi_historical_csv_to_netcdf_v5.py'
CATEGORIES=('precipitation','temperature','humidity')

def converter_command(category,year=2021,resolution='10min',stage='build',root_dir='opendata.chmi.cz',output=None,elements=None,workers=4,compression=4,duplicate_policy='error'):
    """Create an explicit CLI command. Default: build only; no automatic downloads."""
    if category not in CATEGORIES:raise ValueError(f'Choose one of {CATEGORIES}')
    if resolution!='10min':raise ValueError('Notebook example is intentionally restricted to 10min')
    if stage not in ('build','download','all'):raise ValueError('Unknown stage')
    if not (1900<=int(year)<=2100):raise ValueError('Invalid year')
    if not (0<=int(compression)<=9):raise ValueError('Compression must be 0..9')
    if int(workers)<1:raise ValueError('Workers must be positive')
    if duplicate_policy not in ('error','last'):raise ValueError('Invalid duplicate policy')
    root=Path(root_dir).expanduser().resolve()
    cmd=[sys.executable,str(CONVERTER),category,'--resolution',resolution,'--years',str(year),'--root-dir',str(root),'--stage',stage,'--workers',str(workers),'--compression',str(compression),'--duplicate-policy',duplicate_policy]
    if elements:cmd+=['--elements',*map(str,elements)]
    if output:
        p=Path(output).expanduser().resolve()
        if p.suffix.lower()!='.nc':raise ValueError('Output must be .nc')
        if p.exists():raise FileExistsError('Will not overwrite NetCDF: '+str(p))
        cmd+=['--output',str(p)]
    return cmd

def run_converter(command,approved=False):
    """No download or build without explicit approval in notebook cell."""
    if not approved:
        print('DRY RUN:',' '.join(map(str,command)))
        return None
    if not CONVERTER.is_file():raise FileNotFoundError(CONVERTER)
    return subprocess.run(command,check=True)

def inspect_station_netcdf(path):
    """Metadata-only validation. Verify station dimension, units and history."""
    from netCDF4 import Dataset
    p=Path(path).expanduser().resolve()
    if not p.is_file():raise FileNotFoundError(p)
    with Dataset(p) as ds:
        dims={k:len(v) for k,v in ds.dimensions.items()}
        variables={k:{'dimensions':list(v.dimensions),'units':str(getattr(v,'units',''))} for k,v in ds.variables.items()}
        info={'file':str(p),'dimensions':dims,'variables':variables,'history':str(getattr(ds,'history',''))[:1000]}
        info['groups']=list(ds.groups)
        if 'metadata' in ds.groups:
            info['metadata_groups']=list(ds.groups['metadata'].groups)
        if not any(k in dims for k in ('station','WSI')):raise ValueError('Missing station dimension')
        return info
