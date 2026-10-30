"""Check the seven released files (dimensions, entity counts, required variables, time axis) without reading data."""
from __future__ import annotations
from datetime import datetime, timedelta
from netCDF4 import Dataset, date2num, num2date
import numpy as np
from .config import ROLES

START = datetime(2021, 4, 1)
END_EXCLUSIVE = datetime(2021, 11, 1)
CONTRACT = {
    'raw': ('time', 1, 308160, 'cml_id', 2574),
    'qc': ('time', 1, 308160, 'cml_id', 1987),
    'rain': ('time_1min', 1, 308160, 'cml_id', 1987),
    'radar': ('time', 5, 61632, 'cml_id', 1132),
    'gauge': ('time', 10, 30816, 'WSI', 424),
    'temperature': ('time', 10, 30816, 'station', 251),
    'humidity': ('time', 10, 30816, 'station', 247),
}
REQUIRED = {
    'qc': ('trsl',), 'rain': ('R_sublink_1min',),
    'gauge': ('rr',), 'radar': ('rr',),
    'temperature': ('T',), 'humidity': ('H',),
}


def inspect_one(role, filename):
    """Check one file: entity count, required variables and the 1 April to 31 October 2021 time axis."""
    td, step_minutes, expected_slots, entity, nentity = CONTRACT[role]
    required = REQUIRED.get(role, ())
    with Dataset(filename) as ds:
        problems=[]
        if entity not in ds.dimensions or len(ds.dimensions[entity])!=nentity:
            problems.append(f'Expected {entity}={nentity}')
        for v in required:
            if v not in ds.variables: problems.append(f'Missing variable {v}')
        if td not in ds.variables:
            problems.append(f'Missing coordinate {td}')
            return {'role':role,'file':str(filename),'problems':problems,'ok':False}
        tv=ds.variables[td]
        units=tv.units; calendar=getattr(tv,'calendar','standard')
        x=np.asarray(tv[:],dtype=np.float64)
        low=date2num(START,units,calendar=calendar)
        high=date2num(END_EXCLUSIVE,units,calendar=calendar)
        step=date2num(START+timedelta(minutes=step_minutes),units,calendar=calendar)-low
        i0=int(np.searchsorted(x,low));i1=int(np.searchsorted(x,high))
        selected=x[i0:i1];tol=max(abs(step)*1e-6,1e-8)
        exact=(len(selected)==expected_slots and
               abs(selected[0]-low)<tol and abs(selected[-1]-(high-step))<tol and
               np.all(np.abs(np.diff(selected)-step)<tol)) if len(selected)>0 else False
        if not exact:problems.append(f'April--October axis not continuous: {len(selected)}/{expected_slots}')
        if role=='rain':
            kr=str(getattr(ds,'kr_model',''))
            if kr and '2005' not in kr:problems.append(f'Expected ITU2005, file has {kr}')
        return {
            'role':role,'file':str(filename),'entities':nentity,
            'native_minutes':step_minutes,'source_slots':len(x),
            'selected_slots':len(selected),'required_slots':expected_slots,
            'selected_start':str(num2date(float(selected[0]),units,calendar=calendar)) if len(selected) else None,
            'selected_end':str(num2date(float(selected[-1]),units,calendar=calendar)) if len(selected) else None,
            'is_physically_cropped':len(x)==expected_slots,
            'problems':problems,'ok':not problems,
        }


def inspect_all(paths):
    return [inspect_one(role,paths[role]) for role in ROLES]
