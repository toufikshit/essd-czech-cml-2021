"""Auditable entry points for figures, full retrieval and optional validation."""
from __future__ import annotations
from pathlib import Path
import subprocess
import sys
from .rainfall import ROOT

AVAILABILITY = ROOT/'original/code/figures/availability_native_timestep_percentage.py'
PRODUCER = ROOT/'original/code/cml_rainfall/cml_graf_itu2005_rebuild.py'
VALIDATION = ROOT/'original/code/figures/rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py'
UNVERIFIED_V3 = ROOT/'original/notebook_references/wetdry_validation/cml_graf_rstd_hess_radar_gauge_v3.py'


def invoke(argv):
    """Execute the unchanged source with explicit arguments and propagated exit code."""
    print('ORIGINAL COMMAND:', ' '.join(map(str,argv)),flush=True)
    subprocess.run([str(x) for x in argv],check=True)


def availability(paths, output, dpi=300):
    """Native-timestep availability from the five input types; new output files."""
    invoke([sys.executable,AVAILABILITY,
        '--cml-file',paths['qc'],'--radar-file',paths['radar'],
        '--gauge-file',paths['gauge'],'--temperature-file',paths['temperature'],
        '--humidity-file',paths['humidity'],'--out-dir',output,'--dpi',str(dpi)])


def rebuild(paths, output, workers=4, python_executable=None):
    """Run the UNCHANGED rainfall producer using a configurable scientific Python.

    ``python_executable=None`` uses the Jupyter kernel's Python.  HPC users may
    choose a tested scientific environment independently of their notebook
    kernel.  This changes execution plumbing only, not any rainfall formulas.
    """
    output=Path(output)
    target=output/'recomputed_CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN.nc'
    if target.exists():raise FileExistsError(f'Refusing to overwrite {target}')
    science_python=Path(python_executable).expanduser() if python_executable else Path(sys.executable)
    if not science_python.is_file():
        raise FileNotFoundError(f'Scientific Python executable not found: {science_python}')
    print('SCIENTIFIC PYTHON:',science_python,flush=True)
    # Fail BEFORE generating a large NetCDF if dependencies are incompatible.
    invoke([science_python,'-c',
            'import numpy, numba, pycomlink, xarray, netCDF4; print("ESSD science imports: OK")'])
    invoke([science_python,PRODUCER,'--input-nc',paths['qc'],
            '--output-nc',target,'--summary-json',output/'rebuild_summary.json',
            '--workers',str(workers)])


def stage_symlink(src, dest):
    """Small workspace uses links to existing files; does not copy multi-GB data."""
    dest=Path(dest);src=Path(src).resolve()
    if dest.is_symlink() and dest.resolve()==src:return
    if dest.exists() or dest.is_symlink():
        raise FileExistsError(f'Cannot replace existing workspace path: {dest}')
    dest.symlink_to(src)


def validation(paths,output,detector=None,allow_unverified_v3=False,
               skip_event=True,only_check_sources=False):
    """Publication figures; no silently substituted Graf detector.

    Exact source file `cml_graf_rstd_hess_radar_gauge.py` is unrecovered.
    Passing v3 is deliberately an explicit experiment, not an exact replay.
    """
    if bool(detector)==bool(allow_unverified_v3):
        raise ValueError('Specify either exact --detector or --allow-unverified-v3')
    src=UNVERIFIED_V3 if allow_unverified_v3 else Path(detector).expanduser().resolve()
    if not src.is_file():raise FileNotFoundError(src)
    if allow_unverified_v3:print('WARNING: v3 detector is NOT PROVEN the exact publication source.',flush=True)
    workspace=output/'validation_workspace'
    workspace.mkdir(parents=True,exist_ok=True)
    for key in ('qc','gauge','radar','temperature','humidity'):
        stage_symlink(paths[key],workspace/paths[key].name)
    stage_symlink(src,workspace/'cml_graf_rstd_hess_radar_gauge.py')
    command=[sys.executable,VALIDATION,'--root',workspace,'--graf-rain-nc',paths['rain']]
    if skip_event:command.append('--skip-event')
    if only_check_sources:command.append('--only-check-sources')
    invoke(command)
    return workspace/'figures'
