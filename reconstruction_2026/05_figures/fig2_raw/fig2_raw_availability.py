"""Stage 6: Figure 2 (availability) with the CML row from the RAW dataset (2574 CMLs).

Runs the UNMODIFIED producer original/code/figures/availability_native_timestep_percentage.py (sha256 b8d7b983...)
with its own command-line options. The only addition: when the producer opens the released RAW CML file, the dataset
is restricted to the 2574 CML IDs of the approved RAW population (stage2/RAW_APPROVED_POPULATION.csv); the released
RAW file also holds 6 R1 CMLs that are not part of the RAW population. Nothing else is changed:
  * CML rule of the producer: a CML is available at a minute when at least one sublink has both TX and RX
    (paired TX/RX observation, RAW values, no QC mask);
  * radar, gauge, temperature and humidity files are the same as for the current Fig. 2 (stage5b run_fig2.sh), so
    their rows are unchanged.
Usage: python -I fig2_raw_availability.py --cml-file <released RAW nc> --radar-file ... --gauge-file ...
       --temperature-file ... --humidity-file ... --out-dir <dir>      (all options of the producer)
"""
import csv, importlib.util, sys
from pathlib import Path
import xarray as xr

R = Path(r'C:\Users\taouf\ESSD_RELEASE')
PRODUCER = R / 'essd-czech-cml-2021' / 'original' / 'code' / 'figures' / 'availability_native_timestep_percentage.py'
POP = R / 'FINAL_RECONSTRUCTION_20261027' / 'stage2' / 'RAW_APPROVED_POPULATION.csv'

ids = sorted({r['cml_id'] for r in csv.DictReader(open(POP, encoding='utf-8'))})
assert len(ids) == 2574, len(ids)
cml_file = Path(sys.argv[sys.argv.index('--cml-file') + 1]).resolve()

_open = xr.open_dataset


def open_restricted(path, *args, **kw):
    ds = _open(path, *args, **kw)
    if Path(path).resolve() == cml_file:
        missing = sorted(set(ids) - set(ds['cml_id'].values.astype(str)))
        assert not missing, missing[:5]
        ds = ds.sel(cml_id=ids)
        print(f'[RAW population] {ds.sizes["cml_id"]} CMLs selected from {cml_file.name}')
    return ds


xr.open_dataset = open_restricted
spec = importlib.util.spec_from_file_location('availability_producer', PRODUCER)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod      # needed by the producer's dataclass
spec.loader.exec_module(mod)
mod.main()
