#!/usr/bin/env python3
"""Stage 5B / Agent 3b - Fig. 4 runner (reproduction AND regeneration).

Runs the UNCHANGED published producer
  essd-czech-cml-2021/original/code/figures/
  rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py
(functions event_windows, event_plot, style, station_series, hav, ...) by
importing the file from its original location (importlib, no copy, no edit).

Only the input paths differ between the two modes; they are CLI arguments.
The producer's own resolve_paths() builds file names under one --root and
falls back to /home/jovyan; it is bypassed here and the six paths are passed
explicitly (the same six datasets that main() opens). event_windows() and
event_plot() are then called exactly as in main() (only_event branch):
the first window that yields an eligible CML is used, with --event-cml-id
optional (main() passes preferred=None unless --event-cml-id is given).

The hexbin / confusion branch is NOT run (only_event mode of the producer).

Outputs go to --out (the producer writes FINAL_GRAF_event_<t0>_<t1>.png/.pdf
and FINAL_GRAF_event_metadata.json there).

Example (reproduction, original inputs):
  python -I f4_01_event_plot_run.py --mode original --out ..\\reproduction\\fig4
Example (regeneration, corrected inputs):
  python -I f4_01_event_plot_run.py --mode stage4 --out ..\\work\\fig4_new
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

# Never write __pycache__ into the read-only repository when importing the producer.
sys.dont_write_bytecode = True

import pandas as pd
import xarray as xr

ESSD = Path(r'C:\Users\taouf\ESSD_RELEASE')
PRODUCER = (ESSD / 'essd-czech-cml-2021' / 'original' / 'code' / 'figures' /
            'rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py')
FH = ESSD / 'extracted' / 'full' / 'ESSD_CZ_2021_HANDOFF' / 'data' / 'source_full_history'
REC = ESSD / 'FINAL_RECONSTRUCTION_20261027'
ORIGINAL = {
    'cml': FH / 'qc' / 'CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc',
    'rain': FH / 'rain' / 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_20200101_20211231.nc',
    'radar': FH / 'radar' / 'CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc',
    'gauge': FH / 'gauge' / 'CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc',
    'temp': FH / 'temperature' / 'CHMI_10MIN_TEMPERATURE_2021.nc',
    'hum': FH / 'humidity' / 'CHMI_10MIN_HUMIDITY_2021.nc',
}
STAGE4 = dict(ORIGINAL)
STAGE4['cml'] = REC / 'stage3' / 'QC_CORRECTED_CANDIDATE_REV1.nc'
STAGE4['rain'] = REC / 'stage4' / 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc'


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load_producer(path: Path):
    spec = importlib.util.spec_from_file_location('essd_fig456_producer', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', choices=['original', 'stage4'], required=True)
    ap.add_argument('--producer', type=Path, default=PRODUCER)
    for k in ORIGINAL:
        ap.add_argument(f'--{k}', type=Path, default=None, help=f'override {k} path')
    ap.add_argument('--start', default='2021-04-01 00:00:00')
    ap.add_argument('--end', default='2021-10-31 23:59:59')
    ap.add_argument('--event-cml-id', default=None)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    paths = dict(ORIGINAL if a.mode == 'original' else STAGE4)
    for k in paths:
        if getattr(a, k) is not None:
            paths[k] = getattr(a, k)
    for k, p in paths.items():
        if not p.is_file():
            raise FileNotFoundError(f'{k}: {p}')
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    print('PRODUCER', a.producer, 'sha256', sha256(a.producer), flush=True)
    for k, v in paths.items():
        print(f'{k}: {v}', flush=True)
    g = load_producer(a.producer)
    start, end = pd.Timestamp(a.start), pd.Timestamp(a.end)
    result = {'mode': a.mode, 'inputs': {k: str(v) for k, v in paths.items()}}
    # Same open calls as producer main().
    with xr.open_dataset(paths['cml']) as cml, xr.open_dataset(paths['rain']) as rain, \
         xr.open_dataset(paths['radar']) as radar, xr.open_dataset(paths['gauge']) as gauge, \
         xr.open_dataset(paths['temp']) as temp, xr.open_dataset(paths['hum']) as hum:
        windows = g.event_windows(gauge, start, end)
        last = None
        for t0, t1, reason in windows:
            print(f'TRY DIFFERENT EVENT: {t0} -> {t1}; source: {reason}', flush=True)
            try:
                result['event'] = g.event_plot(cml, rain, radar, gauge, temp, hum, out, t0, t1,
                                               preferred=a.event_cml_id,
                                               require_preferred=(a.event_cml_id is not None))
                result['event']['candidate_origin'] = reason
                print(f'SELECTED NEW EVENT: {t0} to {t1}; physical CML {result["event"]["cml_id"]}', flush=True)
                break
            except RuntimeError as exc:
                last = exc
                print(f'NO ELIGIBLE LINK for {reason}: {exc}', flush=True)
        else:
            raise RuntimeError(f'No event found. Last error: {last}')
    (out / 'run_result.json').write_text(json.dumps(result, indent=2, default=str) + '\n')
    print('DONE', flush=True)


if __name__ == '__main__':
    sys.exit(main())
