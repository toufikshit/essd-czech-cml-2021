"""Figure 4: time series of one CML, the radar path rainfall, the nearest rain gauge and the weather station,
13-14 July 2021.

Runs the functions event_windows() and event_plot() of the unchanged producer
original/code/figures/rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py, as its main()
does for one event. The CML of the figure is 45035_45097_2740 (--event-cml-id).
Inputs: the QC file, the 5-minute rainfall of scripts/derive_5min.py, radar, gauge, temperature and humidity files.
Outputs in <output_dir>/figures/fig04: FINAL_GRAF_event_<t0>_<t1>.png/.pdf, FINAL_GRAF_event_metadata.json.
Usage: python scripts/figures/fig04_event_series.py [--config config/user_paths.json]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from essd_cz.config import ORIGINAL, load_config  # noqa: E402
from essd_cz.rain5 import derived_5min_name  # noqa: E402

PRODUCER = ORIGINAL / 'figures' / 'rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
    ap.add_argument('--rain5', type=Path, default=None, help='5-minute rainfall (default: output of derive_5min.py)')
    ap.add_argument('--start', default='2021-04-01 00:00:00')
    ap.add_argument('--end', default='2021-10-31 23:59:59')
    ap.add_argument('--event-cml-id', default='45035_45097_2740')
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    cfg = load_config(a.config)
    ds = cfg['datasets']
    paths = {'cml': ds['qc'], 'rain': a.rain5 or cfg['output_dir'] / 'rainfall_5min' / derived_5min_name(ds['rain'].name),
             'radar': ds['radar'], 'gauge': ds['gauge'], 'temp': ds['temperature'], 'hum': ds['humidity']}
    for k, p in paths.items():
        if not Path(p).is_file():
            raise FileNotFoundError(f'{k}: {p}' + (' (run scripts/derive_5min.py first)' if k == 'rain' else ''))
    out = (a.out or cfg['output_dir'] / 'figures' / 'fig04').resolve()
    out.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location('essd_fig456_producer', PRODUCER)
    g = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g)
    start, end = pd.Timestamp(a.start), pd.Timestamp(a.end)
    result = {'inputs': {k: str(v) for k, v in paths.items()}}
    # Same open calls as producer main().
    with xr.open_dataset(paths['cml']) as cml, xr.open_dataset(paths['rain']) as rain, \
         xr.open_dataset(paths['radar']) as radar, xr.open_dataset(paths['gauge']) as gauge, \
         xr.open_dataset(paths['temp']) as temp, xr.open_dataset(paths['hum']) as hum:
        windows = g.event_windows(gauge, start, end)
        last = None
        for t0, t1, reason in windows:
            print(f'TRY EVENT: {t0} -> {t1}; source: {reason}', flush=True)
            try:
                result['event'] = g.event_plot(cml, rain, radar, gauge, temp, hum, out, t0, t1,
                                               preferred=a.event_cml_id,
                                               require_preferred=(a.event_cml_id is not None))
                result['event']['candidate_origin'] = reason
                print(f'SELECTED EVENT: {t0} to {t1}; CML {result["event"]["cml_id"]}', flush=True)
                break
            except RuntimeError as exc:
                last = exc
                print(f'NO ELIGIBLE LINK for {reason}: {exc}', flush=True)
        else:
            raise RuntimeError(f'No event found. Last error: {last}')
    (out / 'run_result.json').write_text(json.dumps(result, indent=2, default=str) + '\n')
    print('DONE', out, flush=True)


if __name__ == '__main__':
    sys.exit(main())
