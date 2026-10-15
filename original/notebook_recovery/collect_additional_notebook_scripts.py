#!/usr/bin/env python3
"""Read-only, exact-path collection of scripts explicitly referenced by two ESSD notebooks.

Run on the Jupyter server where /home/jovyan/<script.py> is located.
Never executes the referenced scripts and never searches directories recursively.
The original QC producer of the March 2026 cml_preprocessed_v3.nc is NOT identified here.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path('/home/jovyan')
SCRIPTS = {
    'availability_2d_histogram_independent_masked_no_inside_labels.py': ('figures_legacy', 'dataste_papaer.ipynb cell 9'),
    'batch_cml_vs_radar_joint_events_essd_1h.py': ('event_maps_legacy', 'dataste_papaer.ipynb cell 3'),
    'batch_cml_vs_radar_joint_events_essd_multiscale_cartopy_v3.py': ('event_maps_legacy', 'dataste_papaer.ipynb cell 2'),
    'batch_cml_vs_radar_joint_events_essd_vertical.py': ('event_maps_legacy', 'dataste_papaer.ipynb cell 0'),
    'blettner_event_rainfall_maps.py': ('event_maps_legacy', 'dataste_papaer.ipynb cell 7'),
    'cml_compare_single_combo_pycml_itu03_pastorek_hexbin_oneplot_reference.py': ('validation_legacy', 'dataste_papaer.ipynb cell 12'),
    'cml_compare_single_combo_wet_radar_only.py': ('validation_legacy', 'DATASET_PUBLICATION.ipynb cell 9'),
    'cml_dual_wetdry_pastorek_optimized.py': ('rainfall_experiments', 'DATASET_PUBLICATION.ipynb cell 16'),
    'cml_estimate_and_grid_cml_radar_gauge_merged_parallel_timewindow.py': ('rainfall_experiments', 'DATASET_PUBLICATION.ipynb cell 14'),
    'cml_estimate_and_grid_cml_radar_gauge_merged_parallel_timewindow_v2.py': ('rainfall_experiments', 'DATASET_PUBLICATION.ipynb cells 8, 13'),
    'cml_graf_blettner_publishable_pipeline_v2.py': ('rainfall_experiments', 'DATASET_PUBLICATION.ipynb cell 17'),
    'cml_graf_rstd_hess_radar_gauge_v3.py': ('wetdry_validation', 'dataste_papaer.ipynb cell 11'),
    'cml_multi_waa_itu_fullperiod_optimized.py': ('rainfall_experiments', 'DATASET_PUBLICATION.ipynb cell 18'),
    'cml_wet_dry_rsd_contingency_matrix_only.py': ('wetdry_validation', 'dataste_papaer.ipynb cell 8'),
    'compare_two_cml_rain_products_same_samples_fast.py': ('qc_sensitivity', 'DATASET_PUBLICATION.ipynb cell 12'),
    'filter_rainrate_empty_only_if_trsl_empty.py': ('qc_sensitivity', 'dataste_papaer.ipynb cell 13'),
    'inspect_essd_three_netcdf.py': ('publication_diagnostics', 'DATASET_PUBLICATION.ipynb cell 1'),
    'run_filter_strategy_rainfall_experiment.py': ('qc_sensitivity', 'DATASET_PUBLICATION.ipynb cell 13'),
}

README = '''Additional source-code files REFERENCED by the two older dataset-paper notebooks.
Original script bytes are copied unchanged; no NetCDFs and no execution of scripts.
Notebook cell numbers are zero-based and show reference locations, not proof of success.
These older experiments often use ITU2003 and the full 2020–2021 period: they are NOT the
final April–October 2021 Graf/ITU2005/Pastorek rainfall producer.
This collection does not recover the original March 2026 six-step QC producer of
cml_preprocessed_v3.nc. The later strict QC v4_3 code is already in the initial ZIP.
Before any public code release, determine copyright/software licence and review for secrets.
'''

def build(out: Path, root: Path):
    included = []
    absent = []
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        z.writestr('README.txt', README)
        for filename, (category, notebook_cell) in sorted(SCRIPTS.items()):
            p = root / filename
            if not p.is_file():
                absent.append(str(p))
                print('MISSING', p)
                continue
            if p.stat().st_size > 20_000_000:
                absent.append(str(p) + ' (too large: >20 MB)')
                print('SKIPPED large source', p)
                continue
            data = p.read_bytes()
            arch = f'notebook_references/{category}/{filename}'
            z.writestr(arch, data)
            included.append({'source_path': str(p), 'notebook_reference': notebook_cell,
                             'zip_member': arch, 'sha256': hashlib.sha256(data).hexdigest(),
                             'bytes': len(data)})
            print('COPIED', p.name, '->', category, f'({len(data):,} bytes)')
        z.writestr('manifest.json', json.dumps({
            'created_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
            'source': 'notebook-referenced scripts only, exact paths',
            'included': included, 'not_found_or_skipped': absent,
            'original_march_qc_producer': 'NOT IDENTIFIED',
            'note': 'References are not proof of execution or creation of the scientific datasets.'
        }, indent=2) + '\n')
    with zipfile.ZipFile(out) as z:
        bad = z.testzip()
        if bad: raise RuntimeError('Corrupt ZIP member: ' + bad)
    print('\nARCHIVE:', out)
    print('FOUND:', len(included), 'MISSING:', len(absent), 'ZIP size:',out.stat().st_size)
    print('No source files or datasets were modified.')

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, default=ROOT, help='Source folder containing the original Python scripts')
    ap.add_argument('--out', type=Path, default=Path('ESSD_ADDITIONAL_NOTEBOOK_REFERENCED_SCRIPTS.zip'))
    a = ap.parse_args()
    build(a.out,a.root)
