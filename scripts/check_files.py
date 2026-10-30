"""Check the seven released NetCDF files: names, dimensions, entity counts, required variables and time axis.

Reads only the metadata and the time axis. Prints PASS or FAIL for every file.
Usage: python scripts/check_files.py [--config config/user_paths.json]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from essd_cz.catalog import inspect_all  # noqa: E402
from essd_cz.config import load_paths  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
a = ap.parse_args()
results = inspect_all(load_paths(a.config))
for item in results:
    print(f"{item['role']:12} {item.get('entities', '?'):>5} entities, {item.get('selected_slots', '?'):>6} time steps: "
          f"{'PASS' if item['ok'] else 'FAIL'}  {Path(item['file']).name}")
    for issue in item['problems']:
        print('    ERROR:', issue)
sys.exit(0 if all(r['ok'] for r in results) else 1)
