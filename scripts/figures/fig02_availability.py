"""Figure 2: percentage of CMLs, radar paths, rain gauges and weather stations with valid data at every time step.

Runs the unchanged original/code/figures/availability_native_timestep_percentage.py on the seven released files.
A CML is available at a minute when at least one sublink has both TX and RX (raw values, no QC mask); 2574 CMLs.
The published figure is the producer's PNG without its top 125 pixel rows (the year label); this script writes both.
Outputs in <output_dir>/figures/fig02: the producer's PNG/PDF and five CSV files, and
availability_native_timestep_percentage_v1.png (the published figure).
Usage: python scripts/figures/fig02_availability.py [--config config/user_paths.json]
"""
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from essd_cz.config import ORIGINAL, load_config  # noqa: E402

PRODUCER = ORIGINAL / 'figures' / 'availability_native_timestep_percentage.py'
CROP_TOP = 125

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
ap.add_argument('--out-dir', default=None)
a = ap.parse_args()
cfg = load_config(a.config)
ds = cfg['datasets']
out = Path(a.out_dir) if a.out_dir else cfg['output_dir'] / 'figures' / 'fig02'
out.mkdir(parents=True, exist_ok=True)

sys.argv = [str(PRODUCER), '--cml-file', str(ds['raw']), '--radar-file', str(ds['radar']),
            '--gauge-file', str(ds['gauge']), '--temperature-file', str(ds['temperature']),
            '--humidity-file', str(ds['humidity']), '--out-dir', str(out)]
spec = importlib.util.spec_from_file_location('availability_producer', PRODUCER)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod      # needed by the producer's dataclass
spec.loader.exec_module(mod)
mod.main()

src = out / 'availability_native_timestep_percentage.png'
im = Image.open(src)
top = np.asarray(im.convert('RGB'))[:CROP_TOP]
cols = np.where((top < 250).any(axis=(0, 2)))[0]
print(f'[crop] non-white columns in the removed band: {cols.min() if cols.size else None}-'
      f'{cols.max() if cols.size else None} (year label only)')
dst = out / 'availability_native_timestep_percentage_v1.png'
im.crop((0, CROP_TOP, im.size[0], im.size[1])).save(dst, dpi=(600, 600))
print('[crop] saved', dst)
