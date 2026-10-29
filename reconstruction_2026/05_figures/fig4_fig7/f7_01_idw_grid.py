#!/usr/bin/env python3
"""Stage 5B / Agent 3b - Fig. 7 CML IDW grid (unchanged connected-workflow gridder).

Calls the UNCHANGED repository function
  essd-czech-cml-2021/src/essd_cz/connected.py :: grid_final_rain
which itself loads the UNCHANGED original IDW functions build_static_idw and
idw_field from original/code/figures/02_blettner_aproct_dataset_figures_grid.py
(AST extraction, as in the published HPC run).  Nothing is re-implemented.

Parameters exposed exactly as in notebooks/03_Grid_CML_Rainfall.ipynb:
START/END grid times, resolution, nnear, power, max_dist, pad, extent and
max_grid_cells (keyword arguments of grid_final_rain; max_grid_cells only
guards the grid size and does not change any value).

Typical calls:
  # A. exact re-run of the HPC evidence grid (0.02 deg demo) for verification
  python -I f7_01_idw_grid.py --rain <full-history 5-min rain> --qc <full-history QC> \
      --resolution 0.02 --out ..\\work\\grid_original_0p02.nc
  # B. published-resolution grid (0.01 deg, 'original gridding parameter')
  python -I f7_01_idw_grid.py ... --resolution 0.01 --max-grid-cells 1000000 --out ...
  # C. Stage 4 grid on the IDENTICAL node lattice: pass --extent of grid B
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # never write __pycache__ into the read-only repository

REPO = Path(r'C:\Users\taouf\ESSD_RELEASE\essd-czech-cml-2021')


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--repo', type=Path, default=REPO)
    ap.add_argument('--rain', type=Path, required=True)
    ap.add_argument('--qc', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--start', default='2021-07-13 15:05:00')
    ap.add_argument('--end', default='2021-07-13 23:05:00')
    ap.add_argument('--resolution', type=float, default=0.02)
    ap.add_argument('--nnear', type=int, default=15)
    ap.add_argument('--power', type=float, default=2.0)
    ap.add_argument('--max-dist', type=float, default=0.3)
    ap.add_argument('--pad', type=float, default=0.05)
    ap.add_argument('--extent', type=float, nargs=4, default=None)
    ap.add_argument('--max-grid-cells', type=int, default=100000)
    a = ap.parse_args()
    sys.path.insert(0, str(a.repo / 'src'))
    from essd_cz import connected  # noqa: E402
    print('connected.py sha256', sha256(a.repo / 'src' / 'essd_cz' / 'connected.py'))
    print('02_blettner (IDW source) sha256',
          sha256(a.repo / 'original' / 'code' / 'figures' / '02_blettner_aproct_dataset_figures_grid.py'))
    print('rain', a.rain); print('qc', a.qc)
    out = connected.grid_final_rain(a.rain, a.qc, a.out, a.start, a.end,
                                    resolution=a.resolution, pad=a.pad, nnear=a.nnear,
                                    power=a.power, max_dist=a.max_dist, extent=a.extent,
                                    max_grid_cells=a.max_grid_cells)
    print('WROTE', out, 'sha256', sha256(out))


if __name__ == '__main__':
    main()
