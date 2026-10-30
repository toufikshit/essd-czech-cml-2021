"""Water vapour density (WVD) from E-band CMLs: validation against CHMI stations (Table S13) and Figure 8.

Runs the unchanged scripts of original/code/wvd/:
  1. 33_v6_binned_validation_hexbins.py: pairs the V6 WVD series with the nearest CHMI temperature and humidity
     station (<= 5 km), Graf wet/dry branch, quantization bins; helper 21_validate_wvd_v6_vs_original.py;
  2. plot_wvp_support_hexbins_no_titles.py: Figure 8 and its metrics for E-band links with L >= 4 km and
     delta-rho-q <= 0.20 g m-3.
Inputs: the three WVD supporting files of the deposit (extras.wvd_v6_nc, extras.wvd_v6_diagnostics_csv,
extras.wvd_quantization_csv), the RAW, temperature and humidity files.
Outputs in <output_dir>/wvd: WVD_V6_QUANTIZATION_BINNED_VALIDATION/ and WVP_SUPPORT_HEXBINS/
(wvp_support_all_E_band.png = Figure 8, wvp_support_metrics.csv).
The V6 retrieval itself (19_wvd_eband_progressive_baseline_v6.py, run in parallel by 20_run_wvd_v6_parallel_hpc.py,
and the quantization of 31_run_v6_quantization_parallel.py / 30c_v6_quantization_chunk.py) produced the three
supporting files; see docs/METHODS.md.
Usage: python scripts/wvd_validation.py [--config config/user_paths.json]
"""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from essd_cz.config import ORIGINAL, extra, load_config  # noqa: E402

WVD = ORIGINAL / 'wvd'
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
ap.add_argument('--dpi', type=int, default=300)
a = ap.parse_args()
cfg = load_config(a.config)
ds = cfg['datasets']
root = cfg['output_dir'] / 'wvd'
val = root / 'WVD_V6_QUANTIZATION_BINNED_VALIDATION'
val.mkdir(parents=True, exist_ok=True)

cmd = [sys.executable, str(WVD / '33_v6_binned_validation_hexbins.py'), '--root', str(root),
       '--v6-nc', str(extra(a.config, 'wvd_v6_nc')),
       '--v6-diagnostics-csv', str(extra(a.config, 'wvd_v6_diagnostics_csv')),
       '--quant-csv', str(extra(a.config, 'wvd_quantization_csv')),
       '--raw-cml-nc', str(ds['raw']), '--humidity-nc', str(ds['humidity']), '--temperature-nc', str(ds['temperature']),
       '--validation-helper-script', str(WVD / '21_validate_wvd_v6_vs_original.py'),
       '--station-assignments-csv', str(val / 'station_assignments.csv'), '--out-dir', str(val),
       '--branch', 'GRAF', '--max-station-distance-km', '5', '--dpi', str(a.dpi)]
print('RUN:', ' '.join(cmd), flush=True)
subprocess.run(cmd, check=True)

cmd = [sys.executable, str(WVD / 'plot_wvp_support_hexbins_no_titles.py'), '--root', str(root),
       '--out-dir', str(root / 'WVP_SUPPORT_HEXBINS'), '--dpi', str(a.dpi), '--gridsize', '55', '--plot-max', '35']
print('RUN:', ' '.join(cmd), flush=True)
subprocess.run(cmd, check=True)
