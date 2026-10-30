"""Figures 5 (gauge detection, contingency) and 6 (CML against gauge rain rate for gauge-wet intervals).

Uses hexbin() and plot_confusion() of the unchanged producer
original/code/figures/rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py on the joint
validation samples written by scripts/gauge_validation.py. Fig. 6: gauge-wet samples (gauge >= 0.1 mm/h);
Fig. 5: all joint samples.
Outputs in <output_dir>/figures/fig05_fig06: the two figures (PNG/PDF) and figure_values.json.
Usage: python scripts/figures/fig05_fig06_validation.py [--config config/user_paths.json] [--samples NPZ]
"""
import argparse, importlib.util, json, sys
from pathlib import Path
sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from essd_cz.config import ORIGINAL, load_config  # noqa: E402

PRODUCER = ORIGINAL / 'figures' / 'rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py'
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
ap.add_argument('--samples', default=None, help='default: <output_dir>/validation/graf_samples.npz')
ap.add_argument('--out-dir', default=None)
a = ap.parse_args()
cfg = load_config(a.config)
samples = Path(a.samples) if a.samples else cfg['output_dir'] / 'validation' / 'graf_samples.npz'
if not samples.is_file():
    raise SystemExit(f'{samples} not found (run scripts/gauge_validation.py first)')
spec = importlib.util.spec_from_file_location('essd_final_graf_figures', PRODUCER)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
out = Path(a.out_dir) if a.out_dir else cfg['output_dir'] / 'figures' / 'fig05_fig06'
out.mkdir(parents=True, exist_ok=True)
S = np.load(samples)
x = S['x']; y = S['y']; p = S['pred']
wet = x >= 0.1
TP = int((p & wet).sum()); FP = int((p & ~wet).sum()); FN = int((~p & wet).sum()); TN = int((~p & ~wet).sum())
stats = m.hexbin(x[wet], y[wet], out / 'FINAL_GRAF_hexbin_gauge_WETONLY_202104_202110')
m.plot_confusion(TP, FP, FN, TN, out / 'FINAL_GRAF_confusion_GAUGE_same_population_202104_202110')
res = {'samples': str(samples.resolve()), 'fig6_gauge_wet_stats': stats,
       'fig5_confusion': {'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN, 'POD_pct': 100 * TP / (TP + FN),
                          'false_alarm_pct_of_dry': 100 * FP / (FP + TN),
                          'miss_pct': 100 * FN / (TP + FN), 'correct_negative_pct': 100 * TN / (FP + TN)}}
(out / 'figure_values.json').write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
