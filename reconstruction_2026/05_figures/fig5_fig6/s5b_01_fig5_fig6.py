"""Stage 5B: Figs. 5 (contingency) and 6 (gauge-wet hexbin) with the ORIGINAL plotting functions.

`hexbin()`, `plot_confusion()`, `style()` and `save()` are imported unchanged from the published producer
essd-czech-cml-2021/original/code/figures/rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py
(sha256 7358e242...). The plotted samples are the joint validation samples saved by Stage 5A
(stage5a/<run>/<tag>_samples.npz), which reproduce the published validation exactly (control) and give the
approved Case A result (Stage 4, with blackout-reconstructed samples). Fig. 6 = hexbin of gauge-wet samples
(gauge >= 0.1 mm/h); Fig. 5 = confusion of all joint samples. Nothing else is computed.

Usage: python -I s5b_01_fig5_fig6.py --samples NPZ --out-dir DIR
"""
import argparse, importlib.util, json, os, sys
from pathlib import Path
sys.dont_write_bytecode = True  # never write __pycache__ into the repository
import numpy as np

PRODUCER = Path(r'C:\Users\taouf\ESSD_RELEASE\essd-czech-cml-2021\original\code\figures\rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py')
ap = argparse.ArgumentParser(); ap.add_argument('--samples', required=True); ap.add_argument('--out-dir', required=True)
a = ap.parse_args()
spec = importlib.util.spec_from_file_location('essd_final_graf_figures', PRODUCER)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
S = np.load(a.samples)
x = S['x']; y = S['y']; p = S['pred']
wet = x >= 0.1
TP = int((p & wet).sum()); FP = int((p & ~wet).sum()); FN = int((~p & wet).sum()); TN = int((~p & ~wet).sum())
stats = m.hexbin(x[wet], y[wet], out / 'FINAL_GRAF_hexbin_gauge_WETONLY_202104_202110')
m.plot_confusion(TP, FP, FN, TN, out / 'FINAL_GRAF_confusion_GAUGE_same_population_202104_202110')
res = {'samples': str(Path(a.samples).resolve()), 'fig6_gauge_wet_stats': stats,
       'fig5_confusion': {'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN, 'POD_pct': 100 * TP / (TP + FN), 'false_alarm_pct_of_dry': 100 * FP / (FP + TN),
                          'miss_pct': 100 * FN / (TP + FN), 'correct_negative_pct': 100 * TN / (FP + TN)}}
(out / 'figure_values.json').write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
