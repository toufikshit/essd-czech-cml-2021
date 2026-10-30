"""Stage 5A: old vs new gauge validation on identical samples; population vs processing changes; extremes.

Inputs: the per-sample files written by s5a_01 (control_released_samples.npz, stage4_new_samples.npz) and,
for path length, the two CML files. A sample is (cml_id, sublink, gauge, 10-minute gauge timestamp).
Groups:
  common      sample valid in both products (identical CML, sublink, gauge, timestamp): processing change
  old_only    sample only in the released validation, split by why (CML removed / CML kept but sample lost)
  new_only    sample only in the new validation (CML new-only / CML common but sample gained)
Every metric (N, RMSE, MBE, PCC on all samples and on gauge-wet >= 0.1; confusion) is reported with and
without blackout-reconstructed samples (bit 8, new product only). Extremes are described, not filtered.
Outputs: VALIDATION_BEFORE_AFTER.csv, RAIN_DETECTION_METRICS.csv, RAINFALL_VALIDATION_METRICS.csv,
EXTREMES_IN_VALIDATION.csv and a printed summary.

Usage: python -I s5a_02_compare.py [--stage5a DIR]
"""
import argparse, json, math, os
import numpy as np
import pandas as pd
import netCDF4 as nc

ap = argparse.ArgumentParser()
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap.add_argument('--stage5a', default=H)
ap.add_argument('--old-cml', default=r'C:\Users\taouf\ESSD_RELEASE\extracted\full\ESSD_CZ_2021_HANDOFF\data\source_full_history\qc\CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc')
ap.add_argument('--new-cml', default=os.path.join(H, '..', 'stage3', 'QC_CORRECTED_CANDIDATE_REV1.nc'))
ap.add_argument('--stage3', default=os.path.join(H, '..', 'stage3'))
a = ap.parse_args()
O = np.load(os.path.join(a.stage5a, 'control', 'control_released_samples.npz'))
N = np.load(os.path.join(a.stage5a, 'new', 'stage4_new_samples.npz'))
assert np.array_equal(O['gauge_time_ns'], N['gauge_time_ns'])
oc = O['cml_ids'].astype(str); ncid = N['cml_ids'].astype(str)
union = sorted(set(oc) | set(ncid)); code = {c: i for i, c in enumerate(union)}
def keys(D, ids):
    cc = np.array([code[c] for c in ids], np.int64)[D['cml']]
    return ((cc * 2 + D['sub'].astype(np.int64)) * 1024 + D['gauge'].astype(np.int64)) * 65536 + D['tidx'].astype(np.int64)
ko, kn = keys(O, oc), keys(N, ncid)
assert len(np.unique(ko)) == len(ko) and len(np.unique(kn)) == len(kn)
_, io, inn = np.intersect1d(ko, kn, assume_unique=True, return_indices=True)
mo = np.ones(len(ko), bool); mo[io] = False
mn = np.ones(len(kn), bool); mn[inn] = False

def stats(x, y):
    if len(x) == 0:
        return dict(N=0, RMSE=np.nan, MBE=np.nan, PCC=np.nan)
    x = x.astype(float); y = y.astype(float); d = y - x
    return dict(N=int(len(x)), RMSE=float(np.sqrt(np.mean(d * d))), MBE=float(np.mean(d)),
                PCC=float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 and x.std() > 0 and y.std() > 0 else np.nan)
def conf(x, p):
    t = x >= 0.1
    TP = int((p & t).sum()); FP = int((p & ~t).sum()); FN = int((~p & t).sum()); TN = int((~p & ~t).sum())
    den = math.sqrt(float(TP + FP) * float(TP + FN) * float(TN + FP) * float(TN + FN))
    return dict(N=TP + FP + FN + TN, TP=TP, FP=FP, FN=FN, TN=TN, POD=TP / (TP + FN) if TP + FN else np.nan,
                FAR=FP / (TP + FP) if TP + FP else np.nan, CSI=TP / (TP + FP + FN) if TP + FP + FN else np.nan,
                MCC=(TP * TN - FP * FN) / den if den else np.nan, correct_negative=TN / (TN + FP) if TN + FP else np.nan)
rows_r, rows_d = [], []
def add(label, product, D, sel, bit8_note):
    x = D['x'][sel]; y = D['y'][sel]; p = D['pred'][sel]; w = x >= 0.1
    rows_r.append(dict(group=label, product=product, blackout=bit8_note, subset='all_period', **stats(x, y)))
    rows_r.append(dict(group=label, product=product, blackout=bit8_note, subset='gauge_wet_ge_0.1', **stats(x[w], y[w])))
    rows_d.append(dict(group=label, product=product, blackout=bit8_note, **conf(x, p)))

allo = np.ones(len(ko), bool); alln = np.ones(len(kn), bool); nb8 = ~N['bit8']
add('full population', 'released (control)', O, allo, 'n/a')
add('full population', 'Stage 4', N, alln, 'with bit 8')
add('full population', 'Stage 4', N, nb8, 'without bit 8')
co = np.zeros(len(ko), bool); co[io] = True
cn = np.zeros(len(kn), bool); cn[inn] = True
add('common samples (processing change)', 'released (control)', O, co, 'n/a')
add('common samples (processing change)', 'Stage 4', N, cn, 'with bit 8')
# without bit 8 on common samples: drop the same samples from both
b8c = N['bit8'][inn]
coo = np.zeros(len(ko), bool); coo[io[~b8c]] = True
cnn = np.zeros(len(kn), bool); cnn[inn[~b8c]] = True
add('common samples (processing change)', 'released (control)', O, coo, 'without Stage-4 bit-8 samples')
add('common samples (processing change)', 'Stage 4', N, cnn, 'without bit 8')
# population groups
new_set = set(ncid); old_set = set(oc)
o_cml_removed = np.array([oc[i] not in new_set for i in O['cml']])
n_cml_new = np.array([ncid[i] not in old_set for i in N['cml']])
add('released-only samples: CML not in Stage 4', 'released (control)', O, mo & o_cml_removed, 'n/a')
add('released-only samples: CML kept, sample lost', 'released (control)', O, mo & ~o_cml_removed, 'n/a')
add('Stage-4-only samples: new CML', 'Stage 4', N, mn & n_cml_new, 'with bit 8')
add('Stage-4-only samples: new CML', 'Stage 4', N, mn & n_cml_new & nb8, 'without bit 8')
add('Stage-4-only samples: CML kept, sample gained', 'Stage 4', N, mn & ~n_cml_new, 'with bit 8')
add('Stage-4-only samples: CML kept, sample gained', 'Stage 4', N, mn & ~n_cml_new & nb8, 'without bit 8')
R = pd.DataFrame(rows_r); Dd = pd.DataFrame(rows_d)
R.to_csv(os.path.join(a.stage5a, 'RAINFALL_VALIDATION_METRICS.csv'), index=False)
Dd.to_csv(os.path.join(a.stage5a, 'RAIN_DETECTION_METRICS.csv'), index=False)

# before/after table (headline)
pub = dict(N=1354953, RMSE=3.33, MBE=-0.93, PCC=0.59)
w = R[R.subset == 'gauge_wet_ge_0.1'].set_index(['group', 'product', 'blackout'])
ba = pd.DataFrame([
    dict(quantity='published (manuscript)', **pub),
    dict(quantity='control reproduction of the released product', **w.loc[('full population', 'released (control)', 'n/a')].to_dict()),
    dict(quantity='Stage 4, with blackout-reconstructed samples', **w.loc[('full population', 'Stage 4', 'with bit 8')].to_dict()),
    dict(quantity='Stage 4, without blackout-reconstructed samples', **w.loc[('full population', 'Stage 4', 'without bit 8')].to_dict()),
    dict(quantity='common samples: released', **w.loc[('common samples (processing change)', 'released (control)', 'n/a')].to_dict()),
    dict(quantity='common samples: Stage 4 with bit 8', **w.loc[('common samples (processing change)', 'Stage 4', 'with bit 8')].to_dict()),
    dict(quantity='common samples without bit 8: released', **w.loc[('common samples (processing change)', 'released (control)', 'without Stage-4 bit-8 samples')].to_dict()),
    dict(quantity='common samples without bit 8: Stage 4', **w.loc[('common samples (processing change)', 'Stage 4', 'without bit 8')].to_dict()),
])
ba.to_csv(os.path.join(a.stage5a, 'VALIDATION_BEFORE_AFTER.csv'), index=False)

# extremes (description only)
L = {}
for p in (a.old_cml, a.new_cml):
    with nc.Dataset(p) as d:
        for c, l in zip([str(x) for x in d['cml_id'][:]], np.asarray(d['length'][:], float)):
            L[c] = l
ext = []
for prod, D, ids, sel in (('released', O, oc, allo), ('Stage 4', N, ncid, alln)):
    y = D['y'][sel].astype(float); x = D['x'][sel].astype(float); e2 = (y - x) ** 2
    ln = np.array([L[c] for c in ids])[D['cml'][sel]]
    b8 = D['bit8'][sel]
    for thr in (30, 50, 100):
        m = y >= thr
        ext.append(dict(product=prod, cml_10min_rate_ge=thr, n_samples=int(m.sum()),
                        share_of_squared_error=float(e2[m].sum() / e2.sum()),
                        n_with_gauge_lt_1mmh=int((m & (x < 1)).sum()), n_short_link_lt_1km=int((m & (ln < 1)).sum()),
                        n_bit8=int((m & b8).sum()), median_gauge_mmh=float(np.median(x[m])) if m.any() else np.nan))
E = pd.DataFrame(ext); E.to_csv(os.path.join(a.stage5a, 'EXTREMES_IN_VALIDATION.csv'), index=False)

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20)
print('SAMPLES: released', len(ko), 'Stage 4', len(kn), 'common', len(io), 'released-only', int(mo.sum()), 'Stage-4-only', int(mn.sum()))
print('  released-only: CML removed', int((mo & o_cml_removed).sum()), '; CML kept', int((mo & ~o_cml_removed).sum()))
print('  Stage-4-only: new CML', int((mn & n_cml_new).sum()), '; CML kept', int((mn & ~n_cml_new).sum()))
print('  Stage 4 bit-8 samples:', int(N['bit8'].sum()), '; of which gauge-wet:', int((N['bit8'] & (N['x'] >= 0.1)).sum()))
print('\nCONTRIBUTING: released CMLs', len(set(oc[O['cml']])), 'sublinks', len(set(zip(O['cml'], O['sub']))), 'gauges', len(set(O['gauge'])),
      '| Stage 4 CMLs', len(set(ncid[N['cml']])), 'sublinks', len(set(zip(N['cml'], N['sub']))), 'gauges', len(set(N['gauge'])))
print('\nRAINFALL METRICS'); print(R.round(4).to_string(index=False))
print('\nDETECTION'); print(Dd.round(4).to_string(index=False))
print('\nBEFORE / AFTER (gauge-wet >= 0.1)'); print(ba.round(4).to_string(index=False))
print('\nEXTREMES'); print(E.round(4).to_string(index=False))
