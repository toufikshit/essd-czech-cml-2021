"""Stage 3 revision 1: exclusion tables, population reconciliation, comparison with the released QC and consistency checks.

Inputs (read-only): stage3/QC_SUBLINK_RESULTS.csv and QC_CORRECTED_CANDIDATE.nc (from s3_01),
stage2/RAW_APPROVED_POPULATION.csv, stage2/RAW_EXCLUDED_CMLS_ALL.csv, stage2/QC_PREDEFINED_MINUTE_MASKS.csv,
the full-history RAW (ids only) and the released QC (ids and detection limits only).
Outputs: QC_EXCLUDED_CMLS.csv, QC_EXCLUDED_SUBLINKS.csv, QC_BEFORE_AFTER_COMPARISON.csv,
QC_POPULATION_RECONCILIATION.csv, and a printed summary (saved by the caller).

Usage: python -I s3_02_tables_checks.py [--stage3 DIR]
"""
import argparse, os
import netCDF4 as nc
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
here = os.path.dirname(os.path.abspath(__file__))
ap.add_argument('--stage3', default=os.path.dirname(here))
ap.add_argument('--raw-nc', default=r'C:\Users\taouf\ESSD_RELEASE\extracted\full\ESSD_CZ_2021_HANDOFF\data\source_full_history\raw\CZ_CML_RAW_1MIN_20200101_20211231.nc')
ap.add_argument('--old-qc', default=r'C:\Users\taouf\ESSD_RELEASE\zenodo_staging\primary\CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc')
a = ap.parse_args()
d = a.stage3
s2 = os.path.join(d, '..', 'stage2')
res = pd.read_csv(os.path.join(d, 'QC_SUBLINK_RESULTS_REV1.csv'))
pop = pd.read_csv(os.path.join(s2, 'RAW_APPROVED_POPULATION.csv'))
excl_all = pd.read_csv(os.path.join(s2, 'RAW_EXCLUDED_CMLS_ALL.csv'))
masks = pd.read_csv(os.path.join(s2, 'QC_PREDEFINED_MINUTE_MASKS.csv'))

with nc.Dataset(a.raw_nc) as r:
    raw_ids = set(str(x) for x in r['cml_id'][:])
with nc.Dataset(a.old_qc) as o:
    old_ids = [str(x) for x in o['cml_id'][:]]
    old_dl = dict(zip(old_ids, np.ma.filled(o['detection_limit_mmhr'][:], np.nan)))
with nc.Dataset(os.path.join(d, 'QC_CORRECTED_CANDIDATE_REV1.nc')) as q:
    new_ids = [str(x) for x in q['cml_id'][:]]
    st_nc = np.asarray(q['qc_sublink_status'][:])
old_set, new_set = set(old_ids), set(new_ids)
approved = set(pop.cml_id)
nofreq = set(pop[pop.no_frequency_metadata == 1].cml_id)
assessed = set(res.cml_id)
excluded_before = set(excl_all.cml_id)                     # 13 (Stage 1B) + 6 (R1)
empty255 = old_set - raw_ids - excluded_before             # QC-only CMLs that are not the 13

# ---------------- checks ----------------
chk = {}
chk['approved RAW CMLs'] = len(approved)
chk['approved RAW subset of full-history RAW ids'] = approved <= raw_ids
chk['assessed = approved - no-frequency'] = assessed == approved - nofreq and len(assessed) == 2566
chk['new QC subset of approved RAW (strict)'] = new_set < approved
chk['new QC ids unique'] = len(new_ids) == len(new_set)
chk['no 13/R1 CML in new QC'] = not (new_set & excluded_before)
chk['no empty-255 CML in new QC'] = not (new_set & empty255)
chk['no no-frequency CML in new QC'] = not (new_set & nofreq)
ret_cml = set(res[res.cml_retained == 1].cml_id)
chk['NetCDF CMLs = retained CMLs of the table'] = new_set == ret_cml
chk['NetCDF sublink status = table status'] = all(
    st_nc[s, i] == res[(res.cml_id == c) & (res.sublink == s)].status_code.iloc[0] for i, c in enumerate(new_ids) for s in (0, 1))
chk['every retained CML has >= 1 retained sublink'] = bool((res[res.cml_retained == 1].groupby('cml_id').status_code.apply(lambda x: (x == 0).any())).all())
chk['empty-255 count'] = len(empty255)

# ---------------- exclusion tables ----------------
res['excluded'] = (res.status_code != 0).astype(int)
sub_ex = res[res.status_code != 0].copy()
sub_ex['scope'] = np.where(sub_ex.cml_retained == 1, 'sublink (CML retained)', 'sublink of excluded CML')
cols = ['cml_id', 'sublink', 'scope', 'status', 'p5_nonoperational', 'detection_limit_mmhr', 'frac_long', 'n_valid_long',
        'frac_short', 'n_valid_short', 'n_tl_observed_after_masks', 'frequency_MHz', 'length_km']
sub_ex[cols].to_csv(os.path.join(d, 'QC_EXCLUDED_SUBLINKS_REV1.csv'), index=False)

def cml_reason(g):
    return ' | '.join(f's{int(s)}:{st}' for s, st in zip(g.sublink, g.status))
cml_ex = (res[res.cml_retained == 0].groupby('cml_id', sort=False).apply(cml_reason).rename('sublink_reasons').reset_index())
def primary(txt):
    parts = [p.split(':', 1)[1] for p in txt.split(' | ')]
    s = set(parts)
    if s <= {'excluded_no_observed_total_loss_in_period', 'excluded_p5_tx_zero_nonoperational'}:
        return 'no usable observations in the period' if 'excluded_p5_tx_zero_nonoperational' not in s else 'no usable observations (P5 direction and no data)'
    if s & {'excluded_detection_limit', 'excluded_fluctuation_and_detection_limit'} and not s & {'excluded_fluctuation'}:
        return 'detection limit (all sublinks with data)' if 'excluded_fluctuation_and_detection_limit' not in s else 'detection limit and fluctuation'
    if s & {'excluded_fluctuation'} and not s & {'excluded_detection_limit', 'excluded_fluctuation_and_detection_limit'}:
        return 'fluctuation (all sublinks with data)'
    return 'mixed: ' + ', '.join(sorted(s))
cml_ex['qc_stage'] = 'QC rules (Stage 3)'
cml_ex['primary_reason'] = cml_ex.sublink_reasons.map(primary)
pre = pd.DataFrame({'cml_id': sorted(nofreq), 'sublink_reasons': '', 'qc_stage': 'before QC (Stage 2, P3)',
                    'primary_reason': 'no frequency metadata (not rainfall-ready)'})
cml_ex = pd.concat([pre, cml_ex], ignore_index=True)
cml_ex['in_released_qc'] = cml_ex.cml_id.isin(old_set).astype(int)
cml_ex.to_csv(os.path.join(d, 'QC_EXCLUDED_CMLS_REV1.csv'), index=False)

# ---------------- reconciliation ----------------
n_sub_ret = int((res.status_code == 0).sum())
rec = pd.DataFrame([
    ('Approved RAW input', len(approved), 2 * len(approved)),
    ('Excluded before QC: no frequency (P3)', len(nofreq), 2 * len(nofreq)),
    ('Assessed by QC', len(assessed), len(res)),
    ('Excluded by QC (CMLs without a passing sublink)', len(assessed) - len(new_set), int(2 * (len(assessed) - len(new_set)))),
    ('Excluded sublinks inside retained CMLs', 0, int(((res.cml_retained == 1) & (res.status_code != 0)).sum())),
    ('Retained after QC', len(new_set), n_sub_ret),
], columns=['population', 'cmls', 'sublinks'])
rec.to_csv(os.path.join(d, 'QC_POPULATION_RECONCILIATION_REV1.csv'), index=False)
chk['reconciliation CMLs'] = len(approved) == len(nofreq) + (len(assessed) - len(new_set)) + len(new_set)
chk['reconciliation sublinks'] = 2 * len(approved) == 2 * len(nofreq) + int(2 * (len(assessed) - len(new_set))) + int(((res.cml_retained == 1) & (res.status_code != 0)).sum()) + n_sub_ret

# ---------------- comparison with released QC ----------------
allids = sorted(old_set | approved | excluded_before)
dl_new = res.pivot(index='cml_id', columns='sublink', values='detection_limit_mmhr')
st_new = res.pivot(index='cml_id', columns='sublink', values='status')
rows = []
e13 = set(excl_all[excl_all.stage.astype(str) == '1B'].cml_id); r1 = set(excl_all[excl_all.stage.astype(str) == '2'].cml_id)
for c in allids:
    if c in e13: why = 'not in approved RAW: author exclusion of the 13 TX = 0 CMLs (Stage 1B)'
    elif c in r1: why = 'not in approved RAW: R1 non-operational (Stage 2)'
    elif c in empty255: why = 'not in approved RAW: empty metadata-only CML'
    elif c in nofreq: why = 'excluded before QC: no frequency'
    elif c in new_set: why = 'retained'
    else: why = 'excluded by QC: ' + primary(cml_ex.set_index('cml_id').sublink_reasons.get(c, ''))
    rows.append(dict(cml_id=c, in_released_qc=int(c in old_set), in_approved_raw=int(c in approved), in_new_qc=int(c in new_set),
                     new_status=why, released_detection_limit_mmhr_2003_cml=old_dl.get(c, np.nan),
                     new_detection_limit_mmhr_s0=dl_new.loc[c, 0] if c in dl_new.index else np.nan,
                     new_detection_limit_mmhr_s1=dl_new.loc[c, 1] if c in dl_new.index else np.nan,
                     new_status_s0=st_new.loc[c, 0] if c in st_new.index else '', new_status_s1=st_new.loc[c, 1] if c in st_new.index else ''))
cmp_ = pd.DataFrame(rows)
cmp_['change'] = np.select([(cmp_.in_released_qc == 1) & (cmp_.in_new_qc == 1), (cmp_.in_released_qc == 1) & (cmp_.in_new_qc == 0),
                            (cmp_.in_released_qc == 0) & (cmp_.in_new_qc == 1)], ['in both', 'released only', 'new only'], 'in neither')
cmp_.to_csv(os.path.join(d, 'QC_BEFORE_AFTER_COMPARISON_REV1.csv'), index=False)

# ---------------- summary ----------------
pd.set_option('display.width', 220)
print('CHECKS'); [print(f'  {k}: {v}') for k, v in chk.items()]
print('\nRECONCILIATION'); print(rec.to_string(index=False))
print('\nSUBLINK STATUS (assessed)'); print(res.status.value_counts().to_string())
print('\nSUBLINK STATUS in retained CMLs'); print(res[res.cml_retained == 1].status.value_counts().to_string())
print('\nCML EXCLUSION REASONS'); print(cml_ex.primary_reason.value_counts().to_string())
print('\nCOMPARISON WITH RELEASED QC (2219)'); print(cmp_.groupby(['change', 'new_status']).size().to_string())
rr = res[res.cml_retained == 1]
ret = res[res.status_code == 0]
print('\nPROCESSING COUNTS (retained sublinks | all assessed sublinks)')
for k in ['n_p6_masked', 'n_tx_invalid', 'n_rx_invalid', 'n_plateau_masked', 'n_blackout_gaps', 'n_blackout_minutes', 'n_blackout_rejected_gaps', 'n_blackout_rejected_minutes', 'n_interp_gaps', 'n_interp_minutes', 'n_tl_observed_after_masks']:
    print(f'  {k}: {int(ret[k].sum()):,} | {int(res[k].sum()):,}')
print(f'  sublinks with >= 1 reconstructed blackout minute (retained): {int((ret.n_blackout_minutes > 0).sum())}')
print(f'  P6 minutes in window: {int(((masks.time_index >= 656640) & (masks.time_index < 964800)).sum())} of {len(masks)}')
print('\nDETECTION LIMIT (ITU 2005, per sublink, assessed sublinks with frequency):')
print(res.detection_limit_mmhr.describe().round(3).to_string())
both = cmp_[(cmp_.change == 'in both')]
print(f'  sublinks with R_min >= 2: {int((res.detection_limit_mmhr >= 2).sum())}; NaN: {int(res.detection_limit_mmhr.isna().sum())}')
print(f'  CMLs in both: released 2003 CML value vs new per-sublink min: median ratio '
      f'{np.nanmedian(np.fmin(both.new_detection_limit_mmhr_s0, both.new_detection_limit_mmhr_s1) / both.released_detection_limit_mmhr_2003_cml):.3f}')
print('\nE-BAND (>= 70 GHz) CMLs: new QC', int(res[(res.cml_retained == 1) & (res.frequency_MHz >= 70000)].cml_id.nunique()))
print('ALL CHECKS PASS' if all(v for k, v in chk.items() if isinstance(v, bool)) else 'CHECK FAILURE')

# ---------------- revision 1 vs first candidate ----------------
v0 = pd.read_csv(os.path.join(d, 'QC_SUBLINK_RESULTS.csv'))
m = v0.merge(res, on=['cml_id', 'sublink'], suffixes=('_v0', '_r1'))
assert len(m) == len(v0) == len(res)
print('\nREVISION 1 vs FIRST CANDIDATE (5132 sublinks)')
print('  status changed:', int((m.status_code_v0 != m.status_code_r1).sum()))
print('  retained CML set identical:', set(v0[v0.cml_retained == 1].cml_id) == set(res[res.cml_retained == 1].cml_id))
for k in ['n_blackout_gaps', 'n_blackout_minutes', 'n_interp_gaps', 'n_interp_minutes', 'n_plateau_masked', 'n_tl_observed_after_masks']:
    print(f'  {k}: {int(m[k + "_v0"].sum()):,} -> {int(m[k + "_r1"].sum()):,}  (retained sublinks: {int(m[m.status_code_r1 == 0][k + "_v0"].sum()):,} -> {int(m[m.status_code_r1 == 0][k + "_r1"].sum()):,})')
print('  removed blackout minutes = rejected minutes:', int(m.n_blackout_minutes_v0.sum() - m.n_blackout_minutes_r1.sum()) == int(m.n_blackout_rejected_minutes.sum()))
for k in ['frac_long', 'frac_short', 'detection_limit_mmhr']:
    dd = (m[k + '_v0'] - m[k + '_r1']).abs()
    print(f'  max |{k} change|: {np.nanmax(dd.values):.3g}; sublinks changed: {int((dd > 0).sum())}')
