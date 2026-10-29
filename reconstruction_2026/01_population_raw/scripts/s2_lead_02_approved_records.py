"""Stage 2 (lead): reproducible records of the approved Stage 2 decisions (author approval 2026-10-28).

Writes (no data file is modified):
  RAW_APPROVED_POPULATION.csv      2574 CMLs x 2 sublinks: approved raw population with documentation flags
  RAW_EXCLUDED_CMLS_ALL.csv        every CML excluded from the publication RAW/QC, with stage, rule and reason
  QC_PREDEFINED_SUBLINK_EXCLUSIONS.csv  the 15 TX = 0 sublinks (P5) that may not contribute QC observations
  QC_PREDEFINED_MINUTE_MASKS.csv   exact minutes of the P6 mixed-placeholder values (TX and RX masked in QC)
Inputs: RAW_LINK_AUDIT.csv, RAW_PROPOSED_EXCLUSIONS.csv, ../stage1b/STAGE1B_EXCLUDED_13_CMLS.csv and, for the
exact P6 minutes only, the full-history RAW NetCDF (read-only, only the P6 CMLs).

Usage: python -I s2_lead_02_approved_records.py [--stage2 DIR] [--raw-nc PATH]
"""
import argparse, os
import numpy as np
import pandas as pd
import netCDF4 as nc

ap = argparse.ArgumentParser()
ap.add_argument('--stage2', default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ap.add_argument('--raw-nc', default=r'C:\Users\taouf\ESSD_RELEASE\extracted\full\ESSD_CZ_2021_HANDOFF\data\source_full_history\raw\CZ_CML_RAW_1MIN_20200101_20211231.nc')
a = ap.parse_args()
d = a.stage2
au = pd.read_csv(os.path.join(d, 'RAW_LINK_AUDIT.csv'))
ex = pd.read_csv(os.path.join(d, 'RAW_PROPOSED_EXCLUSIONS.csv'))
e13 = pd.read_csv(os.path.join(d, '..', 'stage1b', 'STAGE1B_EXCLUDED_13_CMLS.csv'))

r1 = sorted(ex[ex.proposal.isin(['P2a', 'P2b', 'P2c'])].cml_id.unique())
assert len(r1) == 6
nofreq = sorted(ex[ex.proposal == 'P3'].cml_id.unique())
p5 = ex[(ex.proposal == 'P5') & ~ex.cml_id.isin(r1)][['cml_id', 'sublink']].copy()
p5['sublink'] = p5.sublink.astype(int)
assert len(nofreq) == 8 and len(p5) == 15

# approved population
pop = au[~au.cml_id.isin(r1)].copy()
assert pop.cml_id.nunique() == 2574 and len(pop) == 5148
pop['no_frequency_metadata'] = pop.cml_id.isin(nofreq).astype(int)
key = set(zip(p5.cml_id, p5.sublink))
pop['sublink_tx_zero_nonoperational'] = [int((c, s) in key) for c, s in zip(pop.cml_id, pop.sublink)]
pop['sublink_empty'] = (pop.sublink_class == 5).astype(int)
pop['rainfall_ready_eligible_cml'] = (~pop.cml_id.isin(nofreq)).astype(int)
cols = ['cml_id', 'sublink', 'frequency_MHz', 'length_km', 'polarization', 'in_released_qc', 'sublink_class', 'cml_class',
        'no_frequency_metadata', 'sublink_tx_zero_nonoperational', 'sublink_empty', 'rainfall_ready_eligible_cml']
pop[cols].to_csv(os.path.join(d, 'RAW_APPROVED_POPULATION.csv'), index=False)

# every exclusion from the publication datasets
rows = [dict(cml_id=c, stage='1B', rule='author dataset decision', reason='TX = 0 dBm at every valid minute; non-operational',
             in_released_raw=0, in_released_qc=1) for c in e13.cml_id]
for c in r1:
    p = ex[ex.cml_id == c].iloc[0]
    rows.append(dict(cml_id=c, stage='2', rule='R1 (' + p.proposal + ')', reason=p.reason, in_released_raw=1,
                     in_released_qc=int(au[au.cml_id == c].in_released_qc.iloc[0])))
rec = pd.DataFrame(rows)
rec.to_csv(os.path.join(d, 'RAW_EXCLUDED_CMLS_ALL.csv'), index=False)

p5.assign(rule='P5 (author 2026-10-28)', reason='TX = 0 dBm at every valid minute; direction non-operational; opposite direction retained',
          raw_action='kept, documented', qc_action='no valid observations; sublink excluded').to_csv(
    os.path.join(d, 'QC_PREDEFINED_SUBLINK_EXCLUSIONS.csv'), index=False)

# exact P6 minutes (TX >= 60 dBm) in the approved population
ds = nc.Dataset(a.raw_nc)
ids = np.array([str(x) for x in ds['cml_id'][:]])
t0 = pd.Timestamp('2020-01-01')
p6 = ex[(ex.proposal == 'P6') & ~ex.cml_id.isin(r1)]
out = []
for c, s in zip(p6.cml_id, p6.sublink.astype(int)):
    j = int(np.flatnonzero(ids == c)[0])
    tx = np.ma.filled(ds['tx'][:, s, j].astype('f8'), np.nan)
    rx = np.ma.filled(ds['rx'][:, s, j].astype('f8'), np.nan)
    for i in np.flatnonzero(tx >= 60):
        nb = lambda v: v[i - 1] if i > 0 else np.nan
        out.append(dict(cml_id=c, sublink=s, time_utc=str(t0 + pd.Timedelta(minutes=int(i))), time_index=int(i),
                        tx_raw_dBm=tx[i], rx_raw_dBm=rx[i], tx_prev=nb(tx), rx_prev=nb(rx),
                        tx_next=tx[i + 1], rx_next=rx[i + 1],
                        rule='P6 (author 2026-10-28)', qc_action='mask TX and RX of this minute', raw_action='kept, documented'))
m = pd.DataFrame(out)
m.to_csv(os.path.join(d, 'QC_PREDEFINED_MINUTE_MASKS.csv'), index=False)
print('approved RAW:', pop.cml_id.nunique(), 'CMLs /', len(pop), 'sublinks; no-frequency CMLs:', len(nofreq),
      '; P5 sublinks:', len(p5), '; rainfall-ready eligible CMLs:', int(pop.drop_duplicates('cml_id').rainfall_ready_eligible_cml.sum()))
print('excluded from publication RAW/QC:', len(rec), '(13 at 1B + 6 by R1); of R1 in released QC:', int(rec[rec.stage == '2'].in_released_qc.sum()))
print('P6 minutes in approved population:', len(m), 'in', m.groupby(['cml_id', 'sublink']).ngroups, 'sublinks')
print(m[['cml_id', 'sublink', 'time_utc', 'tx_raw_dBm', 'rx_raw_dBm', 'tx_prev', 'rx_prev', 'tx_next', 'rx_next']].to_string(index=False))
