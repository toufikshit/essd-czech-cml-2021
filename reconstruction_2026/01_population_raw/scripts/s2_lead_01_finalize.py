"""Stage 2 (lead): final Stage 2 tables from the verified agent drafts (no data file is read).

Inputs: agents/RAW_LINK_AUDIT_draft.csv (Agent 1, 5160 rows; Agent 4: 0 mismatches),
        agents/RAW_PROPOSED_EXCLUSIONS_draft.csv (Agent 1, 54 rows).
Outputs: RAW_LINK_AUDIT.csv, RAW_PROPOSED_EXCLUSIONS.csv, RAW_PROPOSED_ADDITIONS.csv.
Changes to the drafts (from Agent 4's review and the lead's reading of the P6 values):
  * audit: class-2 wording corrected; flag for sublinks whose RX minutes are > 50 % without TX;
  * P2 split into P2a (TX placeholder links) and P2b (one-sided records), with the observations they hold;
  * P2c added: 41920_323338_4235 (TX = 0 on both sublinks except one placeholder-mixed minute);
  * P6 recast as mixed-placeholder minutes (TX and RX of that minute), proposed as a QC record rule;
  * every row gets the Agent 4 verdict and the lead recommendation.
Nothing is executed on the data.

Usage: python -I s2_lead_01_finalize.py [--stage2 DIR]
"""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--stage2', default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
a = ap.parse_args()
d = a.stage2
au = pd.read_csv(os.path.join(d, 'agents', 'RAW_LINK_AUDIT_draft.csv'))
ex = pd.read_csv(os.path.join(d, 'agents', 'RAW_PROPOSED_EXCLUSIONS_draft.csv'))
assert len(au) == 5160 and au.cml_id.nunique() == 2580 and not au.cml_in_13.any()

# ---- audit table ----
au['sublink_flags'] = au.sublink_flags.str.replace('2:no both-finite minute in 2021-04-01..10-31',
                                                   '2:no paired TX-RX minute in 2021-04-01..10-31 (data end before 2021-01-20)', regex=False)
au['rx_only_gt50pct_flag'] = ((au.rx_finite > 0) & (au.rx_only / au.rx_finite.where(au.rx_finite > 0) > 0.5)).astype(int)
au['verified_by_agent4'] = 1
au.to_csv(os.path.join(d, 'RAW_LINK_AUDIT.csv'), index=False)
sub = au.set_index(['cml_id', 'sublink'])

# ---- proposals ----
ex['agent4_verdict'] = ''
ex['lead_recommendation'] = ''
p2a = {'10004_11582_3122', '27040_27099_3747'}
m = ex.proposal == 'P2'
for i in ex[m].index:
    c = ex.at[i, 'cml_id']
    g = au[au.cml_id == c]
    ex.at[i, 'proposal'] = 'P2a' if c in p2a else 'P2b'
    ex.at[i, 'obs_minutes_affected'] = int(g.tx_finite.sum() + g.rx_finite.sum())  # single-quantity values lost
    ex.at[i, 'agent4_verdict'] = ('confirmed; TX is a placeholder at every minute (NaN in RAW); consistent with the 13'
                                  if c in p2a else 'confirmed; ordinary one-sided record, could be kept with a flag')
    ex.at[i, 'lead_recommendation'] = 'exclude CML (non-operational rule R1: no paired TX-RX minute)'
ex.loc[ex.proposal == 'P3', 'agent4_verdict'] = 'confirmed; not a raw-validity reason (observations valid), prefer keep and flag'
ex.loc[ex.proposal == 'P3', 'lead_recommendation'] = 'keep in raw with metadata flag; QC excludes (no detection limit without frequency)'
ex.loc[ex.proposal == 'P4', 'agent4_verdict'] = 'confirmed; sublink dimension is fixed at 2, so no data change'
ex.loc[ex.proposal == 'P4', 'lead_recommendation'] = 'no data change; document the empty sublink'
ex.loc[ex.proposal == 'P5', 'agent4_verdict'] = 'confirmed; author decision, not equivalent to the 13 (other sublink operates)'
ex.loc[ex.proposal == 'P5', 'lead_recommendation'] = 'keep in raw; decide the QC rule for TX = 0 sublinks in Stage 3'
r = sub.loc[('41920_323338_4235', 0)]
p2c = dict(proposal='P2c', scope='whole CML', cml_id='41920_323338_4235', sublink='both',
           reason='TX = 0 dBm at every valid minute on both sublinks, except one placeholder-mixed minute (127.5 dBm) on sublink 0',
           evidence=f'[F] sl0: tx_zero {int(r.tx_zero_minutes)} of tx_finite {int(r.tx_finite)}; the other minute is the P6 value 127.5 = (255+0)/2; sl1: TX = 0 at all 41534 minutes',
           decision_type='author decision: same criterion as the 13 (CML level)',
           cmls_before=2575, sublinks_before=5150, cmls_after=2574, sublinks_after=5148,
           obs_minutes_affected=int(au[au.cml_id == '41920_323338_4235'].both_finite.sum()),
           agent4_verdict='raised by Agent 4: meets the 13 criterion once the P6 value is set aside',
           lead_recommendation='exclude CML (rule R1, consistent with the 13)')
ex = pd.concat([ex, pd.DataFrame([p2c])], ignore_index=True)
m6 = ex.proposal == 'P6'
ex.loc[m6, 'scope'] = 'individual minutes (TX and RX of the same minute)'
ex.loc[m6, 'reason'] = ex.loc[m6, 'reason'].str.replace('; see distribution and basis in findings',
    '; value equals the 1-min mean of an operator placeholder (TX 255, RX -99.9) and a valid sample', regex=False)
ex.loc[m6, 'evidence'] = ('[F] value gap 30-78 dBm; e.g. (255+18)/2 = 136.5 with RX (-99.9-48.9)/2 = -74.4; '
                          '[C] readers resample to 1 min by mean before RAW masks only exact sentinels')
ex.loc[m6, 'agent4_verdict'] = 'confirmed; corrupted records (RX drops 20-30 dB at 9 of 13); handle as a record-level rule'
ex.loc[m6, 'lead_recommendation'] = 'keep raw values unchanged and documented; QC masks TX and RX of these minutes'
ex.loc[m6, ['cmls_before', 'cmls_after']] = 2574
order = {'P2a': 1, 'P2b': 2, 'P2c': 3, 'P3': 4, 'P4': 5, 'P5': 6, 'P6': 7}
ex = ex.sort_values(['proposal', 'cml_id'], key=lambda s: s.map(order) if s.name == 'proposal' else s).reset_index(drop=True)
ex.to_csv(os.path.join(d, 'RAW_PROPOSED_EXCLUSIONS.csv'), index=False)

pd.DataFrame(columns=['cml_id', 'reason', 'evidence', 'decision']).to_csv(os.path.join(d, 'RAW_PROPOSED_ADDITIONS.csv'), index=False)

print('audit rows', len(au), '| rx_only>50% sublinks', int(au.rx_only_gt50pct_flag.sum()))
print(ex.groupby('proposal').agg(rows=('cml_id', 'size'), cmls=('cml_id', 'nunique'), minutes=('obs_minutes_affected', 'sum')).to_string())
rec_ex = ex[ex.lead_recommendation.str.startswith('exclude CML')].cml_id.unique()
print('lead-recommended raw CML exclusions:', len(rec_ex), sorted(rec_ex), '-> raw', 2580 - len(rec_ex), 'CMLs /', 5160 - 2 * len(rec_ex), 'sublinks')
