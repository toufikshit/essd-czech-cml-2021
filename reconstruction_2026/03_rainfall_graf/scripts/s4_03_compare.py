"""Stage 4: compare the Stage 4 rainfall with the previously validated (released) Graf-ITU2005-Pastorek product.

Read-only. Common period: 2021-04-01 00:00 .. 2021-10-31 23:59 (both files). Matching by cml_id and sublink.
Inputs: new 1-min rainfall + approved QC REV1 (for trsl, source flags); old released 1-min rainfall + old released
QC (for its trsl). Outputs (stage4/):
  RAINFALL_BEFORE_AFTER.csv       one row per CML x sublink in either product
  GRAF_THRESHOLD_COMPARISON.csv   thresholds, q80 and dry-baseline level per sublink
  RAINFALL_BEFORE_AFTER_CML.csv   accumulation per CML from R_mean_1min
  and a printed summary (saved by the caller).
Definitions: valid minute = graf_valid_1min == 1 (finite rate); wet = graf_wet_1min == 1; accumulation =
sum(R)/60 over valid minutes (mm); dry-baseline level = median trsl over valid dry minutes (the constant baseline
equals trsl at dry minutes and holds the last dry value during wet minutes).

Usage: python -I s4_03_compare.py [--new-rain ...] [--new-qc ...] [--old-rain ...] [--old-qc ...] [--out-dir ...]
"""
import argparse, os
import numpy as np
import pandas as pd
import netCDF4 as nc

H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument('--new-rain', default=os.path.join(H, 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc'))
ap.add_argument('--new-qc', default=os.path.join(H, '..', 'stage3', 'QC_CORRECTED_CANDIDATE_REV1.nc'))
ap.add_argument('--old-rain', default=r'C:\Users\taouf\ESSD_RELEASE\zenodo_staging\primary\CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc')
ap.add_argument('--old-qc', default=r'C:\Users\taouf\ESSD_RELEASE\zenodo_staging\primary\CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc')
ap.add_argument('--out-dir', default=H)
ap.add_argument('--batch', type=int, default=32)
ap.add_argument('--old-only-limit', type=int, default=0, help='test only')
a = ap.parse_args()

NR, NQ, OR, OQ = (nc.Dataset(p) for p in (a.new_rain, a.new_qc, a.old_rain, a.old_qc))
for d in (NR, NQ, OR, OQ):
    d.set_auto_mask(False)
nids = [str(x) for x in NR['cml_id'][:]]; oids = [str(x) for x in OR['cml_id'][:]]
qids = [str(x) for x in NQ['cml_id'][:]]
qpos = {c: i for i, c in enumerate(qids)}
assert set(nids) <= set(qids) and oids == [str(x) for x in OQ['cml_id'][:]]
assert np.array_equal(NR['time_1min'][:], OR['time_1min'][:]), 'time axes differ'
assert len(NR['time_1min'][:]) == 308160
opos = {c: i for i, c in enumerate(oids)}


def runs(idx):
    idx = np.asarray(idx)
    return np.split(np.arange(idx.size), np.flatnonzero(np.diff(idx) != 1) + 1)


def read(ds, var, idx):
    """(time, sublink, k) for the given cml indices (any order)."""
    idx = np.asarray(idx)
    order = np.argsort(idx); s = idx[order]
    parts = [ds[var][:, :, int(s[r[0]]):int(s[r[-1]]) + 1] for r in runs(s)]
    arr = np.concatenate(parts, axis=2)
    out = np.empty_like(arr); out[:, :, order] = arr
    return out


def side(R, W, V, TL, SF=None):
    v = (V == 1) & np.isfinite(R)
    w = v & (W == 1)
    Rv = np.where(v, R, 0.0)
    dry = v & ~w
    base = np.array([np.nanmedian(np.where(dry[:, s], TL[:, s], np.nan)) if dry[:, s].any() else np.nan for s in range(R.shape[1])])
    o = dict(valid_min=v.sum(0), wet_min=w.sum(0), acc_mm=Rv.sum(0) / 60.0, max_R=np.where(v, R, -np.inf).max(0),
             n_R_ge10=(v & (R >= 10)).sum(0), n_R_ge30=(v & (R >= 30)).sum(0), dry_baseline_median_db=base)
    o['max_R'] = np.where(np.isfinite(o['max_R']), o['max_R'], np.nan)
    if SF is not None:
        for bit, nm in ((8, 'blackout'), (16, 'interp')):
            m = v & ((SF & bit) > 0)
            o[f'valid_min_{nm}'] = m.sum(0)
            o[f'acc_mm_{nm}'] = np.where(m, R, 0.0).sum(0) / 60.0
            o[f'n_R_ge10_{nm}'] = (m & (R >= 10)).sum(0)
            o[f'n_R_ge30_{nm}'] = (m & (R >= 30)).sum(0)
    return o, v, w


rows, crow = [], []
pool = dict(n=0, sx=0.0, sy=0.0, sxx=0.0, syy=0.0, sxy=0.0)
for b0 in range(0, len(nids), a.batch):
    b1 = min(b0 + a.batch, len(nids)); ids = nids[b0:b1]
    R = NR['R_sublink_1min'][:, :, b0:b1]; W = NR['graf_wet_1min'][:, :, b0:b1]; V = NR['graf_valid_1min'][:, :, b0:b1]
    SF = NR['trsl_source_flag_1min'][:, :, b0:b1]; TL = read(NQ, 'trsl', [qpos[c] for c in ids]); RM = NR['R_mean_1min'][:, b0:b1]
    thn = NR['graf_threshold_db'][:, b0:b1]; qn = NR['graf_rstd_q80_db'][:, b0:b1]; stn = NQ['qc_sublink_status'][:, [qpos[c] for c in ids]]
    com = [k for k, c in enumerate(ids) if c in opos]
    if com:
        oi = [opos[ids[k]] for k in com]
        Ro = read(OR, 'R_sublink_1min', oi); Wo = read(OR, 'graf_wet_1min', oi); Vo = read(OR, 'graf_valid_1min', oi)
        TLo = read(OQ, 'trsl', oi); RMo = np.stack([OR['R_mean_1min'][:, i] for i in oi], 1)
        tho = OR['graf_threshold_db'][:, oi]; qo = OR['graf_rstd_q80_db'][:, oi]
    for k, c in enumerate(ids):
        sn, vn, wn = side(R[:, :, k], W[:, :, k], V[:, :, k], TL[:, :, k], SF[:, :, k])
        rmn = RM[:, k]; acc_cml_new = np.nansum(rmn) / 60.0
        old = k in com
        if old:
            j = com.index(k)
            so, vo, wo = side(Ro[:, :, j], Wo[:, :, j], Vo[:, :, j], TLo[:, :, j])
            acc_cml_old = np.nansum(RMo[:, j]) / 60.0
        crow.append(dict(cml_id=c, in_old=int(old), in_new=1, acc_cml_old_mm=acc_cml_old if old else np.nan, acc_cml_new_mm=acc_cml_new))
        for s in range(2):
            r = dict(cml_id=c, sublink=s, in_old=int(old), in_new=1, new_qc_status=int(stn[s, k]),
                     thr_new=thn[s, k], q80_new=qn[s, k])
            r.update({f'{key}_new': val[s] for key, val in sn.items()})
            if old:
                r.update(thr_old=tho[s, j], q80_old=qo[s, j])
                r.update({f'{key}_old': val[s] for key, val in so.items()})
                both = vn[:, s] & vo[:, s]
                x = Ro[both, s, j].astype(np.float64); y = R[both, s, k].astype(np.float64)
                r.update(both_valid_min=int(both.sum()), wet_both=int((wn[:, s] & wo[:, s]).sum()),
                         wet_old_dry_new=int((wo[:, s] & vn[:, s] & ~wn[:, s]).sum()),
                         dry_old_wet_new=int((wn[:, s] & vo[:, s] & ~wo[:, s]).sum()),
                         dry_both=int((both & ~wn[:, s] & ~wo[:, s]).sum()),
                         acc_old_common_mm=x.sum() / 60.0, acc_new_common_mm=y.sum() / 60.0,
                         corr_common=np.corrcoef(x, y)[0, 1] if len(x) > 2 and x.std() > 0 and y.std() > 0 else np.nan,
                         max_abs_diff_mmh=np.abs(x - y).max() if len(x) else np.nan,
                         n_min_diff_gt_0p1=int((np.abs(x - y) > 0.1).sum()))
                pool['n'] += len(x); pool['sx'] += x.sum(); pool['sy'] += y.sum(); pool['sxx'] += (x * x).sum()
                pool['syy'] += (y * y).sum(); pool['sxy'] += (x * y).sum()
            rows.append(r)
    print(f'new {b1}/{len(nids)}', flush=True)

oonly = [c for c in oids if c not in set(nids)][:a.old_only_limit or None]
oi_all = [opos[c] for c in oonly]
for b0 in range(0, len(oi_all), a.batch):
    oi = oi_all[b0:b0 + a.batch]
    Ro = read(OR, 'R_sublink_1min', oi); Wo = read(OR, 'graf_wet_1min', oi); Vo = read(OR, 'graf_valid_1min', oi)
    TLo = read(OQ, 'trsl', oi); tho = OR['graf_threshold_db'][:, oi]; qo = OR['graf_rstd_q80_db'][:, oi]
    for j, i in enumerate(oi):
        so, _, _ = side(Ro[:, :, j], Wo[:, :, j], Vo[:, :, j], TLo[:, :, j])
        crow.append(dict(cml_id=oids[i], in_old=1, in_new=0, acc_cml_old_mm=np.nansum(OR['R_mean_1min'][:, i]) / 60.0, acc_cml_new_mm=np.nan))
        for s in range(2):
            r = dict(cml_id=oids[i], sublink=s, in_old=1, in_new=0, thr_old=tho[s, j], q80_old=qo[s, j])
            r.update({f'{key}_old': val[s] for key, val in so.items()})
            rows.append(r)
    print(f'old-only {min(b0 + a.batch, len(oi_all))}/{len(oi_all)}', flush=True)

df = pd.DataFrame(rows); cm = pd.DataFrame(crow)
df.to_csv(os.path.join(a.out_dir, 'RAINFALL_BEFORE_AFTER.csv'), index=False)
cm.to_csv(os.path.join(a.out_dir, 'RAINFALL_BEFORE_AFTER_CML.csv'), index=False)
tcols = ['cml_id', 'sublink', 'in_old', 'in_new', 'new_qc_status', 'thr_old', 'thr_new', 'q80_old', 'q80_new',
         'dry_baseline_median_db_old', 'dry_baseline_median_db_new']
df[tcols].to_csv(os.path.join(a.out_dir, 'GRAF_THRESHOLD_COMPARISON.csv'), index=False)

# ---------------- summary ----------------
pd.set_option('display.width', 220)
def g(col, sel=None):
    x = df if sel is None else df[sel]
    return x[col].sum()
c_old = cm[cm.in_old == 1]; c_new = cm[cm.in_new == 1]; c_both = cm[(cm.in_old == 1) & (cm.in_new == 1)]
s_old = df[df.in_old == 1]; s_new = df[df.in_new == 1]; s_both = df[(df.in_old == 1) & (df.in_new == 1)]
print('POPULATION')
print(f'  CMLs old {len(c_old)}, new {len(c_new)}, common {len(c_both)}; old-only {len(c_old) - len(c_both)}, new-only {len(c_new) - len(c_both)}')
print(f'  CMLs with >= 1 valid rain minute: old {int((c_old.acc_cml_old_mm.notna() & (s_old.groupby("cml_id").valid_min_old.sum().reindex(c_old.cml_id).values > 0)).sum())}, '
      f'new {int((s_new.groupby("cml_id").valid_min_new.sum() > 0).sum())}')
print(f'  sublinks with >= 1 valid minute: old {int((s_old.valid_min_old > 0).sum())} of {len(s_old)}, new {int((s_new.valid_min_new > 0).sum())} of {len(s_new)}')
print('VALID / MISSING RAINFALL (directional sublink-minutes)')
print(f'  valid: old {int(s_old.valid_min_old.sum()):,}  new {int(s_new.valid_min_new.sum()):,}')
print(f'  missing among sublinks with >= 1 valid minute: old {int(((s_old.valid_min_old > 0) * (308160 - s_old.valid_min_old)).sum()):,}  '
      f'new {int(((s_new.valid_min_new > 0) * (308160 - s_new.valid_min_new)).sum()):,}')
sb = s_both[(s_both.valid_min_old > 0) & (s_both.valid_min_new > 0)]
print(f'COMMON SUBLINKS with data in both: {len(sb)}')
dt = sb.thr_new - sb.thr_old; rt = sb.thr_new / sb.thr_old
print(f'  Graf threshold: median diff {dt.median():.4f} dB, median ratio {rt.median():.4f}; |diff| > 5 %: {int((abs(rt - 1) > 0.05).sum())}; identical: {int((dt == 0).sum())}')
print('  threshold ratio quantiles', rt.quantile([0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]).round(4).to_dict())
db = sb.dry_baseline_median_db_new - sb.dry_baseline_median_db_old
print(f'  dry-baseline level: median diff {db.median():.3f} dB; |diff| > 1 dB: {int((abs(db) > 1).sum())}; quantiles', db.quantile([0.01, 0.05, 0.5, 0.95, 0.99]).round(3).to_dict())
ww, wd, dw, dd = (int(sb[k].sum()) for k in ('wet_both', 'wet_old_dry_new', 'dry_old_wet_new', 'dry_both'))
tot = ww + wd + dw + dd
print(f'  wet/dry on minutes valid in both ({tot:,}): wet-wet {ww:,}, old wet/new dry {wd:,}, old dry/new wet {dw:,}, dry-dry {dd:,}; agreement {(ww + dd) / tot:.5f}')
p = pool; n = p['n']
corr = (n * p['sxy'] - p['sx'] * p['sy']) / np.sqrt((n * p['sxx'] - p['sx'] ** 2) * (n * p['syy'] - p['sy'] ** 2))
print(f'  rain rate on common valid minutes: n {n:,}, pooled correlation {corr:.5f}, mean old {p["sx"] / n:.4f} new {p["sy"] / n:.4f} mm/h')
print(f'  accumulation on common valid minutes: old {sb.acc_old_common_mm.sum():,.0f} mm, new {sb.acc_new_common_mm.sum():,.0f} mm (sum over sublinks)')
print(f'  sublinks identical on common minutes (no minute differs by > 0.1 mm/h): {int((sb.n_min_diff_gt_0p1 == 0).sum())}')
sb2 = sb.assign(acc_diff=sb.acc_mm_new - sb.acc_mm_old)
print('  largest accumulation differences (sublink, full period):')
print(sb2.reindex(sb2.acc_diff.abs().sort_values(ascending=False).index)[['cml_id', 'sublink', 'acc_mm_old', 'acc_mm_new', 'valid_min_old', 'valid_min_new', 'thr_old', 'thr_new']].head(12).round(2).to_string(index=False))
print('  availability change (valid minutes new - old), quantiles:', (sb.valid_min_new - sb.valid_min_old).quantile([0.01, 0.05, 0.5, 0.95, 0.99]).to_dict())
print('ACCUMULATED RAINFALL (CML mean series, mm over Apr-Oct)')
print(f'  old population mean {c_old.acc_cml_old_mm.mean():.1f} (n {len(c_old)}); new population mean {c_new.acc_cml_new_mm.mean():.1f} (n {len(c_new)})')
print(f'  common CMLs: old mean {c_both.acc_cml_old_mm.mean():.1f}, new mean {c_both.acc_cml_new_mm.mean():.1f}, median of per-CML ratio {(c_both.acc_cml_new_mm / c_both.acc_cml_old_mm).median():.4f}')
print(f'  new-only CMLs mean {c_new[c_new.in_old == 0].acc_cml_new_mm.mean():.1f}; old-only CMLs mean {c_old[c_old.in_new == 0].acc_cml_old_mm.mean():.1f}')
print('FLAGGED INPUT IN THE NEW RAINFALL (all new sublinks)')
tv = s_new.valid_min_new.sum(); ta = s_new.acc_mm_new.sum()
for nm in ('blackout', 'interp'):
    print(f'  {nm}: valid minutes {int(s_new[f"valid_min_{nm}_new"].sum()):,} ({s_new[f"valid_min_{nm}_new"].sum() / tv:.4%}); '
          f'accumulation {s_new[f"acc_mm_{nm}_new"].sum():,.0f} mm of {ta:,.0f} ({s_new[f"acc_mm_{nm}_new"].sum() / ta:.3%}); '
          f'R >= 10: {int(s_new[f"n_R_ge10_{nm}_new"].sum()):,} of {int(s_new.n_R_ge10_new.sum()):,} ({s_new[f"n_R_ge10_{nm}_new"].sum() / s_new.n_R_ge10_new.sum():.3%}); '
          f'R >= 30: {int(s_new[f"n_R_ge30_{nm}_new"].sum()):,} of {int(s_new.n_R_ge30_new.sum()):,} ({s_new[f"n_R_ge30_{nm}_new"].sum() / s_new.n_R_ge30_new.sum():.3%})')
print(f'  intensity minutes R>=10 / R>=30: old {int(s_old.n_R_ge10_old.sum()):,} / {int(s_old.n_R_ge30_old.sum()):,}; new {int(s_new.n_R_ge10_new.sum()):,} / {int(s_new.n_R_ge30_new.sum()):,}')
print(f'  max 1-min rate: old {s_old.max_R_old.max():.1f}, new {s_new.max_R_new.max():.1f} mm/h')
