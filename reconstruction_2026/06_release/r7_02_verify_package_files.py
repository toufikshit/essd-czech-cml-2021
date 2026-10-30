"""Stage 7: verify the files of the corrected Zenodo package by content (not by name).

* RAW: cml_id set = approved RAW population (2574); every variable equals the previous RAW release restricted to these
  CMLs (NaN-aware, read in time blocks).
* QC and rainfall: every variable equals the approved Stage 3 / Stage 4 file; only the listed global attributes differ.
* Populations: QC 1987 CMLs and 3731 retained sublinks, QC CMLs are a subset of RAW, rainfall cml_id order = QC,
  finite directional 1-min rainfall values = 631,742,984; radar 1132 paths; gauges 424; temperature 251; humidity 247.
Usage: python -I r7_02_verify_package_files.py --package <package dir>
"""
import argparse, csv
from pathlib import Path
import numpy as np
import netCDF4 as nc

R = Path(r'C:\Users\taouf\ESSD_RELEASE')
S = R / 'FINAL_RECONSTRUCTION_20261027'
ap = argparse.ArgumentParser()
ap.add_argument('--package', required=True)
a = ap.parse_args()
P = Path(a.package) / 'primary'
OK = []


def check(label, cond):
    OK.append(bool(cond)); print(('PASS ' if cond else 'FAIL ') + label)


def ids(d):
    return np.array([str(x) for x in d['cml_id'][:]])


def same(x, y):
    x = np.asarray(x); y = np.asarray(y)
    if x.dtype.kind == 'f':
        return x.shape == y.shape and np.array_equal(x, y, equal_nan=True)
    return x.shape == y.shape and np.array_equal(x, y)


def compare_all(new, old, sel=None, tdim=None):
    """Compare every variable of new with old (optionally old restricted to cml indices sel)."""
    bad = []
    for name, v in new.variables.items():
        w = old[name]; v.set_auto_maskandscale(False); w.set_auto_maskandscale(False)
        if v.dimensions != w.dimensions:
            bad.append(name + ' dims'); continue
        ax = v.dimensions.index('cml_id') if 'cml_id' in v.dimensions else None
        if v.ndim >= 2 and v.dimensions[0] == tdim:
            n = len(new.dimensions[tdim])
            for t0 in range(0, n, 10080):
                b_new = v[t0:t0 + 10080]; b_old = w[t0:t0 + 10080]
                if sel is not None and ax is not None:
                    b_old = np.take(b_old, sel, axis=ax)
                if not same(b_new, b_old):
                    bad.append(f'{name} block {t0}'); break
        else:
            b_new = v[...]; b_old = w[...]
            if sel is not None and ax is not None:
                b_old = np.take(np.asarray(b_old), sel, axis=ax)
            if not same(b_new, b_old):
                bad.append(name)
    return bad


# RAW
approved = sorted({r['cml_id'] for r in csv.DictReader(open(S / 'stage2' / 'RAW_APPROVED_POPULATION.csv', encoding='utf-8'))})
raw = nc.Dataset(P / 'CZ_CML_RAW_1MIN_20210401_20211031.nc')
old_raw = nc.Dataset(R / 'zenodo_staging' / 'primary' / 'CZ_CML_RAW_1MIN_20210401_20211031.nc')
rid = ids(raw)
check(f'RAW cml_id = approved RAW population ({len(rid)} CMLs, {len(set(rid))} unique)', sorted(rid) == approved and len(rid) == 2574)
sel = np.flatnonzero(np.isin(ids(old_raw), rid))
check('RAW CML order kept from the previous release', np.array_equal(ids(old_raw)[sel], rid))
bad = compare_all(raw, old_raw, sel, 'time')
check('RAW every variable equal to the previous release restricted to the 2574 CMLs' + (f' {bad}' if bad else ''), not bad)
check('RAW dimensions ' + str({k: len(v) for k, v in raw.dimensions.items()}), len(raw.dimensions['time']) == 308160)

# QC
qc = nc.Dataset(P / 'CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc')
qc0 = nc.Dataset(S / 'stage3' / 'QC_CORRECTED_CANDIDATE_REV1.nc')
bad = compare_all(qc, qc0, None, 'time')
check('QC every variable equal to approved QC REV1' + (f' {bad}' if bad else ''), not bad)
diff_attrs = sorted(k for k in set(qc.ncattrs()) | set(qc0.ncattrs()) if str(getattr(qc, k, None)) != str(getattr(qc0, k, None)))
check(f'QC changed global attributes {diff_attrs}', set(diff_attrs) <= {'title', 'processing_level', 'history'})
qid = ids(qc); st = qc['qc_sublink_status'][:]
check(f'QC {len(qid)} CMLs, {(np.asarray(st) == 0).sum()} retained sublinks, all in RAW', len(qid) == 1987 and (np.asarray(st) == 0).sum() == 3731 and set(qid) <= set(rid))

# Rainfall
rn = nc.Dataset(P / 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc')
rn0 = nc.Dataset(S / 'stage4' / 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc')
bad = compare_all(rn, rn0, None, 'time_1min')
check('Rainfall every variable equal to approved Stage 4 product' + (f' {bad}' if bad else ''), not bad)
diff_attrs = sorted(k for k in set(rn.ncattrs()) | set(rn0.ncattrs()) if str(getattr(rn, k, None)) != str(getattr(rn0, k, None)))
check(f'Rainfall changed global attributes {diff_attrs}', set(diff_attrs) <= {'title', 'source_cml_file', 'source_cml_sha256', 'comment', 'history'})
check('Rainfall cml_id identical to QC (1987, same order)', np.array_equal(ids(rn), qid))
v = rn['R_sublink_1min']; n = 0
for t0 in range(0, len(rn.dimensions['time_1min']), 10080):
    n += int(np.isfinite(np.asarray(v[t0:t0 + 10080], dtype=np.float32)).sum())
check(f'Rainfall finite directional 1-min values {n:,} (approved 631,742,984)', n == 631742984)

# Other products
for f, dim, exp in [('CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc', 'cml_id', 1132), ('CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc', 'WSI', 424),
                    ('CHMI_10MIN_TEMPERATURE_20210401_20211031.nc', 'station', 251), ('CHMI_10MIN_HUMIDITY_20210401_20211031.nc', 'station', 247)]:
    d = nc.Dataset(P / f); check(f'{f}: {dim} = {len(d.dimensions[dim])}', len(d.dimensions[dim]) == exp); d.close()
print('ALL PASS' if all(OK) else 'SOME CHECKS FAILED')
