"""Stage 3 (lead): read-only verification of QC_CORRECTED_CANDIDATE.nc (no QC rerun, no rewrite).

Checks file size, SHA-256, dimensions, variables and attributes, and reads small samples at the start, middle and
end of the time axis.  Usage: python -I s3_03_file_verification.py [--nc PATH]
"""
import argparse, hashlib, os, json
import numpy as np
import netCDF4 as nc

ap = argparse.ArgumentParser()
ap.add_argument('--nc', default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'QC_CORRECTED_CANDIDATE.nc'))
a = ap.parse_args()
h = hashlib.sha256()
with open(a.nc, 'rb') as f:
    for ch in iter(lambda: f.read(1 << 24), b''):
        h.update(ch)
print('file', a.nc); print('size_bytes', os.path.getsize(a.nc)); print('sha256', h.hexdigest())
d = nc.Dataset(a.nc, 'r')
print('format', d.data_model); print('dimensions', {k: len(v) for k, v in d.dimensions.items()})
assert len(d.dimensions['time']) == 308160 and len(d.dimensions['sublink_id']) == 2 and len(d.dimensions['cml_id']) == 1987
need = ['tx', 'rx', 'trsl', 'qc_flag', 'qc_sublink_status', 'qc_detection_limit_mmhr', 'qc_fraction_5h_windows_std_gt_2dB',
        'qc_fraction_1h_windows_std_gt_0p8dB', 'time', 'sublink_id', 'cml_id', 'length', 'frequency', 'polarization',
        'site_a_latitude', 'site_a_longitude', 'site_b_latitude', 'site_b_longitude']
missing = [v for v in need if v not in d.variables]
print('variables', list(d.variables)); print('missing required', missing); assert not missing
for v in ('tx', 'rx', 'trsl', 'qc_flag'):
    x = d[v]; print(f'  {v}: {x.dtype} {x.dimensions} chunks={x.chunking()} filters={ {k: w for k, w in x.filters().items() if w} }')
print('qc_flag flag_masks', list(d['qc_flag'].flag_masks), d['qc_flag'].flag_meanings)
print('qc_sublink_status', d['qc_sublink_status'].flag_meanings)
for k in ('qc_period', 'population', 'n_cml_assessed', 'n_cml_retained', 'source_raw_sha256', 'qc_script_sha256', 'time_coverage_start', 'time_coverage_end'):
    print(f'  attr {k}: {d.getncattr(k)}')
P = json.loads(d.qc_parameters_json); print('  qc_parameters_json keys', len(P), 'long_min_valid', P['long_min_valid'], 'short_min_valid', P['short_min_valid'], 'plateau_min_valid', P['plateau_min_valid'])
t = d['time'][:]; print('time', d['time'].units, int(t[0]), int(t[-1]), 'step1', bool((np.diff(t) == 1).all()))
st = d['qc_sublink_status'][:]; print('status counts', {int(k): int(v) for k, v in zip(*np.unique(st, return_counts=True))})
print('every CML has a retained sublink', bool(((st == 0).any(axis=0)).all()))
ids = [str(x) for x in d['cml_id'][:]]; print('cml ids unique', len(set(ids)) == len(ids), ids[0], ids[-1])
print('polarization values', set(np.asarray(d['polarization'][:]).ravel().tolist()))
for name, sl in (('start', slice(0, 1440)), ('middle', slice(154080 - 720, 154080 + 720)), ('end', slice(308160 - 1440, 308160))):
    for c in (0, 993, 1986):
        tx = np.ma.filled(d['tx'][sl, :, c], np.nan); rx = np.ma.filled(d['rx'][sl, :, c], np.nan)
        tl = np.ma.filled(d['trsl'][sl, :, c], np.nan); fl = d['qc_flag'][sl, :, c]
        obs = np.isfinite(tx) & np.isfinite(rx) & ((fl & 24) == 0)
        diff = np.nanmax(np.abs(tl[obs] - (tx[obs] - rx[obs]))) if obs.any() else 0.0
        print(f'  {name:6s} cml[{c}] {ids[c]}: finite tx {int(np.isfinite(tx).sum())} rx {int(np.isfinite(rx).sum())} trsl {int(np.isfinite(tl).sum())}; '
              f'flags {int((fl > 0).sum())}; max|trsl-(tx-rx)| at unflagged {diff:.2e}')
d.close()
print('FILE VERIFICATION PASS')
