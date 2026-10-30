"""Stage 7: build the three corrected CML release files for the Zenodo package (no data values change).

1. RAW: CZ_CML_RAW_1MIN_20210401_20211031.nc restricted to the 2574 CMLs of the approved RAW population
   (stage2/RAW_APPROVED_POPULATION.csv). Source: the previously released RAW file (2580 CMLs = 2574 + 6 R1 CMLs).
   Same variables, dtypes, chunking, compression and attributes; one line added to `history`.
2. QC: copy of stage3/QC_CORRECTED_CANDIDATE_REV1.nc under the release name; only the global attributes `title`,
   `processing_level` and `history` are changed (the words "candidate, not released" removed).
3. Rainfall: copy of the Stage 4 1-min product under the release name; only the global attributes `title`,
   `source_cml_file`, `source_cml_sha256`, `comment` and `history` are changed.
Usage: python -I r7_01_build_release_files.py --out <package>/primary
"""
import argparse, csv, hashlib, shutil, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import netCDF4 as nc

R = Path(r'C:\Users\taouf\ESSD_RELEASE')
S = R / 'FINAL_RECONSTRUCTION_20261027'
ap = argparse.ArgumentParser()
ap.add_argument('--raw-src', default=str(R / 'zenodo_staging' / 'primary' / 'CZ_CML_RAW_1MIN_20210401_20211031.nc'))
ap.add_argument('--qc-src', default=str(S / 'stage3' / 'QC_CORRECTED_CANDIDATE_REV1.nc'))
ap.add_argument('--rain-src', default=str(S / 'stage4' / 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc'))
ap.add_argument('--population', default=str(S / 'stage2' / 'RAW_APPROVED_POPULATION.csv'))
ap.add_argument('--out', required=True)
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
RAW_NAME = 'CZ_CML_RAW_1MIN_20210401_20211031.nc'
QC_NAME = 'CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc'
RAIN_NAME = 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc'


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------ 1. RAW 2574
ids = sorted({r['cml_id'] for r in csv.DictReader(open(a.population, encoding='utf-8'))})
assert len(ids) == 2574
src = nc.Dataset(a.raw_src)
src_ids = np.array([str(x) for x in src['cml_id'][:]])
keep = np.isin(src_ids, ids)
assert keep.sum() == 2574, keep.sum()
kidx = np.flatnonzero(keep)
removed = sorted(src_ids[~keep].tolist())
print('RAW source', src.dimensions['cml_id'].size, 'CMLs; kept', keep.sum(), '; removed', removed)
dst_path = out / RAW_NAME
dst = nc.Dataset(dst_path, 'w', format=src.file_format)
for name, dim in src.dimensions.items():
    dst.createDimension(name, 2574 if name == 'cml_id' else len(dim))
for v in src.variables.values():
    kw = {}
    ch = v.chunking()
    if ch != 'contiguous':
        kw['chunksizes'] = ch
    flt = v.filters() or {}
    if flt.get('zlib'):
        kw.update(zlib=True, complevel=flt.get('complevel', 4), shuffle=flt.get('shuffle', False))
    fill = v.getncattr('_FillValue') if '_FillValue' in v.ncattrs() else None
    dtype = str if v.dtype == str else v.dtype
    w = dst.createVariable(v.name, dtype, v.dimensions, fill_value=fill, **kw)
    w.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != '_FillValue'})
    v.set_auto_maskandscale(False); w.set_auto_maskandscale(False)
    if 'cml_id' not in v.dimensions:
        w[...] = v[...]
    elif v.dimensions[0] == 'time':
        ax = v.dimensions.index('cml_id')
        for t0 in range(0, len(src.dimensions['time']), 10080):
            blk = v[t0:t0 + 10080]
            w[t0:t0 + 10080] = np.take(blk, kidx, axis=ax)
    else:
        ax = v.dimensions.index('cml_id')
        w[...] = np.take(v[...], kidx, axis=ax)
g = {k: src.getncattr(k) for k in src.ncattrs()}
g['history'] = (g.get('history', '') + f'\n{NOW} restricted to the 2574 CMLs of the approved RAW population; '
                f'6 CMLs without usable paired observations removed ({", ".join(removed)}); values unchanged.')
dst.setncatts(g)
dst.close(); src.close()
print('wrote', dst_path, dst_path.stat().st_size, sha256(dst_path))

# ------------------------------------------------------------------ 2. QC
qc_path = out / QC_NAME
shutil.copyfile(a.qc_src, qc_path)
d = nc.Dataset(qc_path, 'a')
d.setncattr('title', 'Quality-controlled commercial microwave link signal levels, Czech Republic, 1 April - 31 October 2021')
d.setncattr('processing_level', 'Quality-controlled CML signal data')
d.setncattr('history', d.getncattr('history') + f'\n{NOW} released as {QC_NAME}; only global attributes title, '
            'processing_level and history changed.')
d.close()
qc_sha = sha256(qc_path)
print('wrote', qc_path, qc_path.stat().st_size, qc_sha)

# ------------------------------------------------------------------ 3. Rainfall
rain_path = out / RAIN_NAME
shutil.copyfile(a.rain_src, rain_path)
d = nc.Dataset(rain_path, 'a')
d.setncattr('title', '1-minute CML rainfall estimates using Graf RSTD, Pastorek wet antenna correction and '
            'ITU-R P.838-3 (2005), Czech Republic, 1 April - 31 October 2021')
d.setncattr('source_cml_file', QC_NAME)
d.setncattr('source_cml_sha256', qc_sha)
d.setncattr('comment', 'Graf thresholds and the dry baseline are calculated on the period 2021-04-01 to 2021-10-31. '
            'Input total loss includes flagged blackout-reconstructed and interpolated minutes; see '
            'trsl_source_flag_1min. No TX or RX values were generated.')
hist = d.getncattr('history') + '\n' if 'history' in d.ncattrs() else ''
d.setncattr('history', hist + f'{NOW} released as {RAIN_NAME}; only global attributes title, source_cml_file, '
            'source_cml_sha256, comment and history changed.')
d.close()
print('wrote', rain_path, rain_path.stat().st_size, sha256(rain_path))
