"""Stage 4: derive the 5-minute rainfall from the Stage 4 native 1-minute product.

Uses the repository's unchanged aggregation `essd_cz.rain5.derive_5min_dataset` (imported read-only from
essd-czech-cml-2021/src): a directional 5-minute rate exists only where all five 1-minute samples are valid and is
their mean; wet if any valid minute is wet; the CML mean is taken over the finite directional 5-minute rates.
Added (documentation of flag propagation, no effect on the values): `blackout_minutes_5min` and
`interpolated_minutes_5min` = number of the valid 1-minute samples of the interval whose total loss was
blackout-reconstructed (trsl_source_flag bit 8) or interpolated (bit 16).

Usage: python -I s4_02_derive_5min.py --native ../<1min>.nc --output ../<5min>.nc [--repo DIR]
"""
import argparse, hashlib, os, sys
import numpy as np
import xarray as xr

ap = argparse.ArgumentParser()
ap.add_argument('--native', required=True)
ap.add_argument('--output', required=True)
ap.add_argument('--repo', default=r'C:\Users\taouf\ESSD_RELEASE\essd-czech-cml-2021')
a = ap.parse_args()
sys.path.insert(0, os.path.join(a.repo, 'src'))
from essd_cz import rain5  # noqa: E402

if os.path.exists(a.output):
    sys.exit(f'refusing to overwrite {a.output}')
with xr.open_dataset(a.native, decode_times=True, mask_and_scale=False) as rain:
    assert rain5.is_native_1min(rain)
    der = rain5.derive_5min_dataset(rain, link_batch=64)
    nt, ns, nl = rain.sizes['time_1min'], rain.sizes['sublink_id'], rain.sizes['cml_id']
    bo = np.zeros((nt // 5, ns, nl), np.uint8)
    it = np.zeros((nt // 5, ns, nl), np.uint8)
    for b0 in range(0, nl, 64):
        b1 = min(b0 + 64, nl)
        part = rain.isel(cml_id=slice(b0, b1))
        R1 = np.asarray(part['R_sublink_1min'].transpose('time_1min', 'sublink_id', 'cml_id').values, np.float32)
        V1 = np.asarray(part['graf_valid_1min'].transpose('time_1min', 'sublink_id', 'cml_id').values, np.uint8)
        SF = np.asarray(part['trsl_source_flag_1min'].transpose('time_1min', 'sublink_id', 'cml_id').values, np.uint8)
        vv = (V1 == 1) & np.isfinite(R1)
        n5 = nt // 5
        bo[:, :, b0:b1] = (((SF & 8) > 0) & vv).reshape(n5, 5, ns, b1 - b0).sum(1)
        it[:, :, b0:b1] = (((SF & 16) > 0) & vv).reshape(n5, 5, ns, b1 - b0).sum(1)
    der['blackout_minutes_5min'] = (('time_5min', 'sublink_id', 'cml_id'), bo,
                                    {'long_name': 'Valid 1-minute samples in the interval computed from blackout-reconstructed total loss'})
    der['interpolated_minutes_5min'] = (('time_5min', 'sublink_id', 'cml_id'), it,
                                        {'long_name': 'Valid 1-minute samples in the interval computed from interpolated total loss'})
    der.attrs['derived_from'] = os.path.basename(a.native)
    h = hashlib.sha256()
    with open(a.native, 'rb') as fh:
        for ch in iter(lambda: fh.read(1 << 24), b''):
            h.update(ch)
    der.attrs['derived_from_sha256'] = h.hexdigest()
enc = {k: {'zlib': True, 'complevel': 4, 'shuffle': True} for k in
       ('R_sublink_5min', 'R_mean_5min', 'graf_wet_5min', 'graf_valid_5min', 'wet_minutes_5min',
        'valid_minutes_5min', 'blackout_minutes_5min', 'interpolated_minutes_5min')}
enc['time_5min'] = {'units': 'seconds since 1970-01-01 00:00:00 UTC', 'calendar': 'standard', 'dtype': 'int64'}
tmp = a.output + '.partial'
der.to_netcdf(tmp, engine='netcdf4', encoding=enc)
os.replace(tmp, a.output)
print('DONE', a.output, os.path.getsize(a.output))
print('dims', dict(der.sizes))
print('finite directional 5-min rates', int(np.isfinite(der['R_sublink_5min'].values).sum()))
print('5-min intervals with >= 1 blackout minute', int((bo > 0).sum()), '; with >= 1 interpolated minute', int((it > 0).sum()))
