"""Derive the 5-minute CML rainfall from the released 1-minute rainfall file.

Uses essd_cz.rain5.derive_5min_dataset: a directional 5-minute rate exists only where all five 1-minute samples
are valid and is their mean; wet if any valid minute is wet; the CML mean is taken over the finite directional
5-minute rates. Two counts are added: blackout_minutes_5min and interpolated_minutes_5min = number of the valid
1-minute samples of the interval whose total loss was blackout-reconstructed (trsl_source_flag bit 8) or
interpolated (bit 16).
The 5-minute series are used by the gauge validation and Figs. 4 and 7. They are not a release file.
Usage: python scripts/derive_5min.py [--config config/user_paths.json] [--native RAIN_1MIN.nc] [--output OUT.nc]
"""
import argparse, hashlib, os, sys
import numpy as np
import xarray as xr

ap = argparse.ArgumentParser()
ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
ap.add_argument('--native', default=None, help='1-minute rainfall (default: rain of the configuration)')
ap.add_argument('--output', default=None, help='default: <output_dir>/rainfall_5min/<name>')
a = ap.parse_args()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
from essd_cz import rain5  # noqa: E402
from essd_cz.config import load_config  # noqa: E402
cfg = load_config(a.config)
a.native = a.native or str(cfg['datasets']['rain'])
a.output = a.output or str(cfg['output_dir'] / 'rainfall_5min' / rain5.derived_5min_name(os.path.basename(a.native)))
os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)

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
