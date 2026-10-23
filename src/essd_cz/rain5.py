"""Derive the analysis 5-minute rainfall from the released native 1-minute product.

The published CML rainfall NetCDF is native 1-minute
(`CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc`). All
manuscript rainfall benchmarks use 5-minute quantities, which are an exact
function of that release: a directional 5-minute rate exists only where ALL
FIVE 1-minute samples in the wall-clock-aligned block are valid, and is their
arithmetic mean; the 5-minute wet flag is set when any of those five valid
minutes was wet; the physical-link mean is taken over the finite DIRECTIONAL
5-minute rates afterwards (never by averaging `R_mean_1min`, whose validity
mask differs). These are the unchanged semantics of the original producer.
The derivation was verified bit-for-bit against the validated full 5-minute
reference for all 2,219 links (see docs/EVIDENCE.md).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

AGGREGATION_NOTE = ('5-minute values derived from the released native 1-minute product: '
                    'mean of five complete valid 1-minute directional samples; wet if any '
                    'valid wet minute; physical mean over finite directional 5-minute rates')
META_1D = ('site_a_latitude', 'site_a_longitude', 'site_b_latitude', 'site_b_longitude')
META_2D = ('length', 'frequency', 'polarization', 'graf_threshold_db', 'graf_rstd_q80_db')


def is_native_1min(dataset_or_path) -> bool:
    """A rain product is the native release when it carries R_sublink_1min."""
    if isinstance(dataset_or_path, xr.Dataset):
        return 'R_sublink_1min' in dataset_or_path
    with xr.open_dataset(dataset_or_path, decode_times=False) as ds:
        return 'R_sublink_1min' in ds


def derived_5min_name(native_name: str) -> str:
    """CZ_CML_RAINRATE_..._1MIN_a_b.nc -> CZ_CML_RAINRATE_..._5MIN_DERIVED_a_b.nc"""
    if '_1MIN_' not in native_name:
        raise ValueError(f'Not a native 1-minute product name: {native_name!r}')
    return native_name.replace('_1MIN_', '_5MIN_DERIVED_', 1)


def check_native_axis(times: pd.DatetimeIndex):
    """The release axis must be a regular, wall-clock 5-min-aligned 1-min grid."""
    t = times.asi8
    if len(t) % 5:
        raise ValueError(f'1-minute time axis length {len(t)} is not divisible by 5')
    if len(t) > 1 and not np.all(np.diff(t) == 60_000_000_000):
        raise ValueError('Released product time axis is not a regular 1-minute grid')
    if t[0] % 300_000_000_000:
        raise ValueError('1-minute axis does not start on a 5-minute boundary')
    return times[::5]


def _blocks(R1, W1, V1):
    """Exact original aggregation of (time,sublink,link) 1-min arrays to 5-min."""
    nt, ns, nl = R1.shape
    n5 = nt // 5
    vv = (V1 == 1) & np.isfinite(R1)
    rb = R1.reshape(n5, 5, ns, nl)
    vb = vv.reshape(n5, 5, ns, nl)
    wb = ((W1 == 1) & vv).reshape(n5, 5, ns, nl)
    counts = vb.sum(1).astype(np.uint8)
    wc = wb.sum(1).astype(np.uint8)
    complete = counts == 5
    r5 = np.full((n5, ns, nl), np.nan, np.float32)
    r5[complete] = np.mean(np.where(vb, rb, np.nan), axis=1, dtype=np.float64)[complete].astype(np.float32)
    w5 = np.full((n5, ns, nl), 255, np.uint8)
    w5[complete] = (wc[complete] > 0).astype(np.uint8)
    cnt = np.sum(np.isfinite(r5), axis=1)
    sm = np.nansum(r5, axis=1)
    rm5 = np.full((n5, nl), np.nan, np.float32)
    ok = cnt > 0
    rm5[ok] = (sm[ok] / cnt[ok]).astype(np.float32)
    return {'R_sublink_5min': r5, 'R_mean_5min': rm5, 'graf_wet_5min': w5,
            'graf_valid_5min': complete.astype(np.uint8),
            'wet_minutes_5min': wc, 'valid_minutes_5min': counts}


def derive_link_5min(rain: xr.Dataset, link_index: int) -> xr.Dataset:
    """Derive one physical link's 5-minute dataset from the opened 1-min release."""
    return derive_5min_dataset(rain.isel(cml_id=[int(link_index)]))


def derive_5min_dataset(rain: xr.Dataset, link_batch: int = 128) -> xr.Dataset:
    """Derive the full analysis 5-minute dataset (legacy product schema) in memory."""
    t5 = check_native_axis(pd.DatetimeIndex(rain['time_1min'].values))
    ns = rain.sizes['sublink_id']
    nl = rain.sizes['cml_id']
    n5 = len(t5)
    full = {'R_sublink_5min': np.full((n5, ns, nl), np.nan, np.float32),
            'R_mean_5min': np.full((n5, nl), np.nan, np.float32),
            'graf_wet_5min': np.full((n5, ns, nl), 255, np.uint8),
            'graf_valid_5min': np.zeros((n5, ns, nl), np.uint8),
            'wet_minutes_5min': np.zeros((n5, ns, nl), np.uint8),
            'valid_minutes_5min': np.zeros((n5, ns, nl), np.uint8)}
    for b0 in range(0, nl, max(1, int(link_batch))):
        b1 = min(b0 + link_batch, nl)
        part = rain.isel(cml_id=slice(b0, b1))
        R1 = np.asarray(part['R_sublink_1min'].transpose('time_1min', 'sublink_id', 'cml_id').values, np.float32)
        W1 = np.asarray(part['graf_wet_1min'].transpose('time_1min', 'sublink_id', 'cml_id').values, np.uint8)
        V1 = np.asarray(part['graf_valid_1min'].transpose('time_1min', 'sublink_id', 'cml_id').values, np.uint8)
        got = _blocks(R1, W1, V1)
        for key, arr in got.items():
            if key == 'R_mean_5min':
                full[key][:, b0:b1] = arr
            else:
                full[key][:, :, b0:b1] = arr
    data = {'R_sublink_5min': (('time_5min', 'sublink_id', 'cml_id'), full['R_sublink_5min'], {'units': 'mm h-1'}),
            'R_mean_5min': (('time_5min', 'cml_id'), full['R_mean_5min'], {'units': 'mm h-1'}),
            'graf_wet_5min': (('time_5min', 'sublink_id', 'cml_id'), full['graf_wet_5min']),
            'graf_valid_5min': (('time_5min', 'sublink_id', 'cml_id'), full['graf_valid_5min']),
            'wet_minutes_5min': (('time_5min', 'sublink_id', 'cml_id'), full['wet_minutes_5min']),
            'valid_minutes_5min': (('time_5min', 'sublink_id', 'cml_id'), full['valid_minutes_5min'])}
    for name in META_2D:
        if name in rain:
            data[name] = (rain[name].dims, rain[name].values, dict(rain[name].attrs))
    for name in META_1D:
        if name in rain:
            data[name] = (rain[name].dims, rain[name].values, dict(rain[name].attrs))
    attrs = dict(rain.attrs)
    attrs['title'] = 'Analysis 5-minute CML rainfall derived from the released native 1-minute product'
    attrs['derived_from'] = str(rain.encoding.get('source', 'native 1-minute release'))
    attrs['derivation'] = AGGREGATION_NOTE
    return xr.Dataset(data,
                      coords={'time_5min': t5, 'sublink_id': rain['sublink_id'].values,
                              'cml_id': rain['cml_id'].values},
                      attrs=attrs)


def derive_5min_file(native_nc, output_nc, link_batch: int = 128) -> Path:
    """Materialize the analysis 5-minute NetCDF next to a workspace, never in-place."""
    native_nc = Path(native_nc).expanduser().resolve()
    output_nc = Path(output_nc).expanduser().resolve()
    if output_nc.exists():
        raise FileExistsError(f'Refusing to overwrite existing derived file: {output_nc}')
    if output_nc == native_nc:
        raise ValueError('Derived output must not replace the released input')
    output_nc.parent.mkdir(parents=True, exist_ok=True)
    with xr.open_dataset(native_nc, decode_times=True, mask_and_scale=False) as rain:
        if not is_native_1min(rain):
            raise ValueError(f'Not a native 1-minute rainfall product: {native_nc}')
        derived = derive_5min_dataset(rain, link_batch=link_batch)
        derived.attrs['derived_from'] = native_nc.name
    # The legacy product stores time_5min as int64 epoch seconds; original
    # scripts read those raw values, so the derived file keeps that encoding.
    time_enc = {'units': 'seconds since 1970-01-01 00:00:00 UTC',
                'calendar': 'standard', 'dtype': 'int64'}
    try:
        import netCDF4  # noqa: F401
        engine = 'netcdf4'
        enc = {k: {'zlib': True, 'complevel': 4, 'shuffle': True}
               for k in ('R_sublink_5min', 'R_mean_5min', 'graf_wet_5min',
                         'graf_valid_5min', 'wet_minutes_5min', 'valid_minutes_5min')}
        enc['time_5min'] = time_enc
    except ImportError:
        engine = 'scipy'
        # Classic NetCDF has no int64; epoch seconds fit int32 (fixture path only).
        enc = {'time_5min': dict(time_enc, dtype='int32')}
        # Classic NetCDF has no unsigned byte; keep flag VALUES, change dtype only.
        for key in ('graf_wet_5min', 'graf_valid_5min', 'wet_minutes_5min', 'valid_minutes_5min'):
            derived[key] = derived[key].astype(np.int16)
    tmp = output_nc.with_name(output_nc.name + '.partial')
    if tmp.exists():
        raise FileExistsError(f'Partial derived file already exists, inspect before retry: {tmp}')
    try:
        derived.to_netcdf(tmp, engine=engine, encoding=enc)
        tmp.rename(output_nc)
    finally:
        tmp.unlink(missing_ok=True)
    return output_nc
