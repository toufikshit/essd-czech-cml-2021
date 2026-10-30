"""Validation of the 5-minute CML rainfall against the CHMI rain gauges (no figures).

Method: nearest CHMI gauge to the CML midpoint within 5 km; gauge rr (10 min, end-labelled) * 6 = mm/h;
CML rate = mean of the two complete 5-minute directional bins starting at t-10 and t-5; Graf RSTD detection
recomputed on the 1-minute trsl of the CML file (window 60, q 0.8 x 1.12, 60/45 minimum periods), aggregated
'any' over the ten minutes of the gauge interval with all ten valid; one joint mask (rain valid, detection valid,
gauge 0..150 mm/h, CML 0..150 mm/h); gauge wet >= 0.1 mm/h; per directional sublink.
The functions gauge_validation() and compute_stats() are those of the figure producer
original/code/figures/rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py, without the
plots. The Graf helper functions are loaded from original/code/validation/cml_graf_rstd_hess_radar_gauge_v3.py.
Every joint sample is saved (cml, sublink, gauge, gauge-time index, x, y, predicted wet); Figs. 5 and 6 are drawn
from these samples (scripts/figures/fig05_fig06_validation.py). With --blackout-qc, a sample is marked bit8 if any
of the ten 1-minute QC minutes of its interval has qc_flag bit 8; statistics are then reported for all samples and
for samples without bit 8.

Inputs: the QC file, the 5-minute rainfall of scripts/derive_5min.py and the gauge file.
Outputs in <output_dir>/validation: <tag>_summary.json (statistics), <tag>_samples.npz, <tag>_audit.csv.
Usage: python scripts/gauge_validation.py [--config config/user_paths.json]
"""
import argparse, hashlib, importlib.util, json, math, os, re, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import xarray as xr

MAX_R = 150.0
START, END = pd.Timestamp('2021-04-01'), pd.Timestamp('2021-10-31 23:59:59')
ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / 'original' / 'code' / 'validation' / 'cml_graf_rstd_hess_radar_gauge_v3.py'


def key(value):
    if isinstance(value, bytes):
        value = value.decode('utf-8', 'replace')
    if isinstance(value, str):
        result = value.strip()
        if re.fullmatch(r'[+-]?[0-9]+\.0+', result):
            return result.split('.')[0]
        return result
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)) and np.isfinite(value):
        if float(value).is_integer():
            return str(int(value))
    return str(value).strip()


def ids(ds, name='cml_id'):
    return [key(v) for v in ds[name].values]


def positions(ds, name='cml_id'):
    out = dict((v, i) for i, v in enumerate(ids(ds, name)))
    if len(out) != ds.sizes[name]:
        raise ValueError(f'Duplicate {name} values in dataset')
    return out


def select_rain_rate(rain):
    da = rain['R_sublink_5min']
    if set(da.dims) != {'time_5min', 'sublink_id', 'cml_id'}:
        raise ValueError(da.dims)
    return da


def compute_stats(x, y):
    if len(x) == 0:
        return {'N': 0, 'RMSE': None, 'MBE': None, 'PCC': None}
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    diff = y - x
    return {'N': int(len(x)), 'RMSE': float(np.sqrt(np.mean(diff * diff))),
            'MBE': float(np.mean(diff)),
            'PCC': float(np.corrcoef(x, y)[0, 1]) if len(x) >= 2 and x.std() > 0 and y.std() > 0 else None}


def numeric(da, *dims):
    arr = np.asarray(da.values, dtype=float)
    if set(da.dims) != set(dims):
        raise ValueError(f'Expected {dims}, found {da.dims}')
    return np.transpose(arr, [da.dims.index(s) for s in dims])


def skill(TP, FP, FN, TN):
    den = math.sqrt(float(TP + FP) * float(TP + FN) * float(TN + FP) * float(TN + FN))
    return {'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN, 'N': TP + FP + FN + TN,
            'POD': TP / (TP + FN) if TP + FN else None, 'FAR': FP / (TP + FP) if TP + FP else None,
            'CSI': TP / (TP + FP + FN) if TP + FP + FN else None, 'MCC': (TP * TN - FP * FN) / den if den else None,
            'correct_negative_fraction': TN / (TN + FP) if TN + FP else None}


def summarise(x, y, p, threshold=0.1):
    truth = x >= threshold
    TP = int((p & truth).sum()); FP = int((p & ~truth).sum()); FN = int((~p & truth).sum()); TN = int((~p & ~truth).sum())
    return {'all_period': compute_stats(x, y), 'gauge_wet_ge_0_1': compute_stats(x[truth], y[truth]),
            'confusion': skill(TP, FP, FN, TN)}


def gauge_validation(cml, rain, gauge, g, start, end, threshold=0.1, bqc=None):
    rate = select_rain_rate(rain).sel(time_5min=slice(start - pd.Timedelta(minutes=10), end))
    rt = pd.DatetimeIndex(rate.time_5min.values)
    times = pd.DatetimeIndex(gauge.time.values)
    keep = np.flatnonzero((times >= start + pd.Timedelta(minutes=10)) & (times <= end))
    gt = times[keep]
    i1 = rt.get_indexer(gt - pd.Timedelta(minutes=10))
    i2 = rt.get_indexer(gt - pd.Timedelta(minutes=5))
    exact = (i1 >= 0) & (i2 >= 0)
    if not exact.any():
        raise ValueError('No exact 10min gauge / 5min final Graf rainfall matches')
    i1 = np.clip(i1, 0, len(rt) - 1); i2 = np.clip(i2, 0, len(rt) - 1)
    st = pd.DatetimeIndex(cml.time.values)
    imap = g.build_reference_intervals(st, gt, 'end', 10.)
    if int(imap['expected']) != 10:
        raise ValueError('Expected ten 1min Graf samples per gauge interval')
    rainpos = positions(rain); cpos = positions(cml)
    cids = np.asarray(cml['cml_id'].values).astype(str)
    arg = SimpleNamespace(site_a_lon='site_a_longitude', site_a_lat='site_a_latitude',
                          site_b_lon='site_b_longitude', site_b_lat='site_b_latitude',
                          gauge_lon=None, gauge_lat=None, gauge_id=None, max_gauge_dist_km=5)
    info = g.prepare_gauge_mapping(cml=cml, cml_ids=cids, gauge=gauge, gauge_var='rr', args=arg)
    mapping = info['mapping']; mapping = mapping[(mapping.gauge_position >= 0) & (mapping.distance_km <= 5)].copy()
    print(f'Gauge midpoint matched CMLs <=5km: {len(mapping):,}', flush=True)
    gv = info['rr']; gdim = info['station_dim']; tdim = info['time_dim']
    sid_c = ids(cml, 'sublink_id'); sid_r = ids(rate, 'sublink_id')
    if len(set(sid_c)) != len(sid_c) or set(sid_c) != set(sid_r):
        raise ValueError(f'Directional sublinks differ: QC {sid_c}; Graf rainfall {sid_r}')
    ri = {sid: i for i, sid in enumerate(sid_r)}
    source = g.detect_trsl_source(cml, None)
    if bqc is not None:
        bpos = positions(bqc)
        bt = pd.DatetimeIndex(bqc.time.values)
        bmap = g.build_reference_intervals(bt, gt, 'end', 10.)
    cache = {}; audit = []; S = {k: [] for k in ('cml', 'sub', 'gauge', 'tidx', 'x', 'y', 'pred', 'bit8')}
    pre_rain = pre_det = 0; rain_membership = 0
    for k, (raw_cid, station, gid, dist) in enumerate(mapping[['cml_id', 'gauge_position', 'gauge_id', 'distance_km']].itertuples(index=False, name=None), 1):
        cid = key(raw_cid)
        if cid not in rainpos or cid not in cpos:
            continue
        rain_membership += 1
        station = int(station)
        if station not in cache:
            cache[station] = np.asarray(gv.isel({gdim: station, tdim: keep}).values, float).reshape(-1) * 6
        ref = cache[station]
        cp = cpos[cid]
        sig = (cml['tx'].isel(cml_id=cp) - cml['rx'].isel(cml_id=cp)) if source is None else cml[source].isel(cml_id=cp)
        sig = numeric(sig, 'time', 'sublink_id')
        wet, valid, *_ = g.graf_rstd_matrix(sig, window=60, quantile=.8, factor=1.12,
                                            threshold_min_periods=60, detection_min_periods=45)
        pred, predvalid = g.aggregate_to_reference(wet, valid, imap, 1., 'any')
        if pred.shape != (len(gt), len(sid_c)):
            raise ValueError(f'Graf detector shape mismatch: {pred.shape}')
        rainfall = numeric(rate.isel(cml_id=rainpos[cid]), 'sublink_id', 'time_5min')
        base = exact & imap['complete'] & np.isfinite(ref) & (ref >= 0) & (ref <= MAX_R)
        if bqc is not None:
            fl = numeric(bqc['qc_flag'].isel(cml_id=bpos[cid]), 'time', 'sublink_id').astype(np.int64)
            b8 = (fl & 8) > 0
            b8ref, _ = g.aggregate_to_reference(b8, np.ones_like(b8), bmap, 1., 'any')
        for i, sid in enumerate(sid_c):
            j = ri[sid]
            a = rainfall[j, i1]; b = rainfall[j, i2]
            y = (a + b) * .5
            rainvalid = base & np.isfinite(a) & np.isfinite(b) & np.isfinite(y) & (y >= 0) & (y <= MAX_R)
            detvalid = base & predvalid[:, i]
            pre_rain += int(rainvalid.sum()); pre_det += int(detvalid.sum())
            ok = rainvalid & detvalid
            if not ok.any():
                continue
            x = ref[ok].astype(np.float32); yy = y[ok].astype(np.float32)
            predicted = np.asarray(pred[ok, i], bool)
            idx = np.flatnonzero(ok)
            S['cml'].append(np.full(idx.size, cp, np.int32)); S['sub'].append(np.full(idx.size, i, np.int8))
            S['gauge'].append(np.full(idx.size, station, np.int32)); S['tidx'].append(idx.astype(np.int32))
            S['x'].append(x); S['y'].append(yy); S['pred'].append(predicted)
            S['bit8'].append(np.asarray(b8ref[ok, i], bool) if bqc is not None else np.zeros(idx.size, bool))
            truth = x >= threshold
            audit.append({'cml_id': cid, 'sublink_id': sid, 'gauge_id': key(gid), 'distance_km': float(dist), 'N': len(x),
                          'TP': int((predicted & truth).sum()), 'FP': int((predicted & ~truth).sum()),
                          'FN': int((~predicted & truth).sum()), 'TN': int((~predicted & ~truth).sum()),
                          'n_bit8': int(S['bit8'][-1].sum())})
        if k % 150 == 0:
            print(f'Gauge validated {k}/{len(mapping)} midpoint-matched physical CMLs', flush=True)
    S = {k: np.concatenate(v) for k, v in S.items()}
    res = summarise(S['x'], S['y'], S['pred'])
    res['population'] = {'assigned_midpoint_CMLs': len(mapping), 'matched_CMLs_in_final_rainfall': rain_membership,
                         'contributing_physical_CMLs': len({d['cml_id'] for d in audit}),
                         'contributing_directional_sublinks': len(audit),
                         'distinct_gauges': len({d['gauge_id'] for d in audit}),
                         'rainvalid_pre_intersection': pre_rain, 'detvalid_pre_intersection': pre_det}
    if bqc is not None:
        m = ~S['bit8']
        res['without_bit8'] = summarise(S['x'][m], S['y'][m], S['pred'][m])
        res['without_bit8']['n_samples_removed'] = int(S['bit8'].sum())
    return res, S, audit, gt, cids


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for ch in iter(lambda: f.read(1 << 24), b''):
            h.update(ch)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
    ap.add_argument('--cml-nc', default=None, help='QC file (default: qc of the configuration)')
    ap.add_argument('--rain5-nc', default=None, help='5-minute rainfall (default: output of scripts/derive_5min.py)')
    ap.add_argument('--gauge-nc', default=None); ap.add_argument('--out-dir', default=None)
    ap.add_argument('--blackout-qc', default=None, help='QC file for the bit-8 marking (default: the QC file)')
    ap.add_argument('--no-blackout-split', action='store_true'); ap.add_argument('--tag', default='graf')
    a = ap.parse_args()
    import sys
    sys.path.insert(0, str(ROOT / 'src'))
    from essd_cz.config import load_config
    from essd_cz.rain5 import derived_5min_name
    cfg = load_config(a.config)
    ds = cfg['datasets']
    a.cml_nc = a.cml_nc or str(ds['qc'])
    a.rain5_nc = a.rain5_nc or str(cfg['output_dir'] / 'rainfall_5min' / derived_5min_name(ds['rain'].name))
    a.gauge_nc = a.gauge_nc or str(ds['gauge'])
    a.out_dir = a.out_dir or str(cfg['output_dir'] / 'validation')
    if not a.no_blackout_split:
        a.blackout_qc = a.blackout_qc or a.cml_nc
    t0 = time.time()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    if (out / f'{a.tag}_summary.json').exists():
        raise SystemExit('refusing to overwrite')
    spec = importlib.util.spec_from_file_location('essd_exact_graf', HELPER)
    g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
    cml = xr.open_dataset(a.cml_nc); rain = xr.open_dataset(a.rain5_nc); gauge = xr.open_dataset(a.gauge_nc)
    bqc = xr.open_dataset(a.blackout_qc) if a.blackout_qc else None
    res, S, audit, gt, cids = gauge_validation(cml, rain, gauge, g, START, END, bqc=bqc)
    res['inputs'] = {k: {'path': str(Path(p).resolve()), 'sha256': sha256(p)} for k, p in
                     (('cml', a.cml_nc), ('rain5', a.rain5_nc), ('gauge', a.gauge_nc), ('helper', HELPER), ('script', __file__))
                     + ((('blackout_qc', a.blackout_qc),) if a.blackout_qc else ())}
    res['method'] = {'reference': 'nearest CHMI gauge to CML midpoint <=5km; end-labelled 10min *6; wet>=0.1',
                     'rainfall': 'two complete 5min bins at t-10 and t-5', 'detection': 'Graf RSTD on the CML file 1min trsl, Q80*1.12; all 10min valid',
                     'validation_mask': 'identical (link,gauge,sublink,10min timestamp) for confusion and statistics', 'in_range': '0..150 mm/h inclusive'}
    res['elapsed_s'] = time.time() - t0
    (out / f'{a.tag}_summary.json').write_text(json.dumps(res, indent=2))
    pd.DataFrame(audit).to_csv(out / f'{a.tag}_audit.csv', index=False)
    np.savez_compressed(out / f'{a.tag}_samples.npz', cml_ids=np.asarray(cids), gauge_time_ns=gt.asi8, **S)
    print(json.dumps({k: v for k, v in res.items() if k != 'inputs'}, indent=2))


if __name__ == '__main__':
    main()
