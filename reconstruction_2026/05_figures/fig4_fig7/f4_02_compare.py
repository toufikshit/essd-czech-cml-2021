#!/usr/bin/env python3
"""Stage 5B / Agent 3b - Fig. 4 numerical old-vs-new comparison for the event CML.

For the event window 2021-07-13 12:00 .. 2021-07-14 01:00 (producer event_windows) and
the given CML(s) this reports, old (full-history QC + full-history 5-min Graf rain) vs
new (QC REV1 + Stage 4 5-min / 1-min rain):
  * panel (a): TRSL per sublink (finite minutes, max, identical minutes, max |diff|);
  * panel (b): the CML series exactly as event_plot builds it (mean of valid sublinks of
    R_sublink_5min): finite steps, max, event accumulation (sum/12);
  * Graf detection: graf_wet_5min / graf_valid_5min counts, graf_threshold_db;
  * Stage 4 1-min trsl_source_flag_1min bits 8 (blackout reconstructed) and 16 (interpolated),
    QC REV1 qc_flag bits, blackout/interpolated minutes in the 5-min file;
  * selection score of event_plot (min of max CML, radar, gauge rate) for each CML, using the
    producer's own functions (hav, positions, select_rain_rate) imported unchanged.
Usage: python -I f4_02_compare.py [--cml 45035_45097_2740 28038_28079_2496] [--json out.json]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from f4_01_event_plot_run import ORIGINAL, STAGE4, PRODUCER  # noqa: E402

REC = Path(r'C:\Users\taouf\ESSD_RELEASE\FINAL_RECONSTRUCTION_20261027')
RAIN1 = REC / 'stage4' / 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc'
T0, T1 = pd.Timestamp('2021-07-13 12:00'), pd.Timestamp('2021-07-14 01:00')


def load_producer():
    spec = importlib.util.spec_from_file_location('essd_fig456_producer', PRODUCER)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def series(paths, cid, g):
    out = {}
    with xr.open_dataset(paths['cml']) as cml, xr.open_dataset(paths['rain']) as rain, \
         xr.open_dataset(paths['radar']) as radar, xr.open_dataset(paths['gauge']) as gauge:
        cp = g.positions(cml); rp = g.positions(rain); dp = g.positions(radar)
        if cid not in cp or cid not in rp:
            return {'present_in_cml': cid in cp, 'present_in_rain': cid in rp}
        k = cp[cid]
        st = pd.DatetimeIndex(cml.time.values); sel = np.flatnonzero((st >= T0) & (st <= T1))
        trsl = g.numeric(cml['trsl'].isel(cml_id=k, time=sel), 'time', 'sublink_id')
        out['trsl'] = trsl; out['t1'] = st[sel]
        if 'qc_flag' in cml:
            out['qc_flag'] = g.numeric(cml['qc_flag'].isel(cml_id=k, time=sel), 'time', 'sublink_id')
        rate = g.select_rain_rate(rain).sel(time_5min=slice(T0, T1))
        rv = g.numeric(rate.isel(cml_id=rp[cid]), 'sublink_id', 'time_5min')
        valid = np.isfinite(rv); cnt = valid.sum(axis=0)
        cr = np.full(rv.shape[1], np.nan)
        np.divide(np.where(valid, rv, 0).sum(axis=0), cnt, out=cr, where=cnt > 0)
        out['cr'] = cr; out['rv_sub'] = rv; out['t5'] = pd.DatetimeIndex(rate.time_5min.values)
        rr = rain.sel(time_5min=slice(T0, T1)).isel(cml_id=rp[cid])
        for v in ('graf_wet_5min', 'graf_valid_5min', 'blackout_minutes_5min', 'interpolated_minutes_5min'):
            if v in rr:
                out[v] = g.numeric(rr[v], 'time_5min', 'sublink_id') if set(rr[v].dims) == {'time_5min', 'sublink_id'} else None
        out['graf_threshold_db'] = np.asarray(rain['graf_threshold_db'].isel(cml_id=rp[cid]).values, float).tolist()
        # selection score as in event_plot
        t_ra = pd.DatetimeIndex(radar.time.values); rsel = np.flatnonzero((t_ra >= T0) & (t_ra <= T1))
        radv = np.asarray(radar['rr'].isel(cml_id=dp[cid], time=rsel).values, float).reshape(-1) * 12 if cid in dp else np.array([np.nan])
        st_dim = next(d for d in gauge['rr'].dims if d != 'time')
        lat = np.asarray(gauge['lat'].values, float); lon = np.asarray(gauge['lon'].values, float)
        mlat = (float(cml.site_a_latitude[k]) + float(cml.site_b_latitude[k])) / 2
        mlon = (float(cml.site_a_longitude[k]) + float(cml.site_b_longitude[k])) / 2
        d = g.hav(mlat, mlon, lat, lon); d[~np.isfinite(d)] = np.inf; gi = int(d.argmin())
        t_g = pd.DatetimeIndex(gauge.time.values); gsel = np.flatnonzero((t_g >= T0) & (t_g <= T1))
        gv = np.asarray(gauge['rr'].isel({st_dim: gi, 'time': gsel}).values, float).reshape(-1) * 6
        out['score'] = [float(np.nanmax(cr)), float(np.nanmax(radv)), float(np.nanmax(gv))]
        out['gauge_km'] = float(d[gi])
    return out


def summarize(o, n, cid):
    r = {'cml_id': cid}
    if 'trsl' not in o or 'trsl' not in n:
        r['note'] = f'old present {"trsl" in o}, new present {"trsl" in n}'
        return r
    if len(o['t1']) != len(n['t1']) or not (o['t1'] == n['t1']).all():
        raise RuntimeError('time axes differ')
    for s in range(o['trsl'].shape[1]):
        a, b = o['trsl'][:, s], n['trsl'][:, s]
        both = np.isfinite(a) & np.isfinite(b)
        r[f'sublink{s+1}_trsl'] = dict(old_finite=int(np.isfinite(a).sum()), new_finite=int(np.isfinite(b).sum()),
                                       identical_minutes=int((both & (a == b)).sum()),
                                       max_abs_diff_db=float(np.abs(a[both] - b[both]).max()) if both.any() else None,
                                       old_max=float(np.nanmax(a)), new_max=float(np.nanmax(b)))
    if 'qc_flag' in n:
        f = n['qc_flag'].astype(int)
        r['new_qc_flag_minutes'] = {f'bit{b}': int(((f & b) > 0).sum()) for b in (1, 2, 4, 8, 16)}
    for lab, d in (('old', o), ('new', n)):
        cr = d['cr']
        r[f'{lab}_panel_b'] = dict(finite_5min=int(np.isfinite(cr).sum()), n_5min=int(len(cr)), max_mmph=float(np.nanmax(cr)),
                                   acc_mm=float(np.nansum(cr) / 12.0),
                                   sublink_acc_mm=[float(np.nansum(x) / 12.0) for x in d['rv_sub']],
                                   graf_wet_5min_per_sublink=(np.nansum(d['graf_wet_5min'] == 1, axis=0).astype(int).tolist()
                                                              if d.get('graf_wet_5min') is not None else None),
                                   graf_valid_5min_per_sublink=(np.nansum(d['graf_valid_5min'] == 1, axis=0).astype(int).tolist()
                                                                if d.get('graf_valid_5min') is not None else None),
                                   graf_threshold_db=d['graf_threshold_db'], selection_score_inputs=d['score'],
                                   selection_score=min(d['score']), gauge_km=d['gauge_km'])
    both = np.isfinite(o['cr']) & np.isfinite(n['cr'])
    r['panel_b_diff'] = dict(identical_5min=int((both & (o['cr'] == n['cr'])).sum()),
                             max_abs_diff_mmph=float(np.abs(o['cr'][both] - n['cr'][both]).max()),
                             mean_diff_mmph=float((n['cr'][both] - o['cr'][both]).mean()))
    if n.get('blackout_minutes_5min') is not None:
        r['new_5min_blackout_minutes'] = int(np.nansum(n['blackout_minutes_5min']))
        r['new_5min_interpolated_minutes'] = int(np.nansum(n['interpolated_minutes_5min']))
    return r


def flags_1min(cid):
    with xr.open_dataset(RAIN1) as r:
        ids = [str(v) for v in r.cml_id.values]
        if cid not in ids:
            return None
        t = pd.DatetimeIndex(r.time_1min.values); sel = np.flatnonzero((t >= T0) & (t <= T1))
        f = np.asarray(r['trsl_source_flag_1min'].isel(cml_id=ids.index(cid), time_1min=sel).values).astype(int)
        w = np.asarray(r['graf_wet_1min'].isel(cml_id=ids.index(cid), time_1min=sel).values)
        R1 = np.asarray(r['R_sublink_1min'].isel(cml_id=ids.index(cid), time_1min=sel).values, float)
        b8 = (f & 8) > 0
        return dict(minutes=int(len(sel)), bit8_blackout_reconstructed=int(b8.sum()), bit16_interpolated=int(((f & 16) > 0).sum()),
                    graf_wet_1min_per_sublink=np.nansum(w == 1, axis=0).astype(int).tolist(),
                    rain_from_bit8_minutes_mm=float(np.nansum(np.where(b8, R1, 0)) / 60.0),
                    max_rate_at_bit8_mmph=float(np.nanmax(np.where(b8, R1, np.nan))) if b8.any() else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cml', nargs='+', default=['45035_45097_2740', '28038_28079_2496'])
    ap.add_argument('--json', type=Path, default=Path(__file__).resolve().parent.parent / 'reproduction' / 'fig4_old_vs_new.json')
    a = ap.parse_args()
    g = load_producer()
    res = []
    for cid in a.cml:
        o = series(ORIGINAL, cid, g); n = series(STAGE4, cid, g)
        s = summarize(o, n, cid); s['stage4_1min_flags'] = flags_1min(cid)
        res.append(s)
        print(json.dumps(s, indent=1))
    a.json.write_text(json.dumps(res, indent=2) + '\n')


if __name__ == '__main__':
    main()
