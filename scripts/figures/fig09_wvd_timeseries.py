"""Figure 9: time series of CML water vapour density (WVD) and the station reference for one sublink.

Uses the output of scripts/wvd_validation.py (the validation pairs behind Fig. 8 and Table S13) and the same
functions: the V6 WVD of the Graf branch (wvd_internal_anchor), two complete 5-minute values per 10-minute interval
(strict_5_to_10), and the station WVD from relative humidity and temperature of the assigned CHMI station
(rh_t_to_wvd), both from original/code/wvd/21_validate_wvd_v6_vs_original.py. No smoothing, offset or time shift.

Selection (fixed rule, not chosen for good agreement):
  * population: the 45 sublinks of Fig. 8 (select_population of plot_wvp_support_hexbins_no_titles.py, 71-76 and
    81-86 GHz, L >= 4 km, delta-rho-q <= 0.20 g m-3);
  * sublink: the one whose RMSE over all its pairs is closest to the median RMSE of the 45 sublinks;
  * period: the 14-day window (starting at 00:00 UTC) with the most matched pairs of that sublink.
The script checks that the rebuilt pairs equal the stored validation pairs of the sublink.
Output: <output_dir>/figures/fig09/f09_wvd_timeseries.png/.pdf (style of Fig. 4) and f09_selection.json (selection,
window statistics and the cause of every 10-minute step without a pair).
Usage: python scripts/figures/fig09_wvd_timeseries.py [--config config/user_paths.json] [--days 14]
       [--cml-id ID --sublink-id N --start YYYY-MM-DD]   (to plot another sublink or period)
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from essd_cz.config import ORIGINAL, extra, load_config  # noqa: E402

WVD = ORIGINAL / 'wvd'


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(b, float) - np.asarray(a, float)) ** 2)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
    ap.add_argument('--days', type=int, default=14)
    ap.add_argument('--cml-id', default=None)
    ap.add_argument('--sublink-id', type=int, default=None)
    ap.add_argument('--start', default=None, help='first day of the window (default: window with most pairs)')
    a = ap.parse_args()
    cfg = load_config(a.config)
    src = cfg['output_dir'] / 'wvd' / 'WVD_V6_QUANTIZATION_BINNED_VALIDATION'
    if not (src / 'validation_pairs_graf.npz').is_file():
        raise SystemExit(f'{src} has no validation pairs (run scripts/wvd_validation.py first)')
    out = cfg['output_dir'] / 'figures' / 'fig09'
    out.mkdir(parents=True, exist_ok=True)

    plotmod = load(WVD / 'plot_wvp_support_hexbins_no_titles.py', 'wvp_support_plot')
    helper = load(WVD / '21_validate_wvd_v6_vs_original.py', 'v6_validation_helper')
    links = pd.read_csv(src / 'validation_population_links.csv')
    with np.load(src / 'validation_pairs_graf.npz') as z:
        ref, est, link_idx = z['ref'], z['est'], z['link_idx']
    _, _, use, met = plotmod.select_population(links, link_idx.astype(np.int64), ref.astype(float),
                                               est.astype(float), 71.0, 86.0)
    print(f"Fig. 8 population: {met['n_sublinks']} sublinks, {met['n_links']} CMLs, N {met['N']}")
    use = use.copy()
    use['rmse'] = [rmse(ref[link_idx == i], est[link_idx == i]) for i in use.link_idx]
    use['n'] = [int((link_idx == i).sum()) for i in use.link_idx]
    if a.cml_id:
        row = use[(use.cml_id == a.cml_id) & (use.sublink_id == a.sublink_id)].iloc[0]
        rule = 'chosen with --cml-id/--sublink-id'
    else:
        med = float(use.rmse.median())
        row = use.assign(d=(use.rmse - med).abs()).sort_values(['d', 'n'], ascending=[True, False]).iloc[0]
        rule = f'RMSE closest to the median RMSE of the {len(use)} sublinks ({med:.2f} g m-3)'
    li = int(row.link_idx)
    print(f'sublink {row.cml_id} / {row.sublink_id}: {row.frequency_ghz} GHz, {row.length_km:.2f} km, '
          f'station {row.station_id} ({row.station_distance_to_path_km:.1f} km), RMSE {row.rmse:.2f}, n {row.n}; {rule}')

    # rebuild the pairs of this sublink with their times, exactly as script 33 does
    ds = cfg['datasets']
    with xr.open_dataset(extra(a.config, 'wvd_v6_nc')) as v6, xr.open_dataset(ds['humidity']) as hum, \
            xr.open_dataset(ds['temperature']) as temp:
        t5 = helper.decode_retrieval_time(v6)
        mi = [str(x).lower() for x in helper.METHODS].index('graf')
        vi = int(row.series_index)
        w5 = np.asarray(v6['wvd_internal_anchor'].isel(method=mi, series=vi).values, dtype=float)
        s5 = np.asarray(v6['retrieval_status'].isel(method=mi, series=vi).values)
        tb5 = np.asarray(v6['hardware_transition_buffer_5min'].isel(series=vi).values)
        hs5 = np.asarray(v6['hardware_state_id'].isel(series=vi).values)
        t10, w10 = helper.strict_5_to_10(t5, w5, s5)
        sid = str(row.station_id)
        hs = pd.Series(np.asarray(hum['H'].sel(station=sid).values, dtype=float),
                       index=pd.DatetimeIndex(pd.to_datetime(hum['H'].sel(station=sid)['time'].values)))
        ts = pd.Series(np.asarray(temp['T'].sel(station=sid).values, dtype=float),
                       index=pd.DatetimeIndex(pd.to_datetime(temp['T'].sel(station=sid)['time'].values)))
        truth = helper.rh_t_to_wvd(hs.reindex(t10).to_numpy(float), ts.reindex(t10).to_numpy(float))
    common = np.isfinite(truth) & np.isfinite(w10)
    same = (np.array_equal(truth[common].astype(np.float32), ref[link_idx == li])
            and np.array_equal(w10[common].astype(np.float32), est[link_idx == li]))
    print('rebuilt pairs equal the stored validation pairs:', same)
    if not same:
        raise SystemExit('rebuilt pairs differ from the validation pairs')

    paired = pd.DataFrame({'cml': np.where(common, w10, np.nan), 'station': np.where(common, truth, np.nan)},
                          index=t10)
    days = pd.date_range(paired.index[0].normalize(), paired.index[-1].normalize() - pd.Timedelta(days=a.days - 1))
    counts = pd.Series([int(paired.loc[d:d + pd.Timedelta(days=a.days) - pd.Timedelta(minutes=10), 'cml'].notna().sum())
                        for d in days], index=days)
    start = pd.Timestamp(a.start) if a.start else counts.idxmax()
    end = start + pd.Timedelta(days=a.days) - pd.Timedelta(minutes=10)
    w = paired.loc[start:end]
    ok = w.cml.notna()
    sel = {'cml_id': row.cml_id, 'sublink_id': int(row.sublink_id), 'frequency_ghz': float(row.frequency_ghz),
           'length_km': float(row.length_km), 'station_id': sid,
           'station_distance_to_path_km': float(row.station_distance_to_path_km),
           'selection_rule': rule, 'sublink_rmse_all_pairs_gm3': float(row.rmse), 'sublink_pairs': int(row.n),
           'window_start_utc': str(start), 'window_end_utc': str(end), 'window_rule':
           'given with --start' if a.start else f'{a.days}-day window with the most pairs',
           'window_pairs': int(ok.sum()), 'window_rmse_gm3': rmse(w.station[ok], w.cml[ok]),
           'window_bias_gm3': float((w.cml[ok] - w.station[ok]).mean()),
           'window_station_min_max_gm3': [float(w.station[ok].min()), float(w.station[ok].max())]}

    # cause of every 10-min step of the window without a pair. A V6 value exists only for "safe dry" minutes: valid
    # total loss, known Graf state, no Graf-wet minute within the wet buffer (6 h), no hardware transition (12 h).
    import netCDF4
    with netCDF4.Dataset(ds['raw']) as f:
        f.set_auto_mask(False)
        j = [str(x) for x in f['cml_id'][:]].index(str(row.cml_id))
        traw = pd.to_datetime(f['time'][:], unit='m', origin='2020-01-01')
        tx = f['tx'][:, int(row.sublink_id), j].astype(float)
        rx = f['rx'][:, int(row.sublink_id), j].astype(float)
    okraw = pd.Series(np.isfinite(tx) & np.isfinite(rx) & (tx > -99) & (tx != 255) & (rx > -99.9) & (rx < -20),
                      index=traw)
    s5f = pd.DataFrame({'trans': tb5, 'state': hs5}, index=t5)
    st_ok = pd.Series(np.isfinite(truth), index=t10)
    causes = {}
    for t in w.index[~ok]:
        if not np.isfinite(paired.loc[t, 'cml']) and np.isfinite(w10[t10.get_loc(t)]):
            c = 'station value missing'
        else:
            b = s5f.loc[t - pd.Timedelta(minutes=10):t - pd.Timedelta(minutes=5)]
            if (b.trans > 0).any() or (b.state < 0).any():
                c = 'hardware transition'
            elif not okraw.loc[t - pd.Timedelta(minutes=10):t - pd.Timedelta(minutes=1)].all():
                c = 'missing raw signal minutes'
            else:
                c = 'classified as wet by the Graf method, or within 6 h of a wet period'
            if not st_ok[t]:
                c += '; station value also missing'
        causes[c] = causes.get(c, 0) + 1
    sel['window_10min_steps'] = int(len(w))
    sel['window_gap_causes_10min_steps'] = causes
    print(json.dumps(sel, indent=1))

    # style of Fig. 4 (style() and save() of the event-plot producer)
    prod = load(ORIGINAL / 'figures' / 'rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py',
                'essd_fig4_style')
    from matplotlib.ticker import AutoMinorLocator, MaxNLocator
    prod.style()
    fig, ax = plt.subplots(figsize=(8.6, 3.0), constrained_layout=True)
    ax.plot(w.index, w.cml, color='#1F4E79', lw=1.3, label='CML')
    ax.plot(w.index, w.station, color='#C05E25', lw=1.2, ls='--', label='Station')
    ax.set_ylabel('WVD [g m$^{-3}$]')
    ax.set_xlabel('Time [UTC]')
    ax.legend(frameon=False, ncol=2, loc='upper right', fontsize=8.5)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_minor_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_position(('outward', 4))
    ax.set_xlim(start, end); ax.margins(x=0)
    lo, hi = ax.get_ylim()
    ticks = MaxNLocator(nbins=5).tick_values(lo, hi)
    ax.set_ylim(float(ticks[0]), float(ticks[-1])); ax.set_yticks(ticks)
    ax.tick_params(which='both', direction='out', top=False, right=False)
    prod.save(fig, out / 'f09_wvd_timeseries')
    (out / 'f09_selection.json').write_text(json.dumps(sel, indent=1))


if __name__ == '__main__':
    main()
