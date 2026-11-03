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
Output: <output_dir>/figures/fig09/f09_wvd_timeseries.png/.pdf and f09_selection.json.
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
    print(json.dumps(sel, indent=1))

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
    fig, ax = plt.subplots(figsize=(10, 3.6))
    ax.plot(w.index, w.station, color='#1f6fb0', lw=1.1, label='Station (CHMI, from RH and T)')
    ax.plot(w.index, w.cml, color='#d95f02', lw=1.1, label='CML')
    ax.set_ylabel('WVD (g m$^{-3}$)')
    ax.set_xlabel('Time (UTC)')
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
    ax.set_xlim(start, end)
    ax.grid(True, lw=0.4, color='0.85')
    ax.legend(loc='upper left', frameon=False, ncol=2)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / 'f09_wvd_timeseries.png', dpi=300)
    fig.savefig(out / 'f09_wvd_timeseries.pdf')
    (out / 'f09_selection.json').write_text(json.dumps(sel, indent=1))
    print('saved', out / 'f09_wvd_timeseries.png')


if __name__ == '__main__':
    main()
