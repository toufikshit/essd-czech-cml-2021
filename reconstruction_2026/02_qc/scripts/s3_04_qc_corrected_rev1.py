"""Stage 3, revision 1 (author decision 2026-10-28): corrected QC candidate for the ESSD Czech CML 2021 dataset (author-approved settings, 2026-10-28).

Input (read-only): the full-history RAW NetCDF (approved population of the Stage 2 records). Output: a NEW QC
candidate NetCDF plus per-sublink results. No existing product is modified.

Approved settings (stage3/STAGE3_QC_SETTINGS_PROPOSAL.md, author decision in PROCESSING_LOG.md, 2026-10-28):
  * period 2021-04-01 00:00 .. 2021-10-31 23:59 UTC (release window), 1-minute data;
  * population: 2574 approved RAW CMLs; the 8 CMLs without frequency are excluded before QC; 2566 assessed;
  * P5: the 15 TX = 0 directions are non-operational (no valid observations); the opposite direction is kept;
  * P6: the mixed-placeholder minutes are masked (TX and RX);
  * invalid values: TX <= -99 dBm or TX = 255 or |x| >= 1e35; RX <= -99.9 dBm or RX > -20 dBm;
  * low-RSL plateau: centred 3-minute RSL window with >= 2 valid samples, max < -85 dBm and population std
    (ddof = 0) < 0.5 dB; TX and RX masked with +/-2 minutes padding;
  * blackout, option C: an RSL gap (whole run of missing RX) of at most 60 minutes with RSL < -65 dBm on at least
    one edge. TX and RX stay missing; the total loss of the documented method (series max TX minus series min RX,
    standalone v3 / Blettner et al. 2023) is written to trsl only and flagged. REVISION 1: a gap that contains a
    minute removed by the plateau filter or masked by P6 is not reconstructed; it stays missing in trsl and keeps
    its flags (counted as rejected). Reconstructed minutes are excluded from the fluctuation statistics;
  * short gaps: linear interpolation of total loss over runs of <= 5 missing minutes that are missing in RAW (not
    created by a QC mask), with observed (not reconstructed) total loss on both sides; flagged;
  * fluctuation: per sublink, centred rolling windows of 300 and 60 minutes, population std; a window is valid only
    with 300 (5 h) or 60 (1 h) valid minutes; the sublink fails if more than 10 % of valid 5-h windows have
    std > 2 dB or more than 33 % of valid 1-h windows have std > 0.8 dB, or if it has no valid window;
  * detection limit per sublink: R_min = (0.33 dB / (k L))^(1/alpha), ITU-R P.838-3 (2005) k and alpha from
    pycomlink.processing.k_R_relation.a_b(f_GHz, 'V', approx_type='ITU_2005'); fail if R_min >= 2 mm/h or NaN;
  * a sublink is excluded if it is a P5 direction, has no observed total loss in the period, fails the
    fluctuation test or the detection limit; excluded sublinks are set to missing in the product;
  * a CML is retained if at least one sublink passes; CMLs without a passing sublink are not written.

Pass 1 decides every sublink and writes the per-sublink table; pass 2 recomputes the retained CMLs with the same
kernel, checks that the decisions are identical, and writes the NetCDF.

Usage:
  python -I s3_01_qc_corrected.py --out-dir ../  [--raw-nc ...] [--stage2 ...] [--only-cmls id1,id2] [--tag test]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from datetime import datetime, timezone

import netCDF4 as nc
import numpy as np
import pandas as pd
from numba import njit, prange
from pycomlink.processing.k_R_relation import a_b

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RAW_DEFAULT = r'C:\Users\taouf\ESSD_RELEASE\extracted\full\ESSD_CZ_2021_HANDOFF\data\source_full_history\raw\CZ_CML_RAW_1MIN_20200101_20211231.nc'
OLD_QC_DEFAULT = r'C:\Users\taouf\ESSD_RELEASE\zenodo_staging\primary\CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc'

T0 = pd.Timestamp('2020-01-01 00:00')
WIN_START = pd.Timestamp('2021-04-01 00:00')
WIN_END = pd.Timestamp('2021-10-31 23:59')
IA = int((WIN_START - T0) / pd.Timedelta(minutes=1))
IB = int((WIN_END - T0) / pd.Timedelta(minutes=1)) + 1

P = dict(
    tx_min_dbm=-99.0, tx_sentinel=255.0, rx_min_dbm=-99.9, rx_max_dbm=-20.0, huge_abs=1.0e35,
    plateau_max_dbm=-85.0, plateau_std_db=0.5, plateau_min_valid=2, plateau_pad=2,
    blackout_low_dbm=-65.0, blackout_max_gap=60,
    short_gap_max=5,
    long_window=300, long_threshold_db=2.0, long_max_fraction=0.10, long_min_valid=300,
    short_window=60, short_threshold_db=0.8, short_max_fraction=0.33, short_min_valid=60,
    a_quantization_db=0.33, detection_limit_max_mmhr=2.0, itu_approximation='ITU_2005', polarization='V',
)

F_INVALID, F_P6, F_PLATEAU, F_BLACKOUT, F_INTERP = 1, 2, 4, 8, 16
STATUS = {0: 'retained', 1: 'excluded_p5_tx_zero_nonoperational', 2: 'excluded_no_observed_total_loss_in_period',
          3: 'excluded_fluctuation', 4: 'excluded_detection_limit', 5: 'excluded_fluctuation_and_detection_limit'}
STAT_NAMES = ['n_tx_raw', 'n_rx_raw', 'n_tl_raw', 'n_tx_invalid', 'n_rx_invalid', 'n_p6_masked', 'n_plateau_masked',
              'n_tl_observed_after_masks', 'n_blackout_gaps', 'n_blackout_minutes', 'n_interp_gaps', 'n_interp_minutes',
              'frac_long', 'n_valid_long', 'n_exceed_long', 'frac_short', 'n_valid_short', 'n_exceed_short',
              'n_blackout_rejected_gaps', 'n_blackout_rejected_minutes']


@njit(cache=True)
def rolling_exceed(x, window, min_valid, thr):
    """Centred rolling population std (pandas-like bounds); a window is valid with >= min_valid finite values."""
    n = x.size
    m = 0.0
    k = 0
    for i in range(n):
        if math.isfinite(x[i]):
            m += x[i]
            k += 1
    if k == 0:
        return math.nan, 0, 0
    m /= k
    cs = np.zeros(n + 1)
    cs2 = np.zeros(n + 1)
    cc = np.zeros(n + 1, dtype=np.int64)
    for i in range(n):
        cs[i + 1] = cs[i]
        cs2[i + 1] = cs2[i]
        cc[i + 1] = cc[i]
        v = x[i]
        if math.isfinite(v):
            d = v - m
            cs[i + 1] += d
            cs2[i + 1] += d * d
            cc[i + 1] += 1
    left = window // 2
    right = window - left
    ng = 0
    ne = 0
    for i in range(n):
        a = max(i - left, 0)
        b = min(i + right, n)
        cnt = cc[b] - cc[a]
        if cnt >= min_valid:
            mu = (cs[b] - cs[a]) / cnt
            var = (cs2[b] - cs2[a]) / cnt - mu * mu
            if var < 0.0:
                var = 0.0
            ng += 1
            if math.sqrt(var) > thr:
                ne += 1
    if ng == 0:
        return math.nan, 0, 0
    return ne / ng, ng, ne


@njit(cache=True)
def qc_series(tx_in, rx_in, p6):
    n = tx_in.size
    tx = tx_in.copy()
    rx = rx_in.copy()
    flag = np.zeros(n, dtype=np.uint8)
    st = np.full(20, math.nan)
    n_tx = 0; n_rx = 0; n_tl = 0
    for i in range(n):
        ft = math.isfinite(tx[i]); fr = math.isfinite(rx[i])
        if ft:
            n_tx += 1
        if fr:
            n_rx += 1
        if ft and fr:
            n_tl += 1
    st[0] = n_tx; st[1] = n_rx; st[2] = n_tl
    # P6 mixed-placeholder minutes
    n_p6 = 0
    for i in range(n):
        if p6[i]:
            if math.isfinite(tx[i]) or math.isfinite(rx[i]):
                n_p6 += 1
            tx[i] = math.nan; rx[i] = math.nan; flag[i] |= 2
    st[5] = n_p6
    # invalid values
    nbt = 0; nbr = 0
    for i in range(n):
        t = tx[i]
        if math.isfinite(t) and (t <= -99.0 or abs(t - 255.0) < 1e-6 or abs(t) >= 1.0e35):
            tx[i] = math.nan; flag[i] |= 1; nbt += 1
        r = rx[i]
        if math.isfinite(r) and (r <= -99.9 or r > -20.0 or abs(r) >= 1.0e35):
            rx[i] = math.nan; flag[i] |= 1; nbr += 1
    st[3] = nbt; st[4] = nbr
    # low-RSL plateau (>= 2 valid samples in the centred 3-minute window)
    base = np.zeros(n, dtype=np.uint8)
    for i in range(n):
        a = max(i - 1, 0); b = min(i + 2, n)
        c = 0; s = 0.0; s2 = 0.0; vmax = -math.inf
        for j in range(a, b):
            v = rx[j]
            if math.isfinite(v):
                c += 1; s += v; s2 += v * v
                if v > vmax:
                    vmax = v
        if c >= 2:
            mu = s / c
            var = s2 / c - mu * mu
            if var < 0.0:
                var = 0.0
            if vmax < -85.0 and math.sqrt(var) < 0.5:
                base[i] = 1
    n_pl = 0
    for i in range(n):
        a = max(i - 2, 0); b = min(i + 3, n)
        hit = False
        for j in range(a, b):
            if base[j] == 1:
                hit = True
                break
        if hit:
            if math.isfinite(tx[i]) or math.isfinite(rx[i]):
                n_pl += 1
                flag[i] |= 4
            tx[i] = math.nan; rx[i] = math.nan
    st[6] = n_pl
    # total loss of observed values
    tl = np.full(n, math.nan)
    n_obs = 0
    for i in range(n):
        if math.isfinite(tx[i]) and math.isfinite(rx[i]):
            tl[i] = tx[i] - rx[i]
            n_obs += 1
    st[7] = n_obs
    # blackout, option C: reconstructed total loss in trsl only
    txmax = -math.inf; rxmin = math.inf
    for i in range(n):
        if math.isfinite(tx[i]) and tx[i] > txmax:
            txmax = tx[i]
        if math.isfinite(rx[i]) and rx[i] < rxmin:
            rxmin = rx[i]
    nbg = 0; nbm = 0; nrg = 0; nrm = 0
    if math.isfinite(txmax) and math.isfinite(rxmin):
        recon = txmax - rxmin
        i = 0
        while i < n:
            if math.isfinite(rx[i]):
                i += 1
                continue
            j = i + 1
            while j < n and not math.isfinite(rx[j]):
                j += 1
            gap = j - i
            if gap <= 60:
                left_low = i > 0 and math.isfinite(rx[i - 1]) and rx[i - 1] < -65.0
                right_low = j < n and math.isfinite(rx[j]) and rx[j] < -65.0
                has_masked = False  # revision 1: P6-masked (bit 2) or plateau-removed (bit 4) minute in the gap
                for q in range(i, j):
                    if p6[q] or (flag[q] & 4) != 0:
                        has_masked = True
                        break
                if left_low or right_low:
                    if has_masked:
                        nrg += 1; nrm += gap
                    else:
                        for q in range(i, j):
                            tl[q] = recon
                            flag[q] |= 8
                        nbg += 1; nbm += gap
            i = j
    st[8] = nbg; st[9] = nbm; st[18] = nrg; st[19] = nrm
    # short-gap interpolation of RAW gaps only
    nig = 0; nim = 0
    i = 0
    while i < n:
        if math.isfinite(tl[i]):
            i += 1
            continue
        j = i + 1
        while j < n and not math.isfinite(tl[j]):
            j += 1
        gap = j - i
        if gap <= 5 and i > 0 and j < n and (flag[i - 1] & 8) == 0 and (flag[j] & 8) == 0:
            raw_only = True
            for q in range(i, j):
                if (flag[q] & 7) != 0:
                    raw_only = False
                    break
            if raw_only:
                lv = tl[i - 1]; rv = tl[j]
                for q in range(gap):
                    tl[i + q] = lv + (rv - lv) * ((q + 1.0) / (gap + 1.0))
                    flag[i + q] |= 16
                nig += 1; nim += gap
        i = j
    st[10] = nig; st[11] = nim
    # fluctuation statistics without reconstructed blackout minutes
    x = tl.copy()
    for i in range(n):
        if flag[i] & 8:
            x[i] = math.nan
    fl, gl, el = rolling_exceed(x, 300, 300, 2.0)
    fs, gs, es = rolling_exceed(x, 60, 60, 0.8)
    st[12] = fl; st[13] = gl; st[14] = el; st[15] = fs; st[16] = gs; st[17] = es
    return tx, rx, tl, flag, st


@njit(cache=True, parallel=True)
def qc_batch(tx, rx, p6):
    """tx, rx: (series, time) float64; p6: (series, time) bool."""
    ns, n = tx.shape
    otx = np.empty((ns, n)); orx = np.empty((ns, n)); otl = np.empty((ns, n))
    ofl = np.empty((ns, n), dtype=np.uint8); ost = np.empty((ns, 20))
    for s in prange(ns):
        a, b, c, d, e = qc_series(tx[s], rx[s], p6[s])
        otx[s] = a; orx[s] = b; otl[s] = c; ofl[s] = d; ost[s] = e
    return otx, orx, otl, ofl, ost


def sha256(path, block=1 << 24):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for ch in iter(lambda: f.read(block), b''):
            h.update(ch)
    return h.hexdigest()


def decide(st, p5, rmin):
    """Return status code from kernel statistics."""
    if p5:
        return 1
    if st[7] == 0:
        return 2
    fl, gl, fs, gs = st[12], st[13], st[15], st[16]
    pass_f = gl > 0 and gs > 0 and fl <= P['long_max_fraction'] and fs <= P['short_max_fraction']
    pass_d = math.isfinite(rmin) and rmin < P['detection_limit_max_mmhr']
    if pass_f and pass_d:
        return 0
    if not pass_f and not pass_d:
        return 5
    return 3 if not pass_f else 4


def read_block(ds, cols):
    """Read the window for the given RAW cml indices (sorted) -> (series, time) float64 with series = (cml, sublink)."""
    cols = np.asarray(cols)
    runs = np.split(np.arange(cols.size), np.flatnonzero(np.diff(cols) != 1) + 1)  # contiguous column groups
    out = {}
    for v in ('tx', 'rx'):
        parts = [np.ma.filled(ds[v][IA:IB, :, int(cols[r[0]]):int(cols[r[-1]]) + 1].astype(np.float64), np.nan)
                 for r in runs]
        arr = np.concatenate(parts, axis=2)  # (time, sublink, cml)
        out[v] = np.ascontiguousarray(arr.transpose(2, 1, 0).reshape(-1, IB - IA))
    return out['tx'], out['rx']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--raw-nc', default=RAW_DEFAULT)
    ap.add_argument('--stage2', default=os.path.join(ROOT, 'stage2'))
    ap.add_argument('--out-dir', default=os.path.join(ROOT, 'stage3'))
    ap.add_argument('--only-cmls', default='')
    ap.add_argument('--tag', default='')
    ap.add_argument('--batch', type=int, default=48)
    a = ap.parse_args()
    t_start = time.time()
    sfx = ('_' + a.tag) if a.tag else ''
    out_nc = os.path.join(a.out_dir, f'QC_CORRECTED_CANDIDATE_REV1{sfx}.nc')
    out_sub = os.path.join(a.out_dir, f'QC_SUBLINK_RESULTS_REV1{sfx}.csv')
    if os.path.exists(out_nc):
        sys.exit(f'refusing to overwrite {out_nc}')

    pop = pd.read_csv(os.path.join(a.stage2, 'RAW_APPROVED_POPULATION.csv'))
    masks = pd.read_csv(os.path.join(a.stage2, 'QC_PREDEFINED_MINUTE_MASKS.csv'))
    elig = pop[pop.rainfall_ready_eligible_cml == 1]
    cml_list = list(dict.fromkeys(elig.cml_id))
    if a.only_cmls:
        want = a.only_cmls.split(',')
        cml_list = [c for c in cml_list if c in want]
    p5 = set(zip(pop[pop.sublink_tx_zero_nonoperational == 1].cml_id, pop[pop.sublink_tx_zero_nonoperational == 1].sublink))
    masks['w'] = masks.time_index - IA
    masks_in = masks[(masks.w >= 0) & (masks.w < IB - IA)]
    print(f'window index {IA}..{IB} ({IB - IA} min); CMLs assessed {len(cml_list)}; P5 sublinks {len(p5)}; '
          f'P6 minutes {len(masks)} of which in window {len(masks_in)}', flush=True)

    ds = nc.Dataset(a.raw_nc)
    raw_ids = np.array([str(x) for x in ds['cml_id'][:]])
    pos = {c: i for i, c in enumerate(raw_ids)}
    assert all(c in pos for c in cml_list)
    order = np.array(sorted(cml_list, key=lambda c: pos[c]))
    freq = np.ma.filled(ds['frequency'][:].astype(np.float64), np.nan)  # (sublink, cml) MHz
    length = np.ma.filled(ds['length'][:].astype(np.float64), np.nan)

    # detection limit per sublink (ITU-R P.838-3 2005, vertical)
    dl = {}
    for c in order:
        j = pos[c]
        for s in (0, 1):
            f_ghz = freq[s, j] / 1000.0
            L = length[j]
            if math.isfinite(f_ghz) and math.isfinite(L) and L > 0:
                k, al = a_b(f_ghz, 'V', approx_type='ITU_2005')
                k = float(np.ravel(k)[0]); al = float(np.ravel(al)[0])
                r = (P['a_quantization_db'] / (k * L)) ** (1.0 / al)
            else:
                k = al = r = math.nan
            dl[(c, s)] = (k, al, r)

    def p6_array(cmls):
        arr = np.zeros((len(cmls) * 2, IB - IA), dtype=np.bool_)
        for r in masks_in.itertuples():
            if r.cml_id in cmls:
                arr[list(cmls).index(r.cml_id) * 2 + int(r.sublink), int(r.w)] = True
        return arr

    # ---------------- pass 1: decisions ----------------
    rows = []
    for b0 in range(0, len(order), a.batch):
        cmls = order[b0:b0 + a.batch]
        cols = np.array([pos[c] for c in cmls])
        tx, rx = read_block(ds, cols)
        _, _, _, _, st = qc_batch(tx, rx, p6_array(cmls))
        for ci, c in enumerate(cmls):
            for s in (0, 1):
                k, al, r = dl[(c, s)]
                sv = st[ci * 2 + s]
                code = decide(sv, (c, s) in p5, r)
                row = dict(cml_id=c, sublink=s, frequency_MHz=freq[s, pos[c]], length_km=length[pos[c]],
                           itu2005_k=k, itu2005_alpha=al, detection_limit_mmhr=r, p5_nonoperational=int((c, s) in p5))
                row.update({nm: sv[i] for i, nm in enumerate(STAT_NAMES)})
                row.update(status_code=code, status=STATUS[code])
                rows.append(row)
        print(f'  pass 1: {min(b0 + a.batch, len(order))}/{len(order)} CMLs  {time.time() - t_start:.0f}s', flush=True)
    res = pd.DataFrame(rows)
    res['cml_retained'] = res.groupby('cml_id').status_code.transform(lambda x: int((x == 0).any()))
    res.to_csv(out_sub, index=False)
    kept = [c for c in order if res[(res.cml_id == c)].cml_retained.iloc[0] == 1]
    print(f'pass 1 done: retained CMLs {len(kept)} of {len(order)}; retained sublinks {(res.status_code == 0).sum()}', flush=True)

    # ---------------- pass 2: write retained CMLs ----------------
    old = nc.Dataset(OLD_QC_DEFAULT)
    nk = len(kept)
    out = nc.Dataset(out_nc, 'w', format='NETCDF4')
    out.createDimension('time', IB - IA); out.createDimension('sublink_id', 2); out.createDimension('cml_id', nk)
    tv = out.createVariable('time', 'i8', ('time',))
    tv.setncatts({k: old['time'].getncattr(k) for k in old['time'].ncattrs()})
    tv[:] = np.arange(IA, IB, dtype=np.int64)
    sv_ = out.createVariable('sublink_id', 'i8', ('sublink_id',))
    sv_.setncatts({k: old['sublink_id'].getncattr(k) for k in old['sublink_id'].ncattrs()}); sv_[:] = [0, 1]
    cv = out.createVariable('cml_id', str, ('cml_id',))
    cv.setncatts({k: old['cml_id'].getncattr(k) for k in old['cml_id'].ncattrs()})
    cv[:] = np.array(kept, dtype=object)
    kidx = np.array([pos[c] for c in kept])
    for name in ('length', 'site_a_latitude', 'site_b_latitude', 'site_a_longitude', 'site_b_longitude'):
        v = out.createVariable(name, 'f4', ('cml_id',), fill_value=np.float32(np.nan))
        v.setncatts({k: old[name].getncattr(k) for k in old[name].ncattrs() if k != '_FillValue'})
        v[:] = np.ma.filled(ds[name][:].astype(np.float32), np.nan)[kidx]
    fv = out.createVariable('frequency', 'f4', ('sublink_id', 'cml_id'), fill_value=np.float32(np.nan))
    fv.setncatts({k: old['frequency'].getncattr(k) for k in old['frequency'].ncattrs() if k != '_FillValue'})
    fv[:] = freq[:, kidx].astype(np.float32)
    pv = out.createVariable('polarization', str, ('sublink_id', 'cml_id'))
    pv.long_name = 'Microwave polarization'
    pv.comment = 'Vertical for every CML (source metadata, confirmed by the author).'
    pv[:] = np.full((2, nk), 'Vertical', dtype=object)
    coords = 'length frequency polarization site_a_latitude site_a_longitude site_b_latitude site_b_longitude qc_detection_limit_mmhr'
    ch = (10080, 2, min(64, nk))
    def var3(name, dtype, fill, attrs):
        v = out.createVariable(name, dtype, ('time', 'sublink_id', 'cml_id'), zlib=True, complevel=4, shuffle=True,
                               chunksizes=ch, fill_value=fill)
        v.setncatts(attrs)
        return v
    vtx = var3('tx', 'f4', np.float32(np.nan), dict(long_name='Transmitted signal level after quality control', units='dBm',
               coordinates=coords, coverage_content_type='physicalMeasurement',
               comment='Observed TX from the raw product. Missing where QC masked the minute (invalid value, mixed-placeholder minute, low-RSL plateau) or where the sublink is excluded. Never filled.'))
    vrx = var3('rx', 'f4', np.float32(np.nan), dict(long_name='Received signal level after quality control', units='dBm',
               coordinates=coords, coverage_content_type='physicalMeasurement',
               comment='Observed RX from the raw product. Missing where QC masked the minute or where the sublink is excluded. Never filled; blackout minutes stay missing.'))
    vtl = var3('trsl', 'f4', np.float32(np.nan), dict(long_name='Total signal loss (TX - RX) after quality control', units='dB',
               coordinates=coords, coverage_content_type='physicalMeasurement',
               comment='TX - RX where both are observed. In addition, flagged minutes: qc_flag bit 8 = blackout total loss reconstructed (adaptation of Blettner et al., 2023): trsl = max(tx) - min(rx) of the same sublink over its observations in the QC period 2021-04-01 to 2021-10-31 after masking; RX is missing at these minutes, TX may be observed; bit 16 = linear interpolation over a raw gap of at most 5 minutes (TX and/or RX missing). At flagged minutes trsl is not equal to tx - rx. Contains the baseline link loss and weather-related attenuation.'))
    vfl = var3('qc_flag', 'u1', False, dict(long_name='Per-minute quality-control flag', units='1',
               flag_masks=np.array([1, 2, 4, 8, 16], dtype=np.uint8),
               flag_meanings='invalid_value_masked mixed_placeholder_minute_masked low_rsl_plateau_masked blackout_total_loss_reconstructed short_gap_interpolated',
               comment='Bit field per minute and sublink. Bit 4 marks minutes that held a TX or RX value removed by the low-RSL plateau filter. Bit 8: blackout total loss reconstructed in trsl only (RX missing, TX may be observed); blackout gaps containing a bit-2 or bit-4 minute are not reconstructed and stay missing. Values of excluded sublinks are missing; their flags are kept.'))
    def var2(name, dtype, fill, attrs):
        v = out.createVariable(name, dtype, ('sublink_id', 'cml_id'), fill_value=fill)
        v.setncatts(attrs)
        return v
    vdl = var2('qc_detection_limit_mmhr', 'f4', np.float32(np.nan), dict(long_name='Sublink minimum detectable rain rate', units='mm h-1',
               comment='R_min = (0.33 dB / (k L))^(1/alpha), ITU-R P.838-3 (2005) k and alpha (pycomlink a_b, ITU_2005), vertical polarization, sublink frequency, path length. Sublinks with R_min >= 2 mm h-1 are excluded.'))
    vst = var2('qc_sublink_status', 'i1', np.int8(-1), dict(long_name='Sublink QC decision', units='1',
               flag_values=np.array(list(STATUS), dtype=np.int8), flag_meanings=' '.join(STATUS.values())))
    vfa = var2('qc_fraction_5h_windows_std_gt_2dB', 'f4', np.float32(np.nan), dict(long_name='Fraction of valid 5-hour windows with total-loss standard deviation above 2 dB', units='1',
               comment='Valid window: 300 valid minutes; reconstructed blackout minutes excluded. Fail if > 0.10.'))
    vfb = var2('qc_fraction_1h_windows_std_gt_0p8dB', 'f4', np.float32(np.nan), dict(long_name='Fraction of valid 1-hour windows with total-loss standard deviation above 0.8 dB', units='1',
               comment='Valid window: 60 valid minutes; reconstructed blackout minutes excluded. Fail if > 0.33.'))
    for b0 in range(0, nk, 64):
        cmls = np.array(kept[b0:b0 + 64])
        cols = np.array([pos[c] for c in cmls])
        tx, rx = read_block(ds, cols)
        otx, orx, otl, ofl, st = qc_batch(tx, rx, p6_array(cmls))
        nb = len(cmls)
        A = np.full((3, nb * 2, IB - IA), np.nan, dtype=np.float32)
        for ci, c in enumerate(cmls):
            for s in (0, 1):
                si = ci * 2 + s
                k, al, r = dl[(c, s)]
                code = decide(st[si], (c, s) in p5, r)
                ref = res[(res.cml_id == c) & (res.sublink == s)].status_code.iloc[0]
                assert code == ref, (c, s, code, ref)
                if code == 0:
                    A[0, si] = otx[si]; A[1, si] = orx[si]; A[2, si] = otl[si]
                vdl[s, b0 + ci] = r
                vst[s, b0 + ci] = code
                vfa[s, b0 + ci] = st[si][12]; vfb[s, b0 + ci] = st[si][15]
        for v, arr in ((vtx, A[0]), (vrx, A[1]), (vtl, A[2])):
            v[:, :, b0:b0 + nb] = arr.reshape(nb, 2, -1).transpose(2, 1, 0)
        vfl[:, :, b0:b0 + nb] = ofl.reshape(nb, 2, -1).transpose(2, 1, 0)
        print(f'  pass 2: {b0 + nb}/{nk} CMLs written  {time.time() - t_start:.0f}s', flush=True)

    now = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    out.setncatts(dict(
        title='Quality-controlled commercial microwave link signal levels, Czech Republic, 1 April - 31 October 2021 (corrected QC candidate)',
        Conventions='CF-1.13, ACDD-1.3', time_coverage_start='2021-04-01T00:00:00Z', time_coverage_end='2021-10-31T23:59:00Z',
        time_coverage_resolution='PT1M', processing_level='Quality-controlled CML signal data (candidate, not released)',
        source_raw_file=os.path.basename(a.raw_nc), source_raw_sha256=sha256(a.raw_nc),
        qc_script=os.path.basename(__file__), qc_script_sha256=sha256(__file__),
        qc_parameters_json=json.dumps(P), qc_period='2021-04-01T00:00 to 2021-10-31T23:59 UTC',
        population='Approved RAW 2574 CMLs; 8 without frequency excluded before QC; 2566 assessed; ' + f'{nk} retained',
        n_cml_assessed=len(order), n_cml_retained=nk,
        history=f'{now} created by {os.path.basename(__file__)} (Stage 3 revision 1, author decision 2026-10-28: blackout gaps containing plateau-removed or P6-masked minutes are not reconstructed) from the full-history RAW; supersedes candidate QC_CORRECTED_CANDIDATE.nc (sha256 e193b2cf2e7dac3fc12c1e08a12c79a793526679f5e91052df93c735786f7f8d), which is preserved; no existing product modified.',
        qc_blackout_method='Adaptation of Blettner et al. (2023): an RX gap (whole run of missing RX after invalid-value, P6 and plateau masking) of at most 60 minutes with RX < -65 dBm at least at one edge is a blackout. TX and RX are not filled; trsl = max(tx) - min(rx), where the maximum TX and minimum RX are taken over the same sublink in the QC reference period 2021-04-01T00:00 to 2021-10-31T23:59 UTC after masking. Gaps containing a plateau-removed or P6-masked minute are not reconstructed. Reconstructed minutes (qc_flag bit 8) are excluded from the fluctuation statistics and never used as interpolation end points.',
    ))
    out.close(); ds.close(); old.close()
    print(f'DONE {out_nc}  ({os.path.getsize(out_nc) / 1e9:.2f} GB)  total {time.time() - t_start:.0f}s')


if __name__ == '__main__':
    main()
