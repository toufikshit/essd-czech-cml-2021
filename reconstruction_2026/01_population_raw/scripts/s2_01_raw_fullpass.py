"""s2_01_raw_fullpass.py - One chunked, read-only pass over the full-history RAW (2580 CMLs x 2 sublinks).

Per sublink (fresh computation, not reused from Stage 1):
  finite TX, finite RX, both finite (full history and release window 2021-04-01..2021-10-31), TX-only and
  RX-only minutes, first and last minute with both finite, monthly both-finite counts (24 months),
  TX min/max, RX min/max, TX==0 minutes, longest run of consecutive both-finite minutes,
  counts of TX above 30, 40, 60, 100 dBm, TX below -30 dBm, RX above -20 dBm, RX <= -99.9 dBm, TX <= -99, TX == 255,
  and a coarse per-sublink TX histogram (5 dB bins) to describe implausible levels.
Global histograms of TX and RX (0.5 dB bins). Results -> s2_01_raw_fullpass.npz.

Usage: python -I s2_01_raw_fullpass.py [--raw NC] [--chunk N] [--out TXT]
"""
import argparse, os, sys, time
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, netCDF4

HERE = os.path.dirname(os.path.abspath(__file__))
p = argparse.ArgumentParser()
p.add_argument("--raw", default=r"C:\Users\taouf\ESSD_RELEASE\extracted\full\ESSD_CZ_2021_HANDOFF\data\source_full_history\raw\CZ_CML_RAW_1MIN_20200101_20211231.nc")
p.add_argument("--chunk", type=int, default=10080)
p.add_argument("--out", default=os.path.join(HERE, "s2_01_raw_fullpass_output.txt"))
a = p.parse_args()
t0 = time.time()
ds = netCDF4.Dataset(a.raw, "r"); ds.set_auto_mask(False)
T = ds["time"]; times = netCDF4.num2date(T[:], T.units)
nt = len(times); ns = 2; nl = ds.dimensions["cml_id"].size
mon = np.array([(t.year - 2020) * 12 + t.month - 1 for t in times], np.int16)
w0 = next(i for i, t in enumerate(times) if (t.year, t.month, t.day) == (2021, 4, 1)); w1 = w0 + 308160
Z = lambda: np.zeros((ns, nl), np.int64)
S = {k: Z() for k in ("tx", "rx", "both", "tx_win", "rx_win", "both_win", "tx_only", "rx_only", "tx_zero",
                      "tx_gt30", "tx_gt40", "tx_gt60", "tx_gt100", "tx_lt_m30", "rx_gt_m20", "rx_le_m999", "tx_le_m99", "tx_255", "maxrun")}
S["first"] = np.full((ns, nl), -1, np.int64); S["last"] = np.full((ns, nl), -1, np.int64)
S["txmin"] = np.full((ns, nl), np.inf); S["txmax"] = np.full((ns, nl), -np.inf)
S["rxmin"] = np.full((ns, nl), np.inf); S["rxmax"] = np.full((ns, nl), -np.inf)
S["monthly"] = np.zeros((24, ns, nl), np.int64)
tx_edges = np.arange(-150, 305, 5.0); S["tx_hist5"] = np.zeros((ns, nl, len(tx_edges) - 1), np.int64)
gb = np.arange(-150, 300.5, 0.5); S["g_tx"] = np.zeros(len(gb) - 1, np.int64)
rb = np.arange(-130, 30.5, 0.5); S["g_rx"] = np.zeros(len(rb) - 1, np.int64)
run = np.zeros((ns, nl), np.int64)
for c0 in range(0, nt, a.chunk):
    c1 = min(nt, c0 + a.chunk)
    tx = np.asarray(ds["tx"][c0:c1], np.float64); rx = np.asarray(ds["rx"][c0:c1], np.float64)
    ft, fr = np.isfinite(tx), np.isfinite(rx); b = ft & fr
    S["tx"] += ft.sum(0); S["rx"] += fr.sum(0); S["both"] += b.sum(0)
    S["tx_only"] += (ft & ~fr).sum(0); S["rx_only"] += (fr & ~ft).sum(0)
    lo, hi = max(c0, w0), min(c1, w1)
    if hi > lo:
        S["tx_win"] += ft[lo - c0:hi - c0].sum(0); S["rx_win"] += fr[lo - c0:hi - c0].sum(0); S["both_win"] += b[lo - c0:hi - c0].sum(0)
    with np.errstate(invalid="ignore"):
        S["tx_zero"] += (tx == 0).sum(0)
        for k, m in (("tx_gt30", tx > 30), ("tx_gt40", tx > 40), ("tx_gt60", tx > 60), ("tx_gt100", tx > 100), ("tx_lt_m30", tx < -30),
                     ("rx_gt_m20", rx > -20), ("rx_le_m999", rx <= -99.9), ("tx_le_m99", tx <= -99), ("tx_255", tx == 255)):
            S[k] += m.sum(0)
    if b.any():
        anyb = b.any(0); fi = np.argmax(b, axis=0) + c0; la = c1 - 1 - np.argmax(b[::-1], axis=0)
        S["first"] = np.where((S["first"] < 0) & anyb, fi, S["first"]); S["last"] = np.where(anyb, la, S["last"])
        for mm in np.unique(mon[c0:c1]):
            S["monthly"][mm] += b[mon[c0:c1] == mm].sum(0)
    S["txmin"] = np.fmin(S["txmin"], np.nanmin(np.where(ft, tx, np.inf), axis=0)); S["txmax"] = np.fmax(S["txmax"], np.nanmax(np.where(ft, tx, -np.inf), axis=0))
    S["rxmin"] = np.fmin(S["rxmin"], np.nanmin(np.where(fr, rx, np.inf), axis=0)); S["rxmax"] = np.fmax(S["rxmax"], np.nanmax(np.where(fr, rx, -np.inf), axis=0))
    if ft.any():
        idx = np.digitize(np.clip(tx, tx_edges[0], tx_edges[-1] - 1e-6), tx_edges) - 1
        for s in range(ns):
            ii = idx[:, s, :]; ok = ft[:, s, :]
            cols = np.broadcast_to(np.arange(nl), ii.shape)
            np.add.at(S["tx_hist5"][s], (cols[ok], ii[ok]), 1)
        S["g_tx"] += np.histogram(tx[ft], gb)[0]
    if fr.any():
        S["g_rx"] += np.histogram(rx[fr], rb)[0]
    for i in range(c1 - c0):
        run = np.where(b[i], run + 1, 0); S["maxrun"] = np.maximum(S["maxrun"], run)
    print(f"  {c1}/{nt} [{time.time()-t0:.0f}s]", flush=True)
S["cml_id"] = np.array([str(x) for x in ds["cml_id"][:]])
for k in ("length", "frequency", "site_a_latitude", "site_a_longitude", "site_b_latitude", "site_b_longitude"):
    S["meta_" + k] = np.asarray(ds[k][:], float)
pv = ds["polarization"][:]; S["meta_pol"] = np.array(netCDF4.chartostring(pv) if pv.dtype.kind == "S" else pv, dtype=str)
S["time0"] = np.array(str(times[0])); S["w0"] = np.array(w0); S["tx_edges"] = tx_edges; S["g_tx_edges"] = gb; S["g_rx_edges"] = rb
ds.close()
np.savez_compressed(os.path.join(HERE, "s2_01_raw_fullpass.npz"), **S)
lines = [f"[F] {a.raw}: {nl} CMLs x {ns} sublinks, {nt} minutes ({times[0]}..{times[-1]}), window index {w0}:{w1}; done in {time.time()-t0:.0f}s",
         f"[F] sublinks with both finite >0: {(S['both']>0).sum()}; CMLs with any both finite: {(S['both'].sum(0)>0).sum()}; total both-finite minutes {S['both'].sum():,}",
         f"[F] TX<=-99 {S['tx_le_m99'].sum()}, TX==255 {S['tx_255'].sum()}, RX>-20 {S['rx_gt_m20'].sum()}, RX<=-99.9 {S['rx_le_m999'].sum()}",
         f"[F] TX>30 {S['tx_gt30'].sum()} ({(S['tx_gt30']>0).sum()} sublinks), >40 {S['tx_gt40'].sum()} ({(S['tx_gt40']>0).sum()}), >60 {S['tx_gt60'].sum()} ({(S['tx_gt60']>0).sum()}), >100 {S['tx_gt100'].sum()} ({(S['tx_gt100']>0).sum()}); TX<-30 {S['tx_lt_m30'].sum()} ({(S['tx_lt_m30']>0).sum()})"]
open(a.out, "w", encoding="utf-8").write("\n".join(lines) + "\n"); print("\n".join(lines))
