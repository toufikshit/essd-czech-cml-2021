"""s2_02_audit_classify.py - Build RAW_LINK_AUDIT_draft.csv (one row per CML x sublink, CML columns repeated),
classify every sublink and CML into one primary class, and draft RAW_PROPOSED_EXCLUSIONS_draft.csv.
Read-only on all inputs; writes only the two draft tables and this script's output text.

Inputs: s2_01_raw_fullpass.npz (fresh pass over the full-history RAW), operator metadata, Stage 1B completeness
table (information only), the 13-CML decision table (for the TX==0 consistency check), validation membership
(gauge benchmark audit, released radar ids, WVD series; information only, never used for classes).

Primary classes (sublink), in this priority order (first match wins):
  5 completely empty                      : no finite TX and no finite RX in 2020-2021
  3 missing essential source observations : TX entirely missing while RX exists, or RX entirely missing while TX exists,
                                            or TX and RX never finite in the same minute
  4 missing essential metadata            : no frequency for the sublink, or no/zero path length, or incomplete endpoints
  6 other documented issue                : TX = 0 dBm at every valid minute (author decision pending, cf. the 13),
                                            or TX values > 60 dBm present (see basis in findings),
                                            or RAW metadata differ from the operator metadata
  2 suitable for raw, likely unsuitable for rainfall retrieval:
                                            no both-finite minute in 2021-04-01..10-31, or no run of 60 consecutive
                                            both-finite minutes (Graf threshold needs complete 60-sample windows,
                                            src/essd_cz/export_1min.py graf(): thr_min=60)
  1 suitable for raw and QC assessment
CML class = the most usable sublink class (order 1 < 2 < 6 < 4 < 3 < 5).

Usage: python -I s2_02_audit_classify.py [--npz NPZ] [--pkg DIR] [--compl CSV] [--t13 CSV] [--gauge CSV]
                                         [--radar NC] [--wvd NC] [--audit-out CSV] [--prop-out CSV] [--out TXT]
"""
import argparse, csv, os, sys
sys.stdout.reconfigure(encoding="utf-8")
from collections import Counter
import numpy as np, pandas as pd, netCDF4

HERE = os.path.dirname(os.path.abspath(__file__)); S2 = os.path.dirname(HERE); ROOT = os.path.dirname(S2)
E = r"C:\Users\taouf\ESSD_RELEASE"
p = argparse.ArgumentParser()
p.add_argument("--npz", default=os.path.join(HERE, "s2_01_raw_fullpass.npz"))
p.add_argument("--pkg", default=os.path.join(ROOT, "stage1b", "metacentrum_extracted", "ESSD_STAGE1B"))
p.add_argument("--compl", default=os.path.join(ROOT, "stage1b", "scripts", "b1_08_completeness_per_sublink.csv"))
p.add_argument("--t13", default=os.path.join(ROOT, "stage1b", "STAGE1B_EXCLUDED_13_CMLS.csv"))
p.add_argument("--gauge", default=E + r"\essd-czech-cml-2021\docs\evidence\hpc_run_2026-10\02_gauge_validation\validation_workspace\figures\FINAL_GRAF_gauge_shared_population_audit.csv")
p.add_argument("--radar", default=E + r"\zenodo_staging\primary\CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc")
p.add_argument("--qc", default=E + r"\zenodo_staging\primary\CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc")
p.add_argument("--wvd", default=E + r"\zenodo_staging\supplementary_wvd\CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc")
p.add_argument("--audit-out", default=os.path.join(S2, "agents", "RAW_LINK_AUDIT_draft.csv"))
p.add_argument("--prop-out", default=os.path.join(S2, "agents", "RAW_PROPOSED_EXCLUSIONS_draft.csv"))
p.add_argument("--out", default=os.path.join(HERE, "s2_02_audit_classify_output.txt"))
a = p.parse_args()
OUT = []
def say(s=""):
    print(s, flush=True); OUT.append(str(s))
z = np.load(a.npz, allow_pickle=True)
ids = [str(x) for x in z["cml_id"]]; nl = len(ids)
m = pd.read_csv(os.path.join(a.pkg, "metadata", "metadata_paired_2021.csv"), sep=";"); m["idcml"] = m.idcml.astype(str); mi = m.set_index("idcml")
grp = m[~m.isPairingUnique].groupby("id1").idcml.apply(lambda s: ";".join(sorted(s))).to_dict()
cp = pd.read_csv(a.compl); cp = cp[cp.source == "RAW"].set_index(["cml_id", "sublink"])
t13 = set(pd.read_csv(a.t13).cml_id.astype(str))
gauge = Counter()
for r in csv.DictReader(open(a.gauge, newline="")):
    gauge[(r["cml_id"], int(r["sublink_id"]))] = int(r["N"])
def ids_of(path, var="cml_id"):
    d = netCDF4.Dataset(path); x = [str(v) for v in d[var][:]]; d.close(); return x
radar = set(ids_of(a.radar)); qc = set(ids_of(a.qc))
w = netCDF4.Dataset(a.wvd); wv = set(zip([str(x) for x in w["cml_id"][:]], [int(str(x)) for x in w["sublink_id"][:]])); w.close()
t0 = pd.Timestamp(str(z["time0"]))
def tm(i):
    return str(t0 + pd.Timedelta(minutes=int(i))) if i >= 0 else ""
months = [f"{2020 + k // 12}-{k % 12 + 1:02d}" for k in range(24)]
U = {1: 0, 2: 1, 6: 2, 4: 3, 3: 4, 5: 5}
rows = []
for j, c in enumerate(ids):
    L = z["meta_length"][j]; co = [z[k][j] for k in ("meta_site_a_latitude", "meta_site_a_longitude", "meta_site_b_latitude", "meta_site_b_longitude")]
    mr = mi.loc[c] if c in mi.index else None
    for s in (0, 1):
        f = z["meta_frequency"][s, j]
        tx, rx, both = int(z["tx"][s, j]), int(z["rx"][s, j]), int(z["both"][s, j])
        txz = int(z["tx_zero"][s, j])
        meta_diff = []
        if mr is not None:
            if not (np.isclose(f, mr[f"freq{s}"], equal_nan=True)): meta_diff.append("frequency")
            if not (np.isclose(L, mr.length, atol=1e-3, equal_nan=True)): meta_diff.append("length")
            if not np.allclose(co, [mr.lat0, mr.lon0, mr.lat1, mr.lon1], atol=1e-5): meta_diff.append("coordinates")
        flags = []
        if tx == 0 and rx == 0: flags.append("5:empty")
        if (tx == 0 and rx > 0): flags.append("3:TX entirely missing, RX present")
        if (rx == 0 and tx > 0): flags.append("3:RX entirely missing, TX present")
        if tx > 0 and rx > 0 and both == 0: flags.append("3:TX and RX never finite together")
        if not np.isfinite(f): flags.append("4:no frequency")
        if not (np.isfinite(L) and L > 0): flags.append("4:no/zero length")
        if not all(np.isfinite(co)): flags.append("4:incomplete coordinates")
        if tx > 0 and txz == tx: flags.append("6:TX = 0 dBm at every valid minute (author decision, cf. the 13)")
        if z["tx_gt60"][s, j] > 0: flags.append(f"6:TX > 60 dBm in {int(z['tx_gt60'][s, j])} minutes")
        if meta_diff: flags.append("6:RAW metadata differ from operator (" + ",".join(meta_diff) + ")")
        if both > 0 and z["both_win"][s, j] == 0: flags.append("2:no both-finite minute in 2021-04-01..10-31")
        if both > 0 and z["maxrun"][s, j] < 60: flags.append("2:no run of 60 consecutive both-finite minutes")
        cls = min([int(x.split(":")[0]) for x in flags], key=lambda k: -U[k] if False else [5, 3, 4, 6, 2].index(k)) if flags else 1
        info = cp.loc[(c, s)] if (c, s) in cp.index else None
        rows.append({
            "cml_id": c, "sublink": s, "in_metadata": int(mr is not None), "in_released_qc": int(c in qc),
            "frequency_MHz": f, "length_km": L, "site_a_lat": co[0], "site_a_lon": co[1], "site_b_lat": co[2], "site_b_lon": co[3],
            "polarization": z["meta_pol"][s, j], "band_metadata": int(mr.band) if mr is not None else "",
            "isPairingUnique": int(mr.isPairingUnique) if mr is not None else "", "parallel_group": grp.get(mr.id1, "") if mr is not None and not mr.isPairingUnique else "",
            "tx_finite": tx, "rx_finite": rx, "both_finite": both, "tx_only": int(z["tx_only"][s, j]), "rx_only": int(z["rx_only"][s, j]),
            "tx_finite_window": int(z["tx_win"][s, j]), "rx_finite_window": int(z["rx_win"][s, j]), "both_finite_window": int(z["both_win"][s, j]),
            "first_valid": tm(z["first"][s, j]), "last_valid": tm(z["last"][s, j]), "longest_run_both_min": int(z["maxrun"][s, j]),
            "monthly_both_2021": ";".join(f"{months[k][5:]}:{int(z['monthly'][k, s, j])}" for k in range(12, 24) if z["monthly"][k, s, j] > 0),
            "tx_min": z["txmin"][s, j] if tx else np.nan, "tx_max": z["txmax"][s, j] if tx else np.nan,
            "rx_min": z["rxmin"][s, j] if rx else np.nan, "rx_max": z["rxmax"][s, j] if rx else np.nan,
            "tx_constant": int(tx > 0 and z["txmin"][s, j] == z["txmax"][s, j]), "tx_constant_value": z["txmin"][s, j] if (tx > 0 and z["txmin"][s, j] == z["txmax"][s, j]) else np.nan,
            "tx_zero_minutes": txz, "tx_all_zero": int(tx > 0 and txz == tx),
            "tx_gt40": int(z["tx_gt40"][s, j]), "tx_gt60": int(z["tx_gt60"][s, j]), "tx_gt100": int(z["tx_gt100"][s, j]), "tx_lt_m30": int(z["tx_lt_m30"][s, j]),
            "rx_gt_m20": int(z["rx_gt_m20"][s, j]), "rx_le_m99_9": int(z["rx_le_m999"][s, j]),
            "operator_valid_minutes_info": int(info.operator_valid) if info is not None else "", "operator_minus_raw_info": int(info.deficit) if info is not None else "",
            "gauge_benchmark_info": int((c, s) in gauge), "radar_info": int(c in radar), "wvd_info": int((c, s) in wv),
            "sublink_flags": " | ".join(flags), "sublink_class": cls})
df = pd.DataFrame(rows)
cmlc = df.groupby("cml_id").sublink_class.apply(lambda v: min(v, key=lambda k: U[k]))
df["cml_class"] = df.cml_id.map(cmlc)
df["cml_in_13"] = df.cml_id.isin(t13).astype(int)
os.makedirs(os.path.dirname(a.audit_out), exist_ok=True)
df.to_csv(a.audit_out, index=False)
say(f"[F] audit rows {len(df)} (CMLs {df.cml_id.nunique()}); any of the 13 in RAW: {int(df.cml_in_13.sum())}")
say(f"[F] sublink classes: {dict(sorted(Counter(df.sublink_class).items()))}")
say(f"[F] CML classes: {dict(sorted(Counter(cmlc).items()))}")
for k in (2, 3, 4, 5, 6):
    sub = df[df.sublink_class == k]
    say(f"  class {k}: {len(sub)} sublinks in {sub.cml_id.nunique()} CMLs; flags {Counter(f.split(' | ')[0] for f in sub.sublink_flags).most_common(6)}")
say(f"[F] all flags (multi-label) counts: {Counter(x for f in df.sublink_flags if f for x in f.split(' | ')).most_common(20)}")
# TX == 0 consistency list
z0 = df[df.tx_all_zero == 1]
say(f"[F] RAW sublinks with TX = 0 at every valid minute: {len(z0)} in {z0.cml_id.nunique()} CMLs; other sublink of these CMLs: " +
    "; ".join(f"{c} s{s} (TX min..max of other sublink {df[(df.cml_id==c)&(df.sublink!=s)].tx_min.iloc[0]}..{df[(df.cml_id==c)&(df.sublink!=s)].tx_max.iloc[0]}, both finite {df[(df.cml_id==c)&(df.sublink!=s)].both_finite.iloc[0]})" for c, s in zip(z0.cml_id, z0.sublink)))
# proposals
P = []
def cml_rows(cls):
    return sorted(cmlc[cmlc == cls].index)
n0c, n0s = nl, 2 * nl
empty_c = cml_rows(5)
for c in empty_c:
    P.append({"proposal": "P1", "scope": "whole CML", "cml_id": c, "sublink": "both", "reason": "completely empty: no finite TX and no finite RX in 2020-2021",
              "evidence": "[F] s2_01 full pass: tx_finite = rx_finite = 0 in both sublinks", "decision_type": "proposed exclusion (no observations lost)"})
c3 = cml_rows(3)
for c in c3:
    d = df[df.cml_id == c]
    P.append({"proposal": "P2", "scope": "whole CML", "cml_id": c, "sublink": "both", "reason": "no sublink has TX and RX in the same minute (" + "; ".join(f"s{r.sublink}: {r.sublink_flags}" for r in d.itertuples()) + ")",
              "evidence": f"[F] s2_01: tx_finite {list(d.tx_finite)}, rx_finite {list(d.rx_finite)}, both {list(d.both_finite)}", "decision_type": "author decision: exclude CML, or keep RX-only series"})
c4 = cml_rows(4)
for c in c4:
    d = df[df.cml_id == c]
    P.append({"proposal": "P3", "scope": "whole CML", "cml_id": c, "sublink": "both", "reason": "missing essential metadata: " + "; ".join(sorted(set(x for f in d.sublink_flags for x in f.split(' | ') if x.startswith('4:')))),
              "evidence": "[F] RAW metadata and operator metadata (freq0/freq1 NaN)", "decision_type": "author decision: exclude CML, or keep raw observations with a metadata flag"})
for r in df[(df.sublink_class.isin([5, 3])) & (~df.cml_id.isin(empty_c + c3))].itertuples():
    P.append({"proposal": "P4", "scope": "one sublink", "cml_id": r.cml_id, "sublink": r.sublink, "reason": r.sublink_flags,
              "evidence": f"[F] s2_01: tx_finite {r.tx_finite}, rx_finite {r.rx_finite}, both {r.both_finite}", "decision_type": "proposed sublink exclusion (sublink has no usable TX-RX pair)"})
for r in z0.itertuples():
    P.append({"proposal": "P5", "scope": "one sublink", "cml_id": r.cml_id, "sublink": r.sublink, "reason": "TX = 0 dBm at every valid minute (same criterion as the 13 excluded CMLs)",
              "evidence": f"[F] s2_01: tx_finite {r.tx_finite} = tx_zero {r.tx_zero_minutes}; other sublink has non-zero TX", "decision_type": "author decision (consistency with the 13): exclude sublink or keep"})
for r in df[df.tx_gt60 > 0].itertuples():
    P.append({"proposal": "P6", "scope": "individual observations", "cml_id": r.cml_id, "sublink": r.sublink,
              "reason": f"{r.tx_gt60} minutes with TX > 60 dBm (max {r.tx_max:.1f}); see distribution and basis in findings",
              "evidence": "[F] s2_01 TX histogram; [C] NB-CZ cell 80 CMLDataProcessor TX_MAX = 60 dBm", "decision_type": "author decision: set these TX values to missing, or keep"})
pr = pd.DataFrame(P)
# before/after counts per proposal (cumulative in order P1..P5; P6 observation-level)
cum_c, cum_s = n0c, n0s; ba = {}
for k in ("P1", "P2", "P3", "P4", "P5", "P6"):
    g = pr[pr.proposal == k]
    dc = int((g.scope == "whole CML").sum()); dsub = 2 * dc + int((g.scope == "one sublink").sum())
    ba[k] = (cum_c, cum_s, cum_c - dc, cum_s - dsub); cum_c -= dc; cum_s -= dsub
pr["cmls_before"] = pr.proposal.map(lambda k: ba[k][0]); pr["sublinks_before"] = pr.proposal.map(lambda k: ba[k][1])
pr["cmls_after"] = pr.proposal.map(lambda k: ba[k][2]); pr["sublinks_after"] = pr.proposal.map(lambda k: ba[k][3])
pr["obs_minutes_affected"] = [int(df[(df.cml_id == r.cml_id)].both_finite.sum()) if r.scope == "whole CML" else
                              (int(df[(df.cml_id == r.cml_id) & (df.sublink == r.sublink)].both_finite.iloc[0]) if r.scope == "one sublink" else
                               int(df[(df.cml_id == r.cml_id) & (df.sublink == r.sublink)].tx_gt60.iloc[0])) for r in pr.itertuples()]
pr.to_csv(a.prop_out, index=False)
for k, v in ba.items():
    g = pr[pr.proposal == k]
    say(f"[F] {k}: {len(g)} rows ({g.scope.iloc[0] if len(g) else '-'}); CMLs {v[0]}->{v[2]}, sublinks {v[1]}->{v[3]} (cumulative); both-finite minutes affected {int(g.obs_minutes_affected.sum()):,}")
open(a.out, "w", encoding="utf-8").write("\n".join(OUT) + "\n")
