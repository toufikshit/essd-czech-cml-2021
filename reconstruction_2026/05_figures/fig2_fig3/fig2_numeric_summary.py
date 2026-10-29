"""fig2_numeric_summary.py - Stage 5B, Figure 2: numbers behind the old and the new CML row.

Reads the per-timestep CSVs written by the unmodified producer for the original inputs (fig2_original) and for the
approved QC REV1 CML file (fig2_new) and prints, per sensor, n_total, overall min/median/mean/max and monthly means
(Apr-Oct 2021), plus the CML row difference new - old. Writes fig2_monthly_old_vs_new.csv.

Usage: python -I -B fig2_numeric_summary.py [--old DIR] [--new DIR] [--out CSV]
"""
import argparse, os
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
p = argparse.ArgumentParser()
p.add_argument("--old", default=os.path.join(HERE, "..", "reproduction", "fig2_original"))
p.add_argument("--new", default=os.path.join(HERE, "..", "reproduction", "fig2_new"))
p.add_argument("--out", default=os.path.join(HERE, "..", "reproduction", "fig2_monthly_old_vs_new.csv"))
a = p.parse_args()
rows = []
for s in ("cml", "radar", "gauge", "temperature", "humidity"):
    for tag, d in (("old", a.old), ("new", a.new)):
        x = pd.read_csv(os.path.join(d, f"availability_native_timestep_percentage_{s}.csv"), parse_dates=["time"])
        v = x.availability_percent
        print(f"{s:12s} {tag}: n_total={x.n_total.iloc[0]}, n_time={len(x)}, min {v.min():.2f} median {v.median():.2f} mean {v.mean():.2f} max {v.max():.2f}; "
              f"mean n_valid {x.n_valid.mean():.1f}")
        m = x.groupby(x.time.dt.to_period("M")).availability_percent.mean()
        for k, val in m.items():
            rows.append({"sensor": s, "version": tag, "month": str(k), "mean_availability_percent": round(val, 3), "n_total": int(x.n_total.iloc[0])})
df = pd.DataFrame(rows)
df.to_csv(a.out, index=False)
piv = df[df.sensor == "cml"].pivot(index="month", columns="version", values="mean_availability_percent")
piv["new_minus_old"] = piv["new"] - piv["old"]
print("\nCML monthly mean availability (%):\n" + piv.round(2).to_string())
o = pd.read_csv(os.path.join(a.old, "availability_native_timestep_percentage_cml.csv")); n = pd.read_csv(os.path.join(a.new, "availability_native_timestep_percentage_cml.csv"))
assert (o.time == n.time).all()
d = n.availability_percent - o.availability_percent
print(f"\nCML per-minute difference new-old: mean {d.mean():+.2f} pp, median {d.median():+.2f}, p5 {d.quantile(.05):+.2f}, p95 {d.quantile(.95):+.2f}, min {d.min():+.2f}, max {d.max():+.2f}")
for thr in (0.1, 10, 20):
    print(f"   minutes with old < {thr}%: {(o.availability_percent < thr).sum()}, new < {thr}%: {(n.availability_percent < thr).sum()}")
print("[saved]", os.path.abspath(a.out))
