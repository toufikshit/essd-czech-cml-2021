#!/usr/bin/env python3
"""
V6 quantization x path length x E-band frequency analysis.

Required input:
  CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv

Optional input:
  pair-level V6-vs-CHMI validation CSV. If omitted, the script searches
  WVD_V6_VALIDATION recursively for a CSV containing cml_id, sublink_id,
  reference/CHMI WVD and V6 WVD columns.

The V6 retrieval itself is not changed. Quantization thresholds are a
pre-specified diagnostic sensitivity analysis; they are not optimized against
CHMI humidity.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BANDS=[("71-76 GHz",71.,76.),("81-86 GHz",81.,86.)]
LENGTH_THRESHOLDS=[2,3,4,5,6]
LENGTH_BINS=[("2-3",2,3),("3-4",3,4),("4-5",4,5),("5-6",5,6),("6-10",6,10.000001)]
QTHR=[0.20,0.25,0.30,0.40,0.50,0.75,1.00]

REF_NAMES=["chmi_wvd","wvd_chmi","reference_wvd","ref_wvd","truth_wvd","wvd_reference","rho_chmi","rho_ref"]
V6_NAMES=["v6_wvd","wvd_v6","wvd_internal_anchor_v6","cml_wvd","estimated_wvd","wvd_est"]
CID_NAMES=["cml_id","cml","link_id"]
SID_NAMES=["sublink_id","sublink","direction","sub_id"]
BRANCH_NAMES=["branch","wetdry_branch","method","retrieval_branch"]


def args_parse():
    root=Path("/home/jovyan/ESSD_publication_dataset")
    p=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--quant-csv",type=Path,default=root/"CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv")
    p.add_argument("--validation-dir",type=Path,default=root/"WVD_V6_VALIDATION")
    p.add_argument("--pairs-csv",type=Path,default=None)
    p.add_argument("--out-dir",type=Path,default=root/"WVD_V6_QUANTIZATION_ANALYSIS")
    p.add_argument("--branch",choices=["GRAF","GAUGE","ALL"],default="GRAF")
    p.add_argument("--max-wvd",type=float,default=35.)
    p.add_argument("--gridsize",type=int,default=55)
    p.add_argument("--dpi",type=int,default=250)
    return p.parse_args()


def norm_id(x):
    s=str(x).strip()
    if s.endswith(".0"):
        try:return str(int(float(s)))
        except Exception:pass
    return s


def first_present(cols,names):
    m={str(c).lower():c for c in cols}
    for n in names:
        if n.lower() in m:return m[n.lower()]
    return None


def prep_quant(path):
    d=pd.read_csv(path)
    req=["cml_id","sublink_id","frequency_ghz","length_km","tl_quantization_db",
         "specific_attenuation_quantization_db_per_km",
         "wvd_quantization_resolution_ref20c_rho10_gm3",
         "wvd_quantization_resolution_standard_min_gm3",
         "wvd_quantization_resolution_standard_median_gm3",
         "wvd_quantization_resolution_standard_max_gm3"]
    miss=[c for c in req if c not in d]
    if miss:raise KeyError(f"Missing quantization columns: {miss}")
    d["cml_id"]=d.cml_id.map(norm_id); d["sublink_id"]=d.sublink_id.map(norm_id)
    for c in req[2:]:d[c]=pd.to_numeric(d[c],errors="coerce")
    f=d.frequency_ghz
    d["band"]="outside"
    d.loc[(f>=71)&(f<=76),"band"]="71-76 GHz"
    d.loc[(f>=81)&(f<=86),"band"]="81-86 GHz"
    d=d[d.band!="outside"].copy()
    if d.duplicated(["cml_id","sublink_id"]).any():
        raise ValueError("Duplicate directional links in quantization CSV")
    return d


def write_quant_tables(q,out):
    rows=[]
    for band in ["all E-band","71-76 GHz","81-86 GHz"]:
        b=q if band=="all E-band" else q[q.band==band]
        rows.append(dict(population=band,n_links=len(b),median_length_km=b.length_km.median(),
            median_qtl_db=b.tl_quantization_db.median(),
            median_qgamma_db_per_km=b.specific_attenuation_quantization_db_per_km.median(),
            median_delta_rho_ref_gm3=b.wvd_quantization_resolution_ref20c_rho10_gm3.median(),
            median_delta_rho_standard_gm3=b.wvd_quantization_resolution_standard_median_gm3.median(),
            median_standard_max_delta_rho_gm3=b.wvd_quantization_resolution_standard_max_gm3.median()))
    pd.DataFrame(rows).to_csv(out/"quantization_link_summary.csv",index=False)

    rows=[]
    for band in ["all E-band","71-76 GHz","81-86 GHz"]:
        b=q if band=="all E-band" else q[q.band==band]
        for L in LENGTH_THRESHOLDS:
            base=b[b.length_km>=L]
            valid=base[base.wvd_quantization_resolution_ref20c_rho10_gm3.notna()]
            for thr in QTHR:
                keep=valid[valid.wvd_quantization_resolution_ref20c_rho10_gm3<=thr]
                rows.append(dict(band=band,min_length_km=L,delta_rho_threshold_gm3=thr,
                    n_links_population=len(base),n_links_with_quantization=len(valid),n_links_retained=len(keep),
                    retention_population_pct=100*len(keep)/len(base) if len(base) else np.nan))
    pd.DataFrame(rows).to_csv(out/"quantization_threshold_retention.csv",index=False)

    rows=[]
    for band,_,_ in BANDS:
        b=q[q.band==band]
        for label,lo,hi in LENGTH_BINS:
            x=b[(b.length_km>=lo)&(b.length_km<hi)]
            rows.append(dict(band=band,length_bin=label,n_links=len(x),
                median_qtl_db=x.tl_quantization_db.median(),
                median_qgamma_db_per_km=x.specific_attenuation_quantization_db_per_km.median(),
                median_delta_rho_ref_gm3=x.wvd_quantization_resolution_ref20c_rho10_gm3.median(),
                q25_delta_rho_ref_gm3=x.wvd_quantization_resolution_ref20c_rho10_gm3.quantile(.25),
                q75_delta_rho_ref_gm3=x.wvd_quantization_resolution_ref20c_rho10_gm3.quantile(.75),
                median_standard_max_delta_rho_gm3=x.wvd_quantization_resolution_standard_max_gm3.median()))
    pd.DataFrame(rows).to_csv(out/"quantization_length_band_summary.csv",index=False)


def plot_quant(q,figdir,dpi):
    fig,axs=plt.subplots(1,2,figsize=(11,4.3),sharey=True)
    for ax,(band,_,_) in zip(axs,BANDS):
        x=q[q.band==band]
        ax.scatter(x.length_km,x.wvd_quantization_resolution_ref20c_rho10_gm3,s=13,alpha=.5)
        xm=[];ym=[]
        for _,lo,hi in LENGTH_BINS:
            z=x[(x.length_km>=lo)&(x.length_km<hi)]
            if len(z):xm.append(z.length_km.median());ym.append(z.wvd_quantization_resolution_ref20c_rho10_gm3.median())
        if xm:ax.plot(xm,ym,marker="o",linewidth=1.3)
        for y in (.25,.5,1.):ax.axhline(y,linestyle="--",linewidth=.7,alpha=.5)
        ax.set_title(f"{band} (n={len(x)})");ax.set_xlabel("Path length [km]");ax.grid(alpha=.2)
    axs[0].set_ylabel(r"Quantization-limited WVD resolution $\Delta\rho_q$ [g m$^{-3}$]\n(reference: 20 °C, 10 g m$^{-3}$)")
    fig.suptitle("V6 link-derived quantization resolution versus path length",y=1.02)
    fig.tight_layout();fig.savefig(figdir/"01_quantization_resolution_vs_length_by_band.png",dpi=dpi,bbox_inches="tight");plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.4,4.7))
    size=18+10*np.clip(q.length_km.to_numpy(float),0,10)
    ax.scatter(q.frequency_ghz,q.wvd_quantization_resolution_ref20c_rho10_gm3,s=size,alpha=.35)
    ax.set_xlabel("Frequency [GHz]");ax.set_ylabel(r"$\Delta\rho_q$ [g m$^{-3}$] at 20 °C, WVD=10 g m$^{-3}$")
    ax.set_title("V6 link-derived quantization resolution versus frequency");ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(figdir/"02_quantization_resolution_vs_frequency.png",dpi=dpi,bbox_inches="tight");plt.close(fig)

    ret=pd.read_csv(figdir.parent/"quantization_threshold_retention.csv")
    fig,axs=plt.subplots(2,1,figsize=(8.5,7),sharex=True,sharey=True)
    for ax,band in zip(axs,["71-76 GHz","81-86 GHz"]):
        b=ret[ret.band==band]
        for L in LENGTH_THRESHOLDS:
            x=b[b.min_length_km==L]
            ax.plot(x.delta_rho_threshold_gm3,x.retention_population_pct,marker="o",label=f"L ≥ {L} km")
        ax.set_title(band);ax.set_ylabel("Directional links retained [%]");ax.set_ylim(0,103);ax.grid(alpha=.2);ax.legend(ncol=3,fontsize=8)
    axs[-1].set_xlabel(r"Maximum $\Delta\rho_q$ [g m$^{-3}$]")
    fig.suptitle("Quantization-threshold sensitivity (diagnostic, not CHMI-optimized)")
    fig.tight_layout();fig.savefig(figdir/"03_quantization_threshold_retention_by_band.png",dpi=dpi,bbox_inches="tight");plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(11,4.3),sharey=True)
    for ax,(band,_,_) in zip(axs,BANDS):
        x=q[q.band==band].sort_values("length_km")
        ax.scatter(x.length_km,x.wvd_quantization_resolution_standard_median_gm3,s=12,alpha=.5)
        n=len(x); step=max(1,n//150)
        xx=x.iloc[::step]
        for _,r in xx.iterrows():
            if np.isfinite(r.wvd_quantization_resolution_standard_min_gm3) and np.isfinite(r.wvd_quantization_resolution_standard_max_gm3):
                ax.plot([r.length_km,r.length_km],[r.wvd_quantization_resolution_standard_min_gm3,r.wvd_quantization_resolution_standard_max_gm3],linewidth=.5,alpha=.25)
        ax.set_title(f"{band} (n={len(x)})");ax.set_xlabel("Path length [km]");ax.grid(alpha=.2)
    axs[0].set_ylabel(r"$\Delta\rho_q$ [g m$^{-3}$], T=10–30 °C, WVD=2–25 g m$^{-3}$")
    fig.suptitle("P.676 quantization-limited WVD resolution over standardized physical domain",y=1.02)
    fig.tight_layout();fig.savefig(figdir/"04_quantization_standard_range_vs_length_by_band.png",dpi=dpi,bbox_inches="tight");plt.close(fig)


def inspect_pairs(path):
    try:h=pd.read_csv(path,nrows=3)
    except Exception:return None
    ref=first_present(h.columns,REF_NAMES);v6=first_present(h.columns,V6_NAMES);cid=first_present(h.columns,CID_NAMES);sid=first_present(h.columns,SID_NAMES);branch=first_present(h.columns,BRANCH_NAMES)
    if all(x is not None for x in [ref,v6,cid,sid]):return dict(path=path,ref=ref,v6=v6,cid=cid,sid=sid,branch=branch)
    return None


def find_pairs(a):
    if a.pairs_csv is not None:
        z=inspect_pairs(a.pairs_csv)
        if z is None:raise ValueError(f"Not recognized as pair-level validation CSV: {a.pairs_csv}")
        return z
    if not a.validation_dir.exists():return None
    z=[]
    for p in a.validation_dir.rglob("*.csv"):
        info=inspect_pairs(p)
        if info is not None:z.append((p.stat().st_size,info))
    return max(z,key=lambda x:x[0])[1] if z else None


def load_pairs(info,q,branch):
    d=pd.read_csv(info["path"])
    ren={info["cid"]:"cml_id",info["sid"]:"sublink_id",info["ref"]:"chmi_wvd",info["v6"]:"v6_wvd"}
    if info["branch"] is not None:ren[info["branch"]]="branch"
    d=d.rename(columns=ren)
    d["cml_id"]=d.cml_id.map(norm_id);d["sublink_id"]=d.sublink_id.map(norm_id)
    d["chmi_wvd"]=pd.to_numeric(d.chmi_wvd,errors="coerce");d["v6_wvd"]=pd.to_numeric(d.v6_wvd,errors="coerce")
    if "branch" in d and branch!="ALL":d=d[d.branch.astype(str).str.upper().str.contains(branch,regex=False)]
    cols=["cml_id","sublink_id","frequency_ghz","length_km","band","tl_quantization_db","wvd_quantization_resolution_ref20c_rho10_gm3"]
    d=d.merge(q[cols],on=["cml_id","sublink_id"],how="inner",validate="many_to_one")
    return d[np.isfinite(d.chmi_wvd)&np.isfinite(d.v6_wvd)].copy()


def metr(x):
    r=x.chmi_wvd.to_numpy(float);e=x.v6_wvd.to_numpy(float);n=len(r)
    if not n:return dict(N=0,PCC=np.nan,RMSE=np.nan,Bias=np.nan,MAE=np.nan)
    er=e-r;p=np.nan
    if n>1 and np.std(r)>0 and np.std(e)>0:p=float(np.corrcoef(r,e)[0,1])
    return dict(N=n,PCC=p,RMSE=float(np.sqrt(np.mean(er*er))),Bias=float(np.mean(er)),MAE=float(np.mean(np.abs(er))))


def validation_tables(pairs,out):
    rows=[]
    for band,_,_ in BANDS:
        b=pairs[pairs.band==band]
        for L in LENGTH_THRESHOLDS:
            x=b[b.length_km>=L]
            rows.append(dict(band=band,min_length_km=L,n_links=x[["cml_id","sublink_id"]].drop_duplicates().shape[0],**metr(x)))
    pd.DataFrame(rows).to_csv(out/"validation_metrics_by_length_band.csv",index=False)
    rows=[]
    for band,_,_ in BANDS:
        b=pairs[pairs.band==band]
        for t in QTHR:
            x=b[b.wvd_quantization_resolution_ref20c_rho10_gm3<=t]
            rows.append(dict(band=band,max_delta_rho_q_gm3=t,n_links=x[["cml_id","sublink_id"]].drop_duplicates().shape[0],**metr(x)))
    pd.DataFrame(rows).to_csv(out/"validation_metrics_by_quantization_band.csv",index=False)


def panel(ax,x,a,title):
    m=metr(x)
    if len(x):ax.hexbin(x.chmi_wvd,x.v6_wvd,gridsize=a.gridsize,extent=(0,a.max_wvd,0,a.max_wvd),mincnt=1,bins="log")
    ax.plot([0,a.max_wvd],[0,a.max_wvd],linestyle="--",linewidth=.8)
    ax.set_xlim(0,a.max_wvd);ax.set_ylim(0,a.max_wvd);ax.set_aspect("equal",adjustable="box");ax.set_title(title,fontsize=9)
    ax.text(.04,.96,f"N={m['N']:,}\nPCC={m['PCC']:.3f}\nRMSE={m['RMSE']:.2f}\nBias={m['Bias']:+.2f}",transform=ax.transAxes,ha="left",va="top",fontsize=7.2,bbox=dict(boxstyle="round,pad=.25",alpha=.85))


def validation_plots(pairs,figdir,a):
    for band,_,_ in BANDS:
        b=pairs[pairs.band==band]
        fig,axs=plt.subplots(1,5,figsize=(15.5,3.2),sharex=True,sharey=True)
        for ax,L in zip(axs,LENGTH_THRESHOLDS):panel(ax,b[b.length_km>=L],a,f"L ≥ {L} km")
        axs[0].set_ylabel(r"V6 WVD [g m$^{-3}$]")
        for ax in axs:ax.set_xlabel(r"CHMI WVD [g m$^{-3}$]")
        fig.suptitle(f"V6 validation by minimum path length — {band}",y=1.02);fig.tight_layout()
        fig.savefig(figdir/f"10_validation_length_{band.replace(' ','_').replace('-','_')}.png",dpi=a.dpi,bbox_inches="tight");plt.close(fig)

        th=[.20,.25,.30,.40,.50,.75]
        fig,axs=plt.subplots(2,3,figsize=(9.5,6.3),sharex=True,sharey=True);axs=axs.ravel()
        for ax,t in zip(axs,th):panel(ax,b[b.wvd_quantization_resolution_ref20c_rho10_gm3<=t],a,rf"$\Delta\rho_q \leq$ {t:.2f}")
        axs[0].set_ylabel(r"V6 WVD [g m$^{-3}$]");axs[3].set_ylabel(r"V6 WVD [g m$^{-3}$]")
        for ax in axs[3:]:ax.set_xlabel(r"CHMI WVD [g m$^{-3}$]")
        fig.suptitle(f"V6 validation by quantization-limited resolution — {band}",y=1.01);fig.tight_layout()
        fig.savefig(figdir/f"11_validation_quantization_{band.replace(' ','_').replace('-','_')}.png",dpi=a.dpi,bbox_inches="tight");plt.close(fig)


def main():
    a=args_parse();out=a.out_dir.expanduser().resolve();figdir=out/"figures";figdir.mkdir(parents=True,exist_ok=True)
    q=prep_quant(a.quant_csv.expanduser().resolve())
    print("="*100);print("V6 QUANTIZATION × LENGTH × FREQUENCY ANALYSIS");print("="*100)
    print("links:",len(q))
    for band,_,_ in BANDS:
        x=q[q.band==band]
        print(f"{band}: n={len(x)} median L={x.length_km.median():.2f} km median qTL={x.tl_quantization_db.median():.3f} dB median dRho={x.wvd_quantization_resolution_ref20c_rho10_gm3.median():.3f} g/m3")
    write_quant_tables(q,out);plot_quant(q,figdir,a.dpi)
    info=find_pairs(a)
    if info is None:
        print("\nNo pair-level V6-vs-CHMI CSV found. Quantization analysis is complete.")
        print("For validation hexbins rerun with: --pairs-csv /path/to/pairs.csv")
        print("Output:",out);return
    print("\nPair-level validation CSV:",info["path"])
    pairs=load_pairs(info,q,a.branch)
    print("paired observations:",f"{len(pairs):,}"," directional links:",pairs[["cml_id","sublink_id"]].drop_duplicates().shape[0])
    validation_tables(pairs,out);validation_plots(pairs,figdir,a)
    print("Output:",out)

if __name__=="__main__":main()
