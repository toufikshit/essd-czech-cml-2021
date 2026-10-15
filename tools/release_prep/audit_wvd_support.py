"""Audit gate: WVD V6 supporting products — series counts and index consistency.

Usage: python -I audit_wvd_support.py <data_root> <evidence_dir> <report.json>
  data_root:   .../handoff_data (contains wvd_support/)
  evidence_dir: .../repository/results/connected_workflow (for executed CSV artifacts)
"""
import json
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

BS = chr(92)


def wl(p):
    p = os.path.abspath(p)
    return (BS * 2 + "?" + BS) + p if len(p) > 240 else p


def main():
    data_root, evidence_dir, report_out = sys.argv[1], sys.argv[2], sys.argv[3]
    r = {}
    v6 = os.path.join(data_root, "wvd_support/wvd_v6_nc/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc")
    with xr.open_dataset(v6) as ds:
        r["v6_nc"] = {
            "n_series": int(ds.sizes.get("series", -1)),
            "dims": {k: int(v) for k, v in ds.sizes.items()},
            "data_vars": sorted(ds.data_vars),
            "time_min": str(np.asarray(ds["time"].values.min())) if "time" in ds else None,
            "time_max": str(np.asarray(ds["time"].values.max())) if "time" in ds else None,
        }
    diag = pd.read_csv(os.path.join(data_root, "wvd_support/wvd_v6_diagnostics_csv/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv"))
    quant = pd.read_csv(os.path.join(data_root, "wvd_support/wvd_quantization_csv/CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv"))
    r["diagnostics_csv"] = {"rows": len(diag),
                            "unique_series_index": int(diag["series_index"].nunique()) if "series_index" in diag else None,
                            "columns": list(diag.columns)[:20]}
    r["quantization_csv"] = {"rows": len(quant),
                             "unique_series_index": int(quant["series_index"].nunique()) if "series_index" in quant else None,
                             "columns": list(quant.columns)[:20]}
    n = r["v6_nc"]["n_series"]
    for key in ("diagnostics_csv", "quantization_csv"):
        df = diag if key == "diagnostics_csv" else quant
        if "series_index" in df:
            si = pd.to_numeric(df["series_index"], errors="coerce")
            r[key]["indices_within_v6_range"] = bool(si.notna().all() and (si >= 0).all() and (si < n).all())

    # executed HPC artifacts
    pop = pd.read_csv(os.path.join(evidence_dir, "WVD_V6_QUANTIZATION_BINNED_VALIDATION/validation_population_links.csv"))
    r["validation_population_links"] = {"rows": len(pop), "columns": list(pop.columns)[:15]}
    met = pd.read_csv(os.path.join(evidence_dir, "07_wvd_paper_population/wvp_support_metrics.csv"))
    r["paper_population_metrics"] = met.to_dict(orient="records")
    counts = pd.read_csv(os.path.join(evidence_dir, "WVD_V6_QUANTIZATION_BINNED_VALIDATION/population_counts.csv"))
    r["population_counts_total_directional_links"] = int(counts["n_directional_links"].sum())

    checks = {
        "v6_has_1352_series": n == 1352,
        "diagnostics_describe_1352_series": r["diagnostics_csv"]["unique_series_index"] == 1352,
        "quantization_describe_1352_series": r["quantization_csv"]["unique_series_index"] == 1352,
        "paper_all_Eband_45_sublinks_31_links": any(
            rec.get("n_sublinks") == 45 and rec.get("n_links") == 31 for rec in r["paper_population_metrics"]),
    }
    r["checks"] = checks
    with open(report_out, "w", encoding="utf-8") as f:
        json.dump(r, f, indent=2)
    print(json.dumps(checks, indent=2))
    print("v6 n_series:", n, "| diagnostics rows:", len(diag), "| quantization rows:", len(quant))


if __name__ == "__main__":
    main()
