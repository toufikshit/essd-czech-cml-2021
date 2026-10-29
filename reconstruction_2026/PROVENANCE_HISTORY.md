# Provenance history

This is the only page that describes how the released dataset differs from
an earlier, never published staging of the files. Everything else in the
repository describes the released dataset only. Use this page if you need to
understand old numbers in `docs/evidence/`, in `docs/CHANGES_FROM_HANDOFF.md`,
in the stage reports of this folder or in the executed notebook outputs.

## 1. The earlier staging (October 2026, not published)

Before the final processing, a set of files was staged for Zenodo but never
published. It differed from the released dataset as follows.

| Item | Earlier staging | Released dataset |
|---|---:|---:|
| RAW population | 2,580 CMLs | 2,574 CMLs, 5,148 sublink slots |
| QC population | 2,219 CMLs (4,438 sublinks) | 1,987 CMLs, 3,731 retained sublinks |
| Rainfall population | 2,219 CMLs (307 without any valid rain minute) | 1,987 CMLs (all with valid minutes) |
| QC window | whole recorded history | 1 April to 31 October 2021 |
| QC producer | six-step producer of `cml_preprocessed_v3.nc`, not recovered | `02_qc/scripts/s3_04_qc_corrected_rev1.py` |
| Graf thresholds and dry baseline calibrated on | 13-19 Jan and 1 Apr to 5 Nov 2021 | 1 April to 31 October 2021 |
| Rainfall exporter | `src/essd_cz/export_1min.py` (requires 2,219 links and a bit-for-bit match with the old 5-minute reference) | `03_rainfall_graf/scripts/s4_01_rain_1min.py` (same retrieval functions) |
| Gauge validation, gauge-wet | N 1,354,953, RMSE 3.33, MBE -0.93, PCC 0.585; POD 75.8 %, correct negatives 87.9 % | N 1,331,687, RMSE 3.63, MBE -0.91, PCC 0.633; POD 74.6 %, correct negatives 88.6 % |
| CML availability (mean of the 1-minute percentage) | 58.67 % | 66.35 % |

Staged files (in the author's `zenodo_staging/primary/`): RAW 2,257,209,979 B,
sha256 `46b130aa144a4e96d9698ca94c06ff625ad7670d8aa362b6e4c49a467806d012`;
QC 3,351,446,358 B, sha256
`badc3bc800904806c85ccc740c0f28448b63094911666216c038f56c5bc44d19`;
1-minute rainfall 493,609,050 B, sha256
`be2cd93843b6babbfe9d091de35d23db81071c66c7d9e1eee3ae03fee7a0780f`.

## 2. Reconciliation of the populations

**RAW, 2,580 to 2,574.** Six CMLs are non-operational (rule R1: no sublink
with paired valid TX and RX at a reported TX level): `10004_11582_3122`,
`10662_10795_3342`, `27040_27099_3747`, `29002_29004_3810`,
`41920_323338_4235`, `49111_329218_4481`. The 13 CMLs with TX = 0 dBm at
every minute were already outside the 2,580. All 19 are listed in
`01_population_raw/RAW_EXCLUDED_CMLS_ALL.csv`.

**QC, 2,219 to 1,987.** 1,804 CMLs are in both. 415 were only in the earlier
QC and 183 are only in the released QC (2,219 - 415 + 183 = 1,987).

* The 415: 255 metadata-only CMLs without observations, the 13 TX = 0 CMLs,
  4 R1 CMLs, 94 that fail the 5-hour/1-hour fluctuation test (the earlier
  producer used 21- and 5-minute windows and tested only sublink 0), 48
  without usable observations in April to October, and 1 whose only
  direction with data is a TX = 0 direction.
* The 183: short links that the earlier QC removed with a CML-level detection
  limit (ITU 2003 table, lowest frequency). With the per-sublink ITU 2005
  detection limit they pass.

Details: `02_qc/STAGE3_QC_REPORT.md` section 5 and
`02_qc/QC_BEFORE_AFTER_COMPARISON_REV1.csv`.

## 3. Comparison of the products

Rainfall (`03_rainfall_graf/STAGE4_RAINFALL_REPORT.md`): on 573,689,197 minutes
valid in both products, wet/dry agreement 99.67 %, pooled correlation of the
1-minute rates 0.957; per-CML accumulation on the 1,804 common CMLs 207.0 to
208.6 mm (mean); Graf threshold median ratio 1.003. Tables:
`RAINFALL_BEFORE_AFTER.csv`, `RAINFALL_BEFORE_AFTER_CML.csv`,
`GRAF_THRESHOLD_COMPARISON.csv`.

Validation (`04_validation/STAGE5A_REPORT.md`):

* Control: the validation code with the earlier inputs reproduces the
  earlier gauge validation exactly (every count identical; PCC to 1e-14).
  This used the helper
  `original/notebook_references/wetdry_validation/cml_graf_rstd_hess_radar_gauge_v3.py`.
* On the 30,750,879 samples present in both, the released processing is
  slightly better: gauge-wet RMSE 3.357 to 3.296, PCC 0.594 to 0.612.
* The higher full-sample RMSE comes from heavy-rain samples that now pass
  thanks to blackout reconstruction, and from the 183 short links.
* Tables: `RAINFALL_VALIDATION_METRICS.csv`, `RAIN_DETECTION_METRICS.csv`,
  `VALIDATION_BEFORE_AFTER.csv`, `EXTREMES_IN_VALIDATION.csv`.

Figures (`05_figures/STAGE5B_REPORT.md`, `FIGURE_PROVENANCE.csv`,
`FIGURE_NUMERICAL_COMPARISON.csv`): the producers of Figs. 2, 4, 5 and 6
reproduce the earlier figures pixel for pixel on the earlier inputs. The Fig.
3 producer was reproduced partly (its original input is not available). The
Fig. 1 producer was not found. The Fig. 7 grid builder
`05_figures/fig7_grid_builder/cml_grid_essd_publication_safe_singleopen.py`
was identified from the author's notebook archive; on the earlier rainfall it
gives 92 % of the CML-row pixels in the same colour class as the earlier
figure (strongly supported, not proven). The HPC folder name
`cmlonly_graf_nearby_comparison` in the Fig. 7 provenance is the original
workflow folder of the Graf grid.

Water vapour: unchanged (`04_validation/WVD_VERIFICATION_REPORT.md`).

## 4. Author decisions

| Date | Step | Decision |
|---|---|---|
| 2026-10-27 | Stage 1B | 13 TX = 0 CMLs excluded as non-operational; 255 empty CMLs not in QC |
| 2026-10-28 | Stage 2 | R1 approved: RAW = 2,574 CMLs / 5,148 sublink slots. P3: 8 CMLs without frequency kept in RAW, excluded before QC. P5: 15 TX = 0 sublinks kept in RAW, no valid QC observations. P6: 15 mixed-placeholder minutes kept in RAW, masked in QC |
| 2026-10-28 | Stage 3 | QC settings approved; revision 1 (no blackout reconstruction in gaps with plateau or P6 minutes) approved: `QC_CORRECTED_CANDIDATE_REV1.nc` |
| 2026-10-28 | Stage 4 | Rainfall approved as computed; blackout minutes kept and flagged; very high short-link rates documented, not filtered |
| 2026-10-28 | Stage 5A | Headline validation with blackout-reconstructed samples; the values without them reported alongside |
| 2026-10-29 | Stage 6 | Keep the Graf method. Figs. 2, 3 (with the approved axis correction), 5 and 6 regenerated; Fig. 8 unchanged |
| 2026-10-29 | Stage 6 | Figs. 1 and 7 to be regenerated from the released QC network and the released Graf rainfall (Fig. 1 redrawn in the published style; Fig. 7 with the unchanged IDW and renderer). The documentation describes the released dataset only |

The stage reports were written during the work; where a report says
"awaiting approval", the decision is the one in this table.

## 5. Verification of the earlier staging

These checks concern the earlier 2,219-CML staging. Their proof files stay
unchanged in `docs/evidence/hpc_run_2026-10/` and
`docs/evidence/local_audit_2026-10/`.

* **Rainfall rebuild.** The HPC run recomputed all 2,219 links over the full
  2020 to 2021 history with the unchanged original script
  (`original/code/cml_rainfall/cml_graf_itu2005_rebuild.py`; parameters and
  the 525 skipped series in
  `hpc_run_2026-10/01_rainfall_rebuild/rebuild_summary.json`). The local
  audit compared it with the validated 5-minute reference value by value
  (`rainfall_identity_report.json`, `rainfall_static_vars_report.json`):
  `R_sublink_5min` (137,782,581 finite values), `R_mean_5min` (84,656,385),
  flags, counters, thresholds and metadata identical. One-link rerun of
  `10004_311866_3129`: 7,666 shared samples, zero error
  (`one_link_comparison_report.json`).
* **Native 1-minute staging file** (`verify_1min_release_20261023.json`):
  308,160 regular slots; 648,844,578 finite directional rates; the 5-minute
  derivation (`src/essd_cz/rain5.py`) bit-identical to the 5-minute
  reference; geometry, thresholds and metadata bit-identical. The exporter is
  preserved as `docs/evidence/hpc_run_2026-10/cml_graf_export_1min_verified.py`.
* **Entity counts of the staging:** raw 2,580, preprocessed 2,219, rainfall
  2,219, radar 1,132, gauges 424, temperature 251, humidity 247; 4,438
  directional sublinks, 1,804 of them E band.
* **Gauge validation of the staging** (`hpc_run_2026-10/02_gauge_validation/.../FINAL_GRAF_validation_summary.json`):
  1,020 contributing CMLs, 218 gauges, N = 34.6 million; all period RMSE
  0.70 mm/h, correlation 0.60; wet periods RMSE 3.33, correlation 0.59;
  POD 0.758, FAR 0.796, MCC 0.350.
* **Availability of the staging:** CML 58.67 %, radar 99.48 %, gauges
  99.74 %, temperature 99.57 %, humidity 99.56 %
  (`hpc_run_2026-10/00_exploration/`,
  `local_audit_2026-10/availability_native_timestep_percentage.png`).
* **Release package:** an earlier 5.8 GiB archive contained a wrong generic
  rainfall file instead of the ITU 2005/Pastorek product; it was never
  published.

## 6. Comparison scripts

These scripts compare the released products with the earlier staging. They
are kept for the record; the released products do not depend on them.

```bash
W=FINAL_RECONSTRUCTION_20261027
# QC tables and the comparison table (needs the staged 2,219-CML QC)
python -I -B $W/stage3/scripts/s3_05_tables_checks_rev1.py --stage3 $W/stage3 --raw-nc <full-history RAW> --old-qc <staged QC>
# Rainfall comparison
python -I -B $W/stage4/scripts/s4_03_compare.py --old-rain <staged 1-min rainfall> --old-qc <staged QC>
# Control validation with the earlier inputs, then the sample-by-sample comparison
python -I -B $W/stage5a/scripts/s5a_01_gauge_validation.py --tag control_released \
    --cml-nc <full-history QC> --rain5-nc <full-history 5-min rainfall> --gauge-nc <full-history gauge> --out-dir $W/stage5a/control
python -I -B $W/stage5a/scripts/s5a_02_compare.py --stage5a $W/stage5a --old-cml <full-history QC>
```

`s3_04_qc_corrected_rev1.py` reads the staged QC file at a fixed path
(`OLD_QC_DEFAULT`) only to copy variable attributes (`time`, `sublink_id`,
`cml_id`, length, coordinates, frequency); see `README.md` section 6.
