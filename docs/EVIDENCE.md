# Evidence

This page lists what was verified, with the numbers and where the proof files
sit. Two sources of evidence are kept apart on purpose:

* **HPC run, October 2026.** The real data execution of the whole workflow on
  the original HPC. Small outputs are committed unchanged under
  `docs/evidence/hpc_run_2026-10/` (see `EVIDENCE_MANIFEST.json`). Large
  outputs stay in the project archive, not in git.
* **Local audit, October 2026.** An independent check on a second machine
  (Windows, see `ENVIRONMENT.md`) against the verified archive data. Reports
  under `docs/evidence/local_audit_2026-10/`.

## 1. The preserved code is unchanged

`python tools_verify_originals.py` confirms all 44 preserved scripts and
notebooks match their recorded SHA256 hashes. This was checked after every
restructuring step and on a fresh clone. The handoff archive itself was also
verified file by file before any work started (257 of 257 files matched).

## 2. The validated 5-minute reference is reproduced exactly

The HPC run recomputed all 2,219 links over the full 2020 to 2021 history
with the unchanged original script
(`original/code/cml_rainfall/cml_graf_itu2005_rebuild.py`; parameters and the
525 skipped series with under 2 valid observations are in
`hpc_run_2026-10/01_rainfall_rebuild/rebuild_summary.json`). The local audit
then compared that recomputed file against the validated 5-minute reference
(the author's internal benchmark, from which the released native 1-minute
product was in turn verified in section 2b), value by value
(`rainfall_identity_report.json`, `rainfall_static_vars_report.json`):

| Check | Result |
|---|---|
| `R_sublink_5min`, 137,782,581 shared finite values | identical, zero differences, identical masks |
| `R_mean_5min`, 84,656,385 values | identical |
| wet flags, validity flags, minute counters (up to 934,323,264 values each) | identical |
| per link thresholds `graf_threshold_db`, `graf_rstd_q80_db`, geometry, metadata | identical |

The earlier single link gate is contained in this result and was also rerun
from scratch on the second machine (`run.py compare-rainfall`, link
`10004_311866_3129`): 7,666 shared samples, zero error, zero mask mismatches
(`one_link_comparison_report.json`).

Caveat: this identity needs the full 2020 to 2021 history. The released April
to October data alone cannot rebuild the thresholds. The thresholds are
stored inside the rainfall file. There was no calibration against gauges or
radar; the thresholds and dry baselines are derived from the CML data.

## 2b. The native 1-minute release reproduces the 5-minute reference exactly

The published rainfall product is the native 1-minute NetCDF
(`CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc`), written
on the original HPC by the verified exporter (preserved byte-identical as
`docs/evidence/hpc_run_2026-10/cml_graf_export_1min_verified.py`; its
maintained, path-configurable copy is `src/essd_cz/export_1min.py`). The
exporter reruns the unchanged retrieval on the full 2020 to 2021 history and
compares every batch against the validated 5-minute file before keeping
anything.

The local audit (October 23, 2026,
`local_audit_2026-10/verify_1min_release_20261023.json`) then independently
verified the released file against the validated April to October 5-minute
benchmark:

| Check | Result |
|---|---|
| time axis: 308,160 regular 1-minute UTC slots, Apr 1 00:00 to Oct 31 23:59 | exact |
| finite 1-minute directional rates | 648,844,578, matching the HPC run |
| finite rate if and only if `graf_valid_1min` = 1; wet flag 255 exactly on invalid minutes | consistent everywhere |
| derived `R_sublink_5min` (five complete valid minutes per block) | bit-identical, zero differences |
| derived `R_mean_5min` (physical mean over finite directional 5-minute rates) | bit-identical |
| derived wet/valid flags and minute counters | identical |
| `R_mean_1min` versus mean of valid directional 1-minute rates | bit-identical |
| geometry, thresholds, metadata versus the reference | bit-identical |

Every 5-minute quantity used in the manuscript analysis is therefore exactly
recoverable from the released 1-minute file alone (`src/essd_cz/rain5.py`);
no separate 5-minute rainfall product is published.

## 3. Entity counts and coverage

`run.py check` on the real files: raw 2,580 CMLs, preprocessed 2,219,
rainfall 2,219, path radar 1,132, gauges 424, temperature 251, humidity 247.
All seven files hold a complete April to October axis at native resolution
(308,160 slots at 1 min, 61,632 at 5 min, 30,816 at 10 min).

## 4. Water vapour

Local audit of the support files (`wvd_support_audit.json`): the V6 NetCDF
holds 1,352 directional series over April to October; the diagnostics and
quantization CSVs describe exactly those 1,352 series. From the HPC run:
station matched validation within 5 km gives 587 directional links and
2,834,230 strict pairs
(`hpc_run_2026-10/WVD_V6_QUANTIZATION_BINNED_VALIDATION/run_summary.json`);
the paper population (E band, length at least 4 km, quantization at most
0.20 g per cubic metre, GRAF branch) is 45 directional sublinks from 31
physical CMLs with N = 236,793 and correlation 0.598
(`hpc_run_2026-10/07_wvd_paper_population/wvp_support_metrics.csv`). These
counts were measured, not forced.

## 5. Gauge validation (substitute detector, read the caveat)

`hpc_run_2026-10/02_gauge_validation/.../FINAL_GRAF_validation_summary.json`:
validation of the recomputed rainfall against CHMI gauges over April to
October 2021 (1,020 contributing physical CMLs, 218 gauges, N = 34.6
million): all period RMSE 0.70 mm/h and correlation 0.60; wet periods RMSE
3.33 and correlation 0.59; detection POD 0.758, FAR 0.796, MCC 0.350. The
exact historical detector used for the paper is not recovered, so these
numbers come from the documented substitute and must not be labelled
identical to the paper.

## 6. Availability

The unchanged availability script, run on the real files during the local
audit, reproduced every mean from the HPC walkthrough to the last digit:
CML 58.67 percent, path radar 99.48, gauges 99.74, temperature 99.57,
humidity 99.56. Figures and tables: `hpc_run_2026-10/00_exploration/` and
`local_audit_2026-10/availability_native_timestep_percentage.png`. The
network holds 4,438 directional sublinks, of which 1,804 are E band.

## 7. Tests

20 of 20 regression tests pass, and the synthetic seven file chain
(`tests/run_connected_mini.py`) runs end to end. These use synthetic
fixtures, so they prove plumbing, not science; the science evidence is the
sections above.

## 8. Corrected processing (October 2026)

Sections 2 to 6 describe the earlier 2,219-CML release. The corrected
processing of October 2026 keeps the Graf method and changes the population
and the quality control. Every step was approved by the author and
reproduced by an independent second implementation; the reports are in
`reconstruction_2026/`.

| Check | Result |
|---|---|
| RAW population | 2,574 CMLs, 5,148 sublink slots (2,580 minus 6 non-operational CMLs) |
| QC population | 1,987 CMLs, 3,731 retained sublinks; 2,574 = 8 + 579 + 1,987 |
| Old vs new QC | 1,804 CMLs in both, 415 only in the old QC, 183 only in the new |
| QC file | `QC_CORRECTED_CANDIDATE_REV1.nc`, sha256 `ca8afa95…`; second implementation on 70 CMLs bit-identical |
| 1-minute rainfall | 1,987 CMLs, 631,742,984 valid directional values, 0 retrieval errors, sha256 `9014567c…`; second implementation on 37 sublinks bit-identical |
| Validation control | the earlier inputs reproduce the published validation exactly (N 1,354,953, RMSE 3.33, MBE −0.93, PCC 0.585) |
| Validation, Case A (headline) | gauge-wet N 1,331,687, RMSE 3.63 mm/h, MBE −0.91 mm/h, PCC 0.633; POD 74.6 %, correct negatives 88.6 % |
| Validation, Case B (no blackout samples) | gauge-wet N 1,327,998, RMSE 3.38 mm/h, MBE −0.89 mm/h, PCC 0.606 |
| Identical samples (30,750,879) | gauge-wet RMSE 3.357 → 3.296, PCC 0.594 → 0.612 |
| Water vapour | inputs unchanged (RAW only); Table S13 and Fig. 8 unchanged |

The rerun commands, parameters and checksums are in
`reconstruction_2026/README.md`; the release files are listed in
`reconstruction_2026/DATA_RELEASE_MANIFEST.md`.
