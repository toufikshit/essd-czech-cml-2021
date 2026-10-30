# Evidence

This page lists what was verified for the released dataset, with the numbers
and where the proof files sit. The processing reports, tables and scripts are
in `reconstruction_2026/`. Every processing step was checked by an
independent second implementation.

Older verification records (an HPC run and a local audit of an earlier,
unpublished staging, kept unchanged under `docs/evidence/`) are explained in
[reconstruction_2026/PROVENANCE_HISTORY.md](../reconstruction_2026/PROVENANCE_HISTORY.md).
Their numbers do not describe the released dataset.

## 1. The preserved code is unchanged

`python tools_verify_originals.py` confirms all 44 preserved scripts and
notebooks match their recorded SHA256 hashes. This was checked after every
restructuring step and on a fresh clone. The handoff archive itself was also
verified file by file before any work started (257 of 257 files matched).

## 2. Populations and quality control

`reconstruction_2026/02_qc/STAGE3_QC_REPORT.md`,
`STAGE3_REVISION1_REPORT.md`:

| Check | Result |
|---|---|
| RAW population | 2,574 CMLs, 5,148 sublink slots; 2,566 with frequency metadata, 904 E band |
| QC population | 1,987 CMLs, 3,731 retained sublinks, 810 E band CMLs; 2,574 = 8 + 579 + 1,987 and 5,148 = 16 + 1,158 + 243 + 3,731 |
| Every QC CML is in RAW | yes |
| QC file | `QC_CORRECTED_CANDIDATE_REV1.nc`, sha256 `ca8afa95…`; time 308,160 one-minute steps; `trsl` = `tx − rx` at every unflagged minute |
| Independent implementation | 70 CMLs: decisions, statistics, `tx`, `rx`, `trsl` and `qc_flag` bit-identical; detection limits of all 5,132 assessed sublinks identical |

## 3. Rainfall

`reconstruction_2026/03_rainfall_graf/STAGE4_RAINFALL_REPORT.md`:

| Check | Result |
|---|---|
| Code | `s4_01_rain_1min.py`; its scientific block (Graf, baseline, Pastorek, ITU 2005; 122 lines) is byte-identical to `src/essd_cz/export_1min.py` |
| Input | the QC file above, checked by sha256 |
| Output | 1,987 CMLs, 631,742,984 valid directional 1-minute values, 0 retrieval errors, sha256 `9014567c…` |
| Independent implementation | 37 sublinks, 11.4 M rates: bit-identical; wet, valid and source flags identical |
| Whole-file checks | ids equal the QC ids; excluded sublinks missing; finite rate if and only if `graf_valid_1min` = 1; `trsl_source_flag_1min` = `qc_flag & 24` at every minute |
| 5-minute derivation | the whole 5-minute series re-derived independently: all 244.9 M values identical |

The Graf thresholds and the dry baselines are computed over the QC period
(1 April to 31 October 2021) and are stored in the rainfall file. There is
no calibration against gauges or radar.

## 4. Entity counts and coverage

`run.py check` (metadata only) on the QC and rainfall files and the four
reference files: QC 1,987 CMLs, rainfall 1,987, path radar 1,132, gauges
424, temperature 251, humidity 247. All hold a complete April to October axis
at native resolution (308,160 slots at 1 min, 61,632 at 5 min, 30,816 at
10 min). The release RAW file (2,574 CMLs) is listed in
`reconstruction_2026/DATA_RELEASE_MANIFEST.md`.

## 5. Water vapour

Local audit of the support files (`wvd_support_audit.json`): the V6 NetCDF
holds 1,352 directional series over April to October; the diagnostics and
quantization CSVs describe exactly those 1,352 series. From the HPC run:
station matched validation within 5 km gives 587 directional links and
2,834,230 strict pairs
(`hpc_run_2026-10/WVD_V6_QUANTIZATION_BINNED_VALIDATION/run_summary.json`);
the paper population (E band, length at least 4 km, quantization at most
0.20 g per cubic metre, GRAF branch) is 45 directional sublinks from 31
physical CMLs with N = 236,793, RMSE 3.46 g per cubic metre and correlation
0.598 (`hpc_run_2026-10/07_wvd_paper_population/wvp_support_metrics.csv`).
These counts were measured, not forced. The WVD retrieval reads the RAW
signals only, and all WVD links are in the released RAW population
(`reconstruction_2026/04_validation/WVD_VERIFICATION_REPORT.md`).

## 6. Gauge validation

`reconstruction_2026/04_validation/STAGE5A_REPORT.md`: nearest CHMI gauge
within 5 km of the CML midpoint, two complete 5-minute CML bins against the
end-labelled 10-minute gauge amount, Graf detection from the QC `trsl` with
`original/notebook_references/wetdry_validation/cml_graf_rstd_hess_radar_gauge_v3.py`.
1,074 contributing CMLs, 2,015 sublinks, 217 gauges.

| Samples | N | RMSE (mm/h) | MBE (mm/h) | PCC |
|---|---:|---:|---:|---:|
| Gauge wet (≥ 0.1 mm/h), all samples | 1,331,687 | 3.63 | −0.91 | 0.633 |
| Gauge wet, without blackout-reconstructed intervals | 1,327,998 | 3.38 | −0.89 | 0.606 |
| All period, all samples | 34,249,231 | 0.768 | −0.019 | 0.634 |

Detection: POD 74.6 %, false alarms 11.4 % of gauge-dry samples, correct
negatives 88.6 %, FAR 0.791, CSI 0.195, MCC 0.354. The independent check
recomputed every metric from the saved samples and reproduced the samples
bit for bit for 30 CMLs.

## 7. Availability

The unchanged availability script on the QC file and the four reference
files gives these means of the native-step percentage: CML 66.35 percent,
path radar 99.48, gauges 99.74, temperature 99.57, humidity 99.56
(`reconstruction_2026/05_figures/fig2_fig3/fig2_numeric_summary_output.txt`).

## 8. Tests

21 of 21 regression tests pass, and the synthetic seven file chain
(`tests/run_connected_mini.py`) runs end to end. These use synthetic
fixtures, so they prove plumbing, not science; the science evidence is the
sections above.
