# Corrected processing, October 2026 (Graf method)

This folder holds the code, parameters, tables and reports of the corrected
processing of the CML data, done between 27 and 29 October 2026. It replaces
the population, quality control (QC), rainfall and validation numbers of the
earlier release. The rainfall method is unchanged: Graf et al. (2020) RSTD
wet/dry detection, constant baseline, Pastorek wet antenna correction and the
ITU-R P.838-3 (2005) power law.

Every step was approved by the author before it ran, and every result was
checked by an independent second implementation. The earlier scripts in
`original/`, `src/` and `docs/evidence/` stay unchanged for provenance. They
describe how the earlier 2,219-CML release was made.

No NetCDF file is stored in git. The data files are listed, with sizes and
checksums, in [DATA_RELEASE_MANIFEST.md](DATA_RELEASE_MANIFEST.md). They will
be published in the Zenodo record (DOI to be added when the data release is
published).

## 1. What changed compared with the earlier release

| Item | Earlier release | Corrected processing |
|---|---:|---:|
| RAW population | 2,580 CMLs | **2,574 CMLs, 5,148 sublink slots** |
| QC population | 2,219 CMLs | **1,987 CMLs, 3,731 retained sublinks** |
| Rainfall population | 2,219 CMLs (307 without any valid rain minute) | **1,987 CMLs (all with valid minutes)** |
| QC window | whole recorded history | 1 April to 31 October 2021 |
| Graf thresholds and dry baseline calibrated on | 13-19 Jan and 1 Apr to 5 Nov 2021 | 1 April to 31 October 2021 |

**RAW (2,580 to 2,574).** Six CMLs are non-operational (rule R1: no sublink
with paired valid TX and RX at a reported TX level). The 13 CMLs with
TX = 0 dBm at every minute were already outside the 2,580. Together these 19
CMLs are listed in `01_population_raw/RAW_EXCLUDED_CMLS_ALL.csv`.

**QC (2,219 to 1,987).** 1,804 CMLs are in both products. 415 are only in the
earlier QC and 183 are only in the new QC (2,219 - 415 + 183 = 1,987).

* The 415 removed CMLs: 255 metadata-only CMLs without observations, the 13
  TX = 0 CMLs, 4 R1 CMLs, 94 that fail the corrected 5-hour/1-hour
  fluctuation test, 48 without usable observations in April to October, and
  1 whose only direction with data is a TX = 0 direction.
* The 183 added CMLs: short links that the earlier QC removed with a
  CML-level detection limit (ITU 2003 table, lowest frequency). With the
  per-sublink ITU 2005 detection limit they pass.
* Reconciliation: 2,574 = 8 (no frequency) + 579 (excluded by QC) + 1,987;
  5,148 = 16 + 1,158 + 243 + 3,731 sublinks.

Details: `02_qc/STAGE3_QC_REPORT.md` section 5.

## 2. Approved processing parameters

**Population (`01_population_raw/`, author decisions of 2026-10-27 and
2026-10-28).**

* R1: 6 non-operational CMLs removed from RAW and QC.
* P3: 8 CMLs without frequency stay in RAW, are excluded before QC; no
  frequency is invented.
* P5: 15 sublinks with TX = 0 at every minute stay in RAW; in QC they have no
  valid observations; the other direction is kept.
* P6: 15 mixed-placeholder minutes (1-minute means of a placeholder and a
  valid value) stay in RAW; QC masks TX and RX there (14 lie inside the QC
  window).

**QC (`02_qc/QC_PROCESSING_PARAMETERS.md` and `..._REV1.md`).**

| Step | Setting |
|---|---|
| Window | 2021-04-01 00:00 to 2021-10-31 23:59 UTC, 308,160 minutes |
| Invalid values | TX <= -99 dBm, TX = 255, RX <= -99.9 dBm, RX > -20 dBm, abs >= 1e35 -> missing (flag 1) |
| Low-RSL plateau | centred 3-minute window, >= 2 valid samples, max < -85 dBm and std < 0.5 dB; TX and RX masked +/- 2 min (flag 4) |
| Blackout (adaptation of Blettner et al., 2023) | RX gap <= 60 min with RX < -65 dBm on at least one edge; `trsl` = max(TX) - min(RX) of the sublink in the window; TX and RX are not filled; gaps that contain a plateau or P6 minute are not reconstructed (flag 8) |
| Short gaps | linear interpolation of `trsl` over <= 5 minutes missing in RAW, never across a QC mask or next to a blackout (flag 16) |
| Fluctuation | windows 5 h and 1 h, valid only with 300 / 60 valid minutes; fail if > 10 % of 5-h windows have std > 2 dB or > 33 % of 1-h windows have std > 0.8 dB |
| Detection limit | per sublink, R_min = (0.33 dB / (k L))^(1/alpha), ITU-R P.838-3 (2005), vertical; fail if R_min >= 2 mm/h |
| Retention | a sublink is removed if it fails a test (values set to missing); a CML is kept if at least one sublink passes |

**Rainfall (`03_rainfall_graf/RAINFALL_PROCESSING_PARAMETERS.md`).**

| Step | Setting (unchanged from the earlier release) |
|---|---|
| Wet/dry | Graf RSTD: centred 60-minute rolling std of `trsl`; threshold = 1.12 x the 0.8 quantile; 60 / 45 minimum valid minutes |
| Baseline | pycomlink `baseline_constant`, `n_average_last_dry = 1` |
| Wet antenna | Pastorek (2021): A_max 14 dB, zeta 0.55, d 0.1 |
| Rain rate | ITU-R P.838-3 (2005), vertical, sublink frequency and path length; dry valid minutes = 0 mm/h |
| 5-minute series | a directional 5-minute rate exists only if all five 1-minute values are valid; wet if any valid minute is wet; CML mean over finite directional 5-minute rates |

The rainfall uses `trsl` as approved, including blackout-reconstructed
(flag 8) and interpolated (flag 16) minutes. The new variable
`trsl_source_flag_1min` tells which minutes these are.

**Validation (`04_validation/`).** Same method as the earlier paper figures:
nearest CHMI gauge within 5 km of the CML midpoint, two complete 5-minute CML
bins against the end-labelled 10-minute gauge amount, gauge wet >= 0.1 mm/h.
The control run with the earlier inputs reproduced the published validation
exactly (N 1,354,953, RMSE 3.33, MBE -0.93, PCC 0.585).

| Gauge-wet samples | N | RMSE (mm/h) | MBE (mm/h) | PCC |
|---|---:|---:|---:|---:|
| **Case A: with blackout-reconstructed samples (headline)** | 1,331,687 | 3.63 | -0.91 | 0.633 |
| Case B: without blackout-reconstructed samples | 1,327,998 | 3.38 | -0.89 | 0.606 |

Detection (Case A): POD 74.6 %, false alarms 11.4 % of gauge-dry samples,
correct negatives 88.6 %, FAR 0.791, MCC 0.354. On the 30,750,879 samples
that are identical in both products, the corrected processing is slightly
better (gauge-wet RMSE 3.357 -> 3.296, PCC 0.594 -> 0.612).

Known limitation, kept and documented (no filter): 104 sublinks exceed
200 mm/h at 1-minute resolution (maximum 795.7 mm/h). Most are sub-kilometre
32-39 GHz links near the receiver floor. Such rates are not credible as
path-averaged rain.

**Water vapour (WVD).** Unchanged. The WVD code reads RAW only, no excluded
CML is among the WVD links, and RAW values did not change. Table S13 and
Fig. 8 stay as published (`04_validation/WVD_VERIFICATION_REPORT.md`). The
WVD code is in `original/code/wvd/`.

## 3. Folder layout

| Folder | Content | Workspace folder it came from |
|---|---|---|
| `01_population_raw/` | approved RAW population, exclusions, QC masks, Stage 1B/2 decisions, audit scripts | `stage1b/`, `stage2/` |
| `02_qc/` | QC script (revision 1), parameters, reports, per-sublink results, exclusion lists | `stage3/` |
| `03_rainfall_graf/` | 1-minute rainfall script, 5-minute derivation, comparison with the earlier product | `stage4/` |
| `04_validation/` | gauge validation script, metrics tables, WVD dependency check | `stage5a/` |
| `05_figures/` | figure runners, figure provenance, Fig. 7 grid builder | `stage5b/` |

`DATA_RELEASE_MANIFEST.md` (and `.csv`) lists the NetCDF files proposed for
the next data release, with sizes, checksums and the files they replace.

`FILE_PROVENANCE.csv` lists every file here with its source in the
reconstruction workspace (`FINAL_RECONSTRUCTION_20261027/`) and its SHA-256.
All files are byte-identical copies. Script outputs are kept next to each
script as `*_output.txt`.

The reports were written during the work. Where a report says "awaiting
approval", the decision that followed is in section 5 below.

## 4. Data files

Approved products (not in git):

| File (workspace name) | Bytes | SHA-256 |
|---|---:|---|
| `stage3/QC_CORRECTED_CANDIDATE_REV1.nc` | 2,409,368,783 | `ca8afa95aacb50dacbd5f7d96f776c56a78c63d1a4888e9f721a29cdcfafde65` |
| `stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc` | 445,051,457 | `9014567c3201be60441d1d8721c809d13ad5fa193e38a7b65606848b868a607e` |
| `stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc` | 338,126,319 | `a28d920be643d8647f21263a247b61e40496b793a44f453a902bdac8018820ee` |

Input of the QC: the full-history RAW
`CZ_CML_RAW_1MIN_20200101_20211231.nc` (sha256 `6acfeed2a2864b8c0e69807e9834dd03368d012fac4a2cae3d430c27b71be8aa`).

**There is no RAW NetCDF with exactly 2,574 CMLs yet.** The approved RAW
population exists as a table (`01_population_raw/RAW_APPROVED_POPULATION.csv`).
The earlier released RAW file holds 2,580 CMLs, that is the 2,574 plus the 6
R1 CMLs. See the manifest for what is still to be produced.

## 5. Author decisions

| Date | Step | Decision |
|---|---|---|
| 2026-10-27 | Stage 1B | RAW = 2,580 CMLs (13 TX = 0 CMLs excluded as non-operational; 255 empty CMLs not in QC) |
| 2026-10-28 | Stage 2 | R1 approved: RAW = 2,574 CMLs / 5,148 sublink slots; P3, P5, P6 as in section 2 |
| 2026-10-28 | Stage 3 | QC settings approved; revision 1 (no blackout reconstruction in gaps with plateau or P6 minutes) approved: `QC_CORRECTED_CANDIDATE_REV1.nc` |
| 2026-10-28 | Stage 4 | Rainfall approved as computed ("keep, report both"); blackout minutes kept and flagged; very high short-link rates documented, not filtered |
| 2026-10-28 | Stage 5A | Headline validation = Case A (with blackout samples); Case B reported alongside |
| 2026-10-29 | Stage 5B/6 | Keep the Graf method. Figures: 1 original, 2 corrected, 3 corrected (axis correction approved), 4 original, 5 and 6 regenerated (Case A), 7 original (no rerun), 8 unchanged |

## 6. How to rerun

The scripts were written for the author's reconstruction workspace. They
expect this layout of sibling folders, and several defaults are absolute
paths on the author's computer (section 7). To rerun, recreate the layout
from this folder, or pass every path on the command line:

| Workspace folder | Copy from |
|---|---|
| `stage1b/` | `01_population_raw/STAGE1B_EXCLUDED_13_CMLS.csv` |
| `stage2/` | `01_population_raw/` (tables and `scripts/`) |
| `stage3/` | `02_qc/` (tables and `scripts/`) |
| `stage4/` | `03_rainfall_graf/` (tables and `scripts/`) |
| `stage5a/` | `04_validation/` (tables and `scripts/`) |

Work in a copy outside the repository: several scripts write their tables
next to themselves or into the stage folder.

Requirements: Python 3.12 with numpy, pandas, xarray, netCDF4, numba and
pycomlink 0.6.0 (`pip install -e ".[science]"`). Always run with `python -I`
so that no local module is picked up, and add `-B` (or set
`PYTHONDONTWRITEBYTECODE=1`) when a script imports repository code, so that
no `__pycache__` is written into the repository. The run times below are
from the author's Windows computer (12 threads).

```bash
W=FINAL_RECONSTRUCTION_20261027            # your workspace
RAW=.../source_full_history/raw/CZ_CML_RAW_1MIN_20200101_20211231.nc

# 1. Population records (writes the four approved tables into stage2/)
python -I -B $W/stage2/scripts/s2_lead_02_approved_records.py --stage2 $W/stage2 --raw-nc $RAW

# 2. QC, revision 1 (about 30 min; writes QC_CORRECTED_CANDIDATE_REV1.nc and QC_SUBLINK_RESULTS_REV1.csv)
python -I -B $W/stage3/scripts/s3_04_qc_corrected_rev1.py --raw-nc $RAW --stage2 $W/stage2 --out-dir $W/stage3
python -I -B $W/stage3/scripts/s3_05_tables_checks_rev1.py --stage3 $W/stage3 --raw-nc $RAW --old-qc <earlier QC>
python -I -B $W/stage3/scripts/s3_03_file_verification.py --nc $W/stage3/QC_CORRECTED_CANDIDATE_REV1.nc

# 3. Rainfall (about 8 min), 5-minute derivation, comparison with the earlier product
python -I -B $W/stage4/scripts/s4_00_check_verbatim.py --repo <this repository>
python -I -B $W/stage4/scripts/s4_01_rain_1min.py \
    --input-nc $W/stage3/QC_CORRECTED_CANDIDATE_REV1.nc \
    --output-nc $W/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc
python -I -B $W/stage4/scripts/s4_02_derive_5min.py --repo <this repository> \
    --native $W/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc \
    --output $W/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc
python -I -B $W/stage4/scripts/s4_03_compare.py --old-rain <earlier 1-min rainfall> --old-qc <earlier QC>

# 4. Gauge validation (Case A and B in one run) and old/new comparison
python -I -B $W/stage5a/scripts/s5a_01_gauge_validation.py --tag stage4_new \
    --cml-nc $W/stage3/QC_CORRECTED_CANDIDATE_REV1.nc \
    --rain5-nc $W/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc \
    --gauge-nc .../source_full_history/gauge/CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc \
    --blackout-qc $W/stage3/QC_CORRECTED_CANDIDATE_REV1.nc --out-dir $W/stage5a/new
python -I -B $W/stage5a/scripts/s5a_02_compare.py --stage5a $W/stage5a --old-cml <earlier full-history QC>

# 5. Figures 5 and 6 from the saved validation samples (original plotting functions)
python -I -B 05_figures/fig5_fig6/s5b_01_fig5_fig6.py --samples $W/stage5a/new/stage4_new_samples.npz --out-dir <out>
```

`s4_01_rain_1min.py` refuses an input whose SHA-256 differs from the approved
QC file (`--expected-input-sha256`). The control validation (`--tag
control_released`) uses the earlier full-history QC and 5-minute rainfall as
inputs.

Figures 2 and 3: `05_figures/fig2_fig3/run_fig2.sh` runs the unchanged
`original/code/figures/availability_native_timestep_percentage.py` with only
`--cml-file` replaced by the QC REV1 file (it expects a copy of that script
next to it); `fig2_pixel_compare.py --apply-crop` applies the published
125-row top crop. `fig3_cml_panels_abc.py --population qc` redraws Fig. 3
from the QC REV1 file with the verbatim notebook cell
(`fig3_cell22_cml_dataset_runing.py`; changed lines in
`fig3_producer_diff.txt`). The axis correction of Fig. 3 approved in Stage 6
is not part of this folder.

Figure 4 and Figure 7 keep the published files. `05_figures/fig4_fig7/`
holds the scripts that proved this: `f4_01_event_plot_run.py` reproduces
Fig. 4 pixel for pixel and shows that the event CML is retained; `f7_*`
build the Fig. 7 CML grid with the grid builder below and compare it with the
published field.

The Fig. 7 chain is: runner `original/code/figures/06_run_multiple_events_vs_radar.py`
-> renderer `original/code/figures/03_blettner_event_maps_timeseries.py`; the
CML grid comes from `05_figures/fig7_grid_builder/cml_grid_essd_publication_safe_singleopen.py`.
That grid builder was not in this repository before. It is a verbatim copy
from the author's notebook archive (`DISDRODB/extracted/matched_code/matched_code/`,
sha256 `4cfb45975c9798f801339ccf406bc7f0e410c861d554bf93749990714167080d`).
Its attribution to the published grid is strongly supported but not proven
(92 % of CML-row pixels in the same colour class). Its argument defaults name
an older ITU 2003 input; they are labels for its own default input, and the
Stage 5B helper `f7_03_safe_singleopen_grid.py` uses only its three gridding
functions (`build_grid`, `load_midpoints`, `idw_fast_series`) on the Graf
5-minute rainfall. The HPC folder name `cmlonly_graf_nearby_comparison` that
appears in the Fig. 7 provenance is the original workflow folder of the
published Graf grid (it also appears in the original runner).

## 7. Hard-coded paths

The scripts are copied byte for byte, so their defaults still point to the
author's computer. All are plain defaults unless marked.

| Default | Used by | Override |
|---|---|---|
| `C:\Users\taouf\ESSD_RELEASE\extracted\full\ESSD_CZ_2021_HANDOFF\data\...` (full-history RAW, QC, rain, gauge, radar, T, H in `source_full_history\`; Fig. 7 event cache in `event_support\`) | s2_01, s2_lead_02, s3_04, s3_05, s5a_02, f4_01, f7_02, f7_05 | `--raw-nc`, `--raw`, `--old-cml`, `--cache`, or the per-file options of f4_01 |
| `C:\Users\taouf\ESSD_RELEASE\zenodo_staging\primary\...` (earlier staged QC and rainfall) | s2_02, s3_05, s4_03 | `--old-qc`, `--old-rain`, `--qc` |
| `C:\Users\taouf\ESSD_RELEASE\zenodo_staging\primary\CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` in `s3_04_qc_corrected_rev1.py` (`OLD_QC_DEFAULT`) | QC pass 2 copies the attributes of `time`, `sublink_id`, `cml_id`, length, coordinates and frequency from this earlier file | **no option**: the file must exist at this path, or the constant must be edited in a copy of the script |
| `C:\Users\taouf\ESSD_RELEASE\essd-czech-cml-2021\...` (this repository) | s4_00, s4_02, f7_01, f7_02 | `--repo` |
| `...\original\notebook_references\wetdry_validation\cml_graf_rstd_hess_radar_gauge_v3.py` in `s5a_01_gauge_validation.py` (`HELPER`) | Graf detection for the validation | **no option**: edit in a copy |
| `...\original\code\figures\rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py` in `s5b_01_fig5_fig6.py` (`PRODUCER`) | Figs. 5 and 6 plotting functions | **no option**: edit in a copy (f4_01 has `--producer`) |
| `C:\Users\taouf\Downloads\DISDRODB\extracted\matched_code\matched_code\cml_grid_essd_publication_safe_singleopen.py` | f7_03 | `--src 05_figures/fig7_grid_builder/cml_grid_essd_publication_safe_singleopen.py` |
| `C:\Users\taouf\ESSD_RELEASE\FINAL_RECONSTRUCTION_20261027\...` | fig1/fig2/fig3 helpers, f4_02, f7_05 | command-line options where present |
| `/home/jovyan/...` | grid builder argument defaults (HPC) | `--input-nc`, `--coords-nc`, `--output-nc` options of the grid builder |

The `*_output.txt` logs and the two validation summaries also contain these
paths. They are records of the actual runs and are kept unchanged.

## 8. Not included

* NetCDF and `.npz` files (too large for git; see the manifest).
* Independent-verification scripts and agent reports of the reconstruction
  workspace (`verification/`, `agents/`). Their results are summarised in each
  stage report.
* The first QC candidate (`QC_CORRECTED_CANDIDATE.nc`) and its scripts,
  except `QC_SUBLINK_RESULTS.csv`, which `s3_05` compares with.
* Figure image files. Their SHA-256 values are in
  `05_figures/FINAL_FIGURES_SHA256SUMS.txt`; the final images are in the
  manuscript.
* Diagnostic scripts for Fig. 1 (the published Fig. 1 is kept; its producer
  was not found).
* Exploratory method tests that were not adopted.
