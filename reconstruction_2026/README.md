# Processing of the released CML products

This folder holds the code, parameters, tables and reports that produced the
released CML products: the raw population, the quality controlled (QC)
signals, the 1-minute rainfall, the gauge validation and the figures. The work
was done between 27 and 29 October 2026. Every step was approved by the
author before it ran, and every result was checked by an independent second
implementation.

The rainfall method is Graf et al. (2020) RSTD wet/dry detection, a constant
baseline, the Pastorek wet antenna correction and the ITU-R P.838-3 (2005)
power law.

No NetCDF file is stored in git. The data files are listed, with sizes and
checksums, in [DATA_RELEASE_MANIFEST.md](DATA_RELEASE_MANIFEST.md). They will
be published in the Zenodo record (DOI to be added when the data release is
published). The history of the processing is in
[PROVENANCE_HISTORY.md](PROVENANCE_HISTORY.md).

## 1. Released populations

| Product | Population |
|---|---|
| RAW signals | 2,574 CMLs, 5,148 sublink slots; 2,566 CMLs with frequency metadata, 904 of them E band |
| QC signals | 1,987 CMLs, 3,731 retained sublinks; 810 CMLs E band |
| Rainfall, 1 minute | 1,987 CMLs, 631,742,984 valid directional values |
| Radar along CML paths | 1,132 paths |
| CHMI stations | 424 rain gauges, 251 temperature, 247 humidity |

Population accounting: 2,574 = 8 (no frequency, excluded before QC) + 579
(excluded by QC) + 1,987; 5,148 = 16 + 1,158 + 243 (excluded sublinks in
retained CMLs) + 3,731 sublinks. Lists: `01_population_raw/` and
`02_qc/QC_EXCLUDED_CMLS_REV1.csv`, `02_qc/QC_EXCLUDED_SUBLINKS_REV1.csv`.

## 2. Processing parameters

**Population (`01_population_raw/`).**

* The raw population is the table `RAW_APPROVED_POPULATION.csv`.
  Non-operational CMLs (no sublink with paired valid TX and RX at a reported
  TX level) are not part of it; they are listed in
  `RAW_EXCLUDED_CMLS_ALL.csv`.
* 8 CMLs without frequency stay in RAW and are excluded before QC; no
  frequency is invented.
* 15 sublinks with TX = 0 at every minute stay in RAW; in QC they have no
  valid observations; the other direction is kept
  (`QC_PREDEFINED_SUBLINK_EXCLUSIONS.csv`).
* 15 mixed-placeholder minutes (1-minute means of a placeholder and a valid
  value) stay in RAW; QC masks TX and RX there, 14 of them inside the QC
  window (`QC_PREDEFINED_MINUTE_MASKS.csv`).

**QC (`02_qc/QC_PROCESSING_PARAMETERS.md` and `..._REV1.md`).**

| Step | Setting |
|---|---|
| Window | 2021-04-01 00:00 to 2021-10-31 23:59 UTC, 308,160 minutes |
| Invalid values | TX <= -99 dBm, TX = 255, RX <= -99.9 dBm, RX > -20 dBm, abs >= 1e35 -> missing (flag 1) |
| Low-RSL plateau | centred 3-minute window, >= 2 valid samples, max < -85 dBm and std < 0.5 dB; TX and RX masked +/- 2 min (flag 4) |
| Blackout (adaptation of Blettner et al., 2023) | RX gap <= 60 min with RX < -65 dBm on at least one edge; `trsl` = max(TX) - min(RX) of the sublink in the window; TX and RX are not filled; gaps that contain a plateau or mixed-placeholder minute are not reconstructed (flag 8) |
| Short gaps | linear interpolation of `trsl` over <= 5 minutes missing in RAW, never across a QC mask or next to a blackout (flag 16) |
| Fluctuation | windows 5 h and 1 h, valid only with 300 / 60 valid minutes; fail if > 10 % of 5-h windows have std > 2 dB or > 33 % of 1-h windows have std > 0.8 dB |
| Detection limit | per sublink, R_min = (0.33 dB / (k L))^(1/alpha), ITU-R P.838-3 (2005), vertical; fail if R_min >= 2 mm/h |
| Retention | a sublink is removed if it fails a test (values set to missing); a CML is kept if at least one sublink passes |

**Rainfall (`03_rainfall_graf/RAINFALL_PROCESSING_PARAMETERS.md`).**

| Step | Setting |
|---|---|
| Wet/dry | Graf RSTD: centred 60-minute rolling std of `trsl`; threshold = 1.12 x the 0.8 quantile; 60 / 45 minimum valid minutes; computed over 1 April to 31 October 2021 |
| Baseline | pycomlink `baseline_constant`, `n_average_last_dry = 1` |
| Wet antenna | Pastorek (2021): A_max 14 dB, zeta 0.55, d 0.1 |
| Rain rate | ITU-R P.838-3 (2005), vertical, sublink frequency and path length; dry valid minutes = 0 mm/h |
| 5-minute series | a directional 5-minute rate exists only if all five 1-minute values are valid; wet if any valid minute is wet; CML mean over finite directional 5-minute rates |

The rainfall uses the QC `trsl`, including blackout-reconstructed (flag 8)
and interpolated (flag 16) minutes. The variable `trsl_source_flag_1min`
tells which minutes these are. The 5-minute series are not a release file;
they are derived from the 1-minute file (`src/essd_cz/rain5.py`).

**Validation (`04_validation/`).** Nearest CHMI gauge within 5 km of the CML
midpoint; two complete 5-minute CML bins against the end-labelled 10-minute
gauge amount; gauge wet >= 0.1 mm/h; Graf detection from the QC `trsl`.

| Gauge-wet samples | N | RMSE (mm/h) | MBE (mm/h) | PCC |
|---|---:|---:|---:|---:|
| **All samples (headline)** | 1,331,687 | 3.63 | -0.91 | 0.633 |
| Without blackout-reconstructed intervals | 1,327,998 | 3.38 | -0.89 | 0.606 |

Detection (all samples): POD 74.6 %, false alarms 11.4 % of gauge-dry
samples, correct negatives 88.6 %, FAR 0.791, MCC 0.354. Contributing: 1,074
CMLs, 2,015 sublinks, 217 gauges.

Known limitation, documented and not filtered: 104 sublinks exceed 200 mm/h
at 1-minute resolution (maximum 795.7 mm/h). Most are sub-kilometre
32-39 GHz links near the receiver floor. Such rates are not credible as
path-averaged rain.

**Water vapour (WVD).** The WVD retrieval reads the RAW signals only. All WVD
links are in the released RAW population. The validation subset is 31 CMLs
and 45 directional sublinks with 236,793 pairs, RMSE 3.46 g m-3 and PCC 0.598
(`04_validation/WVD_VERIFICATION_REPORT.md`). The WVD code is in
`original/code/wvd/`.

## 3. Folder layout

| Folder | Content | Workspace folder it came from |
|---|---|---|
| `01_population_raw/` | raw population, exclusions, QC masks, population decisions, audit scripts | `stage1b/`, `stage2/` |
| `02_qc/` | QC script (revision 1), parameters, reports, per-sublink results, exclusion lists | `stage3/` |
| `03_rainfall_graf/` | 1-minute rainfall script, 5-minute derivation | `stage4/` |
| `04_validation/` | gauge validation script, metrics tables, WVD check | `stage5a/` |
| `05_figures/` | figure runners, figure provenance, Fig. 3 axis correction, Fig. 7 grid builder | `stage5b/` |

`FILE_PROVENANCE.csv` lists every copied file with its source in the
reconstruction workspace (`FINAL_RECONSTRUCTION_20261027/`) and its SHA-256.
All are byte-identical copies. Script outputs are kept next to each script as
`*_output.txt`. The stage reports and CSV tables are records of the work and
are kept unchanged; the context of their comparison tables is in
[PROVENANCE_HISTORY.md](PROVENANCE_HISTORY.md).

## 4. Data files

Approved products (not in git):

| File (workspace name) | Bytes | SHA-256 |
|---|---:|---|
| `stage3/QC_CORRECTED_CANDIDATE_REV1.nc` | 2,409,368,783 | `ca8afa95aacb50dacbd5f7d96f776c56a78c63d1a4888e9f721a29cdcfafde65` |
| `stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc` | 445,051,457 | `9014567c3201be60441d1d8721c809d13ad5fa193e38a7b65606848b868a607e` |
| `stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc` (not released) | 338,126,319 | `a28d920be643d8647f21263a247b61e40496b793a44f453a902bdac8018820ee` |

Input of the QC: the full-history RAW
`CZ_CML_RAW_1MIN_20200101_20211231.nc` (sha256 `6acfeed2a2864b8c0e69807e9834dd03368d012fac4a2cae3d430c27b71be8aa`),
restricted to the raw population. The release RAW NetCDF with the 2,574 CMLs
is still to be written; see the manifest.

## 5. How to rerun

The scripts were written for the author's reconstruction workspace. They
expect this layout of sibling folders, and several defaults are absolute
paths on the author's computer (section 6). To rerun, recreate the layout
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

# 1. Population records (writes the four population tables into stage2/)
python -I -B $W/stage2/scripts/s2_lead_02_approved_records.py --stage2 $W/stage2 --raw-nc $RAW

# 2. QC, revision 1 (about 30 min; writes QC_CORRECTED_CANDIDATE_REV1.nc and QC_SUBLINK_RESULTS_REV1.csv)
python -I -B $W/stage3/scripts/s3_04_qc_corrected_rev1.py --raw-nc $RAW --stage2 $W/stage2 --out-dir $W/stage3
python -I -B $W/stage3/scripts/s3_03_file_verification.py --nc $W/stage3/QC_CORRECTED_CANDIDATE_REV1.nc

# 3. Rainfall (about 8 min) and the 5-minute derivation
python -I -B $W/stage4/scripts/s4_00_check_verbatim.py --repo <this repository>
python -I -B $W/stage4/scripts/s4_01_rain_1min.py \
    --input-nc $W/stage3/QC_CORRECTED_CANDIDATE_REV1.nc \
    --output-nc $W/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc
python -I -B $W/stage4/scripts/s4_02_derive_5min.py --repo <this repository> \
    --native $W/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc \
    --output $W/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc

# 4. Gauge validation (with and without blackout-reconstructed intervals in one run)
python -I -B $W/stage5a/scripts/s5a_01_gauge_validation.py --tag stage4_new \
    --cml-nc $W/stage3/QC_CORRECTED_CANDIDATE_REV1.nc \
    --rain5-nc $W/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc \
    --gauge-nc .../source_full_history/gauge/CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc \
    --blackout-qc $W/stage3/QC_CORRECTED_CANDIDATE_REV1.nc --out-dir $W/stage5a/new

# 5. Figures 5 and 6 from the saved validation samples (original plotting functions)
python -I -B 05_figures/fig5_fig6/s5b_01_fig5_fig6.py --samples $W/stage5a/new/stage4_new_samples.npz --out-dir <out>
```

`s4_01_rain_1min.py` refuses an input whose SHA-256 differs from the approved
QC file (`--expected-input-sha256`). The QC exclusion lists are written by
`s3_05_tables_checks_rev1.py`, whose arguments are given in
[PROVENANCE_HISTORY.md](PROVENANCE_HISTORY.md) section 6.

**Figures.**

* Fig. 1: `05_figures/fig1_fig7/fig1_study_area.py` draws the 1,987 QC CMLs,
  the 424 rain gauges, the 251 weather stations (temperature and humidity
  positions) and the radar coverage in the published style. The radar coverage
  is the set of RADKLIM YW 2017.002 (corrected 2021 supplement) cells with data
  on at least one day from 1 April to 13 July 2021; the 1,132 radar paths lie
  inside it.
* Fig. 2: `05_figures/fig2_fig3/run_fig2.sh` runs the unchanged
  `original/code/figures/availability_native_timestep_percentage.py` with
  `--cml-file` set to the QC file (it expects a copy of that script next to
  it); `fig2_pixel_compare.py --apply-crop` applies the published 125-row top
  crop.
* Fig. 3: `05_figures/fig3_axis_correction/fig3_axis_corrected.py --population qc`
  draws the 1,987 QC CMLs and their 3,731 sublinks with the verbatim notebook
  cell (`fig2_fig3/fig3_cell22_cml_dataset_runing.py`) and the approved axis
  correction.
* Fig. 4: `05_figures/fig4_fig7/f4_01_event_plot_run.py` runs the original
  event-plot producer on the QC and 5-minute rainfall (event CML
  `45035_45097_2740`, 13-14 July 2021).
* Figs. 5 and 6: command above.
* Fig. 7: the released Graf rainfall gridded with the unchanged IDW (nnear 15,
  power 2, maximum distance 0.3 degree, 0.01 degree grid) and drawn with the
  unchanged renderer. Chain: runner
  `original/code/figures/06_run_multiple_events_vs_radar.py` -> renderer
  `original/code/figures/03_blettner_event_maps_timeseries.py`; grid builder
  `05_figures/fig7_grid_builder/cml_grid_essd_publication_safe_singleopen.py`
  (verbatim copy from the author's notebook archive, sha256
  `4cfb45975c9798f801339ccf406bc7f0e410c861d554bf93749990714167080d`). The
  helper `fig4_fig7/f7_03_safe_singleopen_grid.py` uses only its three
  gridding functions (`build_grid`, `load_midpoints`, `idw_fast_series`) on the
  Graf 5-minute rainfall; the grid builder's own argument defaults are not
  used. `05_figures/fig1_fig7/fig7_run.sh` draws the figure (event 4,
  instance 1, 13 July 2021 16:05-23:05 UTC) with the radar row from the
  corrected RADKLIM YW 2017.002 2021 supplement (`YW_2017.002_20210713.nc`).
  Run outputs and figure checksums are in the same folder.
* Fig. 8: water vapour, `original/code/wvd/plot_wvp_support_hexbins_no_titles.py`.

## 6. Hard-coded paths

The scripts are copied byte for byte, so their defaults still point to the
author's computer. All are plain defaults unless marked.

| Default | Used by | Override |
|---|---|---|
| `C:\Users\taouf\ESSD_RELEASE\extracted\full\ESSD_CZ_2021_HANDOFF\data\...` (full-history RAW, QC, rain, gauge, radar, T, H in `source_full_history\`; Fig. 7 event cache in `event_support\`) | s2_01, s2_lead_02, s3_04, s3_05, s5a_02, f4_01, f7_02, f7_05 | `--raw-nc`, `--raw`, `--old-cml`, `--cache`, or the per-file options of f4_01 |
| `C:\Users\taouf\ESSD_RELEASE\zenodo_staging\primary\...` (staged files) | s2_02, s3_05, s4_03 | `--old-qc`, `--old-rain`, `--qc` |
| `C:\Users\taouf\ESSD_RELEASE\zenodo_staging\primary\CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` in `s3_04_qc_corrected_rev1.py` (`OLD_QC_DEFAULT`) | QC pass 2 copies the attributes of `time`, `sublink_id`, `cml_id`, length, coordinates and frequency from this file | **no option**: the file must exist at this path, or the constant must be edited in a copy of the script |
| `C:\Users\taouf\ESSD_RELEASE\essd-czech-cml-2021\...` (this repository) | s4_00, s4_02, f7_01, f7_02 | `--repo` |
| `...\original\notebook_references\wetdry_validation\cml_graf_rstd_hess_radar_gauge_v3.py` in `s5a_01_gauge_validation.py` (`HELPER`) | Graf detection for the validation | **no option**: edit in a copy |
| `...\original\code\figures\rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py` in `s5b_01_fig5_fig6.py` (`PRODUCER`) | Figs. 5 and 6 plotting functions | **no option**: edit in a copy (f4_01 has `--producer`) |
| `C:\Users\taouf\Downloads\DISDRODB\extracted\matched_code\matched_code\cml_grid_essd_publication_safe_singleopen.py` | f7_03 | `--src 05_figures/fig7_grid_builder/cml_grid_essd_publication_safe_singleopen.py` |
| `C:\Users\taouf\ESSD_RELEASE\FINAL_RECONSTRUCTION_20261027\...` | fig2/fig3 helpers, f4_02, f7_05 | command-line options where present |
| `/home/jovyan/...` | grid builder argument defaults (HPC) | `--input-nc`, `--coords-nc`, `--output-nc` options of the grid builder |
| `C:Users	aoufESSD_RELEASE...` (QC file, released station files; RADKLIM supplement days in `external_downloads
| `C:\Users\taouf\ESSD_RELEASE\...` (QC file, released station files; RADKLIM supplement days in `external_downloads\radklim_2021_supplement`) | fig1_study_area.py, fig7_run.sh | `--qc`, `--gauge`, `--temp`, `--hum`, `--radar-nc` of fig1; the variable `R` in fig7_run.sh |

The `*_output.txt` logs and the two validation summaries also contain these
paths. They are records of the actual runs and are kept unchanged.

## 7. Not included

* NetCDF and `.npz` files (too large for git; see the manifest).
* Independent-verification scripts and agent reports of the reconstruction
  workspace (`verification/`, `agents/`). Their results are summarised in each
  stage report.
* Figure image files. Their SHA-256 values are in
  `05_figures/FINAL_FIGURES_SHA256SUMS.txt` and
  `05_figures/fig3_axis_correction/README.md`; the final images are in the
  manuscript.
* Exploratory method tests that were not adopted.
