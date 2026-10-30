# Provenance map

Everything runnable in this repository falls into one of four classes (plus
the processing scripts of the released products, described at the end):

* **Unchanged original.** A preserved scientific source, identical byte for
  byte, hash checked by `tools_verify_originals.py` against
  `original/CONSOLIDATED_SOURCE_MANIFEST.json`.
* **Adapter.** New wrapper code that calls the originals or applies their
  unchanged mathematics to the released files.
* **Not recovered.** The exact historical producer is missing. No substitute
  is passed off as the original.
* **Demonstration.** Ran successfully, but on reduced or synthetic inputs.
  Not evidence of full scale reproduction.

## Unchanged originals (`original/`, 44 files)

Key entries; the full list is in `original/CONSOLIDATED_SOURCE_MANIFEST.json`.

| File | Role |
|---|---|
| `code/cml_rainfall/cml_graf_itu2005_rebuild.py` | final rainfall producer: Graf wet and dry detection, constant baseline, Pastorek wet antenna correction, ITU R P.838 (2005) coefficients |
| `code/cml_qc/...` and `related_versions/` v3, v4_3, v5 | filters related to quality control. None of them is confirmed as the March 2026 six step producer. |
| `code/chmi/chmi_historical_csv_to_netcdf_v5.py` | CHMI CSV to NetCDF converter, used unchanged by notebook 08 |
| `code/radar/dwd_radolan_processor.py` | DWD RADKLIM replay script. The original corrected run is not proven. |
| `code/figures/02_blettner_aproct_dataset_figures_grid.py` | source of the two inverse distance weighting functions reused by notebook 03 |
| `code/figures/03_blettner_event_maps_timeseries.py` | spatial map renderer, called unchanged by notebook 04 |
| `code/wvd/19_...v6.py`, `13_...`, `20_...` | V6 E band water vapour retrieval |
| `code/wvd/33_...hexbins.py`, `21_validate...py`, `plot_wvp_support...py` | water vapour sensitivity, station matching and paper figures |
| `code/publication/prepare_essd_publication_datasets*.py`, `harmonize...v2.py` | historical packaging scripts |
| `notebooks/` and `notebook_recovery/`, `notebook_references/` | recovered historical notebooks and exports. Do not run them as a pipeline. |

## Adapters (`src/essd_cz/`, `notebooks/`, `tests/`, `run.py`)

Path handling, the pipeline manifest, the CLI, and wrappers that call the
unchanged originals. `src/essd_cz/export_1min.py` is a verified HPC
exporter with its hard-coded paths replaced by required arguments (the
byte-identical original sits in
`docs/evidence/hpc_run_2026-10/cml_graf_export_1min_verified.py`); it belongs
to the provenance history. The released rainfall was produced by
`reconstruction_2026/03_rainfall_graf/scripts/s4_01_rain_1min.py`, which
shares its scientific block byte for byte.
`src/essd_cz/rain5.py` applies the producer's unchanged aggregation
semantics to derive the analysis 5-minute series from the 1-minute release;
an independent re-derivation of the whole released file gave identical
values.
Notebook 03 loads the two original gridding functions directly from the
preserved source file and applies them to the derived 5 min rainfall. The
historical manuscript grid came from an older 1 min workflow, so the grid
output is an adapted product and no pixel identity is claimed. Notebook 06 uses a documented compatibility copy of the original
analysis script that changes only two NetCDF reads; its hash and the original
hash are recorded in the evidence files.

## Not recovered

| Item | State |
|---|---|
| the six step producer of an earlier staged QC file | not recovered (see `reconstruction_2026/PROVENANCE_HISTORY.md`); the released QC producer is known |
| the original file of the gauge wet and dry helper | not recovered; the released validation uses `cml_graf_rstd_hess_radar_gauge_v3.py`, behind an explicit switch in notebook 02 |
| the original corrected RADKLIM replay run | script recovered, run not proven |
| the historical 1 min gridding workflow behind the manuscript map | replaced by the documented adapter |
| the upstream raw telecom extraction log | not verified |

## Demonstrations

* `docs/evidence/hpc_run_2026-10/05_wvd/reviewer_wvd_v6_small*`: a 3 series
  water vapour rerun. Never usable for paper figures; notebook 06 refuses
  mismatched inputs.
* `tests/run_connected_mini.py` outputs: a fully synthetic seven file chain.
* The explicit manual event choice in notebook 04 when the historical event
  cache is absent.

## Processing of the released products (`reconstruction_2026/`)

Byte-identical copies of the scripts, tables and reports that produced the
released CML products (population, QC revision 1, Graf rainfall, validation,
figures).
`reconstruction_2026/FILE_PROVENANCE.csv` gives the source and SHA-256 of
every file. The QC script is new code that follows the standalone v3 filter
with documented changes; the rainfall script is a copy of
`src/essd_cz/export_1min.py` in which only the input checks and the new
source flag differ. The Fig. 7 grid builder
`reconstruction_2026/05_figures/fig7_grid_builder/cml_grid_essd_publication_safe_singleopen.py`
is a verbatim copy from the author's notebook archive. The history of the
processing is in `reconstruction_2026/PROVENANCE_HISTORY.md`. Nothing in
`original/` changed, and the new files are not added to
`original/CONSOLIDATED_SOURCE_MANIFEST.json`.

## Numbers and proof

See [EVIDENCE.md](EVIDENCE.md). The honest gap list above matches
[06_PROVENANCE_GAPS.md](06_PROVENANCE_GAPS.md) from the handoff.
