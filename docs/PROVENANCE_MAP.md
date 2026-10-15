# Source and code provenance map

Classification of everything runnable in this repository. Categories:

* **UNCHANGED ORIGINAL** — byte-identical preserved scientific source
  (hash-checked by `tools_verify_originals.py` against
  `original/CONSOLIDATED_SOURCE_MANIFEST.json`).
* **ADAPTER** — new, clearly separated wrapper code that calls originals or
  re-applies original mathematics to the released products.
* **HISTORICALLY UNRESOLVED** — the exact historical producer is not recovered;
  no substitute is passed off as the original.
* **DEMONSTRATION ONLY** — ran successfully, but on reduced/synthetic inputs;
  not evidence of full-scale reproduction.

## Preserved originals (`original/`, 44 files) — UNCHANGED ORIGINAL

Key entries (full list in `original/CONSOLIDATED_SOURCE_MANIFEST.json`):

| File | Role |
|---|---|
| `code/cml_rainfall/cml_graf_itu2005_rebuild.py` | Final Graf RSTD + constant baseline + Pastorek WAA + ITU-R P.838-3 (2005) rainfall producer |
| `code/cml_qc/cml_full_fast_pipeline_no_availability_filter_parallel.py` + `related_versions/v3,v4_3,v5` | QC-related filters; **none confirmed as the March 2026 six-step QC producer** |
| `code/chmi/chmi_historical_csv_to_netcdf_v5.py` | CHMI CSV→NetCDF4 converter (notebook 08 uses it unchanged) |
| `code/radar/dwd_radolan_processor.py` | DWD RADKLIM replay script (original corrected run not proven) |
| `code/figures/02_blettner_aproct_dataset_figures_grid.py` | source of the two IDW functions reused by the stage-03 adapter |
| `code/figures/03_blettner_event_maps_timeseries.py` | original spatial map renderer (notebook 04 calls it unchanged) |
| `code/wvd/19_wvd_eband_progressive_baseline_v6.py`, `13_...`, `20_...` | V6 E-band WVD retrieval |
| `code/wvd/33_v6_binned_validation_hexbins.py`, `21_validate_wvd_v6_vs_original.py`, `plot_wvp_support_hexbins_no_titles.py` | WVD sensitivity, station matching and paper-population figures |
| `code/publication/prepare_essd_publication_datasets*.py`, `harmonize_essd_cml_radar_metadata_v2.py` | historical publication-packaging scripts |
| `notebooks/DATASET_PUBLICATION.ipynb`, `dataste_papaer.ipynb` + `notebook_recovery/`, `notebook_references/` | recovered historical notebooks and cell exports; **do not run as a pipeline** |

## Adapters (`src/essd_cz/`, `notebooks/00–08`, `tests/`) — ADAPTER

* `src/essd_cz/config.py`, `catalog.py`, `cli.py`, `workflows.py`,
  `compare.py`, `rainfall.py`, `metrics.py`, `plots.py`, `notebook_figures.py`,
  `station_workflow.py`, `connected.py` — path handling, staging manifest,
  CLI, and wrappers that invoke the unchanged originals.
* Stage-03 gridding (`connected.grid_final_rain`) loads the two original IDW
  functions **by AST from the preserved source**, then applies them to the
  final 5-min rainfall. The historical manuscript grid came from an earlier
  1-min workflow, so the grid output is adapted, **not asserted
  byte-identical** to historic maps.
* `results/.../33_v6_binned_validation_hexbins_NETCDF4_COMPAT.py` (evidence
  copy in `docs/evidence/`) — a traceable compatibility copy of the original
  WVD analysis changing only two single-series NetCDF reads; original and
  compat SHA256 recorded in `docs/evidence/hpc_run_2026-10/pipeline_manifest.json`
  and `indexing_compat_provenance.json`.

## Historically unresolved — HISTORICALLY UNRESOLVED

| Item | State |
|---|---|
| March 2026 six-step QC producer of `CZ_CML_PREPROCESSED_1MIN_*` | not recovered; July filter scripts are related but not attributable |
| Exact historical gauge wet/dry Graf detector used for the publication's gauge validation | not recovered; `cml_graf_rstd_hess_radar_gauge_v3.py` is a **labelled-unverified** substitute behind an explicit opt-in |
| Original corrected RADKLIM replay run | script recovered, run not proven |
| Historic 1-min Blettner gridding workflow for the manuscript's spatial figure | replaced by the documented 5-min adapter |
| Upstream raw telecom extraction log | unverified |

## Demonstration-only results — DEMONSTRATION ONLY

* `docs/evidence/hpc_run_2026-10/05_wvd/reviewer_wvd_v6_small*` — a 3-series
  V6 retrieval demonstration; **never** usable for paper figures (notebook 06
  refuses mismatched populations).
* `results/mini_connected_demo/` (regenerated locally) — fully synthetic
  7-file chain test.
* Notebook 04's explicit manual event selection — clearly marked; not the
  historical event discovery.

## Real-data evidence

Numeric results of the author's October 2026 HPC execution are summarised with
file-level provenance in [EVIDENCE.md](EVIDENCE.md).
