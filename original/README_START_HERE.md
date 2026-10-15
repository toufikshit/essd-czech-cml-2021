# ESSD Czech Republic — recovered source scripts and notebooks

This is the **combined SOURCE-CODE archive**, not a data deposit. It contains no NetCDF data. All 44 original files were copied byte-for-byte from the two archives supplied by the researcher, and their existing SHA-256 manifest values were checked before merging.

## Where to start

- `code/` — 24 original scripts, including CHMI v5 conversion, ESSD publication preparation, known QC implementations, DWD processing, final CML Graf/ITU2005/Pastorek rainfall, selected figures and WVD V6 analysis.
- `notebooks/` — 2 original dataset-publication notebooks.
- `notebook_references/` — 18 additional original scripts referenced in notebook cells. These mainly document earlier comparison/validation and exploratory workflows. They are **not automatically final publication producers**.
- `notebook_recovery/` — preserved, exact notebook-cell extracts: 2 standalone code extractions and shell/command snippets. These are derived from the notebooks and **not** newly found original producer scripts. Referenced-script gap listings reflect the state of the EARLIER archive; all 18 referenced filenames have now been added.
- `manifest.json` and `notebook_references/source_manifest.json` — original untouched collection manifests; `CONSOLIDATED_SOURCE_MANIFEST.json` and `SHA256SUMS.txt` cover this merged archive.

## Important provenance constraints

- The original `cml_preprocessed_v3.nc` source dataset was created 2026-03-11 and records six QC steps and 2,219 CML links. The source-only recovery does **not** prove which original March 2026 executable generated it. The standalone QC v3/v4_3/v5 and full-fast scripts currently present are **related later versions**, not proven original producers.
- `DATASET_PUBLICATION.ipynb` records an executed `prepare_essd_publication_datasets_fast.py` command; it operates on `cml_preprocessed_v3.nc` and does not create the March QC dataset.
- A later corrected-radar event-map runner (`code/figures/06_run_multiple_events_vs_radar.py`) is included, but its role in generating the released figure is not proven by filename/path alone.
- The final intended rain retrieval uses Graf RSTD, Pastorek WAA and ITU-R P.838-3 (2005), not the older ITU2003 experimental files in `notebook_references/`.
- The seven data products belong to a **separate release**. A previously created 5.8 GiB TAR selected `CZ_CML_RAINRATE_5MIN_20200101_20211231.nc`, not the designated final Graf/ITU2005/Pastorek rainfall product; that TAR must not be treated as the approved final deposit.
- No original raw telecom extraction log/script or definitive RADKLIM upstream processing run was recovered by merging these archives.

## Distribution/licensing

The datasets and the independently authored source code must be licensed separately. As specified in the project, **CML-derived datasets: CC BY-NC 4.0; CHMI/DWD datasets: CC BY 4.0** (with required attribution). Those data terms do **not** automatically license the code or notebooks. Confirm rights for any public distribution. Review notebook outputs for secrets, private paths, or sensitive data.

## Integrity

All original entries were verified against their source manifests (26 + 18, 44/44); no duplicate original source file bytes across archives. `SHA256SUMS.txt` contains per-file checksums of all archive members other than itself. Check extracted files with `sha256sum -c SHA256SUMS.txt` from the ZIP root. `CONSOLIDATED_SOURCE_MANIFEST.json` describes original location, role and byte-level hash for the 44 originals plus 2 notebook-extracted standalone code cells.

This archive does not execute any workflows, change the science files, or claim that manuscript figures have been regenerated.
