# Changes from the verified HPC handoff tree

> This page is a build log of the repository. Numbers in the passes before
> 2026-10-29 (for example 2,580 raw and 2,219 quality controlled links)
> describe an earlier, unpublished staging of the files, not the released
> dataset. See
> [reconstruction_2026/PROVENANCE_HISTORY.md](../reconstruction_2026/PROVENANCE_HISTORY.md).

This repository was assembled mechanically from the full handoff archive
`ESSD_CZ_2021_FULL_HANDOFF.tar`
(SHA256 `b9b76cf49f62189c2456ed88a5db50ed86c628355fc2ad8d2fdf8efa0e26565f`),
whose 257 files all verified against `handoff/FILE_MANIFEST_SHA256.json`
before any rearrangement. The build is scripted and re-verifiable; the build
report lists every copied pair with hashes.

**`original/` is untouched: all 44 preserved scientific sources remain
byte-identical (`python tools_verify_originals.py`).** No scientific formula,
population filter, threshold or unit was changed anywhere.

## Byte-identical relocations

| Handoff path | New path | Reason |
|---|---|---|
| `NOTEBOOK_STATUS.md`, `TEST_STATUS.md`, `TEST_LOG_NOTEBOOKS.md` | `docs/status/` | keep the repository root lean |
| `README.md`, `RUN_FIRST.txt`, `.gitignore` (handoff versions) | `attic/` | superseded by rewritten root files; preserved for history |
| `NOTEBOOK_ADDED_SHA256.json`, `NOTEBOOK_RELEASE_MANIFEST.txt` | `attic/` | describe an older eight-notebook naming scheme; stale for V2 |
| `upgrade_ESSD_wrapper_python_ONLY.py`, (`ESSD_CZ_2021_CONNECTED_WALKTHROUGH_V2_1.zip` not committed) | `attic/` (script only) | the V2→V2.1 migration was already applied on HPC (verified against the zip payload) |
| `notebooks/06_WVD_Sensitivity_Length_Frequency_REFRESHED_V6.ipynb` | `notebooks/06_WVD_Sensitivity_Length_Frequency.ipynb` | newest of three HPC iterations becomes the canonical 06; the two earlier iterations are preserved in `attic/notebook_versions/` |
| selected small outputs of `results/connected_workflow/` | `docs/evidence/hpc_run_2026-10/` | keep the real HPC execution evidence in git without committing large NetCDF/npz files |

## Generated files

`config/active_profile.json` now defaults to the portable `"user"` profile and
`results/connected_workflow_user` workspace (the handoff HPC version is in
`attic/active_profile_hpc_handoff.json`; the author's HPC interpreter remains
available under `science_python_by_profile.hpc`). New: `README.md`,
`LICENSING.md`, `CITATION.cff`, `.gitignore`, `.github/workflows/ci.yml`,
`docs/PROVENANCE_MAP.md`, `docs/EVIDENCE.md`, `docs/ENVIRONMENT.md`,
`examples/`, `attic/README.md`, this file.

## In-cell patches (the only content edits; wrapper/test level, no science)

Recorded with pre/post SHA256 in the build report.

1. `tests/test_connected_chain.py` — the reordered-IDs regression test wrote
   its temporary file with a forced scipy/NETCDF3 engine, which cannot store
   the uint8 wet/dry flags under current xarray; it now mirrors the fixture's
   own engine selection (netCDF4 when available). Failed before, passes after.
2. `notebooks/02_Rainfall_Validation.ipynb` — removed the hardcoded
   `/home/jovyan/...` interpreter fallback (config-driven instead); portability
   requirement enforced by the repository's own tests.
3. `notebooks/04_Spatial_Radar_Maps.ipynb` — `! pip install cartopy` IPython
   shell magic replaced by a comment (not valid Python; broke the repository's
   own notebook-syntax test).
4. Safe defaults restored where the author's executed HPC walkthrough had left
   expensive or unverified toggles enabled: `01 RECOMPUTE=False`,
   `02 RUN_GAUGE_VALIDATION=False`, `02 ALLOW_UNVERIFIED_V3=False`,
   `04 RUN_MAPS=False`, `06 RUN_SENSITIVITY=False`. Required by
   `tests/test_notebook_adapters.py` (and by the 06 notebook's own markdown).

Notebooks otherwise keep the author's executed HPC outputs embedded, as
evidence of the real run.

## Not committed to git (available in the handoff archive / Zenodo staging)

* The seven full-history source NetCDFs and WVD support files (`data/`).
* `results/connected_workflow/` large outputs: recomputed rainfall NetCDF
  (276 MB), IDW grid, `validation_pairs_graf.npz`, per-sensor availability CSV
  dumps, large PNG/PDF duplicates of committed PDFs.
* `results/V2_before_scientific_python_fix_20261014T103003Z/` (pre-upgrade
  backup; superseded, still in the archive).
* `results/mini_connected_demo/` (regenerated locally by
  `tests/run_connected_mini.py`).

## Simplification pass (2026-10-15, second pass)

Done on author instruction, after the first candidate was verified. No
scientific code and no data values changed; the 44 originals still verify.

* Root reduced to README, LICENSE, CITATION.cff, THIRD_PARTY_NOTICES,
  pyproject.toml, requirements.txt, run.py, tools_verify_originals.py.
  run.py and tools_verify_originals.py stay at the root because the notebooks
  call them there.
* `README_CONNECTED_WORKFLOW.md` moved to `docs/CONNECTED_WORKFLOW.md`;
  `CHANGES_FROM_HANDOFF.md` (this file) moved into `docs/`;
  `LICENSING.md` replaced by `docs/LICENCES.md` plus `THIRD_PARTY_NOTICES.md`.
* Removed from the public tree and preserved in the internal archive next to
  the repository (git history also keeps them): `attic/`, `docs/status/`,
  `tools/release_prep/`, `examples/`,
  `config/downloaded_release.example.json`, the split requirements files, the
  large archive verification report, the empty `figures/` folder.
* Added: MIT `LICENSE` (proposed), `THIRD_PARTY_NOTICES.md`, `pyproject.toml`
  (installable package, extras: science, chmi, notebooks, test),
  `data/README.md`, `docs/NOTEBOOK_GUIDE.md`, one "About this notebook"
  markdown cell at the top of each notebook.
* Edited for accuracy: `docs/01_FIRST_STEPS.md` and
  `docs/05_FIGURES_AND_WVD.md` now point at `config/user_paths.json`;
  `docs/ENVIRONMENT.md` points at the package extras; CI installs the package.
* `docs/EVIDENCE.md` and `docs/PROVENANCE_MAP.md` rewritten in plainer
  language, same facts.

## Release filename correction (2026-10-15, third pass)

The seven public data files were renamed so the names state the real
coverage: `*_20210401_20211031.nc` (and the CHMI files likewise). The staged
files were verified by their actual time coordinates before renaming and
revalidated against the sources after; no data value changed.
`src/essd_cz/config.py` now lists the release names in `EXPECTED_NAMES` and
keeps the author's full history names working through `FULL_HISTORY_NAMES`
(the hpc profile and `ACCEPTED_NAMES`). The user paths example, the test
fixture runner and the documentation use the new names.

## Native 1-minute rainfall release (2026-10-23, fourth pass)

The single published CML rainfall product became the **native 1-minute**
NetCDF `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc`,
produced on the original HPC by the verified exporter (handoff 2026-10-23).
The previously planned 5-minute release file is no longer published; it stays
an internal regression reference. The six other primary products are
byte-identical to the previous staging. No scientific formula, threshold,
population or unit changed; the 44 preserved originals still verify.

* New `src/essd_cz/rain5.py`: exact derivation of the analysis 5-minute
  series from the 1-minute release (five complete valid samples per block,
  any-wet flag, physical mean over finite directional 5-minute rates).
  Verified bit for bit against the validated 5-minute reference on all 2,219
  links (`docs/evidence/local_audit_2026-10/verify_1min_release_20261023.json`).
* New `src/essd_cz/export_1min.py`: the verified HPC exporter with the three
  hard-coded `/home/jovyan` path defaults replaced by required arguments; the
  byte-identical original is preserved under
  `docs/evidence/hpc_run_2026-10/cml_graf_export_1min_verified.py`. New CLI
  commands `derive-5min` and `export-1min-rainfall`
  (`workflows.export_native_1min`).
* `config.py` points the public `rain` role at the 1-minute name (the hpc
  profile keeps the full-history 5-minute reference); `catalog.py` inspects
  the 1-minute contract and the legacy reference separately.
* `connected.Project.source_rainfall()` registers the 1-minute release and
  materializes the derived 5-minute file once per workspace; `compare.py`
  derives the published 5-minute values per link when the rain role is the
  native release; `workflows.validation` derives before calling the original
  5-minute plotter.
* Tests: the synthetic rain fixture is now native 1-minute, with the
  corruption scenarios applied to 1-minute values; a new regression test pins
  the complete-block derivation semantics.
* Notebooks 01/02 and the documentation describe the 1-minute release and the
  derivation; no notebook result or figure population changed, because every
  derived 5-minute value is bit-identical to the previous product.

## Corrected processing (2026-10-29, fifth pass)

Added `reconstruction_2026/` with the author-approved corrected processing
(Graf method kept): approved RAW population of 2,574 CMLs, QC revision 1 with
1,987 CMLs and 3,731 retained sublinks, the 1-minute Graf rainfall from that
QC, the gauge validation (Case A headline: N 1,331,687, RMSE 3.63 mm/h,
MBE −0.91 mm/h, PCC 0.633) and the figure runners. All files there are
byte-identical copies from the reconstruction workspace
(`reconstruction_2026/FILE_PROVENANCE.csv`), including the Fig. 7 grid
builder, which was not in the repository before. The README, data
dictionary and evidence pages now give the corrected numbers. No file in
`original/`, `src/`, `notebooks/` or `docs/evidence/` changed; the 44
originals still verify. No NetCDF is in git; the proposed release files are
listed in `reconstruction_2026/DATA_RELEASE_MANIFEST.md`.
