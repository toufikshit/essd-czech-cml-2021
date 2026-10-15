# Changes from the verified HPC handoff tree

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
