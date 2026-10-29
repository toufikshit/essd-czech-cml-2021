# 6. Source-code provenance and unverified portions

These gaps are tracked rather than silently filled by convenient alternatives. The history of the processing, including an earlier unpublished staging of the files, is in [reconstruction_2026/PROVENANCE_HISTORY.md](../reconstruction_2026/PROVENANCE_HISTORY.md).

| Chain stage | Source known? | Reason for restraint |
|---|---|---|
| Original raw telecom data acquisition | No definitive original extraction log | The RAW product (2,574 CMLs) is taken from the full-history raw file restricted to `reconstruction_2026/01_population_raw/RAW_APPROVED_POPULATION.csv`; upstream acquisition script unverified |
| QC producer | Yes | `reconstruction_2026/02_qc/scripts/s3_04_qc_corrected_rev1.py` (1,987 CMLs); parameters in `reconstruction_2026/02_qc/` |
| Final ITU2005/Pastorek rainfall | Yes, script `reconstruction_2026/03_rainfall_graf/scripts/s4_01_rain_1min.py` | Same retrieval functions as `cml_graf_itu2005_rebuild.py` and `src/essd_cz/export_1min.py`; reproduced bit for bit by an independent implementation on 37 sublinks |
| Gauge-validation Graf helper | Used: `cml_graf_rstd_hess_radar_gauge_v3.py` | The validation of the released rainfall uses this helper (`reconstruction_2026/04_validation/`); its exactness check is described in the provenance history |
| DWD 2021 supplemental radar replay | Recovered DWD script but original corrected run not proven | Released path-radar NetCDF is valid reference; do not claim the older script exactly reproduces corrected RADKLIM |
| WVD V6 end-to-end | Partial | Additional supporting files and steps are needed beyond the seven release products |
| Physically April–October-only public NetCDFs | Verified for QC, rainfall and the four reference files | each holds exactly the April–October axis (308,160 / 61,632 / 30,816 slots); the RAW release file is still to be written |
| Dataset release package | Needs final check | file list, sizes and checksums in `reconstruction_2026/DATA_RELEASE_MANIFEST.md` |

The original source tree and historical notebooks are preserved, including outdated models and hard-coded paths. **Do not run all archived notebooks as a pipeline.** The new CLI is the reviewer-safe path. All 44 originals are verified against their original manifests; this check is separate from validating published numerical values.

Data licences from the project: CC BY-NC 4.0 for CML-derived files; CC BY 4.0 for CHMI and DWD. A software licence must be selected with the rights holder before public code distribution; no software licence is implied here.
