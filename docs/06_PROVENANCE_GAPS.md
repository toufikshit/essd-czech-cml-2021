# 6. Source-code provenance and unverified portions

These gaps are tracked rather than silently filled by convenient alternatives.

| Chain stage | Source known? | Reason for restraint |
|---|---|---|
| Original raw telecom data acquisition | No definitive original extraction log | Published 2,580-link RAW product and `prepare_essd_publication_datasets.py` identified; upstream acquisition script unverified |
| March 11, 2026 six-step QC producer | **Not confirmed** | `cml_preprocessed_v3.nc` (2,219 links) records six steps. July scripts v3/v4_3/v5 implement related filters but cannot be attributed to March simply by version name |
| Final ITU2005/Pastorek rainfall | Yes, script `cml_graf_itu2005_rebuild.py` | Real-HPC component execution passed; whole-history comparison will test numeric identity |
| Final gauge-validation Graf helper | Unverified exact source | A compatible `cml_graf_rstd_hess_radar_gauge_v3.py` exists, but its results were not shown to be identical to the publication helper |
| DWD 2021 supplemental radar replay | Recovered DWD script but original corrected run not proven | Released path-radar NetCDF is valid reference; do not claim the older script exactly reproduces corrected RADKLIM |
| WVD V6 end-to-end | Partial | Additional supporting files and steps are needed beyond the seven release products |
| Physically April–October-only public NetCDFs | Not yet verified | HPC 2026-10-13 test read *longer* source files containing all required Apr--Oct observations |
| Dataset release package | Needs final check | Earlier 5.8 GiB TAR included **wrong generic rainfall file**, not final ITU2005/Pastorek rain |

The original source tree and historical notebooks are preserved, including outdated models and hard-coded paths. **Do not run all archived notebooks as a pipeline.** The new CLI is the reviewer-safe path. All 44 originals are verified against their original manifests; this check is separate from validating published numerical values.

Data licences from the project: CC BY-NC 4.0 for CML-derived files; CC BY 4.0 for CHMI and DWD. A software licence must be selected with the rights holder before public code distribution; no software licence is implied here.
