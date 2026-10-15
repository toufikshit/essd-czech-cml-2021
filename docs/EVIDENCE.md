# Numerical evidence and where it comes from

Two independent sources of evidence, kept apart deliberately:

* **HPC run (October 2026)** — the author's real-data execution of the
  connected workflow on the original HPC. Small outputs are committed
  byte-identically under `docs/evidence/hpc_run_2026-10/` (see
  `EVIDENCE_MANIFEST.json`; large NetCDF/npz outputs live in the handoff
  archive, not in git).
* **Release-candidate local audit (October 2026)** — independent
  re-verification on a second machine (Windows, see `docs/ENVIRONMENT.md`)
  against the verified handoff data; reports under
  `docs/evidence/local_audit_2026-10/`.

## 1. Archive integrity (local audit)

* Both handoff archives match their published SHA256; all **257/257** files in
  the full handoff match `FILE_MANIFEST_SHA256.json`
  (`full_handoff_verification.json`).
* **44/44** preserved originals byte-identical (`tools_verify_originals.py`),
  re-verified after repository assembly and on a fresh `git clone`.

## 2. Final Graf/ITU2005/Pastorek rainfall — full-product identity

The author's HPC run recomputed **all 2,219 links over the full 2020–2021
history** with the unchanged original producer
(`original/code/cml_rainfall/cml_graf_itu2005_rebuild.py`; parameters in
`hpc_run_2026-10/01_rainfall_rebuild/rebuild_summary.json`, including 525
directional series skipped for <2 valid observations, matching real network
gaps). The local audit compared that recomputed file against the published
rainfall product element by element (`rainfall_identity_report.json`,
`rainfall_static_vars_report.json`):

| check | result |
|---|---|
| `R_sublink_5min` (137,782,581 shared finite points) | **identical, zero difference, identical masks** |
| `R_mean_5min` (84,656,385) | identical |
| `graf_wet_5min`, `graf_valid_5min`, `wet_minutes_5min`, `valid_minutes_5min` (up to 934,323,264 points each) | identical |
| per-link calibration `graf_threshold_db`, `graf_rstd_q80_db` + all geometry/metadata variables | identical |

This supersedes (and contains) the earlier single-link gate: link
`10004_311866_3129`, 7,666 shared directional samples, RMSE 0 (HPC,
`docs/status/TEST_STATUS.md`). The local audit additionally **re-executed**
the original science path end-to-end on a second platform via
`run.py compare-rainfall` (numba + pycomlink, full-history baseline, link
index 3 = `10004_311866_3129`): 7,666 shared directional samples, RMSE 0.0,
zero mask mismatches — reproducing the HPC gate exactly
(`one_link_comparison_report.json`). **Caveat:** identity holds when recomputing
from the full 2020–2021 QC history; the public April–October release alone
cannot re-derive the Q80/baseline calibration (the per-link thresholds are
shipped inside the rainfall product).

## 3. Entity counts and coverage (local audit, real data)

`run.py check` against the seven full-history sources: raw 2,580 CMLs ·
QC 2,219 CMLs · rain 2,219 · path radar 1,132 · gauges 424 · temperature 251 ·
humidity 247 stations; complete native April–October axes in all seven
(308,160 / 61,632 / 30,816 slots at 1/5/10 min).

## 4. Water vapour (WVD V6)

Local audit of the supporting products (`wvd_support_audit.json`):

* V6 NetCDF: **1,352 directional series** × 61,632 five-min steps
  (2021-04-01..2021-10-31), GRAF + radar method branches.
* Diagnostics CSV and quantization CSV: 1,352 rows each, all `series_index`
  values valid against the NetCDF.
* HPC-executed station-matched validation (≤5 km): 587 directional links,
  2,834,230 strict pairs (`hpc_run_2026-10/WVD_V6_QUANTIZATION_BINNED_VALIDATION/run_summary.json`).
* HPC-executed paper population (E-band, L ≥ 4 km, Δρ_q ≤ 0.20 g m⁻³, GRAF):
  **45 directional sublinks from 31 physical CMLs**, N = 236,793, PCC 0.598
  (`hpc_run_2026-10/07_wvd_paper_population/wvp_support_metrics.csv`).

These counts were measured, not forced; notebook 06 refuses mismatched
populations rather than filtering silently.

## 5. Gauge validation (HPC; detector caveat applies)

`hpc_run_2026-10/02_gauge_validation/.../FINAL_GRAF_validation_summary.json`:
shared-population validation of the recomputed rainfall against CHMI gauges
(Apr–Oct 2021, 1,020 contributing physical CMLs, 218 gauges, N = 34.6 M):
all-period RMSE 0.70 mm h⁻¹ / PCC 0.60; wet (≥0.1 mm) RMSE 3.33 / PCC 0.59;
confusion POD 0.758, FAR 0.796, MCC 0.350. **These numbers come from the
documented substitute pipeline — the exact historical publication detector is
not recovered — and must not be labelled identical to the paper's gauge
validation without numerical confirmation.**

## 6. Availability

HPC notebook 00 walkthrough (committed CSV/figures under
`hpc_run_2026-10/00_exploration/`): Apr–Oct mean availability — CML 58.67%,
path radar 99.48%, gauges 99.74%, temperature 99.57%, humidity 99.56%;
4,438 directional sublinks, 1,804 E-band directional sublinks. The local
audit re-ran the unchanged availability stage (`run.py availability`) against
the real full-history sources and reproduced every mean to the last digit:
CML 58.67% (median 59.53%, max 73.73%), radar 99.48%, gauge 99.74%,
temperature 99.57%, humidity 99.56%
(`local_audit_2026-10/availability_native_timestep_percentage.png`).

## 7. Synthetic regression suite (local)

20/20 pytest cases and the synthetic seven-file connected chain
(`tests/run_connected_mini.py`) pass on the release candidate; see
`docs/status/` for the HPC-era logs.
