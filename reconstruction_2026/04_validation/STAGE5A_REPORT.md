# Stage 5A report: validation of the reconstructed rainfall

**Date:** 2026-10-28.
**Status:** validation complete. Independent verification is in §7. Awaiting author review; Stage 5B (figures) has
not started.

**Untouched:** no RAW, QC, rainfall, radar, gauge or WVD product was modified, and no new filter was introduced.

## 1. Method and code

* **Validation code:** `scripts\s5a_01_gauge_validation.py`.
  * It is a copy of `gauge_validation()` and `compute_stats()` of the published producer
    `essd-czech-cml-2021\original\code\figures\rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py`.
  * It uses the repository helper `cml_graf_rstd_hess_radar_gauge_v3.py`; the 11 functions used are identical in
    every local version.
* **Changes:** no plots; every joint sample saved; optional marking of blackout-reconstructed samples.
* **Method (unchanged):**
  * nearest CHMI gauge to the CML midpoint within 5 km;
  * gauge `rr` (10 min, end-labelled) × 6;
  * CML rate = mean of the two complete 5-minute directional bins at t−10 and t−5;
  * Graf RSTD detection recomputed on the CML file's 1-minute `trsl` (60-minute window, 0.8 quantile × 1.12,
    60/45 minimum periods), "any" over the ten minutes of the gauge interval, all ten valid;
  * one joint mask: rain valid, detection valid, both 0–150 mm h⁻¹;
  * gauge wet ≥ 0.1 mm h⁻¹;
  * per directional sublink.
* **Blackout marking (bit 8):** a sample is marked if any of its ten 1-minute QC minutes has `qc_flag` bit 8.
  "Without bit 8" removes these samples from every statistic, both rainfall and detection.
* **Comparison:** `scripts\s5a_02_compare.py` matches samples on identical (CML, sublink, gauge, timestamp).

## 2. Control test: reproduction of the published validation

| | Published | Control (released inputs) |
|---|---:|---:|
| Gauge-wet N | 1,354,953 | 1,354,953 |
| RMSE | 3.33 (3.32654) | 3.32654 |
| MBE | −0.93 (−0.92628) | −0.92628 |
| PCC | 0.59 (0.58541) | 0.58541 |
| All-period N / RMSE / MBE / PCC | 34,565,739 / 0.7016 / −0.0191 / 0.6003 | identical |
| TP / FP / FN / TN | 1,026,880 / 4,014,661 / 328,073 / 29,196,125 | identical |
| Contributing CMLs / sublinks / gauges | 1020 / 2020 / 218 | identical |

**Exact reproduction:** every count is identical, and the PCC agrees to 1e-14. The inputs are the full-history QC,
the full-history 5-minute rainfall and the full-history gauge file (`control\control_released_summary.json`).

## 3. Validation of the Stage 4 rainfall (same method)

**Inputs:**
* `QC_CORRECTED_CANDIDATE_REV1.nc` for detection;
* the Stage 4 derived 5-minute rainfall;
* the same gauge file.

**Population:** 1074 gauge-matched CMLs (released: 1144, which included empty QC CMLs). Contributing: **1074 CMLs,
2015 sublinks, 217 gauges** (released: 1020 / 2020 / 218).

| Gauge-wet (≥ 0.1 mm h⁻¹) | N | RMSE | MBE | PCC |
|---|---:|---:|---:|---:|
| Published / control | 1,354,953 | 3.33 | −0.93 | 0.585 |
| **Stage 4, with blackout-reconstructed samples** | **1,331,687** | **3.63** | **−0.91** | **0.633** |
| **Stage 4, without blackout-reconstructed samples** | **1,327,998** | **3.38** | **−0.89** | **0.606** |

| All period | N | RMSE | MBE | PCC |
|---|---:|---:|---:|---:|
| Control | 34,565,739 | 0.702 | −0.019 | 0.600 |
| Stage 4 with bit 8 | 34,249,231 | 0.768 | −0.019 | 0.634 |
| Stage 4 without bit 8 | 34,244,692 | 0.714 | −0.019 | 0.615 |

| Detection | TP | FP | FN | TN | POD | FAR | CSI | MCC | Correct negatives |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Control | 1,026,880 | 4,014,661 | 328,073 | 29,196,125 | 0.758 | 0.796 | 0.191 | 0.350 | 87.9 % |
| Stage 4 with bit 8 | 993,974 | 3,757,603 | 337,713 | 29,159,941 | 0.746 | 0.791 | 0.195 | 0.354 | 88.6 % |
| Stage 4 without bit 8 | 990,285 | 3,756,754 | 337,713 | 29,159,940 | 0.746 | 0.791 | 0.195 | 0.353 | 88.6 % |

**Bit-8 samples:** 4,539 samples (3,689 gauge-wet) contain a blackout-reconstructed minute.

Full tables: `RAINFALL_VALIDATION_METRICS.csv`, `RAIN_DETECTION_METRICS.csv`, `VALIDATION_BEFORE_AFTER.csv`.

## 4. Population change vs processing change

Samples: released 34,565,739; Stage 4 34,249,231; **identical in both 30,750,879**.

| Group (gauge-wet) | N | RMSE | MBE | PCC | Reading |
|---|---:|---:|---:|---:|---|
| Common samples, released | 1,194,055 | 3.357 | −0.946 | 0.594 | Same CMLs, gauges and timestamps |
| Common samples, Stage 4 | 1,194,055 | **3.296** | −0.912 | **0.612** | **Processing change alone improves the agreement slightly** |
| Common, without bit 8: released / Stage 4 | 1,193,662 | 3.338 / 3.279 | −0.946 / −0.913 | 0.592 / 0.609 | Same conclusion |
| Released-only, CML removed | 86,714 | 3.301 | −0.975 | 0.491 | Population leaving: 13 TX = 0 CMLs, R1, fluctuation, no data |
| Released-only, CML kept, sample lost | 74,184 | 2.832 | −0.546 | 0.495 | Mainly old placeholder-TX samples now missing |
| Stage-4-only, new CML | 132,817 | 4.064 | −0.686 | 0.617 | The 183 new (mostly short) CMLs: noisier |
| Stage-4-only, CML kept, sample gained | 4,815 | **22.37** | −6.27 | 0.556 | Mostly blackout-reconstructed samples (3,620 of the 15,315 gained samples carry bit 8) |
| …the same group without bit 8 | 1,550 | 11.97 | −4.20 | 0.578 | – |

**Interpretation:**
* On identical samples the corrected processing performs slightly better than the released product (RMSE −0.06,
  PCC +0.017).
* The higher overall RMSE of the Stage 4 validation (3.63) comes from:
  * a few thousand gained samples whose CML rate is a blackout reconstruction (RMSE 22);
  * the newly admitted short links (RMSE 4.06).
* Without the blackout samples the RMSE is 3.38, close to the published 3.33.
* The PCC rises in every Stage 4 variant (0.606–0.633 vs 0.585).
* **What the blackout samples are (Agent 4, gauge evidence)** — this corrects the Stage 4 "overestimation"
  hypothesis:
  * The gauge-wet bit-8 samples are mostly real heavy rain. The gauge median is 21 mm h⁻¹ (mean 29.2) against a
    CML median of 17.2 (mean 22.9).
  * The CML exceeds the gauge in only 44 % of them, and only 2.7 % of the gained ones have a gauge below 1 mm h⁻¹.
  * The reconstruction captures heavy-rain intervals and on average underestimates them, with very large scatter
    (|CML − gauge| > 20 mm h⁻¹ in 31 %).
  * The RMSE rise is therefore mainly a change of sample composition: heavy-rain intervals that blackouts made
    incomplete before now enter the validation. It is not a degradation of the processing.
  * However, 850 bit-8 samples (19 %) are gauge-dry, with a CML mean of 10.9 mm h⁻¹ (199 at ≥ 10 mm h⁻¹). These
    are spurious events.
* **Robustness:** the gauge-wet RMSE is driven by a few thousand heavy-rain samples and moves by about ±0.3 with
  small population changes. The PCC and the identical-sample comparison are the more robust indicators.

## 5. Extreme estimates (description only; no filter)

| CML 10-min rate | Released: samples / share of squared error | Stage 4: samples / share of squared error | Stage 4: of which bit 8 | Stage 4: of which on links < 1 km |
|---|---|---|---:|---:|
| ≥ 30 mm h⁻¹ | 2,295 / 12.4 % | 3,976 / 21.2 % | 985 | 1,397 |
| ≥ 50 mm h⁻¹ | 581 / 6.7 % | 1,130 / 12.0 % | 346 | 474 |
| ≥ 100 mm h⁻¹ | 45 / 1.5 % | 101 / 3.1 % | 29 | 68 |

* Samples with CML ≥ 30 mm h⁻¹ almost double. They carry a fifth of the squared error.
* About a quarter of them come from blackout reconstruction, and about a third from links shorter than 1 km.
* By method, the validation mask excludes CML or gauge values above 150 mm h⁻¹ in both products.
* Table: `EXTREMES_IN_VALIDATION.csv`.

## 6. WVD

The WVD inputs are unchanged (RAW only; no excluded CML; RAW values unchanged), so the verified WVD results are
retained without a rerun.

One artefact minute reaches only the gauge-method WVD series. It is not among the 236,793 graf-method pairs of
Table S13 and Fig. 8, so those are unaffected. See `WVD_VERIFICATION_REPORT.md`.

## 7. Independent verification (Agent 4, `agents\A4_stage5a_verification.md`, scripts `verification\v5_01`–`v5_04`)

| Check | Result |
|---|---|
| `s5a_01` vs the original `gauge_validation()` | Faithful: shared functions identical; MAX_R, period and helper unchanged. Differences (no plots, per-sublink confusion arithmetic, direct statistics instead of the equality assertions, saved samples, bit-8 marking, simplified `select_rain_rate`) cannot change any metric |
| Headline metrics from the saved samples | Control, Stage 4 with and without bit 8, all-period, gauge-wet and detection: identical |
| Own code from the source files, 30 CMLs (common, new short, most bit-8, released-only), 102 sublink comparisons | Identical sample sets; gauge value, CML rate, wet prediction and bit-8 mark bit-identical (879,089 control and 822,069 Stage 4 samples) |
| Bit-8 marking | Definition confirmed (the ten QC minutes = the two 5-minute bins). 4,539 samples (3,689 gauge-wet): 3,620 gained in kept CMLs, 881 common, 38 in new CMLs |
| Identical-sample comparison and population groups | 30,750,879 common samples with identical gauge values. Group sizes confirmed. Gauge-wet sum of squared errors: common samples 13.45 M → 12.97 M (processing improves); new CMLs +2.19 M; gained samples in kept CMLs +2.41 M |
| WVD | Inputs unaffected (confirmed). The artefact appears only in the gauge-method series (see §6) |

No mismatch.

## 8. For the author

* Previous manuscript numbers (1,354,953 / 3.33 / −0.93 / 0.59; POD 75.8 %, correct negatives 87.9 %) are
  superseded. The new values above, with and without blackout samples, are the candidates for Stages 5B and 6.
* Which variant should the manuscript headline: with blackout samples (the product as released) or without them
  (the measured signal only)? Both will be reported as you decided; the choice of headline is yours.
* Stage 5B (figures) needs your separate approval.

## 9. Author decision (2026-10-28)

**Headline: "With blackout samples".**
* The manuscript headlines the Stage 4 product as released: gauge-wet N 1,331,687, RMSE 3.63 mm h⁻¹,
  MBE −0.91 mm h⁻¹, PCC 0.633.
* The values without blackout-reconstructed samples (N 1,327,998, RMSE 3.38, MBE −0.89, PCC 0.606) are reported
  alongside.
* Detection with blackout samples: POD 74.6 %, correct negatives 88.6 %, FAR 0.791, MCC 0.354.

Stage 5A is closed. Stage 5B (figures) is **not** authorised yet.
