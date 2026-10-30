# Stage 3 proposal: QC settings and population handling (for approval; nothing run)

**Date:** 2026-10-28.
**Status:** proposal only. Stage 3 execution needs your separate approval. No QC was run.

**Basis:**
* Stage 1A `stage1\QC_SCRIPT_PROVENANCE.md`: the producer as executed vs the manuscript vs standalone v3;
* the manuscript QC text (l. 307–311);
* standalone v3 (`essd-czech-cml-2021\original\code\cml_qc\related_versions\cml_qc_filter_standalone_v3.py`,
  sha256 1add0b5c…);
* the rainfall retrieval `src\essd_cz\export_1min.py`, for coefficient consistency;
* your Stage 1B and Stage 2 decisions.

## 1. Population handling (fixed by your decisions)

| Step | CMLs | Sublinks | Record |
|---|---:|---:|---|
| Approved RAW input | 2574 | 5148 | `stage2\RAW_APPROVED_POPULATION.csv` |
| Pre-QC: 8 CMLs without frequency, not rainfall-ready (P3) | −8 | −16 | reason "no frequency metadata" |
| Pre-QC: TX = 0 directions non-operational (P5); opposite direction kept | – | 15 marked invalid | `stage2\QC_PREDEFINED_SUBLINK_EXCLUSIONS.csv` |
| Pre-QC: mixed-placeholder minutes, TX and RX masked (P6) | – | 15 minutes in 13 sublinks | `stage2\QC_PREDEFINED_MINUTE_MASKS.csv` |
| Enters the QC rules | **2566** | 5132 (15 of them already invalid) | – |
| QC result | computed | computed | not forced to any count |

About the P6 count: 16 values were identified in Stage 2. One of them, 41920_323338_4235, left with R1, so 15
remain in the approved population.

**Never in QC:**
* the 13 TX = 0 CMLs;
* the 6 R1 CMLs (4 of them are in the released QC and will disappear);
* the 255 empty CMLs.

These are not counted again as QC exclusions. By construction QC_CML_IDs ⊆ RAW_CML_IDs, and this is checked
after the run.

## 2. Proposed QC settings (1-minute data; one row per item of the Stage 3 plan)

| Item | Proposed setting | Source / reason | Decision? |
|---|---|---|---|
| Time interval | **2021-04-01 00:00 – 2021-10-31 23:59** (the release window) | The QC output covers this window, and every decision becomes reproducible from the published RAW. Outside it there are only 2021-01-13..19 and 11-01..05 (about 10.3 days). The originally executed QC used the whole file | **Yes** (A: window, recommended; B: all available data 2021-01-13..11-05) |
| Invalid TX | TX ≤ −99 dBm or = 255 → missing | Manuscript; v3. None of these values remains in RAW, but the rule is kept and logged | no |
| Invalid RX | RX ≤ −99.9 dBm or > −20 dBm → missing | Manuscript; v3. None remains in RAW | no |
| TX = 0 dBm | Only the 15 approved P5 directions are invalid. TX = 0 elsewhere (part of the normal ATPC range on more than 1500 sublinks) stays valid | Your P5 decision; Stage 1B evidence | no |
| Implausible TX | Only the 15 approved P6 minutes, with TX and RX masked. No other limit (TX < −30 dBm is not masked) | Your P6 decision ("no additional speculative masks") | no |
| Low-RSL plateau | Centred 3-minute window: max RSL < −85 dBm and population std (ddof = 0) < 0.5 dB, **with at least 2 valid samples**; mask TX and RX with ±2 min padding | v3 thresholds. v3 accepts 1 valid sample, which flags any isolated reading below −85 as a "plateau"; the executed producer required 2 | **Yes** (≥ 2 recommended, or v3's ≥ 1) |
| Blackout gaps | Gap definition (v3): a whole RX gap of ≤ 60 min with RSL < −65 dBm on at least one edge. Proposed treatment: **flag only, no fill**; values stay missing, and the number of flagged minutes is reported | v3 fills such gaps with the series maximum TX and minimum RX, which creates synthetic observations and conflicts with the safeguard "never generate synthetic TX, RX". The manuscript says "treated" and cites Blettner (2023), whose method fills them | **Yes**, see §4 (A: flag only; B: v3 fill of TX/RX, flagged; C: Agent 4's option, TX/RX missing but the reconstructed total loss written to `trsl` with its own flag) |
| Short gaps | Linear interpolation of total loss for gaps ≤ 5 min with valid values on both sides; every interpolated sample flagged. **Only gaps that are missing in RAW are bridged.** Samples masked by QC (P6, plateau) and blackout-flagged gaps are never interpolated. At interpolated samples `trsl` ≠ `tx − rx`; this is stated in the attributes | Manuscript (l. 307) and v3. v3 would also re-fill plateau-masked gaps and the single P6 minutes (Agent 4) | confirm |
| Fluctuation screening | Per sublink, on total loss after interpolation: centred rolling windows of 300 min (5 h) and 60 min (1 h), population std. A sublink **fails** if more than 10 % of valid 5-h windows have std > 2 dB, or more than 33 % of valid 1-h windows have std > 0.8 dB | Manuscript and v3 (pass if ≤). This corrects the executed producer: windows of 21/5 samples, ≥, sublink 0 only | no |
| Valid window | ≥ 60 finite minutes in the window (so the 1-h window must be complete). A sublink with **no** valid window fails (it cannot be assessed) | v3 (`ROLLING_MIN_PERIODS = 60`). The manuscript does not define "valid" | confirm |
| Detection limit | **Per sublink:** `R_min = (0.33 dB / (k·L))^(1/α)`. k and α from ITU-R P.838-3 (2005), vertical polarization, the sublink's own frequency, metadata length; exactly the same `pycomlink` k–R implementation (`ITU_2005`) as the rainfall retrieval. Sublink fails if R_min ≥ 2 mm h⁻¹ | This is the rain rate whose path attenuation k·R^α·L equals 0.33 dB, i.e. a true minimum detectable rain rate (manuscript l. 309). v3's `0.33/(k·L)` has no α and is not a rain rate. The executed producer used the 2003 table and the CML minimum frequency | **Yes** (coefficient source: 2005, recommended for consistency with the retrieval; or 2003 as in the stored values) |
| Sublink exclusion | A sublink is excluded if it is a P5 direction, has no valid total loss in the window, fails the fluctuation test, or fails the detection limit. **Excluded sublinks are set to missing** in the QC product; their flags and statistics are kept | Manuscript ("sublinks … were excluded"). v3 only flagged them | no |
| CML retention | A CML is retained if at least one sublink passes. It stays one CML (sublink dimension 2, the failed direction missing) | Manuscript l. 309 | no |
| CMLs with no passing sublink | Not written to the QC NetCDF; listed with reasons in `QC_EXCLUDED_CMLS.csv` | Keeps the QC population meaningful (no empty CMLs, unlike the old product) | confirm |
| Polarization | Vertical for every CML (author confirmation) | – | no |

## 3. Outputs (Stage 3 plan) and checks

Everything is written to `stage3\`; no original file is changed.

**Data:**
* `QC_CORRECTED_CANDIDATE.nc`, with:
  * `tx`, `rx` (screened);
  * `trsl` (recomputed);
  * a per-sample `qc_flag` bitmask (invalid value, P6 mask, plateau, blackout flag, interpolated);
  * per-sublink `qc_detection_limit_mmhr`, the fluctuation fractions, and a `qc_sublink_status` with its reason;
  * the CF attributes, the exact parameters in JSON, and the input sha256.

**Tables and reports:**
* `QC_EXCLUDED_CMLS.csv` and `QC_EXCLUDED_SUBLINKS.csv`, one reason per row;
* `QC_BEFORE_AFTER_COMPARISON.csv` against the released 2219-CML QC;
* `QC_PROCESSING_PARAMETERS.md`;
* `STAGE3_QC_REPORT.md`.

**Mandatory checks:**
* QC ⊆ RAW, with no new id;
* the population table reconciles exactly (raw input = excluded + retained, for CMLs and for sublinks);
* placeholders are treated consistently;
* dimensions and metadata are valid.

The old 2219 is compared, not targeted.

**Execution:**
* the v3 kernel is copied to `stage3\scripts\` and changed only as listed above;
* reads are chunked by CML batch;
* it runs locally, or on MetaCentrum with the same command;
* Agent 1 runs it and Agent 4 independently recomputes a sample of sublink decisions and every population count.

**Not affected:** the WVD product (it reads RAW, not QC); the radar and gauge files.

**Downstream (later stages):** the new QC changes the rainfall (Stage 4), the gauge benchmark and the figures
(Stage 5). The 13 and the 4 R1 CMLs leave the benchmark population.

## 3b. Independent review (Agent 4, `agents\A4_stage3_settings_review.md`; code and documents only)

**Confirmed against the v3 code:**
* the plateau test accepts 1 valid sample;
* blackout filling uses the series maximum TX and minimum RX;
* both windows need 60 valid minutes;
* pass if ≤;
* a sublink without a valid window fails;
* v3's detection limit has no α.

**Confirmed against the Stage 2 records:** all the counts.

**Detection limit.** The retrieval uses pycomlink 0.6.0 `calc_R_from_A(..., a_b_approximation='ITU_2005')`.
`a_b(f, 'V', 'ITU_2005')` gives k and α from the ITU-R P.838-3 table, valid for 1–100 GHz; RAW frequencies
(7.45–85.75 GHz) are all inside. `R_min = (0.33/(k·L))^(1/α)` is in mm h⁻¹: the round trip returns exactly 0.33 dB.
It is a quantisation-based lower bound, not the retrieval's actual detection threshold.

**Amendments taken into the proposal:**
* about 10.3 days outside the window;
* the plateau std is the population std;
* the blackout gap definition is fixed;
* interpolation never bridges QC-masked samples or blackout gaps;
* `trsl` ≠ `tx − rx` at interpolated samples is documented.

**Outcomes the data already fix under the release window** (stated so that nobody mistakes them for QC effects):
* 57 of the 2566 CMLs have no paired minute in April–October, so they necessarily fail;
* 26 sublinks enter QC with no possible valid data: the 15 P5 sublinks plus 11 empty sublinks;
* 1 of the 15 P6 masks (26037_26137_3704, 2021-11-04) lies outside the window, so 14 apply.

**Consequences to document later (Stage 6):**
* with ITU 2005 per sublink, every `detection_limit_mmhr` value changes, and the attribute "ITU-R 2003" and the
  manuscript wording change;
* the blackout wording changes if the gaps are not filled.

**Dissent on blackouts.** Agent 4 disagrees with "flag only". Blackouts happen at the highest attenuation, so
leaving them empty biases heavy rain low. It proposes option C. Option A creates no synthetic value but loses those
minutes; options B and C reproduce the cited method but create synthetic values (B in TX/RX, C only in `trsl`). This
is a scientific choice for you.

## 4. Decisions requested

1. **Time interval:** release window (recommended) or all available data.
2. **Plateau:** at least 2 valid samples (recommended) or 1 (v3).
3. **Blackout:**
   * A: flag only, no synthetic values (lead recommendation);
   * B: v3 fill of TX and RX, flagged;
   * C: total loss reconstructed in `trsl` only, flagged (Agent 4's recommendation).
4. **Detection-limit coefficients:** ITU-R P.838-3 2005 per sublink, with the same pycomlink function as the
   retrieval (recommended), or the 2003 table. This changes the released detection-limit values.
5. **Confirm:**
   * the valid-window definition (≥ 60 finite minutes);
   * dropping CMLs without a passing sublink from the QC file;
   * the ≤ 5-min interpolation of RAW gaps only (never across QC masks or blackout gaps).

Then authorise the Stage 3 QC run.
