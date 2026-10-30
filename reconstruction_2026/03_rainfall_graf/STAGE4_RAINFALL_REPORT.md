# Stage 4 report: rainfall from the corrected QC

**Date:** 2026-10-28.
**Status:** rainfall candidates produced. Independent verification is in §8. Awaiting author approval; validation
and figures (Stage 5) have not started.

**Untouched:** no original RAW, QC, rainfall, radar, gauge or WVD product, and no repository code, was modified.

**Settings:** `RAINFALL_PROCESSING_PARAMETERS.md`.

## 1. Inputs, code and execution

* **Input:** `stage3\QC_CORRECTED_CANDIDATE_REV1.nc` (sha256 ca8afa95…, verified at start). 1987 CMLs, 3731
  retained sublinks, 2021-04-01 00:00 – 2021-10-31 23:59.
* **Code:** `scripts\s4_01_rain_1min.py`.
  * Its scientific block (Graf RSTD, constant baseline, Pastorek, ITU 2005; 122 lines) is byte-identical to
    `essd-czech-cml-2021\src\essd_cz\export_1min.py`, and therefore to the verified MetaCentrum exporter
    (sha256 88013e76…). See `scripts\s4_00_check_verbatim_output.txt`.
  * This is the final configuration that produced the released 1-minute product.
* **Removed:** the 2219-CML requirement and the bit-for-bit comparison with the old product.
* **Added:** checks of the new input:
  * sha256;
  * exact time axis;
  * unique ids equal to the QC count;
  * finite frequency and length for retained sublinks;
  * Vertical polarization;
  * excluded sublinks entirely missing;
  * finite rate ⇔ valid flag.
* **Run:** 8.3 min locally, 0 retrieval errors on retained sublinks (`scripts\s4_01_rain_1min_output.txt`).
* **5-minute product:** derived with the repository's unchanged `rain5.derive_5min_dataset`
  (`scripts\s4_02_derive_5min.py`).

## 2. Outputs

| File | Size | Content |
|---|---:|---|
| `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc` | 445,051,457 B (sha256 9014567c3201be60441d1d8721c809d13ad5fa193e38a7b65606848b868a607e) | 308,160 min × 2 × 1987 |
| `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc` | 338,126,319 B | 61,632 × 2 × 1987; 126,329,587 finite directional 5-min rates |

**Variables in the 1-minute file:**
* `R_sublink_1min`, `R_mean_1min`, `graf_wet_1min`, `graf_valid_1min`, `graf_threshold_db`, `graf_rstd_q80_db`;
* new: `trsl_source_flag_1min` (bit 8 = blackout-reconstructed total loss, bit 16 = interpolated).

The 5-minute file adds `blackout_minutes_5min` and `interpolated_minutes_5min`.

## 3. Population

| | Previous product | Stage 4 |
|---|---:|---:|
| CMLs | 2219 | 1987 (1804 common, 415 previous-only, 183 new-only; see `stage3\STAGE3_QC_REPORT.md` §5) |
| CMLs with ≥ 1 valid rain minute | 1912 | **1987** (no empty CMLs) |
| Sublinks with ≥ 1 valid minute | 3796 of 4438 | **3731 of 3974** (all retained QC sublinks) |
| Valid directional rain minutes | 648,844,578 | 631,742,984 |

## 4. Graf thresholds and dry baseline (3410 common sublinks with data in both)

**Thresholds:**
* median change +0.0014 dB, median ratio 1.0034;
* 101 identical;
* 459 (13 %) change by more than 5 %;
* ratio quantiles: 1 % 0.92, 5 % 0.986, 95 % 1.088, 99 % 1.17.
* Causes: the calibration period is now April–October instead of January plus April–November, and the corrected
  QC changes the input (placeholder TX removed, plateau and blackout handling, interpolation only over raw gaps).

**Dry-baseline level** (median `trsl` at valid dry minutes, which equals the constant baseline there):
* median change 0.000 dB;
* 4 sublinks change by more than 1 dB (99 % quantile 0.05 dB).

Per-sublink values: `GRAF_THRESHOLD_COMPARISON.csv`.

## 5. Wet/dry and rainfall (common sublinks, minutes valid in both: 573,689,197)

| | Minutes |
|---|---:|
| Wet in both | 69,716,437 |
| Wet before, dry now | 1,468,860 |
| Dry before, wet now | 436,298 |
| Dry in both | 502,067,602 |
| **Agreement** | **99.67 %** |

* **Rain rate:** pooled correlation 0.957 between old and new 1-min rates. Mean rate 0.0560 → 0.0580 mm h⁻¹.
* **Accumulation on common minutes:** sum over sublinks 535,823 → 554,995 mm (+3.6 %). 277 sublinks are unchanged
  (no minute differs by more than 0.1 mm h⁻¹).
* **Per-CML accumulation (CML mean series, April–October):**
  * common CMLs: 207.0 → 208.6 mm mean, median per-CML ratio 1.013;
  * whole populations: 182.1 mm (2219, including 307 empty CMLs) → 214.7 mm (1987, none empty);
  * new-only CMLs: 275.0 mm (short links); previous-only CMLs: 73.8 mm (mostly empty or placeholder-based).
* **Availability per common sublink:** the median valid minutes are unchanged (5 % −13, 95 % +121, 99 % +401
  minutes).

**Largest sublink differences, and why** (`RAINFALL_BEFORE_AFTER.csv`):

| Sublink | Old → new (mm) | Cause |
|---|---|---|
| 47097_47918_4431 s0 | 328.1 → 8.3 | The old QC had 177,401 of 306,188 `trsl` minutes computed from placeholder TX; they are now missing |
| 32056_38088_3910 s1 | 370.6 → 153.7 | Same cause (146,174 placeholder-TX minutes) |
| 47097_47918_2217 s0 | 148.7 → 4.3 | Same CML pair (RX-only and placeholder-TX record) |
| 28021_28215_3781 s0 | 160.1 → 478.8 | Of the new total, 283.2 mm comes from blackout-reconstructed minutes |
| 31093_31133_3869 s0 | 744.3 → 1109.6 | Blackout-reconstructed minutes: 296.4 mm |
| 45057_45214_4357 s0/s1 | 293.5 → 578.3; 179.0 → 322.9 | Blackout reconstruction (49.3 / 33.3 mm) plus interpolation and threshold changes |

**What Agent 4 confirmed about these causes:**
* **Decreases:** confirmed. The lost rain sat almost entirely on old placeholder-TX minutes (313 of 320 mm, 189 of
  216 mm, 143 of 144 mm).
* **Increases:** only partly explained by blackout reconstruction (17–89 % of the top increases). The rest comes
  from the other QC changes and the threshold and baseline recalibration.

## 6. Contribution of flagged total loss

Over all 3731 new sublinks:

| Source | Valid minutes | Share of accumulation | Share of minutes ≥ 10 mm h⁻¹ | Share of minutes ≥ 30 mm h⁻¹ |
|---|---:|---:|---:|---:|
| Blackout-reconstructed (bit 8) | 36,716 (0.006 %) | **3.03 %** (19,640 of 648,700 mm) | **5.9 %** | **16.9 %** |
| Interpolated (bit 16) | 697,942 (0.11 %) | 0.25 % | 0.44 % | 0.47 % |

* **High intensity:** minutes ≥ 10 / ≥ 30 mm h⁻¹ rise from 445,216 / 52,702 to 525,748 / 82,837.
* **Maximum 1-minute rate:** 646.9 → 795.7 mm h⁻¹.
* Sublinks with a maximum above 200 mm h⁻¹: 70 → 104.

**Interpretation (revised after Agent 4's analysis, `agents\A4_stage4_verification.md`):**

* **The higher maximum rates are mainly a population effect.**
  * Of the 104 sublinks above 200 mm h⁻¹, 95 reach their maximum at a *measured* minute; 57 of those are at
    RX ≤ −85 dBm, near the receiver floor. Only 9 peak at a blackout minute.
  * The 9 highest maxima, including 795.7 mm h⁻¹, are on new-only CMLs: sub-kilometre 32–39 GHz links that now pass
    the per-sublink ITU-2005 detection limit.
  * Sublinks above 200 mm h⁻¹ have a median length of 0.76 km, against 2.80 km for all retained sublinks.
  * Rates of 400–800 mm h⁻¹ are not credible as path-averaged rain. This is a short-link, noise-floor limitation of
    the retrieval. Handling it would need a separate, approved rule; nothing has been changed.
* **Blackout reconstruction probably overestimates heavy rain** (Agent 4's diagnostics on 1216 gaps, [H]):
  * "max TX − min RX" is a median 19.7 dB above the 99.9th percentile of the sublink's own measured loss
    (10–90 %: 5–35 dB), because the maximum TX and minimum RX were never observed together;
  * the absolute edge rule RX < −65 dBm also triggers on links whose normal level is low: 45 % of gaps have an edge
    RX within 20 dB of the sublink median;
  * the jump to the reconstructed value itself makes the Graf window "wet";
  * in 47 % of rainy gaps the gap rate is more than twice the highest measured rate within ±30 minutes, and in
    33 % the surroundings are below 1 mm h⁻¹.
* **Size of the effect:** small on accumulation (3 %), large on intensity statistics (17 % of minutes ≥ 30 mm h⁻¹).
* **Interpolated minutes** have a negligible effect.
* **All such minutes are identifiable** (`trsl_source_flag_1min`, the 5-minute counts), so any of the options in §9
  can be applied later without losing information.

## 7. What to carry into Stage 5A

* The previous rainfall statistics (1,354,953 gauge-wet pairs, RMSE 3.33, bias −0.93, PCC 0.59) are references
  only. The population and the rates have changed.
* The gauge benchmark population changes:
  * the 13 TX = 0 CMLs and the 4 R1 CMLs leave;
  * new-only CMLs may enter if they lie within 5 km of a gauge.

## 8. Independent verification (Agent 4, `agents\A4_stage4_verification.md`, scripts `verification\v4_01`–`v4_06`)

| Check | Result |
|---|---|
| Scientific code identical to the repo exporter and the verified original | Confirmed: `graf`, `estimate`, `pol_hv`, `call_kr`, `decode`, `meta2` are identical, and the calls in `main()` are unchanged. `main()` and `make_native_file()` differ only as authorised (checks replaced, 5-minute reference part removed, source flag added) |
| Own Graf implementation + pycomlink on 37 sublinks (the named cases, random, short new-only, highest maxima, most blackouts) | Bit-identical: all 11.4 M rates; wet, valid and source flags identical; thresholds and q80 within 2e-7 relative; `R_mean_1min` identical |
| Whole-file product checks | Ids = QC REV1; correct time axis; excluded sublinks missing; `trsl_source_flag_1min` = `qc_flag & 24` at every minute; finite rate ⇔ valid; `R_mean` correct |
| 5-minute product | Whole file re-derived independently: all 244.9 M values identical. Totals 126,329,587 / 12,037 / 648,710 confirmed |
| Comparison numbers | All headline numbers recomputed exactly; 27 sampled common sublinks match every column of `RAINFALL_BEFORE_AFTER.csv` |

No mismatch anywhere. The open points are scientific (§6), not computational.

## 9. Decision requested

1. **Approve the Stage 4 rainfall candidates as computed**, with all reconstructed and interpolated minutes flagged,
   **or** request a controlled revision.
2. **Blackout reconstruction:** Agent 4 judges it likely to overestimate heavy rain (§6). Options, each a separate
   controlled revision if you choose one:
   * (a) keep it, and report the Stage 5 statistics both with and without bit-8 minutes;
   * (b) cap the reconstructed loss at the sublink's observed maximum total loss;
   * (c) make the edge criterion relative to each link's RSL distribution;
   * (d) return to flag-only (option A) for the rainfall product.
3. **Very high rates on short links** (up to 796 mm h⁻¹, mostly sub-kilometre new-only links at measured minutes):
   either accept them and document them, or define an approved rule. No rule exists in the current method.
4. Stage 5A (validation) needs your separate authorisation.

## 10. Author decision (2026-10-28)

**Option (a): "Keep, report both".**
* The Stage 4 rainfall candidates are approved as computed:
  * `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc`, sha256 9014567c…;
  * the derived 5-minute file.
* Blackout-reconstructed minutes stay in the product and stay flagged (`trsl_source_flag_1min` bit 8;
  `blackout_minutes_5min`).
* Stage 5A must report every validation statistic both **with and without** the bit-8 minutes. The without case
  excludes 1-minute values with bit 8, and 5-minute intervals with `blackout_minutes_5min` > 0.

**Documented limitation: very high rates on short links.**
* 104 sublinks reach more than 200 mm h⁻¹ at 1-minute resolution, with a maximum of 795.7 mm h⁻¹.
* 95 of them reach their maximum at a measured minute, 57 of those at RX ≤ −85 dBm. Their median length is
  0.76 km, and the 9 highest are sub-kilometre 32–39 GHz links newly admitted by the per-sublink ITU-2005 detection
  limit.
* Rates of this size are not credible as path-averaged rain. They come from the short-link and noise-floor
  behaviour of the unchanged retrieval.
* They are kept, not filtered: no approved rule exists. They are to be described in the data documentation and the
  manuscript.

Stage 5A (validation) is **not** authorised yet.
