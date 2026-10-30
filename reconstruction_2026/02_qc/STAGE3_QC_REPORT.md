# Stage 3 report: corrected QC candidate

**Date:** 2026-10-28.
**Status:** QC candidate produced and verified. Awaiting author review. Stage 4 has not started.

**Settings:** approved by the author on 2026-10-28; see `QC_PROCESSING_PARAMETERS.md`.

**What was and was not touched:** no existing RAW, QC, rainfall or WVD product was modified. No downstream stage
was run.

## 1. Inputs and execution

* **Input:** the full-history RAW `CZ_CML_RAW_1MIN_20200101_20211231.nc` (sha256 6acfeed2…), with the approved
  population and masks from `stage2\` (`RAW_APPROVED_POPULATION.csv`, `QC_PREDEFINED_SUBLINK_EXCLUSIONS.csv`,
  `QC_PREDEFINED_MINUTE_MASKS.csv`).
* **QC script:** `scripts\s3_01_qc_corrected.py` (sha256 f1b32b1f…, stored in the NetCDF). It ran locally in
  28 min, with log `scripts\s3_01_qc_corrected_output.txt`.
  * Pass 1 decided all 2566 sublink pairs.
  * Pass 2 recomputed the retained CMLs and asserted identical decisions before writing.
* **Tables and checks:** `scripts\s3_02_tables_checks.py`, with output `scripts\s3_02_tables_checks_output.txt`.
* **File verification:** `verification\s3_03_file_verification.py`, with output alongside it.
* **Test run:** an 11-CML subset (`test\`), checked minute by minute against RAW before the full run.

## 2. The QC candidate file

| Property | Value |
|---|---|
| File | `QC_CORRECTED_CANDIDATE.nc` (NetCDF-4, zlib level 4, chunks 10080 × 2 × 64) |
| Size | 2,409,367,288 bytes |
| SHA-256 | `e193b2cf2e7dac3fc12c1e08a12c79a793526679f5e91052df93c735786f7f8d` |
| Dimensions | time 308,160 (2021-04-01 00:00 – 2021-10-31 23:59, 1-min steps); sublink_id 2; cml_id 1987 |
| Variables | `tx`, `rx`, `trsl`, `qc_flag`, `qc_sublink_status`, `qc_detection_limit_mmhr`, the two fluctuation fractions, `length`, `frequency`, `polarization` (all Vertical), endpoint coordinates |
| Global attributes | period, population, `qc_parameters_json` (24 parameters), input and script SHA-256 |

**Readability:** samples at the start, middle and end (3 CMLs each) read correctly. `trsl` = `tx − rx` exactly at
unflagged minutes. Every CML has at least one retained sublink.

## 3. Population

| Population | CMLs | Sublinks |
|---|---:|---:|
| Approved RAW input | 2574 | 5148 |
| Excluded before QC: no frequency (P3) | 8 | 16 |
| Assessed by QC | 2566 | 5132 |
| Excluded by QC (no passing sublink) | 579 | 1158 |
| Excluded sublinks inside retained CMLs | – | 243 |
| **Retained after QC** | **1987** | **3731** |
| CMLs in QC but not in RAW | **0** | – |

**Consistency checks (all pass):**
* QC is a strict subset of the approved RAW;
* ids are unique;
* none of the 13 TX = 0 CMLs, the 6 R1 CMLs, the 255 empty CMLs or the 8 no-frequency CMLs is in QC;
* NetCDF ids and sublink status equal the results table;
* both reconciliations are exact: 2574 = 8 + 579 + 1987, and 5148 = 16 + 1158 + 243 + 3731.

**Sublink decisions (5132 assessed):**

| Decision | Sublinks |
|---|---:|
| Retained | 3731 |
| Detection limit ≥ 2 mm h⁻¹ | 852 |
| Fluctuation | 364 |
| No observed total loss in the period | 165 |
| P5 (TX = 0 direction) | 15 |
| Fluctuation and detection limit | 5 |

**Reasons for excluding whole CMLs (`QC_EXCLUDED_CMLS.csv`):**

| Reason | CMLs |
|---|---:|
| Detection limit on every sublink with data | 414 |
| Fluctuation on every sublink with data | 99 |
| No usable observations in the period | 61 |
| No frequency (before QC) | 8 |
| Detection limit and fluctuation | 3 |
| P5 direction plus no data on the other sublink | 2 |

Per-sublink reasons: `QC_EXCLUDED_SUBLINKS.csv` (1401 rows).

**E-band:** 810 retained CMLs have at least one sublink at ≥ 70 GHz.

## 4. Processing counts (retained sublinks | all assessed sublinks)

| Item | Retained | All assessed |
|---|---:|---:|
| Mixed-placeholder minutes masked (P6) | 9 | 14 (the 15th lies outside the period) |
| Invalid TX / RX values | 0 / 0 | 0 / 0 |
| Plateau-masked minutes | 1,262,668 | 2,830,539 |
| Blackout gaps with reconstructed total loss | 8,194 gaps; 57,351 min; in 1520 sublinks | 14,192 gaps; 116,561 min |
| Interpolated short gaps | 651,631 gaps; 706,377 min | 892,374 gaps; 972,894 min |
| Observed total-loss minutes after masking | 631,446,096 | 858,961,242 |

* Reconstructed blackout minutes make up 0.009 % of the retained total-loss minutes; interpolated minutes make up
  0.11 %.
* TX and RX stay missing at every reconstructed minute (checked).

## 5. Difference from the released 2219-CML QC

| Group | CMLs | Explanation |
|---|---:|---|
| In both | **1804** | Retained by both. Inside them, 198 of 3608 sublinks are now excluded (169 fluctuation, 24 no data in the period, 5 P5) and set to missing |
| **Released only** | **415** | See below |
| ↳ Empty metadata-only CMLs | 255 | Not in RAW, and no observations at the source (Stage 1B) |
| ↳ The 13 TX = 0 CMLs | 13 | Your dataset exclusion (Stage 1B) |
| ↳ R1 non-operational | 4 | 10004_11582_3122, 27040_27099_3747, 10662_10795_3342 and 29002_29004_3810 (Stage 2) |
| ↳ Fluctuation | 94 | All sublinks with data now fail. In 177 of the 183 failing sublinks it is the 5-h test (more than 10 % of valid windows with std > 2 dB; median 0.17); 6 fail both tests. The old producer passed them because it used 21- and 5-minute windows and tested only sublink 0 |
| ↳ No usable observations in the period | 48 | 45 have no paired observation in April–October (data only in January or November, which the old QC kept because the producer let empty series pass). 3 have only low-RSL plateau data (receiver floor). For example, 10685_304691_3344 has 31,316 minutes, all at a plateau below −85 dBm |
| ↳ P5 and no data | 1 | 10082_10484_3157: its only direction with data is a TX = 0 direction |
| **New only** | **183** | All were RAW-only links (the Stage 1 group of 629). The old QC excluded them through its detection limit (CML minimum frequency, 2003 table): their old-formula R_min was 2.15–2.82 mm h⁻¹ (median 2.36; one at 0.61, excluded for another reason). With the approved per-sublink ITU-2005 formula, at least one sublink is below 2 mm h⁻¹ and passes the other tests |
| In neither (approved RAW) | 446 | 414 detection limit, 5 fluctuation, 3 detection limit plus fluctuation, 13 no data, 1 P5 plus no data, 8 no frequency, 2 R1 |

**Net change:** 2219 − 415 + 183 = **1987**.

The main drivers:
1. The removal of populations that never belonged in QC: 272 CMLs (255 empty, 13 by your decision, and the 4 R1
   CMLs that are in the released QC).
2. The corrected fluctuation test at the intended 5-h/1-h scale, applied per sublink: −94 CMLs.
3. Links without data in the release window: −48.
4. The physically consistent per-sublink detection limit: +183 short links regained.

The detection limits changed for every sublink. For CMLs in both products, the smallest sublink value with
ITU 2005 is a median 0.59 times the old CML value (2003 table, minimum frequency).

## 6. Independent verification (Agent 4)

Report: `agents\A4_stage3_verification.md`; scripts and outputs: `verification\v3_01`–`v3_04`. Agent 4 wrote its
own implementation from `QC_PROCESSING_PARAMETERS.md` and your decisions before reading the lead's script.

| Check | Result |
|---|---|
| Own implementation on 66 CMLs / 132 sublinks: random ones; P5 directions; all 12 CMLs holding the 14 in-window P6 minutes; most blackout minutes; fluctuation fractions nearest 0.10 and 0.33; detection limits nearest 2 mm h⁻¹; CMLs without window data | All 132 decisions and 20 statistics per sublink identical. `tx`, `rx`, `trsl` bit-identical over 30.2 M minutes each. `qc_flag` identical once the plateau-flag convention (bit 4 only where a value existed) is adopted |
| Detection limit for all 5132 sublinks (pycomlink 0.6.0, ITU 2005, V) | Identical to the CSV (relative difference 0) and to the NetCDF |
| Population from the files themselves | QC ⊂ approved RAW; none of the 13, the 6 R1, the 255 or the 8 no-frequency CMLs present; 2574 = 8 + 579 + 1987 and 5148 = 16 + 1158 + 243 + 3731; vs 2219: 1804 in both, 415 released-only, 183 new-only |
| Whole-file code review (1987 CMLs) | No deviation from the approved settings: `trsl` = `tx − rx` at every unflagged minute; no masked minute holds a value; flag totals equal the CSV; one reconstructed value per sublink |

**Documentation correction** (applied in `QC_PROCESSING_PARAMETERS.md`):
* bit-8 minutes have no RX observation, but TX may be observed (1893 retained minutes);
* bit-16 minutes have TX and/or RX missing;
* bit 4 marks only plateau minutes that held a value.

The NetCDF `trsl` comment carries the old wording. It is not changed now, because the file is not to be rewritten,
and will be corrected in the final product.

## 7. Open points for the author

* **Blackout reconstruction after the plateau mask (decision).** Blackout gaps are defined on RX *after* the
  plateau mask, as in v3 and Blettner (2023), and your approval excluded only gaps containing a P6 minute.
  * Consequence: a plateau run next to low RSL is filled with the reconstructed maximum loss. This happens in 2115
    of 9549 reconstructed runs, and 21,717 minutes carry both the plateau and the blackout flag.
  * This puts maximum attenuation exactly where the low readings were judged unreliable.
  * **Keep it** (faithful to the documented method), **or** exclude plateau-masked minutes from blackout
    reconstruction (they would stay missing)? The second choice needs a rerun of the QC.
* Review the candidate and its exclusion lists.
* Stage 4 (rainfall) needs your approval. It will:
  * change the rainfall population from 2219 to 1987 CMLs;
  * remove the 13 TX = 0 and the 4 R1 CMLs from the gauge benchmark;
  * use the new detection limits only as QC (the retrieval itself is unchanged).
