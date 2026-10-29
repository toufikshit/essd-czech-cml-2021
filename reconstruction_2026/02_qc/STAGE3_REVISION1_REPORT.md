# Stage 3, revision 1: blackout reconstruction excludes plateau-removed and P6-masked minutes

**Date:** 2026-10-28.
**Status:** revised QC candidate produced and independently verified. Awaiting author approval. Stage 4 has not
started.

**Decision:** the author's, of 2026-10-28 (`PROCESSING_LOG.md`).

**Preserved unchanged** (hashes re-checked, `FIRST_CANDIDATE_SHA256.txt`): `QC_CORRECTED_CANDIDATE.nc`
(sha256 e193b2cf…), `QC_SUBLINK_RESULTS.csv`, `s3_01_qc_corrected.py` and every earlier report.

## 1. The change

The single change is in blackout eligibility. A blackout gap that contains a plateau-removed minute (`qc_flag`
bit 4) or a P6-masked minute (bit 2) is no longer reconstructed; `trsl` stays missing and the minutes keep their
flags. Eligible genuine blackout gaps are still reconstructed with option C.

**Formula and reference period:**
* `trsl = max(tx) − min(rx)` of the same sublink;
* over its observations in 2021-04-01 00:00 – 2021-10-31 23:59 UTC, after masking;
* written to `trsl` only, with bit 8;
* RX is missing at these minutes; TX may be observed.

This is documented as an adaptation of Blettner et al. (2023) in `QC_PROCESSING_PARAMETERS_REV1.md` and in the
NetCDF (`qc_blackout_method` and the `trsl` and `qc_flag` comments).

**Script:** `scripts\s3_04_qc_corrected_rev1.py`, sha256 3765739a…. The diff against s3_01 is in
`scripts\s3_04_diff_vs_s3_01.txt`: the eligibility test, two counters and the metadata text.

## 2. Revised candidate

| Property | Value |
|---|---|
| File | `QC_CORRECTED_CANDIDATE_REV1.nc` |
| Size | 2,409,368,783 bytes |
| SHA-256 | `ca8afa95aacb50dacbd5f7d96f776c56a78c63d1a4888e9f721a29cdcfafde65` |
| Dimensions | time 308,160; sublink_id 2; cml_id 1987 |
| File check | Passed read-only (`verification\s3_06_file_verification_rev1_output.txt`): variables, metadata, start/middle/end samples, `trsl` = `tx − rx` at unflagged minutes |

## 3. Population (unchanged by the revision)

| Population | CMLs | Sublinks |
|---|---:|---:|
| Approved RAW | 2574 | 5148 |
| Excluded before QC (no frequency) | 8 | 16 |
| Excluded by QC | 579 | 1158 |
| Excluded sublinks in retained CMLs | – | 243 |
| **Retained** | **1987** | **3731** |
| In QC but not in RAW | 0 | – |

* No sublink decision changed (0 of 5132). The retained CML set is identical to the first candidate.
* The reason: reconstructed minutes were already excluded from the fluctuation statistics and never used for
  interpolation, and the detection limit does not depend on them.
* All consistency checks pass (`scripts\s3_05_tables_checks_rev1_output.txt`):
  * QC is a strict subset of the approved RAW;
  * none of the 13, the 6 R1, the 255 or the 8 no-frequency CMLs is present;
  * both reconciliations are exact.

The difference from the released 2219 is the same as in `STAGE3_QC_REPORT.md` §5 (2219 − 415 + 183 = 1987).

## 4. Blackout counts

| | First candidate | Revision 1 | Rejected by the revision |
|---|---:|---:|---:|
| Retained sublinks, gaps | 8,194 | **6,410** | 1,784 |
| Retained sublinks, minutes | 57,351 | **37,414** | 19,937 |
| Retained sublinks with ≥ 1 reconstructed minute | 1,520 | **1,136** | – |
| All assessed sublinks, gaps / minutes | 14,192 / 116,561 | 9,234 / 54,618 | 4,958 / 61,943 |

* Reconstructed minutes are now 0.006 % of the retained observed total-loss minutes.
* Interpolation (651,631 gaps; 706,377 minutes), plateau masking, P6 masking (14 minutes in the period) and all
  observed values are identical to the first candidate.

## 5. Independent verification (Agent 4; `agents\A4_stage3_verification.md`, section "Revision 1")

| Check | Result |
|---|---|
| Own implementation of the revised rule: 70 CMLs / 140 sublinks, including 26037_26137_3704 and 89084_89099_2244 | 0 mismatches in decisions, statistics and reconstructed/rejected counts; `tx`, `rx`, `trsl` and `qc_flag` bit-identical over 32.7 M minutes each |
| REV1 vs first candidate, all 1987 CMLs | Same ids, time axis, status, detection limits and fractions. `tx`, `rx`: 0 differences. `trsl`, `qc_flag`: 23,640 minutes changed (19,937 in retained sublinks), every one a former bit-8 minute now missing. They form 2115 whole runs, each containing a bit-2 or bit-4 minute; there are no partial runs and no new bit-8 minutes |
| REV1 bit-8 minutes | None coincides with bit 2 or 4; RX always missing; `trsl` always finite; 1750 retained minutes keep an observed TX |
| Metadata | `qc_blackout_method`, the `trsl` and `qc_flag` comments, `history` and `qc_script_sha256` match the code; the first-candidate hash in `history` matches the preserved file |

**Minor remarks, to address when the final product is written (no rewrite now):**
* `qc_parameters_json` does not encode the rejection rule (it is stated in `qc_blackout_method`);
* the title does not mention revision 1;
* rejected blackout minutes have no flag bit of their own (they are counted in `QC_SUBLINK_RESULTS_REV1.csv`).

## 6. Files

* `QC_CORRECTED_CANDIDATE_REV1.nc`
* `QC_SUBLINK_RESULTS_REV1.csv`
* `QC_EXCLUDED_CMLS_REV1.csv`
* `QC_EXCLUDED_SUBLINKS_REV1.csv`
* `QC_BEFORE_AFTER_COMPARISON_REV1.csv`
* `QC_POPULATION_RECONCILIATION_REV1.csv`
* `QC_PROCESSING_PARAMETERS_REV1.md`
* scripts `s3_04`, `s3_05` and their outputs
* `verification\s3_06_*`
* `agents\A4_stage3_verification.md`

## 7. Approval requested

Approve `QC_CORRECTED_CANDIDATE_REV1.nc` as the verified QC candidate for Stage 4, or request changes. Stage 4
(rainfall retrieval) remains paused until you approve it.
