# QC processing parameters (Stage 3 corrected QC candidate)

**Script:** `scripts\s3_01_qc_corrected.py`. Its SHA-256 is written into the NetCDF attribute `qc_script_sha256`.

**Approved:** by the author on 2026-10-28 (`PROCESSING_LOG.md`), following `STAGE3_QC_SETTINGS_PROPOSAL.md`.

**Input:** the full-history RAW `CZ_CML_RAW_1MIN_20200101_20211231.nc` (sha256 6acfeed2…, read-only). Its window is
identical to the released RAW.

**Kernel origin:** the standalone v3 kernel (`cml_qc_filter_standalone_v3.py`), re-implemented with the changes
marked "≠ v3".

| Item | Setting | Relation to v3 / previous QC |
|---|---|---|
| Period | 2021-04-01 00:00 – 2021-10-31 23:59 UTC (308,160 min; time index 656640–964799 in minutes since 2020-01-01) | ≠ previous QC (whole file) |
| Population | Approved RAW 2574 CMLs; 8 without frequency excluded before QC; 2566 assessed (5132 sublinks) | Stage 2 decisions |
| P5 | 15 TX = 0 directions: non-operational, never valid; the opposite direction is processed normally | Stage 2 decision |
| P6 | Mixed-placeholder minutes: TX and RX masked (flag bit 2). 14 of the 15 lie in the period; 26037_26137_3704 (2021-11-04) does not | Stage 2 decision |
| Invalid values | TX ≤ −99 dBm, TX = 255, \|TX\| ≥ 1e35; RX ≤ −99.9 dBm, RX > −20 dBm, \|RX\| ≥ 1e35 → missing (bit 1) | = v3 |
| Low-RSL plateau | Centred 3-minute RSL window with **≥ 2 valid samples**: max < −85 dBm and population std < 0.5 dB; TX and RX masked ±2 min (bit 4) | ≠ v3 (v3: ≥ 1 sample) |
| Blackout (option C) | Whole RX gap ≤ 60 min with RSL < −65 dBm on at least one edge. TX and RX stay missing; `trsl` = series max TX − series min RX of the period (after masking); bit 8. Gaps containing a P6 minute are not reconstructed. Excluded from the fluctuation statistics | ≠ v3 (v3 writes max TX / min RX into TX and RX) |
| Short gaps | Linear interpolation of `trsl` over ≤ 5 missing minutes that are missing in RAW (not QC-masked), with observed total loss on both sides; never next to a reconstructed blackout minute; bit 16 | ≠ v3 (v3 also bridges QC masks) |
| Fluctuation | Per sublink, on `trsl` after interpolation without blackout minutes. Centred windows of 300 and 60 min, population std. **Valid window = 300 / 60 valid minutes.** Fail if > 10 % of valid 5-h windows have std > 2 dB or > 33 % of valid 1-h windows have std > 0.8 dB, or if there is no valid window | ≠ v3 for the 5-h window (v3: ≥ 60 minutes); thresholds = v3 and manuscript |
| Detection limit | Per sublink, `R_min = (0.33 dB / (k L))^(1/α)`; k and α from ITU-R P.838-3 (2005) via `pycomlink.processing.k_R_relation.a_b(f_GHz, 'V', approx_type='ITU_2005')` (pycomlink 0.6.0, as in the rainfall retrieval); sublink frequency, metadata length, vertical. Fail if R_min ≥ 2 mm h⁻¹ or NaN | ≠ v3 (`0.33/(kL)`, no α) and ≠ previous QC (2003 table, CML minimum frequency) |
| Sublink exclusion | P5, no observed total loss in the period, fluctuation fail, or detection-limit fail. TX, RX and `trsl` set to missing; `qc_flag` and statistics kept | ≠ v3 (flag only) |
| CML retention | Retained if ≥ 1 sublink passes; written to the NetCDF. Others are listed in `QC_EXCLUDED_CMLS.csv` | Manuscript |

**Output variables** (`QC_CORRECTED_CANDIDATE.nc`):
* `tx`, `rx`, `trsl` (float32, time × sublink_id × cml_id);
* `qc_flag` (uint8 bit field: 1 invalid, 2 mixed-placeholder minute, 4 plateau, 8 blackout total loss reconstructed,
  16 interpolated);
* per sublink: `qc_sublink_status`, `qc_detection_limit_mmhr`, `qc_fraction_5h_windows_std_gt_2dB`,
  `qc_fraction_1h_windows_std_gt_0p8dB`;
* metadata: `length`, `frequency`, `polarization`, endpoint coordinates;
* global: `qc_parameters_json`, the input sha256 and the script sha256.

At minutes with bit 8 or 16 set, `trsl` is not equal to `tx − rx`.

**Correction (2026-10-28, from Agent 4's verification):**
* At bit-8 minutes there is **no RX observation**; TX may be observed (1893 retained blackout minutes keep an
  observed TX, because the gap is defined on RX).
* At bit-16 minutes TX and/or RX is missing.
* Bit 4 marks only plateau minutes that held a value before masking.

The `trsl` comment inside the candidate NetCDF still says "no TX or RX observation exists". It is left unchanged
(no rewrite of the file without approval) and is to be corrected when the final product is written.
