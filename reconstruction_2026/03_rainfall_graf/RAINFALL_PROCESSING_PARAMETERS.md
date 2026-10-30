# Stage 4 rainfall processing parameters

**Authorised:** by the author on 2026-10-28 (`PROCESSING_LOG.md`).

**Input:** the approved QC candidate `stage3\QC_CORRECTED_CANDIDATE_REV1.nc` (sha256 ca8afa95…, checked at start).

## Code

| Script | Role |
|---|---|
| `scripts\s4_01_rain_1min.py` | Native 1-minute rainfall. Copy of the final maintained exporter `essd-czech-cml-2021\src\essd_cz\export_1min.py` (sha256 abf56a31…). That exporter is the verified MetaCentrum exporter (sha256 88013e76…) with only its data paths made arguments |
| `scripts\s4_00_check_verbatim.py` | Proves the scientific block (imports, `pol_hv`, `call_kr`, `estimate`, `graf`, `decode`, `meta2`; 122 lines) is byte-identical to the maintained exporter and to the verified original |
| `scripts\s4_02_derive_5min.py` | 5-minute product via the repository's unchanged `essd_cz.rain5.derive_5min_dataset` |
| `scripts\s4_03_compare.py` | Comparison with the previously validated product |

The repository code is not modified.

## Retrieval (unchanged)

| Step | Setting |
|---|---|
| Wet/dry | Graf et al. (2020) RSTD: centred 60-minute rolling sample std of `trsl`. Threshold = 1.12 × the 0.8 quantile of the RSTD series (windows with ≥ 60 valid minutes). Classification uses windows with ≥ 45 valid minutes; wet if RSTD > threshold |
| Baseline | `pycomlink.processing.baseline.baseline_constant`, `n_average_last_dry = 1` (baseline = `trsl` at dry minutes, last dry value during wet minutes); A_obs = max(trsl − baseline, 0) |
| Wet antenna | `pycomlink.processing.wet_antenna.waa_pastorek_2021_from_A_obs`, A_max 14 dB, ζ 0.55, d 0.1, clipped to [0, A_obs] |
| Rain rate | `pycomlink.processing.k_R_relation.calc_R_from_A`, ITU_2005 (ITU-R P.838-3), vertical polarization, sublink frequency, path length, R_min 0. Dry valid minutes = 0 mm h⁻¹ |
| CML mean | `R_mean_1min` = mean of the available valid directional rates at each minute |
| 5-minute | A directional rate exists only if all five 1-minute samples are valid, and is their mean. Wet if any valid minute is wet. The CML mean is over finite directional 5-minute rates |

## Changes authorised for Stage 4

1. **Checks replaced.** The 2219-CML requirement and the bit-for-bit comparison with the old 5-minute product are
   removed. They are replaced by checks of the new input:
   * the input sha256;
   * time 308,160 one-minute steps, exactly 2021-04-01 00:00 – 2021-10-31 23:59;
   * 2 sublinks;
   * unique ids, with the count equal to the QC attribute;
   * every CML has a retained sublink;
   * finite frequency and length for retained sublinks;
   * polarization Vertical;
   * excluded sublinks entirely missing;
   * finite rate ⇔ `graf_valid_1min` = 1.
2. **Calibration period.** The Graf thresholds (RSTD quantile) and the dry baseline are computed on April–October
   2021, the period of the QC candidate. Previously: the full recorded history (2021-01-13..19 and
   2021-04-01..11-05).
3. **Flagged input.**
   * `trsl` is used as approved, including blackout-reconstructed (QC bit 8) and interpolated (QC bit 16) minutes.
   * Plateau-removed minutes stay missing, and no TX or RX is generated.
   * The new variable `trsl_source_flag_1min` carries bits 8 and 16 for every minute; the 5-minute file carries
     `blackout_minutes_5min` and `interpolated_minutes_5min`.

**How flags propagate:**
* the rate at a flagged minute comes from the reconstructed or interpolated total loss;
* the wet/dry decision at a minute uses a centred 60-minute window, which may contain flagged minutes;
* the baseline during a wet spell holds the last dry value, which may itself be an interpolated minute.

## Outputs (stage4\)

* `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc`, with `.summary.json`;
* `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc`;
* `RAINFALL_BEFORE_AFTER.csv`, `RAINFALL_BEFORE_AFTER_CML.csv`, `GRAF_THRESHOLD_COMPARISON.csv`;
* `STAGE4_RAINFALL_REPORT.md`.
