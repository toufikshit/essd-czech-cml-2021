# Data dictionary

| Role | File basename | Time step | Entities | Example variable |
|---|---|---|---:|---|
| Raw telecom | `CZ_CML_RAW_1MIN_20210401_20211031.nc` | 1 min | 2,574 CML links (5,148 sublink slots) | link signal |
| QC telecom | `CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` | 1 min | 1,987 CML links, two directional sublinks (3,731 retained) | `trsl`, `qc_flag` |
| Final rainfall | `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc` | 1 min | 1,987 CML links, two directions | `R_sublink_1min` (mm h−1), `R_mean_1min`, wet/valid flags, `trsl_source_flag_1min`, per-link thresholds |
| CHMI rain gauges | `CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc` | 10 min | 424 stations | `rr` (**mm per 10 min**, not mm h−1) |
| Radar along links | `CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc` | 5 min | 1,132 paths | `rr` (**mm per 5 min**, not mm h−1) |
| Air temperature | `CHMI_10MIN_TEMPERATURE_20210401_20211031.nc` | 10 min | 251 stations | `T` |
| Relative humidity | `CHMI_10MIN_HUMIDITY_20210401_20211031.nc` | 10 min | 247 stations | `H` |

The files cover **1 April 00:00 to 31 October 23:59 UTC, 2021**. Native slots in that interval are **308,160** (1-min), **61,632** (5-min), **30,816** (10-min). The released filenames state the coverage window (`20210401_20211031`).

**The rainfall release is native 1-minute.** The analysis 5-minute series used by the manuscript (validation, event series, gridding) are derived from it: a directional 5-minute rate exists only where **all five** 1-minute samples of the wall-clock block are valid and is their arithmetic mean; the 5-minute wet flag is set when any of those valid minutes was wet; the physical-link 5-minute mean is computed from the finite directional 5-minute rates afterwards. `scripts/derive_5min.py` does this derivation (functions in `src/essd_cz/rain5.py`).

**Alignment:** CML `cml_id` is structural and can contain underscores; never coerce it to a floating-point number. Released 1-minute rainfall is labelled by `time_1min`; derived rainfall at the start of each 5-minute window follows `time_5min`. The gauge's 10-minute record is an *end-labelled* accumulation; the validation script aligns two complete CML 5-minute bins to it. For radar data, `rr × 12` converts 5-min mm to mm/h; for gauge data, `rr × 6` converts 10-min mm to mm/h. The validation script does the necessary shift/alignment explicitly; do not compare raw arrays by positional index.

**Populations.** Raw: 2,574 CMLs (5,148 sublink slots), 2,566 with frequency metadata, 904 of them E band; the list is `config/population/RAW_APPROVED_POPULATION.csv`. QC: 1,987 CMLs with 3,731 retained sublinks, 810 E band CMLs; 8 CMLs without frequency are excluded before QC and 579 by QC. Rainfall: the same 1,987 CMLs.

## QC variables (`CZ_CML_PREPROCESSED_1MIN_...nc`)

| Variable | Dimensions | Meaning |
|---|---|---|
| `tx`, `rx` | time × sublink × CML | observed TX and RX in dBm after QC. Never filled; missing where QC masked a minute or the sublink is excluded |
| `trsl` | time × sublink × CML | total loss in dB: `tx − rx` where both exist, plus the flagged values of bits 8 and 16 |
| `qc_flag` | time × sublink × CML | uint8 bit field, see below |
| `qc_sublink_status` | sublink × CML | 0 retained; 1 TX = 0 direction (non-operational); 2 no observed total loss in the period; 3 fluctuation; 4 detection limit; 5 fluctuation and detection limit |
| `qc_detection_limit_mmhr` | sublink × CML | minimum detectable rain rate, R_min = (0.33 dB / (k L))^(1/α), ITU-R P.838-3 (2005), vertical; sublinks with R_min ≥ 2 mm h−1 are excluded |
| `qc_fraction_5h_windows_std_gt_2dB`, `qc_fraction_1h_windows_std_gt_0p8dB` | sublink × CML | fluctuation statistics (fail above 0.10 and 0.33) |

`qc_flag` bits (a minute can carry several bits):

| Bit | Meaning |
|---:|---|
| 1 | invalid value masked (TX ≤ −99 dBm, TX = 255, RX ≤ −99.9 dBm, RX > −20 dBm) |
| 2 | mixed-placeholder minute masked (TX and RX) |
| 4 | minute held a TX or RX value removed by the low-RSL plateau filter |
| 8 | blackout: total loss reconstructed in `trsl` only as max(TX) − min(RX) of the sublink over April–October 2021; RX is missing, TX may be observed |
| 16 | short gap (≤ 5 min) of `trsl` filled by linear interpolation; TX and/or RX missing |

At minutes with bit 8 or 16, `trsl` is not equal to `tx − rx`.

## Rainfall source flag

`trsl_source_flag_1min` (time × sublink × CML, uint8) copies bits 8 and 16 of `qc_flag` at the same minute: 0 = the rain rate comes from measured total loss (or is missing); 8 = from blackout-reconstructed total loss; 16 = from interpolated total loss. Bit-8 minutes are 0.006 % of the valid minutes but 3 % of the accumulation. Leave them out (`trsl_source_flag_1min & 8`) to use measured signal only.

The 5-minute series are not a release file; `scripts/derive_5min.py` derives them from the 1-minute file. It also counts, per 5-minute interval, the valid 1-minute samples that come from bit-8 or bit-16 total loss (`blackout_minutes_5min`, `interpolated_minutes_5min`); the validation without blackout-reconstructed intervals leaves out intervals with `blackout_minutes_5min` > 0.
