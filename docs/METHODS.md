# Processing methods and settings

This page lists the settings of the scripts that produced the released files. The
scripts are in `scripts/`; the unchanged original scientific scripts they call are
in `original/code/`.

## Populations

| Product | Population |
|---|---|
| RAW signals | 2,574 CMLs, 5,148 sublink slots; 2,566 CMLs with frequency metadata, 904 of them E band |
| QC signals | 1,987 CMLs, 3,731 retained sublinks; 810 E band CMLs |
| Rainfall, 1 minute | 1,987 CMLs, 631,742,984 valid directional values |
| Radar along CML paths | 1,132 paths |
| CHMI stations | 424 rain gauges, 251 temperature, 247 humidity |

2,574 CMLs = 8 without frequency (excluded before QC) + 579 excluded by QC +
1,987 retained. The table `config/population/RAW_APPROVED_POPULATION.csv` lists
every sublink of the raw dataset with its frequency, length, the CMLs without
frequency and the 15 sublinks with TX = 0 dBm at every minute.
`config/population/QC_PREDEFINED_MINUTE_MASKS.csv` lists 15 minutes that hold the
1-minute mean of a placeholder and a valid value; QC masks TX and RX there.

## Quality control (`scripts/qc_from_raw.py`)

| Step | Setting |
|---|---|
| Period | 2021-04-01 00:00 to 2021-10-31 23:59 UTC, 308,160 minutes |
| Invalid values | TX <= -99 dBm, TX = 255, RX <= -99.9 dBm, RX > -20 dBm, abs >= 1e35 -> missing (flag 1) |
| Mixed-placeholder minutes | TX and RX masked (flag 2) |
| Low-RSL plateau | centred 3-minute window, >= 2 valid samples, max < -85 dBm and std < 0.5 dB; TX and RX masked +/- 2 min (flag 4) |
| Blackout (adaptation of Blettner et al., 2023) | RX gap <= 60 min with RX < -65 dBm on at least one edge; `trsl` = max(TX) - min(RX) of the sublink in the period; TX and RX are not filled; gaps that contain a plateau or mixed-placeholder minute are not reconstructed (flag 8) |
| Short gaps | linear interpolation of `trsl` over <= 5 minutes missing in RAW, never across a QC mask or next to a blackout (flag 16) |
| Fluctuation | windows 5 h and 1 h, valid only with 300 / 60 valid minutes; fail if > 10 % of 5-h windows have std > 2 dB or > 33 % of 1-h windows have std > 0.8 dB |
| Detection limit | per sublink, R_min = (0.33 dB / (k L))^(1/alpha), ITU-R P.838-3 (2005), vertical; fail if R_min >= 2 mm/h |
| Retention | a sublink is removed if it fails a test (values set to missing); a CML is kept if at least one sublink passes |

## Rainfall (`scripts/rainfall_1min.py`, `scripts/derive_5min.py`)

| Step | Setting |
|---|---|
| Wet/dry | Graf et al. (2020) RSTD: centred 60-minute rolling std of `trsl`; threshold = 1.12 x the 0.8 quantile per sublink; 60 / 45 minimum valid minutes; computed over 1 April to 31 October 2021 |
| Baseline | pycomlink `baseline_constant`, `n_average_last_dry = 1` |
| Wet antenna | Pastorek et al. (2021): A_max 14 dB, zeta 0.55, d 0.1 |
| Rain rate | ITU-R P.838-3 (2005), vertical, sublink frequency and path length; dry valid minutes = 0 mm/h |
| 5-minute series | a directional 5-minute rate exists only if all five 1-minute values are valid; wet if any valid minute is wet; CML mean over the finite directional 5-minute rates |

The rainfall uses the QC `trsl`, including blackout-reconstructed (flag 8) and
interpolated (flag 16) minutes. The variable `trsl_source_flag_1min` tells which
minutes these are.

## Gauge validation (`scripts/gauge_validation.py`)

Nearest CHMI gauge within 5 km of the CML midpoint; two complete 5-minute CML
bins against the end-labelled 10-minute gauge amount; gauge wet >= 0.1 mm/h; Graf
detection from the QC `trsl`.

| Gauge-wet samples | N | RMSE (mm/h) | MBE (mm/h) | PCC |
|---|---:|---:|---:|---:|
| All samples | 1,331,687 | 3.63 | -0.91 | 0.633 |
| Without blackout-reconstructed intervals | 1,327,998 | 3.38 | -0.89 | 0.606 |

## Radar along the CML paths (`scripts/radar_cml_paths.py`)

RADKLIM YW 2017.002 (DWD, 5 min; for 2021 the supplement files published by DWD); path
averages with poligrain intersection weights (path length weights). A path
average is missing when one of its radar pixels is missing.

## Water vapour (`scripts/wvd_validation.py`)

The V6 retrieval (`original/code/wvd/19_wvd_eband_progressive_baseline_v6.py`,
run for every series by `20_run_wvd_v6_parallel_hpc.py`, helper
`13_wvd_raw_60_100_independent.py`) reads the raw signals of the E-band CMLs and
produced the WVD supporting file of the deposit; the quantization per link comes
from `31_run_v6_quantization_parallel.py` / `30c_v6_quantization_chunk.py`. The
WVD file was computed from the full raw record of 2021, because the causal
baseline uses the completed days before each day, including days before April. A
rerun on the released RAW file, which starts on 1 April, therefore starts its
baseline history later. `scripts/wvd_validation.py` uses the released WVD files. The
validation subset is 31 CMLs and 45 directional sublinks with 236,793 pairs,
RMSE 3.46 g m-3 and PCC 0.598.
