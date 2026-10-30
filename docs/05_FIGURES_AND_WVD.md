# 5. Figures, reference validation and supporting water-vapour extraction

## Available from seven released files

`python run.py --config config/user_paths.json availability`

The original `availability_native_timestep_percentage.py` produces CML, radar-path, gauge, temperature and relative-humidity percentage time series and saves a PNG/PDF plus five CSVs. On the released files the means are **66.35%, 99.48%, 99.74%, 99.57%, 99.56%**, respectively. Those are time-mean values for the five different native cadences, not interchangeable rainfall metrics.

`python run.py --config ... validate-gauge --detector /path/to/exact/cml_graf_rstd_hess_radar_gauge.py`

The final paper hexbin and confusion plotting script is preserved unchanged in `original/code/figures/rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py`. The released validation numbers (N 1,331,687, RMSE 3.63 mm/h, MBE −0.91 mm/h, PCC 0.633 for gauge-wet intervals) were computed with a copy of its validation function and the Graf helper `original/notebook_references/wetdry_validation/cml_graf_rstd_hess_radar_gauge_v3.py` (`reconstruction_2026/04_validation/`). The CLI loads that helper only when you pass `--allow-unverified-v3` explicitly; the original helper file itself was not recovered. The plotting script matches the nearest gauge to each CML midpoint (≤5 km), aggregates two complete CML 5-min windows against end-labelled 10-min gauge amounts, and requires the gauge confusion matrix to use the same sample population as the hexbin.

The original `original/code/figures/06_run_multiple_events_vs_radar.py` reads native RADKLIM files and a precomputed CML IDW grid, and is not runnable from the seven published NetCDF files alone. Avoid confusing this with the final Graf/gauge validation figure.

## WVD: a supporting demonstration, not an eighth core release NetCDF

Supporting scripts are under `original/code/wvd/`, notably `13_wvd_raw_60_100_independent.py`, `19_wvd_eband_progressive_baseline_v6.py`, `20_run_wvd_v6_parallel_hpc.py`, `21_validate_wvd_v6_vs_original.py` and `33_v6_binned_validation_hexbins.py`. These include prior V6 intermediate inputs, E-band selection, exclusion masks and temperature/relative-humidity reference handling. They were not proven to run end-to-end on *only* the seven NetCDF files in a clean environment. Therefore this repository **does not claim one-command WVD reproduction**. Read source metadata and `docs/06_PROVENANCE_GAPS.md` before adding a WVD entry point.

## Upstream CHMI / DWD

`original/code/chmi/chmi_historical_csv_to_netcdf_v5.py` is the confirmed converter from CHMI original CSV. `original/code/radar/dwd_radolan_processor.py` is recovered, but the exact corrected 2021 RADKLIM replay is not established; neither upstream script is necessary just to read the released NetCDF products.

## Paper figures

All figures show the released dataset. The figure runners and the figure provenance table are in `reconstruction_2026/05_figures/`.

* Fig. 1: map of the raw dataset (2,574 CMLs), drawn in the published style.
* Fig. 2: availability, with the unchanged original script; CML row from the raw TX and RX of the 2,574 CMLs (mean 67.83%).
* Fig. 3: frequencies and path lengths of the 2,574 raw CMLs (5,132 sublinks with frequency; 8 CMLs without frequency appear only in the path-length panel).
* Fig. 4: event time series (13–14 July 2021), with the original producer.
* Figs. 5 and 6: gauge detection and gauge-wet hexbin with the original plotting functions, all samples: N = 1,331,687, RMSE 3.63 mm/h, bias −0.91 mm/h, correlation 0.633, POD 74.6 %, correct negatives 88.6 %.
* Fig. 7: rainfall fields of the released Graf rainfall, gridded with the unchanged IDW (`reconstruction_2026/05_figures/fig7_grid_builder/cml_grid_essd_publication_safe_singleopen.py`) and drawn with the unchanged renderer `original/code/figures/03_blettner_event_maps_timeseries.py`.
* Fig. 8: water vapour (31 CMLs, 45 sublinks, 236,793 pairs; RMSE 3.46 g m−3, correlation 0.598), with `original/code/wvd/plot_wvp_support_hexbins_no_titles.py`.
