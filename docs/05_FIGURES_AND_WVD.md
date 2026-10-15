# 5. Figures, reference validation and supporting water-vapour extraction

## Available from seven released files

`python run.py --config config/user_paths.json availability`

The original `availability_native_timestep_percentage.py` produces CML, radar-path, gauge, temperature and relative-humidity percentage time series and saves a PNG/PDF plus five CSVs. On the research HPC this completed on real data with means of **58.67%, 99.48%, 99.74%, 99.57%, 99.56%**, respectively. Those are time-mean values for the five different native cadences, not interchangeable rainfall metrics.

`python run.py --config ... validate-gauge --detector /path/to/exact/cml_graf_rstd_hess_radar_gauge.py`

The final paper hexbin and confusion plotting script is preserved unchanged in `original/code/figures/rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py`. It needs an **exact historically executed Graf helper** that was *not* recovered. If you want to test code compatibility only, specify `--allow-unverified-v3` explicitly; this is NOT proof of the original numbers. The plotting script matches the nearest gauge to each CML midpoint (≤5 km), aggregates two complete CML 5-min windows against end-labelled 10-min gauge amounts, and requires the gauge confusion matrix to use the same sample population as the hexbin.

The historical `original/code/figures/06_run_multiple_events_vs_radar.py` reads native RADKLIM files and an earlier event grid, and is not runnable from the seven published NetCDF files alone. Avoid confusing this with the final Graf/gauge validation figure.

## WVD: a supporting demonstration, not an eighth core release NetCDF

Supporting scripts are under `original/code/wvd/`, notably `13_wvd_raw_60_100_independent.py`, `19_wvd_eband_progressive_baseline_v6.py`, `20_run_wvd_v6_parallel_hpc.py`, `21_validate_wvd_v6_vs_original.py` and `33_v6_binned_validation_hexbins.py`. These include prior V6 intermediate inputs, E-band selection, exclusion masks and temperature/relative-humidity reference handling. They were not proven to run end-to-end on *only* the seven NetCDF files in a clean environment. Therefore this repository **does not claim one-command WVD reproduction**. Read source metadata and `docs/06_PROVENANCE_GAPS.md` before adding a WVD entry point.

## Upstream CHMI / DWD

`original/code/chmi/chmi_historical_csv_to_netcdf_v5.py` is the confirmed converter from CHMI original CSV. `original/code/radar/dwd_radolan_processor.py` is recovered, but the exact corrected 2021 RADKLIM replay is not established; neither upstream script is necessary just to read the released NetCDF products.
