# 4. Numeric comparison to the released Graf rainfall NetCDF

The released product is the **native 1-minute** rainfall NetCDF. The
comparator recomputes a link from the QC file with the original functions and
compares the aggregated 5-minute fields against the published values —
deriving them per link from the 1-minute release with the exact
complete-block rule (`src/essd_cz/rain5.py`). The released rainfall was
computed from the released QC file over its whole period (1 April to
31 October 2021), so the QC file is the right input.

## First run

```bash
python run.py --config config/user_paths.json compare-rainfall --link-index 3 --days 14 --threads 2 \
    2>&1 | tee results/compare_link_3.log
```

This calculates **all source minutes for only link index 3** (both directional sublinks). The 14 days control **only the dates compared and saved**, not the duration used to calculate wet/dry and baseline. Outputs: `results/rainfall_comparison/comparison_report.json`, sublink-specific aligned CSV files, and `comparison_link_3.png`.

## What the report measures

- `cml_id` is matched *exactly* across QC and final rainfall; order differences cannot masquerade as agreement.
- `time_5min` (derived from `time_1min` for the native release) is matched *exactly*. A missing timestamp is an error, not silently shifted.
- `R_sublink_5min` and `R_mean_5min`: finite-sample RMSE/bias/max difference, number outside tolerance (`atol=rtol=1e-4` by default), **and missing-value-mask disagreements**. An excellent RMSE with many extra NaNs does not pass.
- `graf_wet_5min`, `graf_valid_5min`, `wet_minutes_5min`, `valid_minutes_5min`: exact flag/count equality where available.
- Graf Q80 and threshold: numeric agreement per direction, if included in the published file.
- Overall `pass` is true only if all present comparable fields pass. If the output has no paired finite rainfall values it cannot pass.

**Interpretation:** If this test fails, inspect `comparison_report.json` and the per-sublink CSV first. Potential reasons include the original detection history window, a different pycomlink version, incorrect source signal units, or a scientific implementation mismatch. Do not automatically change the method to force agreement. The released file was written by `reconstruction_2026/03_rainfall_graf/scripts/s4_01_rain_1min.py` (version `ITU2005_REBUILD_1.0`, same retrieval functions); the comparator is a stored-source replay, not a substitute formula.

## Additional links

You can test several links without rebuilding all 1,987:

```bash
python run.py --config config/user_paths.json compare-rainfall \
    --link-index 3 --link-index 25 --start 2021-06-01 --days 7 --threads 2
```

These indices are examples, not specially selected validation sites. Running more links takes longer. Recompute-before-subsetting is essential; do not compare a short-window smoke test directly to publication values.
