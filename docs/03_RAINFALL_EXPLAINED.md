# 3. Understanding the rainfall algorithm (one directional sublink)

The **original computation is not refactored**. We wrap it in short, named functions. See `src/essd_cz/rainfall.py` for those steps, and the unmodified code under `original/code/cml_rainfall/cml_graf_itu2005_rebuild.py` for the exact implementation.

1. Read the QC TRSL signal (`trsl`) and link geometry, frequency and polarization. A CML is a physical microwave path; it has **two directional sublinks**.
2. `graf()` computes rolling **60-minute standard deviation** of the signal. Its **80th percentile** across the time record, multiplied by **1.12**, is the threshold for wet/dry classification. At least 60 observations are required for estimation of each threshold-window standard deviation, and at least 45 for wet/dry detection. These are the actual source-code constants, not retuned defaults.
3. `estimate()` calls `pycomlink.processing.baseline.baseline_constant` (`n_average_last_dry=1`) to construct a dry reference. It uses the wet/dry mask to exclude wet observations from the baseline.
4. Compute nonnegative observed excess attenuation. Subtract **Pastorek 2021 wet-antenna attenuation** (`A_max=14`, `zeta=0.55`, `d=0.1`). Clamp attenuation to physically nonnegative values.
5. Convert rainfall attenuation into mm/h using **ITU-R P.838-3 (2005)** with the source link's frequency, path length and polarization. Non-wet valid minutes are assigned **0 mm/h**; invalid observations remain missing.
6. Produce a 5-minute rainfall value only if **all five** 1-minute values are valid. Average those five rainfall rates. Wet flag is true when any of the five valid minutes was classified as wet. `R_mean_5min` is the mean over whichever directional sublinks have a finite 5-minute retrieval.

Important: these are per-sublink functions, but some values are **global in time**. If you run only on 14 days, the 80th percentile can change, as can the constant dry baseline. The new numerical test therefore recomputes using the whole available *original QC time record*, THEN selects the requested comparison timestamps.

Do not confuse (a) the original **full-source run** (2020–2021 upstream file) with (b) the **Apr–Oct 2021 dataset segment** or (c) the old **ITU2003** and radar-assisted experiments. Their thresholds, baselines or physics may differ. `src/essd_cz/rainfall.py` deliberately imports the exact original algorithm functions rather than rewriting them.
