# Notebook guide

Run the notebooks in order. Each one registers its output in the pipeline
manifest and the next one reads it from there. You never type a result path
by hand. Configure the data once in `config/user_paths.json`.

"Exact" below means the output was shown to match the released product or a
prior verified run number for number. "Adapted" means the workflow uses the
unchanged original functions on the released files, but the manuscript
figure was made with a different script, so pixel identity is not claimed.
The proof for every "exact" claim is in [EVIDENCE.md](EVIDENCE.md).

## 00 Explore Data

* Needs: all seven NetCDF files.
* Parameters: `RUN_NATIVE_AVAILABILITY` (on, about 2 minutes).
* Makes: network map, CML frequency and length figures, monthly and native
  availability tables and figures.
* Cost: minutes.
* Fidelity: availability uses the unchanged original script that made the
  paper figure (CML 66.35 percent mean and the four reference sensors above
  99.4 percent). The network and characteristics plots are
  reviewer figures.

## 01 Graf Rainfall

* Needs: the released native 1-minute rainfall file (the preprocessed CML
  file only for the optional recomputation).
* Parameters: `RECOMPUTE` (off by default), `WORKERS`.
* Makes: registers the released 1-minute product and materializes its exact
  analysis 5-minute derivation (five complete valid minutes per block) for the
  next notebooks. With `RECOMPUTE=True` it reruns the full original retrieval.
* Cost: about a minute for the one-time derivation. The full recomputation of
  1,987 links takes hours.
* Fidelity: exact for the 5-minute derivation (an independent re-derivation
  of the whole released file gave identical values). The recomputation uses
  the same retrieval functions as the released file.

## 02 Rainfall Validation

* Needs: the registered analysis 5-minute rainfall (derived from the 1-minute
  release in notebook 01), gauges; optional full history comparison needs the
  science extras installed.
* Parameters: `RUN_FULL_HISTORY_COMPARISON`, `RUN_GAUGE_VALIDATION`,
  `ALLOW_UNVERIFIED_V3` (all off by default).
* Makes: rainfall time series figure; optional one link numerical regression;
  optional gauge confusion matrix and hexbins.
* Cost: seconds by default; the optional steps are heavy.
* Fidelity: the one link regression compares a recomputation with the
  released values. The gauge validation uses the Graf helper
  `cml_graf_rstd_hess_radar_gauge_v3.py`, the same helper used for the
  released validation numbers; it runs only with `ALLOW_UNVERIFIED_V3`
  because the original helper file was not recovered.

## 03 Grid CML Rainfall

* Needs: stage 01 rainfall and the preprocessed file for link geometry.
* Parameters: grid resolution, neighbour count, power, maximum distance,
  time window.
* Makes: a 5 min rainfall grid NetCDF built with the unchanged original
  inverse distance weighting functions.
* Cost: minutes.
* Fidelity: adapted. The paper's Fig. 7 grid is built with the grid builder
  in `reconstruction_2026/05_figures/fig7_grid_builder/` (same IDW settings);
  this notebook uses the two original IDW functions directly.

## 04 Spatial Radar Maps

* Needs: the stage 03 grid, a native RADKLIM YW radar folder, cartopy.
* Parameters: `RUN_MAPS` (off by default), panel time, extent.
* Makes: side by side CML and radar map panels with the unchanged original
  renderer.
* Cost: minutes to tens of minutes.
* Fidelity: adapted. Without the historical event cache the event choice is an
  explicit demonstration, not the paper's event discovery.

## 05 E band Water Vapour

* Needs: the V6 support NetCDF from the deposit, or raw CML, gauge and
  temperature files for a small rerun.
* Parameters: `recompute`, `max_series`, pressure.
* Makes: registers the V6 water vapour product (1,352 directional series).
* Cost: seconds to register the full product. A small rerun takes minutes and
  is a demonstration only.
* Fidelity: the registered full product is the original. Small reruns are
  demonstrations and must not feed the paper figures.

## 06 WVD Sensitivity

* Needs: the registered V6 product plus the diagnostics and quantization CSVs
  and the CHMI temperature and humidity files.
* Parameters: `RUN_SENSITIVITY` (off until the input check passes), branch,
  station distance.
* Makes: station matched validation pairs, binned sensitivity metrics and
  figures, with the unchanged original analysis script (plus a documented
  NetCDF read compatibility copy).
* Cost: the full 1,352 series pass is expensive.
* Fidelity: exact original analysis. The notebook refuses mismatched inputs
  instead of filtering silently.

## 07 WVD Paper Population

* Needs: the stage 06 pair cache.
* Parameters: figure dpi, hexbin grid size, plot range.
* Makes: the paper population metrics and hexbins. The selection is 45
  directional sublinks from 31 physical CMLs.
* Cost: minutes.
* Fidelity: exact original script with the paper filters (E band, length at
  least 4 km, quantization at most 0.20 g per cubic metre, GRAF branch).

## 08 CHMI CSV to NetCDF (optional, upstream)

* Needs: internet access to CHMI open data. Not needed to use the released
  files.
* Parameters: `DRY_RUN` (on by default), year, element list.
* Makes: CHMI station NetCDF files with the unchanged original converter.
* Cost: download heavy.
* Fidelity: exact original converter.
