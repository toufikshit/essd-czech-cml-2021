#!/usr/bin/env bash
# Stage 6: Figure 7 (rainfall fields, 13 July 2021) from the released Graf rainfall product, run locally.
# Renderer: essd-czech-cml-2021/original/code/figures/03_blettner_event_maps_timeseries.py, unchanged, with the exact
# arguments of 06_run_multiple_events_vs_radar.py for event 4 / instance 1 (published Fig. 7). Only the inputs differ:
#   --gridded-nc : IDW grid of the Graf Stage 4 5-minute rainfall (1987 QC CMLs), built in Stage 5B with the unchanged
#                  IDW functions (nnear 15, p 2, max_dist 0.3 deg, 0.01 deg) on the lattice of the published grid
#                  (stage5b/fig47/work/grid_safe_stage4_0p01.nc, log stage5b/fig47/scripts/f7_03_safe_grid_stage4_output.txt);
#   --radar-root : corrected DWD RADKLIM YW 2017.002 2021 netCDF supplement (only YW_2017.002_20210713.nc extracted);
#   --selection-cache-csv : the handoff copy of selected_instances_cache.csv (event 4: 2021-07-13 19:05, same extent).
set -euo pipefail
R=/c/Users/taouf/ESSD_RELEASE
OUT=$R/FINAL_RECONSTRUCTION_20261027/stage6_manuscript/figure_work/fig7_out
mkdir -p "$OUT"
python -I $R/essd-czech-cml-2021/original/code/figures/03_blettner_event_maps_timeseries.py \
  --gridded-nc $R/FINAL_RECONSTRUCTION_20261027/stage5b/fig47/work/grid_safe_stage4_0p01.nc \
  --cml-var R --cml-row-label GRAF \
  --cml-geom-nc $R/FINAL_RECONSTRUCTION_20261027/stage3/QC_CORRECTED_CANDIDATE_REV1.nc \
  --radar-root $R/external_downloads/radklim_2021_supplement \
  --selection-cache-csv $R/extracted/full/ESSD_CZ_2021_HANDOFF/data/event_support/historical_event_cache/selected_instances_cache.csv \
  --event-rank 4 --instance-rank 1 --start-time "2021-07-13 16:05:00" \
  --ncols 8 --panel-step-minutes 60 --dpi 600 --skip-timeseries \
  --out-dir "$OUT" 2>&1 | tee "$OUT/run.log"
