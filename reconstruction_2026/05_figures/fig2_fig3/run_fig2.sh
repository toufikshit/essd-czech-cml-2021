#!/bin/bash
# Runs the UNMODIFIED producer availability_native_timestep_percentage.py (sha256 b8d7b983...) twice.
# Only CLI paths differ. Run 1 = original inputs (reproduction); run 2 = CML file replaced by approved QC REV1.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
F=${F:-/c/Users/taouf/ESSD_RELEASE/extracted/full/ESSD_CZ_2021_HANDOFF/data/source_full_history}
QCNEW=${QCNEW:-/c/Users/taouf/ESSD_RELEASE/FINAL_RECONSTRUCTION_20261027/stage3/QC_CORRECTED_CANDIDATE_REV1.nc}
OUT=${OUT:-$HERE/../reproduction}
python -I -B "$HERE/availability_native_timestep_percentage.py" --cml-file "$F/qc/CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc" --radar-file "$F/radar/CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc" --gauge-file "$F/gauge/CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc" --temperature-file "$F/temperature/CHMI_10MIN_TEMPERATURE_2021.nc" --humidity-file "$F/humidity/CHMI_10MIN_HUMIDITY_2021.nc" --out-dir "$OUT/fig2_original"
python -I -B "$HERE/availability_native_timestep_percentage.py" --cml-file "$QCNEW" --radar-file "$F/radar/CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc" --gauge-file "$F/gauge/CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc" --temperature-file "$F/temperature/CHMI_10MIN_TEMPERATURE_2021.nc" --humidity-file "$F/humidity/CHMI_10MIN_HUMIDITY_2021.nc" --out-dir "$OUT/fig2_new"
