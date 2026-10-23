# 2. Seven data products — quick dictionary

| Role | File basename | Time step | Entities | Example variable |
|---|---|---|---:|---|
| Raw telecom | `CZ_CML_RAW_1MIN_20210401_20211031.nc` | 1 min | 2,580 CML links | link signal |
| QC telecom | `CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` | 1 min | 2,219 CML links, two directional sublinks | `trsl` |
| Final rainfall | `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc` | 1 min | 2,219 CML links, two directions | `R_sublink_1min` (mm h−1), `R_mean_1min`, wet/valid flags, per-link thresholds |
| CHMI rain gauges | `CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc` | 10 min | 424 stations | `rr` (**mm per 10 min**, not mm h−1) |
| Radar along links | `CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc` | 5 min | 1,132 paths | `rr` (**mm per 5 min**, not mm h−1) |
| Air temperature | `CHMI_10MIN_TEMPERATURE_20210401_20211031.nc` | 10 min | 251 stations | `T` |
| Relative humidity | `CHMI_10MIN_HUMIDITY_20210401_20211031.nc` | 10 min | 247 stations | `H` |

The release analysis covers **April 1 00:00 through October 31 23:59 UTC, 2021**. Native slots in that interval are **308,160** (1-min), **61,632** (5-min), **30,816** (10-min). The released filenames state the coverage window (`20210401_20211031`). On the research HPC the full history source files keep their longer `20200101_20211231` names; `check` selects the exact intended window from them without changing them (the HPC rain role is the legacy full-history 5-minute reference).

**The rainfall release is native 1-minute.** The analysis 5-minute series used by the manuscript (validation, event series, gridding) are derived from it: a directional 5-minute rate exists only where **all five** 1-minute samples of the wall-clock block are valid and is their arithmetic mean; the 5-minute wet flag is set when any of those valid minutes was wet; the physical-link 5-minute mean is computed from the finite directional 5-minute rates afterwards. Notebook 01 materializes this derivation (`src/essd_cz/rain5.py`, or `python run.py --config ... derive-5min`); it was verified bit for bit against the historical validated 5-minute reference on all 2,219 links.

**Alignment:** CML `cml_id` is structural and can contain underscores; never coerce it to a floating-point number. Released 1-minute rainfall is labelled by `time_1min`; derived rainfall at the start of each 5-minute window follows `time_5min`. The gauge's 10-minute record is an *end-labelled* accumulation; the original final validation script aligns two complete CML 5-minute bins to it. For radar data, `rr × 12` converts 5-min mm to mm/h; for gauge data, `rr × 6` converts 10-min mm to mm/h. The validation script does the necessary shift/alignment explicitly; do not compare raw arrays by positional index.
