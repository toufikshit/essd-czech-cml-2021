# 2. Seven data products — quick dictionary

| Role | File basename | Time step | Entities | Example variable |
|---|---|---|---:|---|
| Raw telecom | `CZ_CML_RAW_1MIN_20200101_20211231.nc` | 1 min | 2,580 CML links | link signal |
| QC telecom | `CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc` | 1 min | 2,219 CML links, two directional sublinks | `trsl` |
| Final rainfall | `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_20200101_20211231.nc` | 5 min | 2,219 CML links, two directions | `R_sublink_5min` (mm h−1), `R_mean_5min` |
| CHMI rain gauges | `CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc` | 10 min | 424 stations | `rr` (**mm per 10 min**, not mm h−1) |
| Radar along links | `CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc` | 5 min | 1,132 paths | `rr` (**mm per 5 min**, not mm h−1) |
| Air temperature | `CHMI_10MIN_TEMPERATURE_2021.nc` | 10 min | 251 stations | `T` |
| Relative humidity | `CHMI_10MIN_HUMIDITY_2021.nc` | 10 min | 247 stations | `H` |

The release analysis covers **April 1 00:00 through October 31 23:59 UTC, 2021**. Native slots in that interval are **308,160** (1-min), **61,632** (5-min), **30,816** (10-min). The filenames retain `20200101_20211231` from upstream processing, so **do not infer actual coverage from filenames**. On the research HPC the files were still longer source products; `check` selects the exact intended window without changing them. **This is not yet proof the final seven physically cropped files are deposited.**

**Alignment:** CML `cml_id` is structural and can contain underscores; never coerce it to a floating-point number. Rainfall at the start of each 5-minute window follows `time_5min`. The gauge's 10-minute record is an *end-labelled* accumulation; the original final validation script aligns two complete CML 5-minute bins to it. For radar data, `rr × 12` converts 5-min mm to mm/h; for gauge data, `rr × 6` converts 10-min mm to mm/h. The validation script does the necessary shift/alignment explicitly; do not compare raw arrays by positional index.
