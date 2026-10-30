# Operational microwave links for hydrometeorological sensing: a multisource dataset for the Czech Republic

DOI (reserved, active when the record is published): https://doi.org/10.5281/zenodo.23398309

Authors: Taoufiq Shit¹, Martin Fencl¹, Martin Mudroch², Vojtěch Bareš¹

¹ Department of Hydraulics and Hydrology, Czech Technical University in
Prague, Prague, Czech Republic\
² Department of Electromagnetic Field, Czech Technical University in Prague,
Prague, Czech Republic

Funded by the Czech Science Foundation (GAČR), grant 24-13677L (MERGOSAT).

## What is in the deposit

All data files are packaged in one archive,
`ESSD_CZ_CML_2021_dataset.zip` (stored without recompression; the NetCDF
files are already internally compressed). Extract it to get the files
described below. `SHA256SUMS.txt` lists the checksum of every extracted
file.

### `primary/` — the seven main products, April to October 2021

Every file covers 2021-04-01 00:00 to 2021-10-31 23:5x UTC at its native
resolution.

| File | Step | Content | Licence |
|---|---|---|---|
| `CZ_CML_RAW_1MIN_20210401_20211031.nc` | 1 min | raw transmitted and received signal levels, 2,574 CMLs, 2 sublinks each | CC BY-NC 4.0 |
| `CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` | 1 min | quality controlled signals and total signal loss, 1,987 CMLs (3,731 sublinks pass quality control), with per-minute quality flags | CC BY-NC 4.0 |
| `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc` | 1 min | CML rainfall in mm/h for the 1,987 quality controlled CMLs, with wet/valid flags, source flags and per sublink thresholds | CC BY-NC 4.0 |
| `CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc` | 10 min | 424 CHMI gauges; `rr` is the amount in mm per 10 min (the variable has no units attribute) | CC BY 4.0 |
| `CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc` | 5 min | RADKLIM radar rainfall averaged along 1,132 CML paths, amount in mm per 5 min | CC BY 4.0 |
| `CHMI_10MIN_TEMPERATURE_20210401_20211031.nc` | 10 min | 251 CHMI temperature stations | CC BY 4.0 |
| `CHMI_10MIN_HUMIDITY_20210401_20211031.nc` | 10 min | 247 CHMI humidity stations | CC BY 4.0 |

The raw file contains 2,574 CMLs. After quality control, 1,987 CMLs were
retained for rainfall processing. The quality controlled file keeps the
quality control parameters in the global attribute `qc_parameters_json`;
`qc_flag` marks masked, reconstructed (bit 8) and interpolated (bit 16)
minutes, and `qc_sublink_status` gives the decision for each sublink.

The rainfall uses Graf RSTD wet/dry classification, a constant dry baseline,
the Pastorek wet antenna correction and ITU-R P.838-3 (2005) coefficients.
Thresholds and baselines are computed on 1 April to 31 October 2021. Rain
rates calculated from reconstructed or interpolated total signal loss are
marked in `trsl_source_flag_1min` (bits 8 and 16).

The rainfall release is **native 1-minute**. The 5-minute CML rainfall used
in the manuscript is derived from it: a directional 5-minute rate exists only
where all five 1-minute samples of the block are valid (their mean; wet when
any valid minute is wet), and the physical-link mean is computed from the
finite directional 5-minute rates. The accompanying software performs this
derivation (`src/essd_cz/rain5.py`). No separate 5-minute CML rainfall file
is published.

The file names state the coverage window (`20210401_20211031`). The
accompanying software expects exactly these names.

### `supplementary_wvd/` — optional water vapour support products

Needed only for the water vapour notebooks. All April to October 2021.
Licence CC BY-NC 4.0.

* `CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc` — E band water vapour retrieval,
  1,352 directional series at 5 min.
* `CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv` — one diagnostics row per
  series.
* `CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv` — one quantization resolution
  row per series.

The selection used in the paper is 45 directional sublinks from 31 physical
CMLs (E band, length at least 4 km, quantization at most 0.20 g per cubic
metre).

Native RADKLIM YW radar fields are not in this deposit. Get them from the DWD
Climate Data Center under its own terms.

### Checksums and dictionary

* `SHA256SUMS.txt` — one checksum per file. Verify after download.
* `DATA_DICTIONARY.md` — the variables of every file, generated from the
  files themselves.

## Conventions

Czech Republic, WGS84 coordinates. Time axes are UTC without a timezone.
Rainfall 5 min bins are rates in mm/h labelled at the bin start. Gauge 10 min
values are amounts labelled at the bin end. CML path lengths in km.

## Known limitations

1. There was no calibration of the CML rainfall against gauges or radar. The
   thresholds and dry baselines are derived from the CML data.
2. Rain rates from total signal loss reconstructed during blackouts are
   uncertain; use `trsl_source_flag_1min` to exclude them if needed. Some very
   short links give very high rain rates; these values are kept.
3. CML availability (at least one sublink with both transmitted and received
   signal level) averages 67.8 percent of 1 min slots over the 2,574 raw
   CMLs. Gaps reflect interruptions of the research data acquisition, not
   shutdown of the links. The reference sensors are above 99 percent.
4. Missing values must not be read as zero rainfall or dry conditions.

## Attribution

* CML signal data provided by T-Mobile Czech Republic a.s.
* Gauge, temperature and humidity data: Czech Hydrometeorological Institute
  (CHMI) open data, CC BY 4.0.
* Path radar rainfall derived from Deutscher Wetterdienst (DWD) RADKLIM YW,
  CC BY 4.0.

## Use with the software

The software repository is https://github.com/toufikshit/essd-czech-cml-2021.
It walks from these files to the paper figures in notebooks. Setup is one
config file; see its README.

## Record structure

One Zenodo record with per file licences, as listed in the tables above. The
record licence field is set to other (open) because the terms differ per
file.
