# Third party notices and attribution

This file maps everything in and around this project to its licence and its
required attribution. The MIT licence in `LICENSE` covers only the project
software. It never covers the data files.

## File to licence mapping

| What | Licence | Notes |
|---|---|---|
| `src/`, `tests/`, `notebooks/`, `run.py`, `tools_verify_originals.py`, `docs/` | MIT (proposed) | project software and documentation written for this repository |
| `original/` (44 preserved scripts and notebooks) | MIT (proposed, pending confirmation) | preserved project research code. Some scripts adapt methods from published work (see below). Confirm ownership before making the repository public. Do not relicense any part that turns out to be third party code. |
| The seven NetCDF data files | CC BY-NC 4.0 or CC BY 4.0 per file | distributed through Zenodo, never through this repository. See the table below and `docs/LICENCES.md`. |
| Python dependencies (numpy, pandas, xarray, netCDF4, scipy, matplotlib, numba, pycomlink, cartopy, requests, beautifulsoup4, tqdm) | each project's own licence | installed from PyPI, not redistributed here |

## Data files and attribution

| Data file | Licence | Attribution |
|---|---|---|
| `CZ_CML_RAW_1MIN_*.nc` | CC BY-NC 4.0 | CML signal data provided by T-Mobile Czech Republic a.s. Noncommercial use only. |
| `CZ_CML_PREPROCESSED_1MIN_*.nc` | CC BY-NC 4.0 | derived from T-Mobile Czech Republic CML data |
| `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_*.nc` | CC BY-NC 4.0 | derived from T-Mobile Czech Republic CML data |
| `CZ_RAIN_GAUGE_10MIN_*.nc` | CC BY 4.0 | Czech Hydrometeorological Institute (CHMI) open data, https://opendata.chmi.cz |
| `CHMI_10MIN_TEMPERATURE_2021.nc` | CC BY 4.0 | CHMI open data |
| `CHMI_10MIN_HUMIDITY_2021.nc` | CC BY 4.0 | CHMI open data |
| `CZ_RADAR_CMLPATH_5MIN_*.nc` | CC BY 4.0 | derived from Deutscher Wetterdienst (DWD) RADKLIM YW radar data, DWD Climate Data Center OpenData |
| WVD V6 support files (optional Zenodo supplement) | CC BY-NC 4.0 | derived from T-Mobile Czech Republic CML data and CHMI station data |

Native RADKLIM YW radar fields are not redistributed by this project. Get
them from the DWD Climate Data Center under its own terms.

## Scientific methods used by the preserved code

The preserved scripts implement or adapt published methods. The methods
themselves are cited in the paper: the Graf wet and dry detection scheme, the
Pastorek wet antenna correction, ITU R P.838 (2005) rain attenuation
coefficients, and processing steps built on the open source pycomlink
library. The pycomlink library is used as an installed dependency under its
own licence; it is not copied into this repository.

## Funding

Czech Science Foundation (GAČR), grant 24-13677L (MERGOSAT).
