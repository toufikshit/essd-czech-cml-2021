# ESSD Czech CML 2021

Code and notebooks for the dataset paper *Operational microwave links for
hydrometeorological sensing: a multisource dataset for the Czech Republic*
(Shit, Fencl, Mudroch and Bareš; Earth System Science Data, in preparation).

The repository contains the Python scripts and notebooks used to process the data
and produce the figures. The README explains how to install the required packages,
set the data paths, and run the code.

> Dataset: Zenodo, https://doi.org/10.5281/zenodo.23398309. The DOI is reserved;
> the files can be downloaded when the record is published. The article DOI will
> be added once the paper has one.

## 1. Download the data

Download `ESSD_CZ_CML_2021_dataset.zip` from the Zenodo record and extract it, for
example into the `data/` folder of this repository. It contains:

| Folder / file | Content |
|---|---|
| `primary/CZ_CML_RAW_1MIN_20210401_20211031.nc` | raw TX and RX, 2,574 CMLs, 1 min |
| `primary/CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` | quality controlled signals, 1,987 CMLs, 1 min |
| `primary/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc` | CML rainfall, 1,987 CMLs, 1 min |
| `primary/CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc` | 424 CHMI rain gauges, 10 min |
| `primary/CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc` | RADKLIM radar along 1,132 CML paths, 5 min |
| `primary/CHMI_10MIN_TEMPERATURE_20210401_20211031.nc` | 251 CHMI stations, 10 min |
| `primary/CHMI_10MIN_HUMIDITY_20210401_20211031.nc` | 247 CHMI stations, 10 min |
| `supplementary_wvd/` | three water vapour (WVD) supporting files |
| `README_DATA.md`, `DATA_DICTIONARY.md`, `SHA256SUMS.txt` | documentation and checksums |

Details and checksums: [data/README.md](data/README.md). Variables:
[docs/DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md).

Two steps need data that are not in the deposit: the radar processing and the
radar parts of Figs. 1 and 7 need daily RADKLIM YW 2017.002 files from the DWD
Climate Data Center (`scripts/radar_cml_paths.py --download` gets them).

## 2. Install

Python 3.10 or newer.

```bash
git clone https://github.com/toufikshit/essd-czech-cml-2021
cd essd-czech-cml-2021
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[all]"
```

`pip install -e .` installs only the packages needed to read the data. The groups
`science` (numba, pycomlink), `figures` (cartopy, pillow), `radar` (poligrain),
`wvd` (itur) and `notebooks` can be installed separately, for example
`pip install -e ".[science,notebooks]"`. Tested with Python 3.12, numpy 2.5,
xarray 2026.9, netCDF4 1.7, numba 0.68, pycomlink 0.6.0, poligrain 0.3.1 and
cartopy 0.26.

## 3. Set the data paths

```bash
cp config/user_paths.example.json config/user_paths.json
```

Edit `config/user_paths.json`: the paths of the seven files, the three WVD files
(`extras`), a folder with RADKLIM files (`extras.radklim_dir`, optional) and the
output folder (`output_dir`, default `results/`). Relative paths are relative to
the `config/` folder; absolute paths also work. The example file already matches
the zip extracted into `data/`. Every script reads this file; another file can be
given with `--config`. Check the paths and files:

```bash
python scripts/check_files.py
```

All seven lines must say PASS.

## 4. Run the code

Each script starts from the released files and writes into `output_dir`. Run them
from the repository folder. `--help` shows the options of every script.

| Order | Script | Input | Output (in `output_dir`) |
|---:|---|---|---|
| 1 | `scripts/qc_from_raw.py` | RAW file | `qc/CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc`, `qc/QC_SUBLINK_RESULTS.csv` |
| 2 | `scripts/rainfall_1min.py` | QC file | `rainfall/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc` |
| 3 | `scripts/derive_5min.py` | 1-minute rainfall | `rainfall_5min/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031.nc` |
| 4 | `scripts/gauge_validation.py` | QC file, 5-minute rainfall, gauges | `validation/graf_summary.json` (statistics), `graf_samples.npz` |
| 5 | `scripts/radar_cml_paths.py` | radar file (path geometry), RADKLIM files | `radar/cml_rr_RADKLIM_2021_CZ.nc` |
| 6 | `scripts/wvd_validation.py` | WVD files, RAW, temperature, humidity | `wvd/WVP_SUPPORT_HEXBINS/` (Table S13, Fig. 8) |
| 7 | `scripts/figures/fig01_study_area.py` | RAW, stations, RADKLIM files | `figures/f01_study_area.png` |
|  | `scripts/figures/fig02_availability.py` | seven files | `figures/fig02/availability_native_timestep_percentage_v1.png` |
|  | `scripts/figures/fig03_cml_characteristics.py` | RAW file | `figures/f03_cml_characteristics.png` |
|  | `scripts/figures/fig04_event_series.py` | QC, 5-minute rainfall, radar, stations | `figures/fig04/FINAL_GRAF_event_20210713_20210714.png` |
|  | `scripts/figures/fig05_fig06_validation.py` | validation samples (step 4) | `figures/fig05_fig06/` |
|  | `scripts/figures/fig07_event_maps.py` | 5-minute rainfall, QC, RADKLIM file of 13 July 2021 | `figures/fig07/hourly_1h/event04_instance01_multitime_hourly1h.png` |

With 12 threads the QC takes about 30 minutes, the rainfall 8 minutes, the gauge
validation 1 hour and Fig. 1 15 minutes; the other steps take minutes, except the
radar processing, which depends on the download. By default, step 2 uses the released QC
file and step 3 the released 1-minute rainfall. To run the whole chain from the
raw data, pass the output of the previous step:

```bash
python scripts/qc_from_raw.py
python scripts/rainfall_1min.py --input-nc results/qc/CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc
python scripts/derive_5min.py --native results/rainfall/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc
python scripts/gauge_validation.py
```

`scripts/qc_from_raw.py --only-cmls id1,id2` and `scripts/rainfall_1min.py
--only-first N` process a few CMLs only, for a quick test.

### Notebooks

The notebooks in `notebooks/` run the same scripts and show their results.

| Notebook | Content |
|---|---|
| `00_Read_Data.ipynb` | open the seven files, plot signals, rainfall and the network |
| `01_QC_and_Rainfall.ipynb` | QC and rainfall for a few CMLs, compared with the released files |
| `02_Rainfall_Validation.ipynb` | 5-minute rainfall, gauge validation, Figs. 5 and 6 |
| `03_Figures.ipynb` | Figs. 1 to 4 and 7 |
| `04_Radar_Processing.ipynb` | radar rainfall along the CML paths |
| `05_Water_Vapour.ipynb` | WVD validation and Fig. 8 |

## 5. Repository content

| Folder | Content |
|---|---|
| `scripts/` | processing, validation and figure scripts |
| `src/essd_cz/` | data paths (`config.py`), file check (`catalog.py`), 5-minute derivation (`rain5.py`) |
| `original/code/` | unchanged scientific scripts called by the scripts: figure producers, radar processing, validation helper, WVD retrieval, CHMI converter |
| `config/` | path template, population tables used by the QC, inputs of Fig. 7 |
| `notebooks/` | notebooks 00 to 05 |
| `docs/` | methods and settings ([docs/METHODS.md](docs/METHODS.md)), data dictionary, licences |
| `tests/` | tests on small artificial files (`python -m pytest tests`) |

`python tools_verify_originals.py` checks that the scripts in `original/` are
unchanged.

## Licences

Software: MIT ([LICENSE](LICENSE)). Data: CML files CC BY-NC 4.0 (noncommercial
use only; provided under an agreement with T-Mobile Czech Republic a.s.), CHMI and
DWD derived files CC BY 4.0. Details: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)
and [docs/LICENCES.md](docs/LICENCES.md).

## Citation

See [CITATION.cff](CITATION.cff). Software by Taoufiq Shit. Dataset and paper by
Taoufiq Shit, Martin Fencl, Martin Mudroch and Vojtěch Bareš (Czech Technical
University in Prague). Please cite both the paper and the dataset. Funded by the
Czech Science Foundation (GAČR), grant 24-13677L (MERGOSAT).
