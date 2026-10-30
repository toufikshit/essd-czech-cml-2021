# Getting the data

The data are not in this repository. They are in the Zenodo record
https://doi.org/10.5281/zenodo.23398309. The DOI is reserved; the files can be
downloaded when the record is published.

Download `ESSD_CZ_CML_2021_dataset.zip` and extract it, for example into this
`data/` folder. `config/user_paths.example.json` already points to
`data/primary/` and `data/supplementary_wvd/`; with the data elsewhere, set the
paths in `config/user_paths.json`. Do not rename the files; the code checks the
names.

The release contains raw signals for 2,574 CMLs, quality controlled signals and
rainfall for 1,987 CMLs, and reference data from CHMI and DWD, 1 April to 31
October 2021. The CML rainfall is at 1 minute; the 5-minute series of the paper
are derived from it with `scripts/derive_5min.py`.

## Files and checksums (SHA-256)

| File | SHA-256 |
|---|---|
| `primary/CZ_CML_RAW_1MIN_20210401_20211031.nc` | `7f7845458d0486622517ae3b2ddfe6a3c6a3ceb7d37f8d179e1224cf674e8d4a` |
| `primary/CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` | `bc2e0c1428f99db9826938df5b953204b9351225d208cd7640acfc8929f4584e` |
| `primary/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc` | `a09251c1582057d697a12319bfbc8e6787e17088a00fd1f461a4422dfb10a5f6` |
| `primary/CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc` | `38dae054ec232f41653a095a5a19bc29c8eb4aa844a237cbd8615c5e0f14262a` |
| `primary/CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc` | `e82d446b8970e6cc20e0448e8875346ef2d0c84564600f71a039f102ed9c1e2c` |
| `primary/CHMI_10MIN_TEMPERATURE_20210401_20211031.nc` | `6a0bb785a9569c7f5d4ec9ab34e7a89622c9fc8d67a93bc3dd4a8553e1a3d3d4` |
| `primary/CHMI_10MIN_HUMIDITY_20210401_20211031.nc` | `221a14c19c5be0470be3d298129700d77cad085d853eda8776163ed8b8537e5b` |
| `supplementary_wvd/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc` | `35bc67ed1b905416bdceb4d44666f706a22ee36e8ba7457f551a628c05b4c4bf` |
| `supplementary_wvd/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv` | `cab58ca361e506b5e45b55ec95851d14a42c5c266bba34db2f29f0107d97fbf7` |
| `supplementary_wvd/CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv` | `67012dca9ded4451723ed2a902cfbbee404b1c5c3aa81e326453549cfd3e498c` |

The zip has 6,368,237,019 bytes and SHA-256
`74445851ced5e5db730f8034d6dc5dd969bf8b11944acd5f5904709a88eabdb1`. After
extracting, check the files with `sha256sum -c SHA256SUMS.txt`, then run
`python scripts/check_files.py`.

## Data not in the deposit

The radar processing (`scripts/radar_cml_paths.py`) and the radar parts of Figs. 1
and 7 need daily RADKLIM YW 2017.002 netCDF files (`YW_2017.002_YYYYMMDD.nc`) of
Deutscher Wetterdienst, Climate Data Center OpenData. Set their folder as
`extras.radklim_dir`, or download them with `scripts/radar_cml_paths.py --download`.

## Licences

CML files are CC BY-NC 4.0 (noncommercial use only). CHMI and DWD derived files
are CC BY 4.0. Details: `THIRD_PARTY_NOTICES.md` in the repository root.
