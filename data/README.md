# Getting the data

The data is not in this repository and never will be. It is distributed as a
Zenodo record (DOI to be added when the data release is published).

The next release contains the corrected products of October 2026: raw
signals for 2,574 CMLs, quality controlled signals and rainfall for 1,987
CMLs. The planned files, sizes and checksums are listed in
`reconstruction_2026/DATA_RELEASE_MANIFEST.md`.

## The seven main files

Download all seven and put them in one folder. Do not rename them. The code
checks the names.

```
CZ_CML_RAW_1MIN_20210401_20211031.nc
CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc
CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc
CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc
CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc
CHMI_10MIN_TEMPERATURE_20210401_20211031.nc
CHMI_10MIN_HUMIDITY_20210401_20211031.nc
```

The file names state the coverage: 1 April to 31 October 2021.

The CML rainfall product is **native 1-minute**. Every 5-minute rainfall
series used in the notebooks (validation, event series, gridding) is derived
from it on first use — mean of five complete valid 1-minute samples per
directional sublink — so a separate 5-minute rainfall download is not needed.
The next release may also include this derived 5-minute file for
convenience; it holds the same values. See `docs/03_RAINFALL_EXPLAINED.md`.

Then copy `config/user_paths.example.json` to `config/user_paths.json` and
set the seven paths. For example, with the data in `data/` inside this
repository:

```json
"raw": "../data/CZ_CML_RAW_1MIN_20210401_20211031.nc"
```

Paths are relative to the `config/` folder. Absolute paths also work. You can
keep the data anywhere, for example on a network drive.

## Optional files

Notebooks 05 to 07 (water vapour) can use three support files from the same
deposit. Set their paths in the `extras` block of `config/user_paths.json`:

```
CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc        -> extras.wvd_v6_nc
CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv -> extras.wvd_v6_diagnostics_csv
CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv      -> extras.wvd_quantization_csv
```

Notebook 04 (radar maps) needs native RADKLIM YW radar files. They are not in
the deposit. Download them from the DWD Climate Data Center and set
`extras.native_radklim_root`.

## Check your download

After configuring, run:

```bash
python run.py --config config/user_paths.json check
```

All seven lines must say PASS. Verify the file checksums against the
`SHA256SUMS.txt` published in the Zenodo deposit:

```bash
sha256sum -c SHA256SUMS.txt
```

## Licences

CML files are CC BY-NC 4.0 (noncommercial use only). CHMI and DWD derived
files are CC BY 4.0. Details: `THIRD_PARTY_NOTICES.md` in the repository root.
