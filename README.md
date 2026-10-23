# ESSD Czech CML 2021

Code and notebooks for the dataset paper *Operational microwave links for
hydrometeorological sensing: a multisource dataset for the Czech Republic*
(Shit, Fencl and Bareš; Earth System Science Data, in preparation).

The dataset contains observations from 2,219 commercial microwave links (CML)
and reference sensors in the Czech Republic for April to October 2021. Nine
notebooks walk from the released NetCDF files to rainfall retrieval, gauge
validation, rainfall maps and E band water vapour analysis. The original
scientific scripts are preserved unchanged under `original/` and the notebooks
call them directly.

> Status: the dataset is not yet on Zenodo. The Zenodo link and the DOIs will
> be added here after the deposit is published.

## Quick start

1. **Get the code.** Clone this repository.

2. **Get the data.** Download the seven NetCDF files (Zenodo link will appear
   here). Put them in one folder. Do not rename them. Details: [data/README.md](data/README.md).

3. **Install.** Python 3.10 or newer.

   ```bash
   pip install -e .
   ```

   This installs the core packages. Extras: `.[science]` for rainfall
   recomputation (numba, pycomlink), `.[chmi]` for the CHMI download notebook,
   `.[test]` for the test suite. Notebook 04 also needs cartopy.

4. **Point the code at your data.** Copy the example and edit the seven paths:

   ```bash
   cp config/user_paths.example.json config/user_paths.json
   ```

   Paths can be absolute or relative to `config/`. This is the only file you
   edit. No notebook contains data paths.

5. **Check before running.**

   ```bash
   python tools_verify_originals.py
   python run.py --config config/user_paths.json check
   ```

   The first command confirms the 44 preserved original scripts are unchanged.
   The second confirms your seven files have the right roles, entity counts
   and time coverage.

6. **Run the notebooks in order.** Start with
   `notebooks/00_Explore_Data.ipynb` and continue 01 to 07. Notebook 08 is an
   optional upstream step that rebuilds the CHMI station files from CHMI open
   data. Each notebook registers its output in a pipeline manifest and the
   next notebook reads it from there. See
   [docs/NOTEBOOK_GUIDE.md](docs/NOTEBOOK_GUIDE.md) for inputs, outputs and
   cost per notebook, and [docs/CONNECTED_WORKFLOW.md](docs/CONNECTED_WORKFLOW.md)
   for the full stage table.

Expensive steps are off by default and marked in each notebook (`RECOMPUTE`,
`RUN_GAUGE_VALIDATION`, `RUN_MAPS`, `RUN_SENSITIVITY`). The default first pass
of every notebook runs in seconds to minutes.

## What is in the data

| File | Step | Content |
|---|---|---|
| `CZ_CML_RAW_1MIN_20210401_20211031.nc` | 1 min | raw signal levels, 2,580 CMLs |
| `CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` | 1 min | quality controlled signals, 2,219 CMLs |
| `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc` | 1 min | final CML rainfall in mm/h, native 1-minute |
| `CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc` | 10 min | 424 CHMI rain gauges; `rr` is the amount in mm per 10 min |
| `CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc` | 5 min | RADKLIM radar rainfall averaged along 1,132 CML paths, amount in mm per 5 min |
| `CHMI_10MIN_TEMPERATURE_20210401_20211031.nc` | 10 min | 251 CHMI temperature stations |
| `CHMI_10MIN_HUMIDITY_20210401_20211031.nc` | 10 min | 247 CHMI humidity stations |

The rainfall retrieval uses CML data only, without gauge or radar
calibration. The method is described in the manuscript and in
[docs/CONNECTED_WORKFLOW.md](docs/CONNECTED_WORKFLOW.md).

The CML rainfall product is released at its **native 1-minute** resolution.
Every 5-minute rainfall series used in the manuscript analysis (validation,
event series, gridding) is derived from it by the notebooks: a directional
5-minute rate exists only where all five 1-minute samples are valid, and the
physical-link mean is computed from the directional 5-minute rates. See
[docs/03_RAINFALL_EXPLAINED.md](docs/03_RAINFALL_EXPLAINED.md).

## Reproducibility

The released rainfall was reproduced exactly: a full recomputation of all
2,219 links with the preserved original script matched the validated
reference in every value and every mask, and the released native 1-minute
file reproduces that same 5-minute reference bit for bit under the documented
five-complete-samples aggregation. The verification numbers and proof files
are in [docs/EVIDENCE.md](docs/EVIDENCE.md), and the provenance of every
piece of code is in [docs/PROVENANCE_MAP.md](docs/PROVENANCE_MAP.md).

## Troubleshooting

* `wrong filename` from `run.py check`: fix the paths in
  `config/user_paths.json`. Do not rename the data files.
* `No output registered for 'X'`: run the previous notebook first.
* `This output manifest was created with another configuration`: change
  `workspace` in `config/active_profile.json` instead of reusing another
  output folder.
* Notebook 04 needs cartopy and a native RADKLIM radar folder. Both are
  optional.

## Licences

The software is MIT. The CML data files are CC BY-NC 4.0, provided under an
agreement with T-Mobile Czech Republic (available on request). The CHMI and
DWD derived data files are CC BY 4.0. Per file details:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Citation

See [CITATION.cff](CITATION.cff). Software by Taoufiq Shit. Dataset and
paper by Taoufiq Shit, Martin Fencl and Vojtěch Bareš (Faculty of Civil
Engineering, Czech Technical University in Prague). Please cite both the
paper and the dataset once their DOIs exist. Funded by the Czech Science
Foundation (GAČR), grant 24-13677L (MERGOSAT).
