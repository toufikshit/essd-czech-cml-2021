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

> Status: release candidate. The dataset is not yet on Zenodo and this
> repository is not yet public. See `docs/LICENCES.md` for the licence state.

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
| `CZ_CML_RAW_1MIN_*.nc` | 1 min | raw signal levels, 2,580 CMLs |
| `CZ_CML_PREPROCESSED_1MIN_*.nc` | 1 min | quality controlled signals, 2,219 CMLs |
| `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_*.nc` | 5 min | final CML rainfall in mm/h |
| `CZ_RAIN_GAUGE_10MIN_*.nc` | 10 min | 424 CHMI rain gauges; `rr` is the amount in mm per 10 min |
| `CZ_RADAR_CMLPATH_5MIN_*.nc` | 5 min | RADKLIM radar rainfall averaged along 1,132 CML paths, amount in mm per 5 min |
| `CHMI_10MIN_TEMPERATURE_2021.nc` | 10 min | 251 CHMI temperature stations |
| `CHMI_10MIN_HUMIDITY_2021.nc` | 10 min | 247 CHMI humidity stations |

The rainfall method is Graf wet and dry detection with a constant baseline,
Pastorek wet antenna correction and ITU R P.838 (2005) coefficients. The
thresholds and baselines are derived from the CML data themselves. There was
no calibration of the CML rainfall against gauges or radar.

## Reproducibility in one paragraph

The released rainfall file was reproduced exactly: a full recomputation of
all 2,219 links with the preserved original script matched the released file
in every value and every mask (137.8 million rain values compared, zero
differences). The published rainfall was computed from the full 2020 to 2021
signal history. The released April to October data alone cannot rebuild the
wet and dry thresholds, so the per link thresholds are stored inside the
rainfall file. The water vapour analysis uses 1,352 directional series, and
the selection used in the paper is 45 directional sublinks from 31 physical
CMLs. Some historical steps are not recovered; see
[docs/PROVENANCE_MAP.md](docs/PROVENANCE_MAP.md) for the honest list and
[docs/EVIDENCE.md](docs/EVIDENCE.md) for the numbers and proof files.

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

Proposed and partly confirmed, see [docs/LICENCES.md](docs/LICENCES.md) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md): MIT for the project code
(pending ownership confirmation), CC BY-NC 4.0 for the CML data files,
CC BY 4.0 for the CHMI and DWD derived data files.

## Citation

See [CITATION.cff](CITATION.cff). Please cite both the paper and the dataset
once their DOIs exist.
