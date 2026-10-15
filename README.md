# ESSD Czech CML 2021 — connected reviewer workflow

Code, notebooks and preserved original scripts for the Czech 2021 commercial
microwave link (CML) rainfall and water-vapour dataset described in the
accompanying ESSD manuscript. Nine connected notebooks walk from the seven
released NetCDF products to rainfall reconstruction, gauge validation, spatial
gridding, radar map panels and E-band water-vapour (WVD) analysis, reusing the
**byte-identical original scientific scripts** preserved under [`original/`](original/).

> **Status:** release candidate prepared for author review. The software
> licence, data licences and Zenodo deposit still require explicit author
> sign-off; see [LICENSING.md](LICENSING.md) and
> [docs/PROVENANCE_MAP.md](docs/PROVENANCE_MAP.md) before distributing.

## 1. Download

* Clone this repository.
* Download the seven released NetCDF data products separately (Zenodo record
  pending; **do not rename the files**):
  `CZ_CML_RAW_1MIN_*.nc`, `CZ_CML_PREPROCESSED_1MIN_*.nc`,
  `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_*.nc`, `CZ_RAIN_GAUGE_10MIN_*.nc`,
  `CZ_RADAR_CMLPATH_5MIN_*.nc`, `CHMI_10MIN_TEMPERATURE_2021.nc`,
  `CHMI_10MIN_HUMIDITY_2021.nc`.
* Optional supporting products (WVD V6 NetCDF and diagnostics/quantization
  CSVs, native RADKLIM YW radar) are separate downloads; notebooks 04–07 say
  exactly which they need.

## 2. Configure (one file)

```bash
cp config/user_paths.example.json config/user_paths.json
```

Edit the seven dataset paths in `config/user_paths.json` (absolute, or relative
to `config/`). `config/active_profile.json` already selects the portable
`"user"` profile. The author's HPC layout remains available as an example in
`config/hpc_source.example.json`. Worked example: [`examples/`](examples/).

## 3. Install

Python ≥ 3.10 with `requirements.txt` (core) plus `requirements-notebooks.txt`
(notebook kernel). `numba`/`pycomlink` are needed only for actual rainfall
recomputation; `cartopy` only for notebook 04 map rendering. Versions this
candidate was re-tested with are pinned in
[docs/ENVIRONMENT.md](docs/ENVIRONMENT.md).

```bash
python -m venv .venv && . .venv/bin/activate   # or your conda env
pip install -r requirements.txt -r requirements-notebooks.txt
```

## 4. Inspect before running anything

```bash
python tools_verify_originals.py            # 44 preserved originals, byte-identical
python -m pytest -q tests                   # 20 synthetic regression tests
python run.py --config config/user_paths.json check   # file roles + time coverage (metadata only)
```

## 5. Run stage 00, then the stages in order

Open `notebooks/00_Explore_Data.ipynb` and run it top to bottom, then proceed
`01 → 07` (08 is optional upstream CHMI preparation). Each notebook registers
its real output in `results/<workspace>/pipeline_manifest.json`; the next
notebook reads that registered artifact — never a hand-typed path. The stage
table, inputs/outputs and caveats are in
[README_CONNECTED_WORKFLOW.md](README_CONNECTED_WORKFLOW.md).

Expensive or external-dependency steps are **off by default** and clearly
labelled in each notebook (`RECOMPUTE`, `RUN_GAUGE_VALIDATION`, `RUN_MAPS`,
`RUN_SENSITIVITY`). Approximate costs on the author's HPC: full 2,219-link
rainfall recompute ≈ hours (4 workers); gauge validation and the 1,352-series
WVD sensitivity pass are also expensive; everything else runs in seconds to
minutes.

## 6. Troubleshooting

* `wrong filename` from `run.py check` — the seven files must keep their
  published names; fix the paths, do not rename data.
* `No output registered for 'X'` — run the previous notebook first; stages
  hand artifacts forward through the manifest.
* `This output manifest was created with another configuration` — change
  `workspace` in `config/active_profile.json` instead of reusing another
  profile's results directory.
* NetCDF/engine errors on machines without `netCDF4` — the tests fall back to
  NETCDF3 fixtures; the real data requires `netCDF4`.
* `cartopy` missing — needed only for notebook 04; install it in the notebook
  environment.

## 7. Data licences

Intended terms (pending verification and author sign-off, see
[LICENSING.md](LICENSING.md)): CML-derived products CC BY-NC 4.0; CHMI station
products and DWD-derived radar products CC BY 4.0. The native RADKLIM YW radar
archive is **not** redistributed here.

## 8. Reproducibility status (honest summary)

* Preserved originals: 44/44 byte-identical (`tools_verify_originals.py`).
* Real-data Graf rainfall: full 2,219-link recompute executed on the author's
  HPC (October 2026); one-link exact comparison passed with zero numerical
  error. See [docs/EVIDENCE.md](docs/EVIDENCE.md).
* The 5-min IDW grid is an **adapted** reviewer demonstration using unchanged
  original IDW functions — not asserted byte-identical to the historic 1-min
  manuscript grid.
* The exact historical gauge wet/dry detector and the March 2026 QC producer
  remain **unrecovered**; substitutes are explicitly labelled unverified.
* Full gap list: [docs/06_PROVENANCE_GAPS.md](docs/06_PROVENANCE_GAPS.md) and
  [docs/PROVENANCE_MAP.md](docs/PROVENANCE_MAP.md).

Large NetCDF data, native radar trees and heavy derived outputs are never
committed to git (see `.gitignore`); they belong in the Zenodo deposit.
