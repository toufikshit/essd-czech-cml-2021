# Licence state

Short version: MIT is proposed for the project software. The CML data files
are CC BY-NC 4.0. The CHMI and DWD derived data files are CC BY 4.0. The full
file to licence mapping and the attribution texts are in
`THIRD_PARTY_NOTICES.md`.

## Software: MIT (proposed)

The `LICENSE` file carries the MIT licence for the project software: `src/`,
`tests/`, `notebooks/`, `run.py`, `tools_verify_originals.py`, the docs, and
the preserved research code under `original/`.

Two confirmations are still needed before the repository goes public:

1. Code ownership. Confirm with the Czech Technical University in Prague that
   the authors may release the project code under MIT, and confirm who the
   copyright holder line should name.
2. Third party content inside `original/`. The preserved scripts adapt
   published methods and use pycomlink. Nothing in them is a copy of someone
   else's source file as far as the provenance record shows, but confirm this
   before publishing. Never relicense third party code under MIT.

MIT never applies to the NetCDF data files.

## Data

| Group | Licence | Basis |
|---|---|---|
| CML files (raw, preprocessed, rainfall, WVD supplement) | CC BY-NC 4.0 | agreement with the network operator, T-Mobile Czech Republic a.s. Noncommercial use only. Final wording of the attribution line needs operator confirmation. |
| CHMI files (gauges, temperature, humidity) | CC BY 4.0 | CHMI official open data policy, https://opendata.chmi.cz |
| DWD derived path rainfall | CC BY 4.0 | DWD Climate Data Center OpenData terms for RADKLIM |

A single Zenodo record must state the licence per file, because CC BY and
CC BY-NC files sit in the same deposit. The data README and the data
dictionary carry the per file licence.

## Still open before publication

- [ ] CTU confirmation of MIT for the project code and the copyright line.
- [ ] T-Mobile attribution wording for the CML files.
- [ ] Final check of the `original/` scripts for third party source code.
- [ ] Decision: one Zenodo record with per file licences, or separate records.
