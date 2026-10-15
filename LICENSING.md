# Licensing status — ACTION REQUIRED BEFORE PUBLICATION

This repository deliberately ships **without** a `LICENSE` file. Nothing here
may be treated as openly licensed until the rights holders decide and a licence
file is committed.

## Software (this repository)

| Component | Status |
|---|---|
| `src/essd_cz/`, `tests/`, `notebooks/` (new wrappers and notebooks) | **No licence selected.** The author/rights holder must choose one (e.g. MIT, BSD-3-Clause, Apache-2.0, EUPL) before the repository is made public. |
| `original/` (44 preserved scientific scripts and notebooks) | Same rights holder decision required; preserved byte-identical for provenance. |
| Third-party dependencies (`pycomlink`, `numba`, …) | Their own licences; not redistributed here. |

## Data (Zenodo deposit — separate from this repository)

Intended terms recorded in the project's own provenance notes
(`original/README_PROVENANCE_AND_LICENSING.txt`), **not yet verified**:

| Product group | Intended licence | Verification needed before deposit |
|---|---|---|
| CML-derived products (raw, QC, rainfall, WVD support) | CC BY-NC 4.0 | Operator agreement permits publication under these terms; confirm attribution text. |
| CHMI station products (gauges, temperature, humidity) | CC BY 4.0 | CHMI open-data terms and required attribution wording. |
| DWD-derived path-averaged radar product | CC BY 4.0 | DWD RADKLIM attribution and derived-product terms. |
| Native RADKLIM YW fields | **not redistributed** | Any future redistribution needs an explicit DWD rights check. |

A single Zenodo record must not blend CC BY and CC BY-NC assets without
per-file licence metadata; the prepared deposit keeps the licence per product
in its data dictionary and README.

## Sign-off checklist (written author approval required for each)

- [ ] Software licence chosen and `LICENSE` committed.
- [ ] CML CC BY-NC 4.0 confirmed with the network operator agreement.
- [ ] CHMI CC BY 4.0 attribution wording confirmed.
- [ ] DWD-derived radar product terms and attribution confirmed.
- [ ] Contributor names/affiliations and funding statements confirmed
      (`CITATION.cff` placeholders filled).
- [ ] Manuscript/dataset cross-citation (DOIs) inserted once known.
