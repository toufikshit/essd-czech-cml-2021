# tools/release_prep/

Scripts used to assemble and audit this release candidate from the verified
`ESSD_CZ_2021_FULL_HANDOFF.tar` archive (October 2026). They are preparation
tooling, not part of the scientific workflow:

* `verify_manifest.py` — verify an extracted handoff against
  `handoff/FILE_MANIFEST_SHA256.json` (257/257 passed).
* `build_github_repo.py` — mechanical, hash-verified assembly of this
  repository from the handoff tree (see `CHANGES_FROM_HANDOFF.md`).
* `audit_rainfall_identity.py` — full element-wise comparison of the published
  vs HPC-recomputed rainfall product (result: identical).
* `audit_wvd_support.py` — WVD V6 series-count and index-consistency audit.
* `zenodo_crop.py` — physical Apr–Oct 2021 crop of one product at native
  resolution, preserving values/encoding/metadata.
* `zenodo_validate.py` — systematic source-vs-cropped validation
  (bit-identity inside the window).
* `zenodo_inventory.py` — SHA256SUMS + data dictionary generator for the
  staging tree.
* `dump_headers.py` — compact NetCDF header inspector.
