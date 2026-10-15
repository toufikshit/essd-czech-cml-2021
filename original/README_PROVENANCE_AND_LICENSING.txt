ESSD Czech Republic publication source-and-notebook archive
===========================================================

Scope
-----
This ZIP collects original Python scripts and the two dataset-paper notebooks
from explicit known filesystem paths. The source files are copied byte-for-byte.
The archive contains NO NetCDF measurements, and this is NOT the seven-product
April-October 2021 publication DATA package. It does not reprocess datasets.

Scientific scope
----------------
Final rainfall: cml_graf_itu2005_rebuild.py, Graf RSTD / Pastorek / ITU2005,
2219 physical links. DO NOT substitute older ITU2003 rainfall products.
V6 WVD: production history identifies 19_wvd_eband_progressive_baseline_v6.py.
This is supporting analysis, not the eighth release NetCDF product.

Notebook provenance
-------------------
DATASET_PUBLICATION.ipynb documents older preparation workflows, including an
executed prepare_essd_publication_datasets_fast.py command.
dataste_papaer.ipynb contains earlier experiments and event maps; not all its
intermediate methods reflect the final released rainfall product.

Unresolved provenance
---------------------
The March 11, 2026 six-step QC NetCDF cml_preprocessed_v3.nc is real and is the
source of the publication QC dataset. The EXACT executable that originally
created it is not established. Later standalone July 2026 scripts and the full
fast pipeline are included as RELATED implementations, not mislabelled as the
original March producer.
Original raw-telecom extraction and corrected-RADKLIM upstream run also need
separate provenance if required for full from-source reconstruction.

Licences
--------
The seven DATA products are handled separately: CML-containing files:
CC BY-NC 4.0 (as specified by the authors); CHMI and DWD-derived files:
CC BY 4.0 with appropriate source attribution. Those licences do not
implicitly license independently authored software or notebooks.
A separate explicit software/notebook licence must be supplied by rights
holders before public redistribution. Review notebook outputs for embedded
credentials, personal data and private file paths before publication.

Integrity
---------
manifest.json records exact input paths, sizes, SHA-256 hashes, UTC modified
times, inclusion reasons and differences from known fingerprints.
No recursive search is executed. A file missing from a known path is reported;
it is not silently replaced by a different version.
