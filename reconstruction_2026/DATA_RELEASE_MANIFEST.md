# Data release manifest (proposal for the next Zenodo version)

Prepared 2026-10-29. Nothing has been uploaded. The data will be published in
the Zenodo record (DOI to be added when the data release is published). The
same content as a table: [DATA_RELEASE_MANIFEST.csv](DATA_RELEASE_MANIFEST.csv).

Checksums were computed on 2026-10-29 by reading the files on disk (all files
of 3 GB or less). They agree with every value recorded earlier (the NetCDF
attributes, the stage reports and `zenodo_staging/SHA256SUMS.txt`). The earlier
QC file (3.35 GB) was not re-hashed; its value comes from
`zenodo_staging/SHA256SUMS.txt`.

Paths are relative to the author's release folder (`ESSD_RELEASE/`).

## 1. Files that change

| Role | Proposed release file | Source on disk now | Bytes | SHA-256 | Replaces |
|---|---|---|---:|---|---|
| RAW | `CZ_CML_RAW_1MIN_20210401_20211031.nc` (2,574 CMLs) | **does not exist** | – | – | staged RAW, 2,580 CMLs, 2,257,209,979 B, `46b130aa144a4e96d9698ca94c06ff625ad7670d8aa362b6e4c49a467806d012` |
| QC | `CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` (1,987 CMLs) | `FINAL_RECONSTRUCTION_20261027/stage3/QC_CORRECTED_CANDIDATE_REV1.nc` | 2,409,368,783 | `ca8afa95aacb50dacbd5f7d96f776c56a78c63d1a4888e9f721a29cdcfafde65` | staged QC, 2,219 CMLs, 3,351,446,358 B, `badc3bc800904806c85ccc740c0f28448b63094911666216c038f56c5bc44d19` |
| Rainfall, 1 min | `CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc` (1,987 CMLs) | `FINAL_RECONSTRUCTION_20261027/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031_STAGE4.nc` | 445,051,457 | `9014567c3201be60441d1d8721c809d13ad5fa193e38a7b65606848b868a607e` | staged rainfall, 2,219 CMLs, 493,609,050 B, `be2cd93843b6babbfe9d091de35d23db81071c66c7d9e1eee3ae03fee7a0780f` |
| Rainfall, 5 min (derived, **not released**) | not a release file; reproducible from the 1-minute file | `FINAL_RECONSTRUCTION_20261027/stage4/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_DERIVED_20210401_20211031_STAGE4.nc` | 338,126,319 | `a28d920be643d8647f21263a247b61e40496b793a44f453a902bdac8018820ee` | nothing (new file) |

### RAW: no 2,574-CML file exists

The approved RAW population (2,574 CMLs, 5,148 sublink slots) is defined by
the table `01_population_raw/RAW_APPROVED_POPULATION.csv`. No NetCDF with
exactly these CMLs has been written. The staged RAW file holds 2,580 CMLs: the
approved 2,574 plus the 6 non-operational R1 CMLs (checked on the CML ids):

`10004_11582_3122`, `10662_10795_3342`, `27040_27099_3747`,
`29002_29004_3810`, `41920_323338_4235`, `49111_329218_4481`.

The QC was computed from the full-history RAW
(`extracted/full/ESSD_CZ_2021_HANDOFF/data/source_full_history/raw/CZ_CML_RAW_1MIN_20200101_20211231.nc`,
2,496,052,589 B, sha256 `6acfeed2a2864b8c0e69807e9834dd03368d012fac4a2cae3d430c27b71be8aa`)
restricted to the approved population. This file is not a release file.

To do (needs author approval): write the release RAW by removing the 6 R1
CMLs from the staged 2,580-CML file, with no other change, and record its
SHA-256. Keeping the 2,580-CML file instead would still contain every QC CML,
but it would contradict the approved RAW population of 2,574.

### QC and rainfall: before release

* The approved files carry workspace names. The repository code expects the
  release names above (`src/essd_cz/config.py`, `EXPECTED_NAMES`), so they
  must be renamed.
* Metadata still says "candidate": the QC title and `processing_level`, and the
  rainfall title ("Stage 4 candidate"). The independent check also noted that
  `qc_parameters_json` does not encode the revision 1 blackout rejection rule
  (it is stated in `qc_blackout_method`) and that the QC title does not mention
  revision 1. Any metadata rewrite changes the SHA-256, so the values above are
  for the files as approved and must be recomputed after a rewrite.
* The 5-minute file is not released, as before and as stated in the manuscript:
  the 5-minute series are derived from the 1-minute file exactly
  (`src/essd_cz/rain5.py`). Its hash is listed only so that a local derivation
  can be checked against the verified one.

## 2. Files that stay the same

| Release file | Bytes | SHA-256 |
|---|---:|---|
| `CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc` | 734,441,295 | `38dae054ec232f41653a095a5a19bc29c8eb4aa844a237cbd8615c5e0f14262a` |
| `CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc` | 55,750,031 | `e82d446b8970e6cc20e0448e8875346ef2d0c84564600f71a039f102ed9c1e2c` |
| `CHMI_10MIN_TEMPERATURE_20210401_20211031.nc` | 53,757,072 | `6a0bb785a9569c7f5d4ec9ab34e7a89622c9fc8d67a93bc3dd4a8553e1a3d3d4` |
| `CHMI_10MIN_HUMIDITY_20210401_20211031.nc` | 5,244,310 | `221a14c19c5be0470be3d298129700d77cad085d853eda8776163ed8b8537e5b` |
| `CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc` (supplement) | 415,709,686 | `35bc67ed1b905416bdceb4d44666f706a22ee36e8ba7457f551a628c05b4c4bf` |
| `CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv` (supplement) | 605,081 | `cab58ca361e506b5e45b55ec95851d14a42c5c266bba34db2f29f0107d97fbf7` |
| `CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv` (supplement) | 308,690 | `67012dca9ded4451723ed2a902cfbbee404b1c5c3aa81e326453549cfd3e498c` |

All are in `zenodo_staging/primary/` or `zenodo_staging/supplementary_wvd/`.

Note on the radar file: it holds 1,132 CML paths. 1,033 of them are in the
corrected QC and 99 are not. Two of the 99 (`10004_11582_3122`,
`27040_27099_3747`) are R1 CMLs outside the approved RAW. The radar file is
unchanged by the approved decisions; whether to drop these two paths is for
the author to decide.

## 3. Documentation in the deposit

`README_DATA.md`, `DATA_DICTIONARY.md` and `SHA256SUMS.txt` in
`zenodo_staging/` still describe the 2,580 / 2,219 release. They need the new
counts, the new variables (`qc_flag`, `qc_sublink_status`,
`qc_detection_limit_mmhr`, `trsl_source_flag_1min`, the 5-minute flag counts)
and the new checksums. The dataset zip `ESSD_CZ_CML_2021_dataset.zip` must be
rebuilt. None of this was done here.

## 4. Repository code to update when the files are swapped

`src/essd_cz/catalog.py` (`CONTRACT`) expects 2,580 RAW, 2,219 QC and 2,219
rainfall CMLs. With the new files, `python run.py ... check` will report
"Expected cml_id=2219" until these counts are changed to 2,574 / 1,987 /
1,987. This was not changed on this branch, because the 2,574-CML RAW file
does not exist yet and the change must match the files that are finally
published.
