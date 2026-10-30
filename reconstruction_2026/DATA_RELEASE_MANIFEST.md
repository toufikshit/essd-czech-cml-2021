# Data release manifest

The released data are in the Zenodo record https://doi.org/10.5281/zenodo.23398309
(DOI reserved; the record becomes public when it is published). The same
content as a table: [DATA_RELEASE_MANIFEST.csv](DATA_RELEASE_MANIFEST.csv).

The archive `ESSD_CZ_CML_2021_dataset.zip` (stored without recompression) has
6,368,237,019 bytes and SHA-256
`74445851ced5e5db730f8034d6dc5dd969bf8b11944acd5f5904709a88eabdb1`
(MD5 `b807b0fb672784e63117312a0d37dddd`). It was built on 2026-10-30 with
`06_release/r7_04_build_verify_zip.py`, which checks every member against
`SHA256SUMS.txt` after reading it back from the archive.

On 2026-10-30 the archive, `README_DATA.md`, `DATA_DICTIONARY.md` and
`SHA256SUMS.txt` replaced the earlier files of the unpublished Zenodo draft
(`06_release/r7_05_zenodo_replace_draft.py`, resumed after a dropped connection
with `r7_05b_zenodo_resume_replace.py`). Zenodo reported the same MD5 for each
file as the local copy. Only the record description was changed; the draft was
not published.

## Files in the archive

| File | Content | Bytes | SHA-256 |
|---|---|---:|---|
| `primary/CZ_CML_RAW_1MIN_20210401_20211031.nc` | raw TX and RX, 2,574 CMLs | 2,247,974,230 | `7f7845458d0486622517ae3b2ddfe6a3c6a3ceb7d37f8d179e1224cf674e8d4a` |
| `primary/CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc` | quality controlled signals, 1,987 CMLs, 3,731 sublinks pass QC | 2,409,369,807 | `bc2e0c1428f99db9826938df5b953204b9351225d208cd7640acfc8929f4584e` |
| `primary/CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc` | Graf rainfall at 1 min, 1,987 CMLs | 445,052,481 | `a09251c1582057d697a12319bfbc8e6787e17088a00fd1f461a4422dfb10a5f6` |
| `primary/CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc` | 424 CHMI gauges | 734,441,295 | `38dae054ec232f41653a095a5a19bc29c8eb4aa844a237cbd8615c5e0f14262a` |
| `primary/CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc` | RADKLIM along 1,132 CML paths | 55,750,031 | `e82d446b8970e6cc20e0448e8875346ef2d0c84564600f71a039f102ed9c1e2c` |
| `primary/CHMI_10MIN_TEMPERATURE_20210401_20211031.nc` | 251 CHMI stations | 53,757,072 | `6a0bb785a9569c7f5d4ec9ab34e7a89622c9fc8d67a93bc3dd4a8553e1a3d3d4` |
| `primary/CHMI_10MIN_HUMIDITY_20210401_20211031.nc` | 247 CHMI stations | 5,244,310 | `221a14c19c5be0470be3d298129700d77cad085d853eda8776163ed8b8537e5b` |
| `supplementary_wvd/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_2021.nc` | WVD retrieval, 1,352 series | 415,709,686 | `35bc67ed1b905416bdceb4d44666f706a22ee36e8ba7457f551a628c05b4c4bf` |
| `supplementary_wvd/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv` | WVD diagnostics | 605,081 | `cab58ca361e506b5e45b55ec95851d14a42c5c266bba34db2f29f0107d97fbf7` |
| `supplementary_wvd/CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv` | WVD quantization | 308,690 | `67012dca9ded4451723ed2a902cfbbee404b1c5c3aa81e326453549cfd3e498c` |
| `README_DATA.md` | data README | 6,022 | `85b7d139d5822aa48332a83fe744fa147d77e2597543f3dbcfd03da423a728d8` |
| `DATA_DICTIONARY.md` | variables of every file | 14,463 | `c19dc9fdbc4dfd1129b28b0de4c7ff37383da4b096e9bd3a3d6d1eaadbce6d22` |
| `SHA256SUMS.txt` | checksums of the files above | 1,383 | `e56b92698f773d50d3f8039bf9b3c0bb54a54ca4aa5566bb906431e8a378fae2` |

## How the three CML files were made

`06_release/r7_01_build_release_files.py`:

* RAW: the RAW 1-minute file restricted to the 2,574 CMLs of
  `01_population_raw/RAW_APPROVED_POPULATION.csv`. Values, variables, chunking
  and compression are unchanged; one line is added to `history`.
* QC: the quality controlled file of `02_qc/` (sha256 of the processed file
  `ca8afa95aacb50dacbd5f7d96f776c56a78c63d1a4888e9f721a29cdcfafde65`) under the
  release name; only the global attributes `title`, `processing_level` and
  `history` are changed.
* Rainfall: the 1-minute file of `03_rainfall_graf/` (sha256 of the processed
  file `9014567c3201be60441d1d8721c809d13ad5fa193e38a7b65606848b868a607e`)
  under the release name; only the global attributes `title`,
  `source_cml_file`, `source_cml_sha256`, `comment` and `history` are changed.

`06_release/r7_02_verify_package_files.py` compares every variable of the
three files with their sources (NaN-aware), checks the populations (RAW
2,574 = approved list; QC 1,987 CMLs and 3,731 retained sublinks, all in RAW;
rainfall CMLs identical to QC; 631,742,984 finite directional 1-minute rain
rates) and the entity counts of the other four products. All checks pass
(`06_release/r7_02_output.txt`). `python run.py check` also passes on the
seven files (`06_release/run_check_output.txt`).

The 5-minute CML rainfall is not a release file. It is derived from the
1-minute file with `src/essd_cz/rain5.py`.

The radar path file is unchanged. It holds 1,132 CML paths; 1,033 of them are
quality controlled CMLs.
