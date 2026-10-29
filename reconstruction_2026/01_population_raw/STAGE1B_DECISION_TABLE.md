# Stage 1B decision table

| # | Decision | Options | Recommendation | Evidence | Affects |
|---|---|---|---|---|---|
| 1 | Raw population | A: 2580 as released · B: 2580 + 13 = 2593 · C: B ± other changes | **B (+13, −0)** | The 13 were excluded only by an undocumented non-zero-TX export filter; their pre-QC observations are in hand and plausible | RAW file (13 CMLs added), QC, rainfall, gauge benchmark (10 CMLs) |
| 2 | TX = 0 dBm | Keep as observed · mask TX as missing · exclude CMLs | **Keep, document the uncertainty** | Operator values (no reader fill); a common ATPC level on this network; link budget gives TX of about −8 to 0 dBm; unreported power control cannot be ruled out | Data documentation; absolute `trsl` flagged as uncertain |
| 3 | 255 empty metadata-only CMLs in QC/rainfall | Keep · remove | **Remove** (they are not in RAW) | No observations in any file; the operator records them as missing | QC, rainfall (no data lost) |
| 4 | 5 RAW CMLs without any finite TX&RX minute | Keep in RAW, let QC exclude · remove from RAW | **Keep in RAW**; QC excludes with a recorded reason | A raw product should not drop observations; 2 of them carry real RX | QC exclusion list |
| 5 | 25 metadata CMLs with operator data but in neither file | Add · leave out and document | **Leave out, document** | Their observations are not in the existing files | Documentation only |
| 6 | Late-April completeness (63.5 % of operator minutes) | Document as a limitation · investigate separately later | **Document now**; any separate investigation is your choice | Operator metadata vs RAW counts; cause unknown | README, availability statements, Fig. 2 text |
| 7 | QC method and parameters | – | Stage 2/3 decision (not part of Stage 1B) | Stage 1A `QC_SCRIPT_PROVENANCE.md` | QC, rainfall, validation |

Details are in `STAGE1B_FINAL_ASSESSMENT.md`.

## Author decisions (2026-10-27)

| # | Decision taken |
|---|---|
| 1 | **RAW population: 2580 CMLs (Option A).** |
| 2 | The 13 TX = 0 CMLs are excluded from RAW and QC as non-operational: an explicit dataset exclusion decision, recorded in `STAGE1B_EXCLUDED_13_CMLS.csv`. |
| 3 | The 255 empty CMLs are excluded from QC; they are disjoint from the 13 and not double-counted. |
| 4–7 | Carried into Stage 2 (raw link audit) and Stage 3 (QC method). |
