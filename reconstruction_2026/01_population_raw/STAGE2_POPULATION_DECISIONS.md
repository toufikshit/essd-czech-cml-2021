# Stage 2: raw population and exclusion decisions (proposals only)

**Date:** 2026-10-27.
**Status:** completed; awaiting author approval. Nothing was added, removed or masked. No NetCDF, QC or rainfall
was run or written.

**Input:** the approved raw population of 2580 CMLs (5160 sublinks), from the full-history RAW
`extracted\full\...\raw\CZ_CML_RAW_1MIN_20200101_20211231.nc` (sha256 6acfeed2…, unchanged).

**Fixed by earlier decisions:**
* the 13 TX = 0 CMLs are excluded (none of them is in the 2580);
* the 255 empty CMLs are not in raw.

Neither group is counted again here.

**Who did what:**
* Agent 1 did a full read of the raw file (`scripts\s2_01`–`s2_03`).
* Agent 4 did an independent full read (`verification\v2_01`–`v2_06`): **0 mismatches in 5160 rows**; every
  proposal, id and count was reproduced.
* The lead wrote the final tables (`scripts\s2_lead_01_finalize.py`; no data read).

## 1. Audit summary (`RAW_LINK_AUDIT.csv`, one row per CML × sublink)

* **Metadata:** frequency, length and endpoints agree with the operator metadata in all 5160 rows; polarization is
  Vertical.
* **No placeholder values remain** (TX ≤ −99 or 255; RX ≤ −99.9 or > −20).
* **Observations:** 918,482,728 minutes with both TX and RX; 2575 CMLs and 5139 sublinks have data.

| Class (primary) | Sublinks | CMLs* | Content |
|---|---:|---:|---|
| 1 Suitable for raw and QC | 4943 | 2507 | – |
| 2 Raw, likely unsuitable for rainfall | 150 | 59 | 149: data end before 2021-01-20 (none in the release window); 1 (86097_86098_4937 sl0): only 8 paired minutes on 2021-06-18, no 60-minute run. **A QC matter, not proposed for raw exclusion** |
| 3 Missing essential observations | 7 | 5 | 6: TX entirely missing, RX present; 1: RX entirely missing |
| 4 Missing essential metadata | 16 | 8 | No frequency on either sublink (operator band 10) |
| 5 Completely empty | 14 | 0 | Empty sublinks; no CML is empty on both sublinks |
| 6 Other documented issue | 30 | 1 | 16: TX = 0 at every minute; 14: 1–2 impossible TX values (§3) |

\* A CML takes the class of its most usable sublink.

Information only, not used for any decision:
* operator-vs-raw completeness (2181 sublinks more than 5 % short; Stage 1B limitation);
* 22 sublinks with more than 50 % of their RX minutes without TX (3.34 M RX-only minutes);
* 21 groups of parallel records (not duplicates);
* gauge, radar and WVD membership.

## 2. Proposals (`RAW_PROPOSED_EXCLUSIONS.csv`)

| ID | Scope | CMLs / sublinks | What it removes | Agent 4 | Lead recommendation |
|---|---|---|---|---|---|
| P1 | Whole CML, completely empty | 0 | – | confirmed | nothing to do |
| **P2a** | Whole CML | 2: 10004_11582_3122, 27040_27099_3747 | TX is a placeholder at every minute (NaN in raw); 1,038,881 RX-only values | consistent with the 13 | **exclude (rule R1)** |
| **P2b** | Whole CML | 3: 10662_10795_3342 (TX only, January), 29002_29004_3810 (RX only, January), 49111_329218_4481 (RX only) | 337,432 single-quantity values; no paired minute | could be kept with a flag | **exclude (rule R1)**; alternative: keep and flag |
| **P2c** | Whole CML | 1: 41920_323338_4235 | TX = 0 on both sublinks except one placeholder-mixed minute; 362,288 paired minutes | raised by Agent 4 | **exclude (rule R1, same criterion as the 13)** |
| P3 | Whole CML | 8 without frequency | 3,049,403 paired minutes | keep and flag | **keep in raw with a metadata flag**; QC excludes them (no detection limit without frequency) |
| P4 | One sublink | 11 empty sublinks in usable CMLs | nothing (the sublink dimension stays) | documentation only | **no data change**; documented |
| P5 | One sublink | 16 with TX = 0 at every minute, while the other sublink works | 1,806,630 paired minutes | author decision; not equivalent to the 13 | **keep in raw**; decide the QC rule in Stage 3 |
| P6 | Individual minutes | 16 minutes in 14 sublinks with TX 78–137 dBm | TX and RX of those minutes | corrupted records; record-level rule | **keep raw values, documented**; QC masks TX and RX of these minutes |

**Rule R1 (proposed, generalising your decision on the 13).** A CML is non-operational, and is excluded from the
publication RAW and QC, when it has no sublink with paired valid TX and RX observations at a reported TX level.
This covers:
* TX constant 0 dBm;
* TX only a placeholder;
* TX or RX entirely missing.

R1 selects exactly the 13 (already excluded) plus the 6 CMLs of P2a, P2b and P2c.

## 3. Mixed-placeholder minutes (basis of P6)

Every impossible TX value equals the 1-minute mean of an operator placeholder and a valid sample in the same
minute:
* (255 + 18)/2 = 136.5, with RX (−99.9 − 48.9)/2 = −74.4 at the same minute;
* (255 + 1)/2 = 128, with RX (−99.9 − 58.6)/2 ≈ −79.2;
* (255 − 10)/2 = 122.5;
* (255 − 99)/2 = 78.

There are no TX values between 30 and 78 dBm [F].

**Mechanism [C].** The readers resample to 1 minute by the mean. Only format_b RX −99.9 is masked before that. The
raw masking of TX ≤ −99, TX = 255 and RX ≤ −99.9 (notebook cell 113) runs on the 1-minute means. A minute in which a
placeholder was averaged with a valid sample therefore survives with a corrupted value.

The same mechanism probably explains part of the 84 TX values below −30 dBm (for example (−99 − 10)/2 = −54.5) [H].

**Limitation.** Minutes where only RX was mixed with a placeholder cannot be identified reliably from 1-minute data,
and no further diagnostic was run (author instruction). Only the impossible TX values are detectable. The number of
affected minutes is very small (16 detectable above 60 dBm, up to 84 below −30 dBm, out of 918 M).

Where 60 dBm comes from: the 30–78 dBm value gap supports it, not the code. The notebook constant `TX_MAX = 60`
(cell 80) is a diagnostic flag that is never applied.

## 4. Population before and after

| Step | CMLs | Sublinks |
|---|---:|---:|
| Approved raw population (your Stage 1B decision) | 2580 | 5160 |
| Rule R1: P2a −2, P2b −3, P2c −1 (lead recommendation) | **2574** | **5148** |
| P3, P4, P5, P6 as recommended (kept in raw; handled in QC or documentation) | 2574 | 5148 |

Alternatives:
* **No raw exclusions:** 2580 / 5160. Everything goes to QC.
* **R1 without P2b:** 2577 / 5154.
* **All of Agent 1's draft plus P2c:** 2566 / 5107 (also removes P3 and P5).

No additions are proposed (`RAW_PROPOSED_ADDITIONS.csv` is empty): the 13 are excluded by your decision, and the
25 metadata candidates have no observations in the existing files.

The 6 R1 CMLs in the released QC and rainfall:
* 10004_11582_3122 and 27040_27099_3747 are in QC with placeholder-based `trsl`;
* 10662_10795_3342 and 29002_29004_3810 are in QC but empty in the release window;
* 41920_323338_4235 and 49111_329218_4481 are not in QC.

## 5. Approval requested

1. **Rule R1:** exclude the 6 CMLs (2580 → 2574), all six or without P2b.
2. **P3:** keep the 8 CMLs without frequency in raw with a flag (recommended), or exclude them.
3. **P5:** keep the 16 TX = 0 sublinks in raw and decide their QC treatment in Stage 3 (recommended), or exclude
   them now.
4. **P6:** keep the raw values and mask the 16 mixed-placeholder minutes (TX and RX) in QC (recommended), or mask
   them in raw.

Nothing is executed until you approve. Stage 3 (QC settings, then the QC run) needs its own approval.

## 6. Author decisions (2026-10-28)

| # | Decision |
|---|---|
| R1 | **Approved.** The 6 non-operational CMLs (P2a, P2b, P2c) are excluded from RAW and QC. **Approved RAW population: 2574 CMLs, 5148 directional sublink slots.** |
| P3 | **Approved:** keep the 8 CMLs without frequency in RAW, with the missing frequency documented in the metadata. Exclude them from the rainfall-ready QC population, because frequency-dependent processing cannot be performed. No frequency values are invented. |
| P5 | **Approved:** keep the TX = 0 sublinks documented in RAW. These directions are non-operational and must not contribute valid observations to the QC rainfall product; the opposite direction is retained where it has valid observations. After R1, **15** sublinks remain (41920_323338_4235 sl1 left with R1). |
| P6 | **Approved:** keep the 16 mixed-placeholder minute values in RAW, with the artefact documented. QC masks both TX and RX of these minutes. No additional speculative masks (for example, TX below −30 dBm is not masked). |
| P4 | Documentation only (11 empty sublinks); no data change. |
| Consistency | QC_CML_IDs ⊆ RAW_CML_IDs. No CML excluded from RAW may remain in QC: when Stage 3 is authorised, the 4 R1 CMLs now in the released QC (10004_11582_3122, 27040_27099_3747, 10662_10795_3342, 29002_29004_3810) are removed. Original full-history files are preserved; every exclusion and QC mask is recorded reproducibly. |

Stage 3 is **not** authorised. The QC settings and the population handling are presented separately for approval
(`stage3\STAGE3_QC_SETTINGS_PROPOSAL.md`).
