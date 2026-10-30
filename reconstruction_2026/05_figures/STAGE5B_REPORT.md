# Stage 5B report: publication figures

**Date:** 2026-10-29.
**Status:** Figs 2, 3, 5 and 6 regenerated; Figs 4 and 8 unchanged; Figs 1 and 7 stopped (author decisions or
MetaCentrum needed). Everything was verified independently by Agent 4. Awaiting author approval; Stage 6 has not
started.

**Untouched:** no original figure, data product, manuscript or repository file was modified. The bytecode caches
written by imports were removed; see `PROCESSING_LOG.md`.

**Files:**
* final candidates: `FINAL_FIGURES\` (+ `SHA256SUMS.txt`);
* originals: `original_figures\`;
* tables: `FIGURE_PROVENANCE.csv`, `FIGURE_STYLE_COMPARISON.md`, `FIGURE_NUMERICAL_COMPARISON.csv`;
* agent reports: `agents\A3a_figures_1_3.md`, `agents\A3b_figures_4_7.md`, `agents\A4_stage5b_verification.md`.

## Per figure

| Fig. | Producer (proof) | Changed by the reconstruction? | Result |
|---|---|---|---|
| 1 Study area | **Not found locally** (thorough search, confirmed) | Yes. The published map draws the old 2219 QC network, including 268 CMLs that are not in RAW and 4 of the 6 R1 CMLs (27040_27099_3747 visibly). 415 paths would disappear and 183 appear | **Stopped** |
| 2 Availability | `original\code\figures\availability_native_timestep_percentage.py`, unmodified (pixel-identical rows) | Yes, the CML row only: denominator 2219 → 1987; mean 58.67 → 66.35 %; May still lowest (49.7 → 57.9 %) | Regenerated |
| 3 CML characteristics | Notebook `cml dataset ()runing.ipynb`, cell 22 (partial reproduction; the original input is unavailable) | Yes: now the QC population, panels (a)/(c) 3731 retained sublinks and (b) 1987 CMLs | Regenerated; tick question open |
| 4 Event time series | Rebuild script `event_plot` (pixel-identical) | No visible change: CML 45035_45097_2740 retained; event total 54.873 → 54.872 mm | Published figure kept |
| 5 Contingency | Rebuild script `plot_confusion` (pixel-identical) | Yes (Case A): hits 74.6 %, false alarms 11.4 %, correct negatives 88.6 % | Regenerated |
| 6 Gauge hexbin | Rebuild script `hexbin` (pixel-identical) | Yes (Case A): N 1,331,687, RMSE 3.63, MBE −0.91, PCC 0.633 | Regenerated |
| 7 IDW fields | Runner `06_run_multiple_events_vs_radar.py` → renderer `03_blettner_event_maps_timeseries.py`; grid builder `cml_grid_essd_publication_safe_singleopen.py` [H, strongly supported] (CML row 92 % same colour class) | Yes: CMLs near the event 1105 → 994. The CML row is markedly wetter at 18:05–20:05, partly from blackout-reconstructed minutes; at 19:05, nodes ≥ 20 mm 166 → 483 | **Stopped**: radar row needs MetaCentrum. Stage 4 grid ready; script `fig47\metacentrum_commands\fig7_metacentrum.sh` |
| 8 WVD | `plot_wvp_support_hexbins_no_titles.py` | No: WVD inputs unchanged (Stage 5A) | Preserved (byte-identical) |

The Case B validation (N 1,327,998, RMSE 3.38, MBE −0.89, PCC 0.606) is not shown in any figure. It belongs in the
manuscript text and the Supplement.

## Independent verification (Agent 4)

Every claim checked holds, with no mismatches:
* the producer copies and their changed lines;
* the pixel identities;
* Fig 2: 1987 CMLs, 66.35 %, +7.68 pp;
* Fig 3: counts, short-path shares and tick behaviour;
* Fig 4: CML retained; the automatic reselection;
* Figs 5/6: Case A values, fonts, sizes and dpi;
* Fig 7: the field numbers recomputed from the grids;
* Fig 8: byte-identical.

No unintended style change.

## Decisions needed

1. **Fig 1:** supply the producer from MetaCentrum, or authorise a newly drawn map (a style change). Also choose the
   population to draw: QC 1987 or RAW 2574.
2. **Fig 2:**
   * Confirm the gauge row: the current gauge file gives 99.5–99.8 %; the published figure shows 98.0–98.4 %, from a
     gauge file that is no longer available.
   * The caption should say "at native time steps", not "monthly".
3. **Fig 3:** accept the producer's automatic, irregular ticks, or authorise a tick-only fix.
4. **Fig 4:** keep the published figure with CML 45035_45097_2740. With the new data, the producer's automatic
   choice would switch to 28038_28079_2496: a heavy-rain case newly captured through blackout minutes, whose nearest
   gauge is 3.9 km away. Recommended: keep the published figure.
5. **Fig 7:** run the prepared MetaCentrum script:
   * Part A is read-only and confirms the grid producer;
   * Part B renders the figure with the identical runner arguments.

   Then return the figure for comparison.

## Text changes to carry into Stage 6 (not made now)

* **Lines 296–299 (Fig. 3 text):**
  * short paths: QC 1987 gives 15.05 % < 1 km and 0 % < 0.1 km; RAW 2574 gives 27.39 % and 0.97 %;
  * E-band: 810 of 1987 CMLs.
* **Line 311 (QC counts):** 2574 raw, 1987 retained.
* **Line 533:** POD 74.6 %, false alarms 11.4 % of gauge-dry, correct negatives 88.6 %. "About four times" still
  holds (false alarms are 3.78 × hits).
* **Line 549:** N 1,331,687, RMSE 3.63, bias −0.91, PCC 0.63, plus the Case B values.
* **Line ~556 (Fig. 7 text):** re-check "lower peak accumulation" once the MetaCentrum figure exists.
* **Captions:** Fig 2 ("native time step"; CML share of the 1987 QC CMLs); Fig 3 (population statement); Fig 6
  (numbers via the figure itself).
