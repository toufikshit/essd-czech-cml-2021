# Figure style comparison (Stage 5B)

"Pixel-identical" means the producer, run on the original inputs, reproduces the published PNG with 0 differing
pixels. That proves both the producer and its style. Verified independently by Agent 4
(`agents\A4_stage5b_verification.md`).

| Fig. | Size / dpi | Producer style changes | Evidence | Status |
|---|---|---|---|---|
| 1 | – | – | Producer not found locally (Agent 3a; confirmed by Agent 4's search: 1773 files, 244 zips) | **Stopped**: author decision |
| 2 | 6848 × 1508 px, same crop of 125 rows | None: the script is byte-identical, only `--cml-file` changes | CML, radar, temperature and humidity rows of the reproduction are pixel-identical to the published figure; the new figure differs from the reproduction only in the CML row | Regenerated |
| 3 | 8574 × 2398 px (published 8571 × 2395) | None in style (DejaVu Serif, colours, bins, layout). Adaptation: loading, and masking of excluded sublinks | Partial reproduction: the original ~2855-link input is not available; bars agree within 2–5 counts. The automatic ticks become irregular on the new data: (a) 10/30/60/90 and 0/200/500/700; (b) y-axis ends at 300 while bars reach 414 | Regenerated; **ticks: author decision** |
| 4 | 3877 × 3517 px, 450 dpi | None | Reproduction pixel-identical. The regenerated image (CML fixed to the published one) is byte-identical to the reproduction | **Published file kept** (no visible change) |
| 5 | 3311 × 3217 px, 450 dpi | None: `plot_confusion` imported unchanged | Control reproduction pixel-identical | Regenerated (Case A) |
| 6 | 3033 × 2561 px, 450 dpi | None: `hexbin` imported unchanged (grid 40, extent 0–150, jet, log norm) | Control reproduction pixel-identical | Regenerated (Case A) |
| 7 | 600 dpi, 0.01° grid | None planned: renderer and grid functions unchanged; lattice forced to the original grid | CML row reproduced to 92 % of pixels in the same colour class; radar row needs MetaCentrum (RADKLIM, cartopy) | **Stopped** at the grid stage; MetaCentrum script ready |
| 8 | unchanged | – | Byte-identical file | Preserved |

**Comparison images (original | new):**
* `fig568\comparison\` (Figs. 5, 6);
* `fig123\comparison\` (Figs. 2, 3);
* `fig47\comparison\` (Fig. 4; Fig. 7 CML row diagnostic, not a candidate).
