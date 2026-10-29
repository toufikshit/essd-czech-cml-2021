# Figure 3: approved axis correction (Stage 6)

`fig3_axis_corrected.py` is the Stage 5B Figure 3 script (`../fig2_fig3/`) with only the tick helper
`nice_ticks()` and the y-limit lines of panels (a) and (b) changed (`fig3_axis_corrected.diff.txt`).
The earlier automatic ticks were irregular and the path-length axes ended at 30.5 km, which cut off the longest
links (33.8 km). The new ticks are regular and cover the full data range. Data, bins, colours, fonts and layout
are unchanged. Population: the 1987 quality controlled CMLs and their 3731 retained sublinks.

Run (paths as in the Stage 5B script):

    python -I -B fig3_axis_corrected.py --population qc --out-png f03_cml_characteristics.png --out-pdf f03_cml_characteristics.pdf

Output used in the manuscript: `f03_cml_characteristics.png`, sha256
2e09978067350daf6c3e6abd7fc226c3dc8f53b3c65ff740c661e98d70f8cdf7 (it replaces the Stage 5B Figure 3 in
`../FINAL_FIGURES_SHA256SUMS.txt`).
