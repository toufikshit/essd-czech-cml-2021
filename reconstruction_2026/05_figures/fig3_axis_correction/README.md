# Figure 3: approved axis correction (Stage 6)

`fig3_axis_corrected.py` is the Stage 5B Figure 3 script (`../fig2_fig3/`) with only the tick helper
`nice_ticks()` and the y-limit lines of panels (a) and (b) changed (`fig3_axis_corrected.diff.txt`).
The earlier automatic ticks were irregular and the path-length axes ended at 30.5 km, which cut off the longest
links (33.8 km). The new ticks are regular and cover the full data range. Data, bins, colours, fonts and layout
are unchanged. Population: the 2574 CMLs of the raw dataset (panel b) and their 5132 sublinks with frequency
metadata (panels a and c).

Run (paths as in the Stage 5B script):

    python -I -B fig3_axis_corrected.py --population raw2574 --raw-full <CZ_CML_RAW_1MIN_20210401_20211031.nc> --out-png f03_cml_characteristics.png --out-pdf f03_cml_characteristics.pdf

Output used in the manuscript: `f03_cml_characteristics.png`, sha256
553d1ffe32d62718f7906f9719bd941abd56796e277b7e09d97d3faaf1f74b8f (it replaces the Stage 5B Figure 3 in
`../FINAL_FIGURES_SHA256SUMS.txt`).
