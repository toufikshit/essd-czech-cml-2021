# Verification of the released rainfall (Supplement Table S7)

| Check | Result |
|---|---|
| Time axis of the 1-minute rainfall | 308,160 UTC slots at 1 min, 1 April 00:00 to 31 October 23:59, 2021 |
| CML population | 1,987 physical links; 3,731 directional sublinks with rainfall |
| Finite directional 1-minute rainfall values | 631,742,984 |
| Independent recomputation of the 1-minute rainfall, Graf flags and thresholds | Identical. A second implementation of the Graf detection with pycomlink, on 37 sublinks (named cases, random sublinks, short links, highest maxima, most blackouts), gave the same 11.4 million rain rates and wet, valid and source flags; thresholds agree within 2e-7 (relative) |
| Derived 5-minute directional rainfall and CML mean from the 1-minute file | Identical. The whole file was derived a second time independently: 244.9 million values identical |
| Gauge-wet comparison | N = 1,331,687; RMSE 3.63 mm/h; MBE -0.91 mm/h; PCC 0.63 |
| Same comparison without intervals with reconstructed total signal loss | N = 1,327,998; RMSE 3.38 mm/h; MBE -0.89 mm/h; PCC 0.61 |

## Checks with the scripts of this repository

The scripts were run on the released files for a few CMLs (45035_45097_2740,
10001_10066_3107, 10001_10385_3108, 10706_10757_3347; notebook
`01_QC_and_Rainfall.ipynb` repeats this):

* `scripts/qc_from_raw.py` on the released RAW file: every variable equals the
  released QC file;
* `scripts/rainfall_1min.py` on that output: every variable equals the released
  rainfall file;
* `scripts/derive_5min.py`: every variable equals the 5-minute series used in
  the paper;
* `scripts/gauge_validation.py`: the validation samples of these CMLs equal those
  of the paper;
* `scripts/radar_cml_paths.py` on the RADKLIM file of 13 July 2021: the path
  averages equal the released radar file (float32).

The figure scripts in `scripts/figures/` reproduce Figs. 1 to 7 of the paper pixel for
pixel. `scripts/wvd_validation.py` reproduces the WVD statistics (31 CMLs, 45 sublinks,
236,793 pairs, RMSE 3.46 g m-3, PCC 0.598) and the content of Fig. 8; the image size
differs by a few pixels from the published file, which was drawn on another computer.
