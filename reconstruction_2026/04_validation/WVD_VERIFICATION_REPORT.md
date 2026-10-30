# WVD dependency check (Stage 5A)

**Question:** do the approved RAW, QC and rainfall changes affect the WVD input? If they do not, the verified WVD
results are retained:
* 31 CMLs, 45 sublinks, 236,793 pairs;
* RMSE 3.46 g m⁻³, PCC 0.598.

The WVD workflow was not rerun.

## Findings

| Check | Result | Evidence |
|---|---|---|
| Which product the WVD code reads | **RAW only.** `19_wvd_eband_progressive_baseline_v6.py`, `30c_v6_quantization_chunk.py` and `33_v6_binned_validation_hexbins.py` take `CZ_CML_RAW_1MIN_20200101_20211231.nc`; none reads the QC or rainfall product | [C] `--cml-nc` defaults, lines 88/978, 369, 197 |
| RAW values | Unchanged. The approved decisions remove whole CMLs from the publication RAW (R1) but change no value. P6 is a QC-only mask; the RAW values are preserved | Stage 2 decisions |
| Excluded CMLs (the 13 TX = 0 and the 6 R1) among the 676 WVD CMLs | **None** | `RAW_EXCLUDED_CMLS_ALL.csv` vs WVD `cml_id` |
| WVD CMLs within the approved RAW (2574) | All 676 | – |
| P5 (TX = 0) directions in WVD CMLs | None | – |
| P6 mixed-placeholder minutes in WVD CMLs | 3 CMLs: 47097_47918_2217, 60205_60206_2218, 89084_89099_2244 | `QC_PREDEFINED_MINUTE_MASKS.csv` |
| …of which in the 45-sublink validation subset (E-band, 4–10 km, Δρ_q ≤ 0.20, H+T station ≤ 5 km) | 60205_60206_2218 s0 (its P6 minute is on s1, so not affected) and **89084_89099_2244 s1** (P6 minutes on s1) | `validation_population_links.csv` |
| Effect inside the published WVD file for 89084_89099_2244 s1 | The 5-minute interval starting 2021-07-12 07:00 has `total_loss_5min` 74.28 dB (neighbours 60.38 / 60.29) and `wvd_internal_anchor` 45.3 g m⁻³ (neighbours 20.6 / 20.4). **Correction (Agent 4):** this appears only in the **gauge** wet/dry method of the WVD file. In the **graf** method, which Table S13 and Fig. 8 use, the spike itself triggers the wet buffer, so the bin is not safe-dry. The 2021-08-11 11:41 minute, and the P6 minutes of the other WVD CMLs, fall in bins with no total loss in either method | WVD NetCDF [F]; `agents\A4_stage5a_verification.md` |
| Effect on the validation statistics | **None for Table S13 and Fig. 8:** the artefact is not among the 236,793 graf-method pairs. If gauge-method statistics were used anywhere, the effect would be about 2×10⁻⁴ g m⁻³ in RMSE | Agent 4 [F] |

## Conclusion

* The WVD inputs are unchanged by the reconstruction, so the verified WVD results (Table S13, Fig. 8) are retained.
  No WVD rerun is needed.
* One artefact minute, a mixed placeholder in RAW, enters one 5-minute value of the gauge-method WVD series of one
  sublink. It is not among the graf-method validation pairs. It is documented here, not filtered (author rule: no
  new filters).
* If you want it removed in the future, that would be a separate, approved WVD revision.
