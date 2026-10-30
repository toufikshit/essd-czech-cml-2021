"""Stage 5B (lead): merge the figure provenance rows (Agents 3a, 3b and the lead rows for Figs 5, 6, 8).
Usage: python -I s5b_03_provenance.py"""
import os, pandas as pd
H = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
a = pd.read_csv(os.path.join(H, 'agents', 'A3a_provenance.csv')); b = pd.read_csv(os.path.join(H, 'agents', 'A3b_provenance.csv'))
P = r'essd-czech-cml-2021\original\code\figures\rebuild_ESSD_FINAL_GRAF_shared_population_NEW_EVENT_FINAL_HEXBIN_SPACING.py'
S = '7358e242925dbe5c2b317d2d2e67213cc54a52d71173adc814ecc66a09dcc6fb'
lead = pd.DataFrame([
    dict(figure=5, published_file='FINAL_GRAF_confusion_GAUGE_same_population_202104_202110.png', producer_file=P, producer_sha256=S,
         producer_function_or_cell='gauge_validation() -> plot_confusion(); figure made by stage5b/fig568/scripts/s5b_01_fig5_fig6.py importing plot_confusion/style/save unchanged',
         original_inputs='Stage 5A control samples (stage5a/control/control_released_samples.npz) = published validation (full-history QC, 5-min rain, gauge)',
         new_inputs='Stage 5A Case A samples (stage5a/new/stage4_new_samples.npz): QC REV1 + Stage 4 5-min rain + same gauge',
         reproduced_original='yes - pixel-identical to the published PNG (0 pixels differ; 3311x3217 at 450 dpi)',
         changed_by_reconstruction='yes - TP 1,026,880->993,974; FP 4,014,661->3,757,603; FN 328,073->337,713; TN 29,196,125->29,159,941; hit 75.8->74.6 % of gauge wet; false alarm 12.1->11.4 % of gauge dry; correct negative 87.9->88.6 %',
         regenerated_file='stage5b/FINAL_FIGURES/FINAL_GRAF_confusion_GAUGE_same_population_202104_202110.png (+pdf)', lines_changed_in_producer='0 (functions imported unchanged)',
         notes='Case A (with blackout-reconstructed samples) only.'),
    dict(figure=6, published_file='FINAL_GRAF_hexbin_gauge_WETONLY_202104_202110.png', producer_file=P, producer_sha256=S,
         producer_function_or_cell='gauge_validation() -> hexbin(x[wet], y[wet]); figure made by s5b_01_fig5_fig6.py importing hexbin/style/save unchanged',
         original_inputs='as Fig. 5', new_inputs='as Fig. 5',
         reproduced_original='yes - pixel-identical to the published PNG (0 pixels differ; 3033x2561 at 450 dpi)',
         changed_by_reconstruction='yes - N 1,354,953->1,331,687; RMSE 3.33->3.63; MBE -0.93->-0.91; PCC 0.585->0.633 mm/h (Case A)',
         regenerated_file='stage5b/FINAL_FIGURES/FINAL_GRAF_hexbin_gauge_WETONLY_202104_202110.png (+pdf)', lines_changed_in_producer='0 (functions imported unchanged)',
         notes='Case A only; Case B (N 1,327,998, RMSE 3.38, MBE -0.89, PCC 0.606) reported in text, not in the figure.'),
    dict(figure=8, published_file='wvp_support_all_E_band.png', producer_file=r'essd-czech-cml-2021\original\code\wvd\plot_wvp_support_hexbins_no_titles.py',
         producer_sha256='see repository original/SHA256SUMS.txt', producer_function_or_cell='WVD V6 support hexbins (verified in the 2026-10-27 verification)',
         original_inputs='WVD V6 NetCDF + CSVs (read RAW only)', new_inputs='unchanged',
         reproduced_original='yes (2026-10-27 data verification: Table S13 and Fig. 8 metrics recomputed)',
         changed_by_reconstruction='no - WVD reads RAW only; RAW values unchanged; no excluded CML among the 676 WVD CMLs; the one P6 artefact minute is only in the gauge-method series, not in the 236,793 graf-method pairs (Stage 5A, Agent 4)',
         regenerated_file='none - published file copied unchanged to stage5b/FINAL_FIGURES (byte-identical)', lines_changed_in_producer='n/a', notes='Preserved by author instruction.'),
])
f = pd.concat([a, b, lead]).sort_values('figure')
f.to_csv(os.path.join(H, 'FIGURE_PROVENANCE.csv'), index=False)
print(f[['figure', 'published_file', 'reproduced_original', 'regenerated_file']].to_string(index=False, max_colwidth=70))
