# Notebook code extraction and remaining ESSD source gaps

Source: uploaded ESSD_CZ_2021_SCRIPTS_AND_NOTEBOOKS.zip.
Two notebooks were parsed as JSON. Saved cell source is copied without scientific changes.
Individual notebook cells in notebook_cells/ are exact source-text extractions (.txt),
NOT automatically runnable shell commands. Many commands overwrite outputs or use legacy data;
do NOT blindly re-run them.

## Executable Python code recovered

- python_recovered/monthly_availability_CML_RADAR_GAUGE_from_dataste_papaer_cell10.py
  is the complete 526-line standalone availability heatmap implementation from the
  dataset-paper notebook's code cell 10 (zero-based). The notebook output records
  monthly availability plot PNG and PDF. Uses CML, radar, gauges (NOT temperature/RH).
  It defaults to calendar year 2021 and has older publication paths. Do not confuse
  with the final April–October five-sensor availability figure generator.
- python_recovered/e_band_percentage_from_DATASET_PUBLICATION_cell05.py
  is the 85-line standalone link frequency summary. This is a notebook diagnostic,
  not a 2021 data product producer.

## QC lineage

The original six-step QC source cml_preprocessed_v3.nc was created 2026-03-11.
Both notebooks READ or use it as an input. They contain no code cell that
produces cml_preprocessed_v3.nc. DATASET_PUBLICATION cell 6 shows execution of
cml_qc_filter_standalone_v4_3.py, which creates a DIFFERENT strict QC NetCDF.
That v4_3 script is already inside the originally uploaded ZIP.
Therefore the original March QC executable is still unproven.

## Remaining scripts actually named in notebook cells

See referenced_scripts_not_in_original_zip.tsv; these names are not present as source
files in the uploaded ZIP. Some are older experiments, diagnostics, or validations,
not required producers of final publication files. The included collection tool
collect_additional_notebook_scripts.py can copy the original .py files from exact
known /home/jovyan paths into another source ZIP when run on Jupyter.
It does not search recursively, load NetCDFs, or execute scientific processing.

## Provenance cautions

- Notebook code/output indicates use or execution only, not independently
  reproducibility of the final validated data products.
- DATASET_PUBLICATION.ipynb fast preparation example used full 2020–2021 and
  prior ITU2003 rainfall; do not replace the final 2021 April–October Graf/ITU2005
  /Pastorek rain product with that older product.
- Event-map notebooks contain older RADKLIM paths; the final corrected script is
  already packaged in code/figures/06_run_multiple_events_vs_radar.py.
- No data arrays are included. Source-code licence must be decided separately.
