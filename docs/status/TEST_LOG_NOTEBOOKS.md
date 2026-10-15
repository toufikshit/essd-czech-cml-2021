# Notebook release candidate: testing evidence (2026-10-13)

**Verified in local sandbox:**

- `python tools_verify_originals.py`: 44/44 original source scripts and notebooks byte-identical to recorded hashes.
- `PYTHONPATH=src:tests python -m pytest -q tests/test_notebook_adapters.py tests/test_mini_structure.py tests/test_comparison_alignment.py`: **15 passed**.
- Every code cell in eight new `.ipynb` files passed Python syntax compilation and every notebook passed nbformat validation.
- Two-day synthetic 7-file fixture: generated network geometry PNG, carrier-frequency/path-length PNG, original five-sensor availability CSV/PNG/PDF (unchanged original source) and monthly reconstructed availability CSV/PNG. Station coordinate metadata were added **only to the synthetic fixture** to exercise static station plots. No science inputs or original files modified.
- Synthetic guards: MHz→GHz without accidental double conversion; lengths remain km; missing station geometry raises; deliberate value/mask/flag/time errors caught by existing tests; CHMI converter dry-run and overwrite guards tested.
- `netCDF4` and `pycomlink` are not installed in this local sandbox. **Full CHMI conversion was not run here**, and this notebook extension is not a new scientific Graf retrieval run.

**Previously verified by the user on HPC (against prior repository):**

- Five-sensor original availability on real Apr–Oct 2021: CML 58.67%, radar 99.48%, gauge 99.74%, temperature 99.57%, humidity 99.56% mean native-timestep availability.
- Original Graf + Pastorek + ITU2005 scientific functions ran on synthetic fixture and real 14-day sample.
- Full-history real link `10004_311866_3129`, 1–14 Apr 2021: 7,666 finite directional values compared, RMSE 0.0, zero finite-mask disagreements and zero Graf-flag disagreements. **One link only**; do not extrapolate to entire population.

**HPC validation still required for this notebook extension:**

1. Confirm seven HPC JSON paths; run notebooks 00, 02, 03 against actual April–October data and compare figures to manuscript.
2. Test standalone CHMI v5 converter against a limited, already-downloaded CHMI source directory (it can fetch metadata if missing), then inspect NetCDF4 values and metadata history. Never overwrite released NetCDFs.
3. Recover exact Graf event/confusion/hexbin detector before claiming exact plotted figure reproduction, or run explicit v3 compatibility experiment and compare output populations.
4. Supply external corrected native RADKLIM + CML IDW event grid and V6 WVD support files before testing notebooks 05 and 06.
5. Keep original files immutable, and re-run all tests after any notebook/interface modification.
