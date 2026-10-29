# 1. Five commands to begin

You need Python 3.10+, `numpy`, `pandas`, `xarray`, `netCDF4`, `matplotlib`, `numba`, and `pycomlink`. The HPC `tf-env` already had these dependencies during the October 13 tests. New users can review `requirements.txt` and install it in an isolated environment.

1. Unzip repository. Download the seven NetCDFs separately. Do not rename them. Put them under `data/` (or edit the paths in `config/user_paths.json`).
2. Run `python tools_verify_originals.py`. This verifies all 44 original files, without opening data.
3. Run `python run.py --config config/user_paths.json check`. It checks native resolutions, entity counts, exact native time coverage and the ITU2005 product label. It is not a validation of observation values.
4. Run `python run.py --config config/user_paths.json availability`. Open the PNG and CSV files under `results/availability/`.
5. Compare one link's recomputed rainfall with the released file using `python run.py --config config/user_paths.json compare-rainfall --link-index 3 --days 14` (needs the science extras).

**Avoid** `rebuild-all-rainfall` until you understand runtime and time-window effects. It recomputes all 1,987 physical links. Start with the comparator for one physical link. Read `docs/04_REPRODUCE_AND_COMPARE.md`.

If you see `sh: [CML]: command not found`, terminal output was accidentally pasted as shell commands. Copy only the commands inside code fences.
