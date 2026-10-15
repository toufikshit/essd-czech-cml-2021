# Tested environments

## Author's HPC (original scientific runs, October 2026)

* Interpreter: `/home/jovyan/my-conda-envs/tf-env/bin/python` (referenced via
  `config/active_profile.json → science_python_by_profile.hpc`).
* Exact package versions of that environment were not exported into the
  handoff; the HPC runs are evidenced by the executed artifacts in
  `docs/evidence/hpc_run_2026-10/`.

## Release-candidate re-test (Windows 11, October 2026)

All 20 tests of `tests/`, `tools_verify_originals.py` and the synthetic
connected mini chain (`tests/run_connected_mini.py`) passed with:

```
Python 3.12.10
numpy==2.5.4          pandas==3.0.6         xarray==2026.9.0
netCDF4==1.7.4        scipy==1.18.1         matplotlib==3.11.2
numba==0.68.0         llvmlite==0.50.0      pycomlink==0.6.0
nbformat==5.11.1      pytest==9.1.1
requests==2.34.2      beautifulsoup4==4.15.0  tqdm==4.70.1
```

These are "known-good" pins, not strict requirements; `requirements.txt`
deliberately stays unpinned minimums. Two environments are distinguished on
purpose:

* **Notebook kernel** — needs `requirements.txt` + `requirements-notebooks.txt`.
* **Science producer interpreter** — the interpreter used for rainfall
  recomputation and WVD retrieval (`numba`, `pycomlink`); set it per profile in
  `config/active_profile.json → science_python_by_profile` if it differs from
  the notebook kernel.

## Continuous integration

`.github/workflows/ci.yml` runs the original-hash check and the full synthetic
test suite on Linux. CI uses synthetic fixtures only — no real data products
and no multi-hour computations.
