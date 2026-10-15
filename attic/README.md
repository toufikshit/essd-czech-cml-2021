# attic/ — superseded repository-management artifacts

Byte-identical copies of handoff files that are no longer the live versions.
Nothing here is needed to run the workflow; everything is kept for
traceability. The preserved **scientific** history lives in `original/`, not
here.

| File | What it is / why it is here |
|---|---|
| `README_handoff_root.md` | the handoff tree's root README, replaced by the rewritten `README.md` |
| `RUN_FIRST_handoff.txt` | handoff quick-start (HPC-first wording), folded into `README.md` |
| `README_NOTEBOOKS_handoff.md` | one-paragraph pointer to `README_CONNECTED_WORKFLOW.md`, superseded by `README.md` |
| `gitignore_handoff.txt` | handoff `.gitignore`, replaced by the extended one |
| `active_profile_hpc_handoff.json` | the HPC-default profile the handoff shipped with (`"profile": "hpc"`); the live `config/active_profile.json` defaults to the portable `"user"` profile |
| `NOTEBOOK_ADDED_SHA256.json`, `NOTEBOOK_RELEASE_MANIFEST.txt` | hash list / manifest of an **older eight-notebook naming scheme** predating the V2 connected notebooks; stale relative to `notebooks/` |
| `upgrade_ESSD_wrapper_python_ONLY.py` | the V2→V2.1 migration script. The migration was already applied on the HPC: `src/essd_cz/workflows.py`, `src/essd_cz/connected.py` and the profile config were verified identical to the V2.1 payload before packaging. The V2.1 zip itself is not committed (binary archive); it remains in the handoff archive. |
| `notebook_versions/06_...__V2_original.ipynb`, `06_...__FIXED_NETCDF4.ipynb` | the two earlier HPC iterations of notebook 06. The shipped `notebooks/06_WVD_Sensitivity_Length_Frequency.ipynb` is the newest iteration (`_REFRESHED_V6`): it adds the NetCDF4 single-series read compatibility, the 1,352-series input audit, and manifest reloading. |
