# Original scientific scripts

The scripts in `code/` are kept unchanged. The scripts in `../scripts/` call them
with the data paths of `config/user_paths.json`. `SHA256SUMS.txt` lists their
checksums; `python tools_verify_originals.py` checks them.

| Folder | Scripts | Used by |
|---|---|---|
| `figures/` | availability (Fig. 2), event series and validation plots (Figs. 4 to 6), event maps (Fig. 7) and its runner, IDW grid builder (Fig. 7) | `scripts/figures/` |
| `validation/` | Graf helper functions of the gauge validation | `scripts/gauge_validation.py` |
| `radar/` | RADKLIM processing along the CML paths | `scripts/radar_cml_paths.py` |
| `wvd/` | WVD retrieval (13, 19, 20), quantization (30c, 31, 32), validation (21, 33) and Fig. 8 | `scripts/wvd_validation.py` |
| `chmi/` | conversion of the CHMI station CSV files to NetCDF | run directly (`--help`) |
