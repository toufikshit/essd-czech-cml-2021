# ESSD Czech 2021 — connected reviewer walkthrough (V2)

**Start with `notebooks/00_Explore_Data.ipynb` and proceed in order.**

Change the paths **once** in `config/user_paths.json` (copied from example), and select `"profile":"user"` **once** in `config/active_profile.json`. On the original HPC, leave the default `"profile":"hpc"`. `workspace` is also controlled in that one file. The same scientific original scripts and notebooks are included byte-for-byte, and `tools_verify_originals.py` checks the 44 originals.

## Paths: demo on the author's HPC OR portable external review

* Author HPC: leave `config/active_profile.json` at `"profile":"hpc"`, using `config/hpc_source.example.json` pointing at the **existing** `/home/jovyan` paths. The code does not store HPC paths in scientific functions.
* Outside user: `cp config/user_paths.example.json config/user_paths.json`; edit the seven paths there and any optional radar/WVD paths. Switch `config/active_profile.json` to `"profile":"user"` **once** (ideally also change `workspace` to `results/connected_workflow_user`). Relative paths in JSON are resolved relative to the configuration file, absolute paths also work. Downloaded data can live anywhere; no source modification is needed.
* Each stage writes into `results/connected_workflow/` and registers its output in `pipeline_manifest.json`. **The next notebook reads that actual artifact**. Input datasets are never modified.

## The connected stages

| Notebook | Previous input → written output | Important condition |
|---|---|---|
| 00 | 7 NetCDF → network map, characteristics, availability monthly CSV/PNG | five-sensor source availability script unmodified |
| 01 | released native 1-min rain → `rainfall_1min` + derived 5-min `rainfall` manifest entries (optional full-history retrieval) | default registers the published 1-min release and its exact 5-min derivation; opt-in expensive original Graf producer |
| 02 | staged rainfall → review rain time series; optional one-link comparison and gauge validation | gauge validation uses the Graf helper v3 (as for the released numbers) behind an explicit switch |
| 03 | stage-01 **5-min Graf rainfall** + QC geometry → `R(time,y,x)` IDW grid NetCDF | reuses unchanged original IDW mathematical functions; **adapter, not the Fig. 7 grid builder** |
| 04 | stage-03 grid + separately configured *native* RADKLIM daily YW → native/hourly map PNG/PDF | `cartopy` + native YW external directory required; if exact selected-event cache absent, opt-in explicit demonstration selection from grid extent is clearly marked (not historical event discovery) |
| 05 | raw CML + gauge/temp → V6 WVD NC, or register original precomputed V6 | full 2020–21 history affects calibration; original V6 supporting product is not in seven main files |
| 06 | stage-05 V6 + **extra quantization and diagnostics CSVs** + CHMI humidity/temp → binned sensitivity metrics and station-matched pairs | original `33_v6_binned_validation_hexbins.py`, no invented data |
| 07 | stage-06 pairs/links → original paper-style filtered WVD hexbins and metrics | GRAF only; L≥4km, quantization Δρ≤0.20g/m³, disjoint E-band |
| 08 | optional CHMI original v5 CSV downloader/converter → upstream CHMI NetCDF4 | not needed just to plot the downloaded seven-product release |

## Scientific and availability caveats

- **Do not confuse longer 2020–21 source files on the author's HPC with the April–October 2021 seven-file release**. The released rainfall thresholds and baselines were computed over the released QC period (1 April to 31 October 2021), so recomputation should use the released QC file.
- **Final rainfall is Graf RSTD + constant baseline + Pastorek WAA + ITU-R 2005**, not ITU2003 or radar-assisted wet-dry.
- **Spatial adapter** follows original IDW code and applies it to the derived 5-min series of the released rainfall. The paper's Figure 7 grid is built with `reconstruction_2026/05_figures/fig7_grid_builder/`; compare outputs scientifically before asserting an exact Figure 7 match.
- The native DWD RADKLIM YW radar fields are not the path-averaged radar NetCDF among the seven. The user must supply their own native YW radar directory, licensed appropriately.
- Water vapour V6 quantization CSVs/validation cache are separate supporting inputs. The QC producer is `reconstruction_2026/02_qc/scripts/s3_04_qc_corrected_rev1.py`; preserve provenance.
- Original code and source checksums are preserved, while the **new wrappers are clearly separate** in `src/essd_cz/connected.py`. Data licences: CML CC BY-NC 4.0; CHMI/DWD CC BY 4.0. Software licence still requires approval.

## Fast synthetic checks

```bash
python tools_verify_originals.py
python -m pytest -q tests/test_mini_structure.py tests/test_comparison_alignment.py tests/test_notebook_adapters.py tests/test_connected_chain.py
```

The local tests confirm plumbing and synthetic NetCDF compatibility only. Run heavyweight original map/V6 and gauge validation on HPC in the provided example environment before declaring full paper reproduction.
