# Independent test status — 13 October 2026

| Gate | Executed here? | Outcome | Interpretation |
|---|---|---|---|
| Freeze & hash recovered originals | Yes | **PASS: 44/44** | Byte-for-byte original code/notebooks preserved |
| Mini seven-file schema / alignment / numerical guardrails | Yes | **PASS: original 9 cases; 15 cases including notebook adapters** | Structural IDs, axes, units, masks, reordered links/directions and negative controls checked |
| Unchanged five-sensor availability (2-day fixture) | Yes | **PASS** | Produced five CSVs, PNG, PDF; all seven synthetic input hashes unchanged |
| Graf + pycomlink + Pastorek + ITU2005 synthetic numerical identity | **Yes, on user HPC** | **PASS** | Two directions, missing-data cases, reordering and deliberate errors tested |
| Real published-rainfall identity (2020–21 QC history) | **Yes, on user HPC** | **PASS: one physical link** | Link `10004_311866_3129`: 7,666 shared directional samples, RMSE 0, zero mask/flag mismatches; not population-wide |
| QC provenance | N/A | **UNCONFIRMED** | Exact March 2026 producer script has not been recovered |

The local runtime does not have `netCDF4` or `pycomlink`. No approximate
scientific implementation has been substituted for these libraries.

**Regression found and fixed before packaging:** in the first wrapper version,
one all-missing directional series raised an exception and discarded the
opposite direction's comparison. The unchanged original producer instead
records the error and writes NaN for that direction while processing the
other. The wrapper was changed *only* to copy that original error-handling
policy and expose `estimate_errors` in its result/report. The corresponding
new regression test failed before the fix and passes after it.

The tiny fixture is **not the actual publication dataset**; it is synthetic
and contains 3 QC links over 2 days, versus the real 2,219 links over the
April–October release period (source inputs may extend to 2020). The offline
NETCDF3 fallback stores flags with 16-bit integers; the real HPC fixture will
use NETCDF4 unsigned byte flags. Passing these tests does **not** certify the
historical Q80 threshold or rainfall product.

## Commands, in order

```bash
python tools_verify_originals.py
python -m pytest -q tests/test_mini_structure.py tests/test_comparison_alignment.py
python tests/run_mini_availability.py --out-dir results/mini_availability
```

Once these pass on your HPC:

```bash
PY=/home/jovyan/my-conda-envs/tf-env/bin/python
$PY tests/run_mini_science.py --out-dir results/mini_hpc_science
```

The scientific synthetic gate and one real physical-link comparison have already passed on the user HPC. Re-run these checks when wrappers change. See `docs/07_SAFE_INCREMENTAL_TESTS.md` and `TEST_LOG_NOTEBOOKS.md`.

**Notebook extension:** Eight notebooks and two new additive modules passed 15 local pytest cases; see `TEST_LOG_NOTEBOOKS.md`. Full CHMI v5 converter and external-asset figures still require HPC tests.
