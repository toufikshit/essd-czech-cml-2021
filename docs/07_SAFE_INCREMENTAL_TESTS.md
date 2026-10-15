# 7. Safe incremental tests — do not refactor scientific algorithms blindly

## What is fixed, and what may change?

**Frozen:** all 44 recovered original scripts/notebooks under `original/`, including
`original/code/cml_rainfall/cml_graf_itu2005_rebuild.py`. Their checksums must
continue to match the original recovery manifest. Do not edit these files to
make a comparison pass; if they differ, record and investigate the discrepancy.

**Allowed new code:** `src/essd_cz/` is the thin layer for input paths, format
checking, ID/timestamp alignment, plotting entry points and metrics. These
modules call the original functions rather than reimplementing Graf, baseline,
Pastorek or ITU2005. Test code lives in `tests/` and writes to new `results/`
folders only. The original QC producer from March 2026 remains unconfirmed.

## Test fixture: comparable structure, NOT a claim of scientific identity

`tests/mini_fixture.py` creates **seven tiny synthetic NetCDF files** with the
original filenames. The simulated input features: physical structural CML IDs
with underscores; three CML links and two directions (raw has an extra link);
1-min TRSL/TX/RX; 5-min time-start-labelled rainfall and radar; 10-min gauge
amounts and CHMI temperature/RH; lengths in km, frequencies in **MHz**;
wet episodes, finite-detection gaps, a wholly missing direction, 255 absent
wet/dry flags, float NaNs, and shuffled-ID scenarios. Only **2 days** of
artificial signals are generated; the seven products have the **same schema
and units**, not the full published entity counts or time coverage.

The fixture is named `SYNTHETIC`; it must never be uploaded as scientific
measurements. On HPC, it is stored in **NETCDF4**; if offline `netCDF4` is
unavailable, a NETCDF3 fallback preserves the values/dimensions but stores
flag arrays as signed 16-bit because the older container has no unsigned byte.
This fallback is for testing structure, NOT the exact published on-disk type.

## Run in four gates — proceed only when the previous gate passed

### Gate 0 — original integrity (seconds)

```bash
python tools_verify_originals.py
```

Expected: `Verified 44/44` and `PASS`. Do not proceed if any original hash
changes unexpectedly.

### Gate 1 — schema, ID, time, missing-data and comparator regression (seconds)

```bash
python -m pytest -q tests/test_mini_structure.py tests/test_comparison_alignment.py
```

Offline tests use synthetic placeholder rain values to isolate *data plumbing*:
we require structural IDs and time labels to match, 255 to stay distinct from
0/1, NaNs to remain meaningful, threshold/mean fields to align, a shuffled
link/direction order to pass, and **deliberately corrupted** values, masks,
flags and timestamps to be rejected. No `pycomlink` is needed at this gate;
**these tests do NOT verify rainfall physics.**

### Gate 2 — original availability plot on synthetic data (seconds)

```bash
python tests/run_mini_availability.py --out-dir results/mini_availability
```

This runs **unchanged** `availability_native_timestep_percentage.py` on all
five small sensor inputs, producing five CSVs, a PNG and a PDF. Seven input
checksums must remain identical. Every test output goes to a new folder; the
runner refuses overwriting a previous test.

### Gate 3 — actual original scientific functions on synthetic data (HPC only)

Activate the scientific environment with `netCDF4`, `numba`, `pycomlink`,
`numpy`, `xarray`, and `pandas`. Then:

```bash
PY=/home/jovyan/my-conda-envs/tf-env/bin/python
$PY tests/run_mini_science.py --out-dir results/mini_hpc_science
```

This builds a 2-day synthetic reference using the **original `graf()`,
`estimate()` and `meta2()`**, directly following the original producer's
five-minute aggregation. It then recomputes two links via the wrapper, compares
both directions and mean rainfall, checks all rain masks, wet flags, valid/wet
minute counts, Graf Q80 and thresholds. It also checks deliberately shuffled
CML IDs/directional order, rejected numerical corruption, rejected flag
corruption and unmodified input file hashes. No release data are opened.

It writes `MINI_SCIENCE_RESULTS.json` only after all gates pass. If it fails,
stop, inspect the report, and fix **the wrapper only** when that is the cause.
Do not edit the original formulas to suppress differences. **Gate 3 was executed successfully by the user on the HPC. It cannot be rerun in this
local build environment, which lacks netCDF4 and pycomlink. Re-run it after
any new scientific-wrapper changes.**

### Gate 4 — real-file numerical comparison (HPC only; after Gate 3 passes)

```bash
PY=/home/jovyan/my-conda-envs/tf-env/bin/python
$PY run.py --config config/hpc_source.example.json check
$PY run.py --config config/hpc_source.example.json \
    compare-rainfall --link-index 3 --days 14 --threads 2 \
    --no-plot
```

The Q80 Graf threshold and constant baseline are computed from **the entire
available 2020–2021 QC history** for that one link, *then* 1–14 April 2021
5-minute results are compared. Both directional series, wet flags and masks
must be checked. No source arrays are written. The first command requires
full 2021-04-01 to 2021-10-31 coverage but accepts longer source files.

The earlier April-only smoke test shows the algorithm runs; it is **not**
proof of numerical identity to the actual source output. Even a passing
2-day synthetic test does not prove historical numerical reproduction.
If the full-period comparison differs, record the exact disagreement before
making a scientific or version change. Many published scripts require
additional intermediate data, and the exact March 2026 QC producer remains
unconfirmed.

## After each intended wrapper change

Repeat **Gate 0 → 1 → 2 → 3**. Do not re-run a full 2219-link rainfall rebuild
just to test a path change. Keep the synthetic fixture files outside the
published data directories. Commit/record the source hash, test report and
Python/environment versions. Only then run Gate 4 again as needed.
