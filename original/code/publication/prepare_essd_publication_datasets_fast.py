#!/usr/bin/env python3
"""
prepare_essd_publication_datasets.py
====================================

Prepare publication-ready Czech precipitation datasets using the temporal
extent of the preprocessed CML dataset as the MASTER time window.

The lower and upper time limits are read directly from cml_preprocessed_v3.nc.
No date range is hard-coded.

Outputs
-------
CZ_CML_RAW_1MIN_<START>_<END>.nc
CZ_CML_PREPROCESSED_1MIN_<START>_<END>.nc
CZ_CML_RAINRATE_5MIN_<START>_<END>.nc
CZ_RADAR_CMLPATH_5MIN_<START>_<END>.nc
CZ_RAIN_GAUGE_10MIN_<START>_<END>.nc

Rules
-----
1. All datasets are cropped to the inclusive MASTER CML time limits.
2. Native temporal resolutions are preserved; no resampling is done here.
3. detection_limit_mmhr is moved to coordinates in the preprocessed CML file.
4. length, frequency and polarization are moved to coordinates in the CML
   rainfall file.
5. CMLs with valid rain rate are defined as links having at least one finite
   R_mean_5min value inside the MASTER time window.
6. The published radar path-average file is restricted to those valid-rain-rate
   CML IDs.
7. The published CML rainfall file is also restricted to the same valid-rain-rate
   CML IDs.
8. Raw and preprocessed CML files preserve their own CML populations and are
   only standardized temporally.
9. Saving is optimized for large publication files: Dask streaming, bounded I/O
   workers, larger CML chunks, compression level 1 by default, and no compression
   overhead for tiny variables.
"""

from __future__ import annotations

import argparse
import json
import time
import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import xarray as xr


def require_dask() -> None:
    try:
        import dask.array  # noqa: F401
    except Exception as exc:
        raise RuntimeError(
            "This script requires dask for memory-safe writing of the large NetCDF files. "
            "Install it with: python -m pip install dask"
        ) from exc


def open_lazy(path: str, time_dim: Optional[str], time_chunk: int, cml_chunk: int) -> xr.Dataset:
    with xr.open_dataset(path) as tmp:
        chunks: Dict[str, int] = {}
        if time_dim and time_dim in tmp.dims:
            chunks[time_dim] = min(int(time_chunk), int(tmp.sizes[time_dim]))
        if "cml_id" in tmp.dims:
            chunks["cml_id"] = min(int(cml_chunk), int(tmp.sizes["cml_id"]))
    return xr.open_dataset(path, chunks=chunks)


def normalize_cml_id(ds: xr.Dataset) -> xr.Dataset:
    if "cml_id" in ds.coords:
        return ds.assign_coords(cml_id=ds["cml_id"].astype(str))
    return ds


def get_time_bounds(path: str, time_dim: str) -> Tuple[pd.Timestamp, pd.Timestamp]:
    with xr.open_dataset(path) as ds:
        if time_dim not in ds.coords and time_dim not in ds.variables:
            raise KeyError(f"{path}: {time_dim!r} not found.")
        t = pd.DatetimeIndex(ds[time_dim].values)

    if t.size == 0:
        raise ValueError(f"{path}: empty time coordinate.")
    t = t[~pd.isna(t)]
    if t.size == 0:
        raise ValueError(f"{path}: all timestamps are missing.")
    return pd.Timestamp(t.min()), pd.Timestamp(t.max())


def crop_time(ds: xr.Dataset, time_dim: str, tmin: pd.Timestamp, tmax: pd.Timestamp) -> xr.Dataset:
    if time_dim not in ds.coords and time_dim not in ds.variables:
        raise KeyError(f"{time_dim!r} not found. Available: {list(ds.variables)}")
    return ds.sel({time_dim: slice(np.datetime64(tmin), np.datetime64(tmax))})


def move_to_coords(ds: xr.Dataset, names: Sequence[str]) -> xr.Dataset:
    missing = [name for name in names if name not in ds.variables]
    if missing:
        raise KeyError(f"Missing variables to convert to coordinates: {missing}")
    return ds.set_coords(list(names))


def clean_source_encoding(ds: xr.Dataset) -> xr.Dataset:
    bad = {
        "chunksizes", "preferred_chunks", "source", "original_shape",
        "contiguous", "compression", "compression_opts", "zlib",
        "complevel", "shuffle", "fletcher32"
    }
    for name in ds.variables:
        enc = dict(ds[name].encoding)
        for key in bad:
            enc.pop(key, None)
        ds[name].encoding = enc
    return ds


def build_encoding(
    ds: xr.Dataset,
    time_dims: Sequence[str],
    compression_level: int,
    time_chunk: int,
    cml_chunk: int,
    compress_min_elements: int,
) -> Dict[str, Dict]:
    """Fast publication encoding: compress only large numeric arrays."""
    encoding: Dict[str, Dict] = {}

    for name, var in ds.variables.items():
        if var.dtype.kind in {"O", "U", "S"}:
            continue

        e: Dict = {}
        n_elem = int(var.size)
        if var.ndim > 0 and compression_level > 0 and n_elem >= int(compress_min_elements):
            # Level 1 is usually much faster than level 4 while still saving substantial space.
            e.update(zlib=True, complevel=int(compression_level), shuffle=True)

        chunksizes = []
        for dim in var.dims:
            size = int(ds.sizes[dim])
            if dim in time_dims:
                chunksizes.append(max(1, min(int(time_chunk), size)))
            elif dim == "cml_id":
                chunksizes.append(max(1, min(int(cml_chunk), size)))
            else:
                chunksizes.append(size)

        if chunksizes:
            e["chunksizes"] = tuple(chunksizes)

        encoding[name] = e

    return encoding


def add_publication_attrs(
    ds: xr.Dataset,
    role: str,
    source_file: str,
    tmin: pd.Timestamp,
    tmax: pd.Timestamp,
) -> xr.Dataset:
    attrs = dict(ds.attrs)
    attrs.update(
        {
            "publication_dataset_role": role,
            "publication_source_file": str(source_file),
            "publication_time_start": tmin.isoformat(),
            "publication_time_end": tmax.isoformat(),
            "publication_time_rule": (
                "Inclusive crop to the lower and upper limits of the "
                "preprocessed CML master time coordinate; native temporal "
                "resolution preserved."
            ),
            "publication_preparation": "prepare_essd_publication_datasets.py",
        }
    )
    ds.attrs = attrs
    return ds


def write_dataset(
    ds: xr.Dataset,
    out_path: Path,
    time_dims: Sequence[str],
    compression_level: int,
    time_chunk: int,
    cml_chunk: int,
    overwrite: bool,
    io_workers: int,
    compress_min_elements: int,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        if overwrite:
            out_path.unlink()
        else:
            raise FileExistsError(f"Output exists: {out_path}. Use --overwrite.")

    ds = clean_source_encoding(ds)
    encoding = build_encoding(
        ds, time_dims, compression_level, time_chunk, cml_chunk, compress_min_elements
    )

    tmp = out_path.with_suffix(out_path.suffix + ".partial")
    if tmp.exists():
        tmp.unlink()

    try:
        print(f"\nWriting: {out_path}", flush=True)
        print(
            f"  save settings: compression={compression_level}, "
            f"time_chunk={time_chunk}, cml_chunk={cml_chunk}, io_workers={io_workers}",
            flush=True,
        )

        # Build the write graph first, then compute with a bounded thread pool.
        # NetCDF writes are internally locked, but parallel reads/decompression can still overlap.
        import dask
        from dask.diagnostics import ProgressBar

        delayed = ds.to_netcdf(
            tmp,
            engine="netcdf4",
            format="NETCDF4",
            encoding=encoding,
            compute=False,
        )
        with dask.config.set(scheduler="threads", num_workers=max(1, int(io_workers))):
            with ProgressBar(minimum=1.0):
                delayed.compute()

        tmp.replace(out_path)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise


def finite_rainrate_cml_ids(
    rain_ds: xr.Dataset,
    rain_var: str,
    cml_dim: str,
    time_dim: str,
    tmin: pd.Timestamp,
    tmax: pd.Timestamp,
) -> np.ndarray:
    if rain_var not in rain_ds.data_vars:
        raise KeyError(
            f"{rain_var!r} not found. Available data variables: {list(rain_ds.data_vars)}"
        )

    da = rain_ds[rain_var].sel(
        {time_dim: slice(np.datetime64(tmin), np.datetime64(tmax))}
    )
    if cml_dim not in da.dims:
        raise ValueError(f"{rain_var} does not contain {cml_dim}. Dims: {da.dims}")

    reduce_dims = [d for d in da.dims if d != cml_dim]
    print(
        f"\nScanning {rain_var} for CMLs with at least one finite rain-rate value ...",
        flush=True,
    )
    valid = np.isfinite(da).any(dim=reduce_dims).compute()
    mask = np.asarray(valid.values, dtype=bool)

    ids = np.asarray(rain_ds[cml_dim].values).astype(str)
    if mask.shape != ids.shape:
        raise RuntimeError(f"Valid-mask shape {mask.shape} != CML ID shape {ids.shape}")
    return ids[mask]


def ordered_intersection(reference_ids: Sequence[str], available_ids: Sequence[str]) -> List[str]:
    aset = set(str(x) for x in available_ids)
    return [str(x) for x in reference_ids if str(x) in aset]


def verify_coord_status(path: Path, expected_coords: Sequence[str]) -> None:
    with xr.open_dataset(path, decode_coords="all") as ds:
        missing = [x for x in expected_coords if x not in ds.coords]
        if missing:
            raise RuntimeError(
                f"{path}: expected coordinates not recognized by xarray: {missing}. "
                f"coords={list(ds.coords)}"
            )


def dataset_summary(path: Path, time_dim: str, cml_dim: Optional[str]) -> Dict:
    with xr.open_dataset(path, decode_coords="all") as ds:
        times = pd.DatetimeIndex(ds[time_dim].values)
        return {
            "file": str(path),
            "time_dimension": time_dim,
            "first_timestamp": str(times[0]) if len(times) else "",
            "last_timestamp": str(times[-1]) if len(times) else "",
            "n_time": int(len(times)),
            "n_cml": int(ds.sizes[cml_dim]) if cml_dim and cml_dim in ds.dims else np.nan,
            "coordinates": ",".join(ds.coords),
            "data_variables": ",".join(ds.data_vars),
            "size_gib": path.stat().st_size / (1024 ** 3),
        }


def prepare(args: argparse.Namespace) -> None:
    require_dask()
    t0 = time.perf_counter()

    master_start, master_end = get_time_bounds(
        args.cml_preprocessed, args.cml_time_dim
    )
    period = f"{master_start:%Y%m%d}_{master_end:%Y%m%d}"

    print("=" * 100)
    print("ESSD PUBLICATION DATASET PREPARATION")
    print("=" * 100)
    print(f"Master CML file    : {args.cml_preprocessed}")
    print(f"Master lower limit : {master_start}")
    print(f"Master upper limit : {master_end}")
    print(f"Period tag         : {period}")
    print("=" * 100)

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "raw_cml": out_dir / f"CZ_CML_RAW_1MIN_{period}.nc",
        "preprocessed_cml": out_dir / f"CZ_CML_PREPROCESSED_1MIN_{period}.nc",
        "cml_rainrate": out_dir / f"CZ_CML_RAINRATE_5MIN_{period}.nc",
        "radar": out_dir / f"CZ_RADAR_CMLPATH_5MIN_{period}.nc",
        "gauge": out_dir / f"CZ_RAIN_GAUGE_10MIN_{period}.nc",
    }

    # Rain-rate file first: determine valid rain-rate CML IDs.
    rain = open_lazy(
        args.cml_rainrate, args.rain_time_dim,
        args.time_chunk_5min, args.cml_chunk
    )
    rain = normalize_cml_id(rain)

    valid_rain_ids = finite_rainrate_cml_ids(
        rain, args.rain_valid_var, args.cml_dim, args.rain_time_dim,
        master_start, master_end
    )
    print(f"Rain-rate CML IDs in source : {rain.sizes[args.cml_dim]:,}")
    print(f"CML IDs with valid rain rate: {len(valid_rain_ids):,}")

    # 1. Raw CML: temporal crop only.
    raw = open_lazy(
        args.cml_raw, args.cml_time_dim,
        args.time_chunk_1min, args.cml_chunk
    )
    raw = normalize_cml_id(raw)
    raw = crop_time(raw, args.cml_time_dim, master_start, master_end)
    raw = add_publication_attrs(
        raw,
        "Raw commercial microwave link observations",
        args.cml_raw,
        master_start,
        master_end,
    )
    write_dataset(
        raw, outputs["raw_cml"], [args.cml_time_dim],
        args.compression_level, args.time_chunk_1min,
        args.cml_chunk, args.overwrite, args.io_workers, args.compress_min_elements
    )
    raw.close()

    # 2. Preprocessed CML: temporal crop + detection limit as coordinate.
    pre = open_lazy(
        args.cml_preprocessed, args.cml_time_dim,
        args.time_chunk_1min, args.cml_chunk
    )
    pre = normalize_cml_id(pre)
    pre = crop_time(pre, args.cml_time_dim, master_start, master_end)
    pre = move_to_coords(pre, [args.detection_limit_var])
    pre = add_publication_attrs(
        pre,
        "Quality-controlled commercial microwave link observations",
        args.cml_preprocessed,
        master_start,
        master_end,
    )
    pre.attrs["publication_auxiliary_coordinates"] = args.detection_limit_var
    write_dataset(
        pre, outputs["preprocessed_cml"], [args.cml_time_dim],
        args.compression_level, args.time_chunk_1min,
        args.cml_chunk, args.overwrite, args.io_workers, args.compress_min_elements
    )
    pre.close()

    # 3. CML rain rate: temporal crop + valid rain IDs + metadata as coordinates.
    rain = crop_time(rain, args.rain_time_dim, master_start, master_end)
    valid_rain_ids = np.asarray(
        ordered_intersection(
            valid_rain_ids,
            np.asarray(rain[args.cml_dim].values).astype(str),
        ),
        dtype=str,
    )
    rain = rain.sel({args.cml_dim: valid_rain_ids})
    rain = move_to_coords(
        rain,
        [args.length_var, args.frequency_var, args.polarization_var],
    )
    rain = add_publication_attrs(
        rain,
        "Derived CML rain-rate estimates with co-located path-averaged radar reference",
        args.cml_rainrate,
        master_start,
        master_end,
    )
    rain.attrs["publication_cml_selection"] = (
        f"Links with at least one finite {args.rain_valid_var} in the publication time window."
    )
    rain.attrs["publication_auxiliary_coordinates"] = (
        f"{args.length_var},{args.frequency_var},{args.polarization_var}"
    )
    write_dataset(
        rain, outputs["cml_rainrate"], [args.rain_time_dim],
        args.compression_level, args.time_chunk_5min,
        args.cml_chunk, args.overwrite, args.io_workers, args.compress_min_elements
    )
    rain.close()

    # 4. Radar: temporal crop + only valid rain-rate CML IDs.
    radar = open_lazy(
        args.radar, args.radar_time_dim,
        args.time_chunk_5min, args.cml_chunk
    )
    radar = normalize_cml_id(radar)
    radar = crop_time(radar, args.radar_time_dim, master_start, master_end)

    radar_ids = np.asarray(radar[args.cml_dim].values).astype(str)
    radar_keep = ordered_intersection(valid_rain_ids, radar_ids)
    missing = sorted(set(valid_rain_ids.tolist()) - set(radar_keep))
    if missing:
        print(
            f"WARNING: {len(missing)} valid-rain-rate CML IDs are absent from radar.",
            flush=True,
        )
    radar = radar.sel({args.cml_dim: radar_keep})
    radar = add_publication_attrs(
        radar,
        "Path-averaged radar precipitation reference for CMLs with valid derived rain rate",
        args.radar,
        master_start,
        master_end,
    )
    radar.attrs["publication_cml_selection"] = (
        f"cml_id restricted to links with at least one finite {args.rain_valid_var}."
    )
    write_dataset(
        radar, outputs["radar"], [args.radar_time_dim],
        args.compression_level, args.time_chunk_5min,
        args.cml_chunk, args.overwrite, args.io_workers, args.compress_min_elements
    )
    radar.close()

    # 5. Gauge: temporal crop only.
    gauge = open_lazy(
        args.gauge, args.gauge_time_dim,
        args.time_chunk_10min, args.cml_chunk
    )
    gauge = crop_time(gauge, args.gauge_time_dim, master_start, master_end)
    gauge = add_publication_attrs(
        gauge,
        "CHMI rain-gauge precipitation observations",
        args.gauge,
        master_start,
        master_end,
    )
    write_dataset(
        gauge, outputs["gauge"], [args.gauge_time_dim],
        args.compression_level, args.time_chunk_10min,
        args.cml_chunk, args.overwrite, args.io_workers, args.compress_min_elements
    )
    gauge.close()

    # Verify requested coordinate roles.
    verify_coord_status(
        outputs["preprocessed_cml"],
        [args.detection_limit_var],
    )
    verify_coord_status(
        outputs["cml_rainrate"],
        [args.length_var, args.frequency_var, args.polarization_var],
    )

    # Inventory + manifest.
    rows = []
    rows.append({"dataset": "raw_cml", **dataset_summary(outputs["raw_cml"], args.cml_time_dim, args.cml_dim)})
    rows.append({"dataset": "preprocessed_cml", **dataset_summary(outputs["preprocessed_cml"], args.cml_time_dim, args.cml_dim)})
    rows.append({"dataset": "cml_rainrate", **dataset_summary(outputs["cml_rainrate"], args.rain_time_dim, args.cml_dim)})
    rows.append({"dataset": "radar_cmlpath", **dataset_summary(outputs["radar"], args.radar_time_dim, args.cml_dim)})
    rows.append({"dataset": "rain_gauge", **dataset_summary(outputs["gauge"], args.gauge_time_dim, None)})

    inventory = out_dir / f"CZ_ESSD_DATASET_INVENTORY_{period}.csv"
    pd.DataFrame(rows).to_csv(inventory, index=False)

    manifest = {
        "master_time_source": str(args.cml_preprocessed),
        "master_time_start": master_start.isoformat(),
        "master_time_end": master_end.isoformat(),
        "period_tag": period,
        "valid_rainrate_cml_count": int(len(valid_rain_ids)),
        "files": {k: str(v) for k, v in outputs.items()},
        "inventory_csv": str(inventory),
    }
    manifest_path = out_dir / f"CZ_ESSD_DATASET_MANIFEST_{period}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("\n" + "=" * 100)
    print("FINAL DATASET INVENTORY")
    print("=" * 100)
    print(pd.DataFrame(rows)[["dataset", "first_timestamp", "last_timestamp", "n_time", "n_cml"]].to_string(index=False))
    print("=" * 100)
    print(f"Inventory: {inventory}")
    print(f"Manifest : {manifest_path}")
    print(f"Elapsed  : {(time.perf_counter() - t0) / 60.0:.1f} min")
    print("=" * 100)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare publication-ready ESSD Czech CML/radar/gauge/rain-rate NetCDF files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("--cml-preprocessed", default="/home/jovyan/cml_preprocessed_v3.nc")
    p.add_argument("--cml-raw", default="/home/jovyan/CML_S2020E2021_1MIN_CZ.nc")
    p.add_argument("--radar", default="/home/jovyan/RADAR_S2021E2021_5MIN_CZ.nc")
    p.add_argument("--gauge", default="/home/jovyan/CHMI_10MIN_v1.nc")
    p.add_argument("--cml-rainrate", default="/home/jovyan/cml_rainrate_S2020_E2021_valid_cml_radar.nc")

    p.add_argument("--out-dir", default="/home/jovyan/ESSD_publication_dataset")

    p.add_argument("--cml-time-dim", default="time")
    p.add_argument("--radar-time-dim", default="time")
    p.add_argument("--gauge-time-dim", default="time")
    p.add_argument("--rain-time-dim", default="time_5min")
    p.add_argument("--cml-dim", default="cml_id")

    p.add_argument("--detection-limit-var", default="detection_limit_mmhr")
    p.add_argument("--rain-valid-var", default="R_mean_5min")
    p.add_argument("--length-var", default="length")
    p.add_argument("--frequency-var", default="frequency")
    p.add_argument("--polarization-var", default="polarization")

    p.add_argument("--time-chunk-1min", type=int, default=10080)
    p.add_argument("--time-chunk-5min", type=int, default=2016)
    p.add_argument("--time-chunk-10min", type=int, default=1008)
    p.add_argument("--cml-chunk", type=int, default=64, help="Larger link chunks reduce NetCDF chunk-count and speed sequential publication writes.")
    p.add_argument("--compression-level", type=int, choices=range(0, 10), default=1, help="Level 1 is the fast default. Use 0 for maximum write speed or 4 for smaller files.")
    p.add_argument("--io-workers", type=int, default=4, help="Bounded Dask thread pool used during chunked writes. 2-4 is usually best on shared storage.")
    p.add_argument("--compress-min-elements", type=int, default=100000, help="Do not spend compression CPU on tiny numeric variables/coordinates.")
    p.add_argument("--overwrite", action="store_true")

    return p.parse_args()


if __name__ == "__main__":
    prepare(parse_args())
