#!/usr/bin/env python3
"""
cml_grid_essd_publication_safe_singleopen.py
============================================

Safe single-open IDW gridding for the filtered ESSD CML publication product.

Why this script exists
----------------------
The merged gridding script opens the same estimation NetCDF inside every grid
job. On this HPC setup the first chunk succeeds, but a later repeated netCDF4/
HDF5 open can terminate Python with SIGSEGV ("NetCDF: HDF error"), even with
--grid-workers 1.

This script preserves the SAME scientific gridding choices but removes that
failure mode:

  * input publication rain-rate NetCDF is opened ONCE;
  * coordinate NetCDF is opened ONCE and closed after coordinates are loaded;
  * output NetCDF is opened ONCE;
  * time chunks are read sequentially from the existing input handle;
  * no ProcessPoolExecutor;
  * no repeated xr.open_dataset() calls;
  * scipy cKDTree query uses workers=1 by default.

Scientific method
-----------------
Input:
    R_mean_5min
    baseline_method = pycml_constant
    kr_model        = ITU_2003
    waa_method      = pastorek_2021

Spatial support:
    CML link midpoint

IDW:
    nnear       = 15
    p           = 2
    max_distance= 0.3 degree
    resolution  = 0.01 degree

NaN handling matches the supplied fast gridding code:
    for each unique finite-link mask, only finite links are used to build the
    KD tree; the nnear closest valid links are weighted by 1/distance**p;
    weights are normalized; exact link/grid coincidences use the link value;
    cells farther than max_distance from the closest valid link are NaN.

Output:
    R(time, y, x), units mm h-1

No radar and no gauge are gridded here, by design.
"""

from __future__ import annotations

import os
# Important on some shared HPC filesystems. Must be set before netCDF4/HDF5 loads.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import argparse
import gc
import math
import time
from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import xarray as xr

try:
    from tqdm import tqdm
except Exception:
    def tqdm(x, **kwargs):
        return x


def build_grid(lon: np.ndarray, lat: np.ndarray, resolution: float) -> Tuple[np.ndarray, np.ndarray]:
    xc = np.arange(float(np.nanmin(lon)), float(np.nanmax(lon)), float(resolution))
    yc = np.arange(float(np.nanmin(lat)), float(np.nanmax(lat)), float(resolution))
    return np.meshgrid(xc, yc)


def choose_input_engine(requested: str) -> Optional[str]:
    if requested == "netcdf4":
        return "netcdf4"
    if requested == "h5netcdf":
        try:
            import h5netcdf  # noqa: F401
            return "h5netcdf"
        except Exception as exc:
            raise RuntimeError("--input-engine h5netcdf requested but h5netcdf is unavailable.") from exc
    # auto: prefer h5netcdf to avoid netCDF4 C-extension file-manager instability
    try:
        import h5netcdf  # noqa: F401
        return "h5netcdf"
    except Exception:
        return "netcdf4"


def open_input(path: Path, engine: Optional[str]) -> xr.Dataset:
    print(f"  Opening input with engine={engine}")
    return xr.open_dataset(
        path,
        engine=engine,
        cache=False,
        decode_times=True,
    )


def as_string_values(da: xr.DataArray) -> list[str]:
    out = []
    for v in np.asarray(da.values).reshape(-1):
        if isinstance(v, (bytes, np.bytes_)):
            out.append(v.decode("utf-8", errors="ignore"))
        else:
            out.append(str(v))
    return out


def select_combo(
    ds: xr.Dataset,
    rain_var: str,
    baseline_method: str,
    kr_model: str,
    waa_method: str,
) -> xr.DataArray:
    if rain_var not in ds:
        raise KeyError(f"{rain_var!r} not found. Available: {list(ds.data_vars)}")

    da = ds[rain_var]
    selections = {
        "baseline_method": baseline_method,
        "kr_model": kr_model,
        "waa_method": waa_method,
    }

    for dim, value in selections.items():
        if dim not in da.dims:
            continue
        values = as_string_values(ds[dim])
        if value not in values:
            raise ValueError(f"{value!r} not found on {dim!r}. Available: {values}")
        # isel avoids any byte/string indexing surprises.
        da = da.isel({dim: values.index(value)})

    required = {"cml_id", "time_5min"}
    missing = required.difference(da.dims)
    if missing:
        raise ValueError(f"{rain_var!r} is missing dimensions {sorted(missing)}; dims={da.dims}")

    extra = [d for d in da.dims if d not in required]
    if extra:
        raise ValueError(
            f"After selecting the requested method combination, unexpected dimensions remain: {extra}"
        )

    return da


def load_midpoints(coords_nc: Path, target_ids: np.ndarray, args) -> Tuple[np.ndarray, np.ndarray]:
    with xr.open_dataset(coords_nc, cache=False) as cds:
        if "cml_id" not in cds.coords:
            raise KeyError("Coordinate file must contain cml_id coordinate.")

        source_ids = np.asarray(cds["cml_id"].values).astype(str)
        pos = {cid: i for i, cid in enumerate(source_ids)}

        required = [
            args.site_a_lon, args.site_a_lat,
            args.site_b_lon, args.site_b_lat,
        ]
        missing = [name for name in required if name not in cds]
        if missing:
            raise KeyError(f"Coordinate file is missing variables: {missing}")

        lon_a_all = np.asarray(cds[args.site_a_lon].values, dtype=np.float64).reshape(-1)
        lat_a_all = np.asarray(cds[args.site_a_lat].values, dtype=np.float64).reshape(-1)
        lon_b_all = np.asarray(cds[args.site_b_lon].values, dtype=np.float64).reshape(-1)
        lat_b_all = np.asarray(cds[args.site_b_lat].values, dtype=np.float64).reshape(-1)

        lon = np.full(len(target_ids), np.nan, dtype=np.float64)
        lat = np.full(len(target_ids), np.nan, dtype=np.float64)

        for i, cid in enumerate(target_ids):
            j = pos.get(str(cid))
            if j is None:
                continue
            lon[i] = 0.5 * (lon_a_all[j] + lon_b_all[j])
            lat[i] = 0.5 * (lat_a_all[j] + lat_b_all[j])

    return lon, lat


def resolve_time_window(
    times: pd.DatetimeIndex,
    tstart: Optional[str],
    tend: Optional[str],
) -> Tuple[int, int, pd.DatetimeIndex]:
    start = pd.Timestamp(tstart) if tstart else times[0]
    end = pd.Timestamp(tend) if tend else times[-1]
    if end < start:
        raise ValueError(f"End time {end} is before start time {start}.")

    i0 = int(np.searchsorted(times.values, np.datetime64(start), side="left"))
    i1 = int(np.searchsorted(times.values, np.datetime64(end), side="right"))
    if i1 <= i0:
        raise ValueError(
            f"No input times overlap {start} -> {end}. "
            f"Available: {times[0]} -> {times[-1]}"
        )
    return i0, i1, times[i0:i1]


def idw_fast_series(
    lon: np.ndarray,
    lat: np.ndarray,
    Z: np.ndarray,
    xg2d: np.ndarray,
    yg2d: np.ndarray,
    nnear: int,
    p: int,
    max_distance: Optional[float],
    batch: int = 32,
    cell_chunk: int = 20000,
    workers: int = 1,
) -> np.ndarray:
    """
    Same numerical IDW implementation as the supplied merged gridding script.
    """
    from scipy.spatial import cKDTree

    Z = np.asarray(Z, dtype="float32")
    n_t, _ = Z.shape
    ny, nx = xg2d.shape
    ncell = ny * nx
    q = np.c_[xg2d.ravel(), yg2d.ravel()]
    out = np.full((n_t, ny, nx), np.nan, dtype="float32")

    finite = np.isfinite(Z)
    keys = np.packbits(finite, axis=1)
    groups: dict[bytes, list[int]] = {}
    for i in range(n_t):
        groups.setdefault(keys[i].tobytes(), []).append(i)

    cell_chunk = max(1, int(cell_chunk))
    batch = max(1, int(batch))

    for idx_list in groups.values():
        idxs = np.asarray(idx_list, dtype=np.int64)
        valid = finite[idxs[0]]
        nv = int(valid.sum())
        if nv == 0:
            continue

        tree = cKDTree(
            np.c_[lon[valid], lat[valid]],
            balanced_tree=False,
            compact_nodes=False,
        )
        k = min(int(nnear), nv)
        d, ix = tree.query(q, k=k, workers=int(workers))
        if k == 1:
            d = d[:, None]
            ix = ix[:, None]

        with np.errstate(divide="ignore", invalid="ignore"):
            w = 1.0 / d ** int(p)
            wn = w / w.sum(axis=1, keepdims=True)

        near0 = d[:, 0] < 1e-10
        far = (d[:, 0] > max_distance) if max_distance is not None else None
        ix0 = ix[:, 0]

        for b0 in range(0, len(idxs), batch):
            bb = idxs[b0:b0 + batch]
            zb = Z[bb][:, valid]
            vals = np.full((len(bb), ncell), np.nan, dtype="float32")

            for c0 in range(0, ncell, cell_chunk):
                c1 = min(c0 + cell_chunk, ncell)
                ixc = ix[c0:c1]
                wnc = wn[c0:c1].astype("float32", copy=False)
                vals[:, c0:c1] = np.einsum(
                    "bck,ck->bc",
                    zb[:, ixc],
                    wnc,
                    optimize=True,
                ).astype("float32", copy=False)

            if near0.any():
                vals[:, near0] = zb[:, ix0[near0]]
            if far is not None:
                vals[:, far] = np.nan

            out[bb] = vals.reshape((len(bb), ny, nx))

        del tree, d, ix, w, wn, near0, far, ix0
        gc.collect()

    return out


class GridWriter:
    def __init__(
        self,
        path: Path,
        times: pd.DatetimeIndex,
        xgrid: np.ndarray,
        ygrid: np.ndarray,
        args,
    ):
        try:
            import netCDF4
        except Exception as exc:
            raise RuntimeError("netCDF4 is required to stream the output file.") from exc

        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            if args.overwrite:
                path.unlink()
            else:
                raise FileExistsError(f"Output exists: {path}. Use --overwrite.")

        root = netCDF4.Dataset(path, "w", format="NETCDF4")
        self.root = root
        root.createDimension("time", len(times))
        root.createDimension("y", len(ygrid))
        root.createDimension("x", len(xgrid))

        tv = root.createVariable("time", "i8", ("time",))
        tv[:] = (times.view("int64") // 10**9).astype("int64")
        tv.units = "seconds since 1970-01-01 00:00:00"
        tv.calendar = "proleptic_gregorian"
        tv.long_name = "CML 5-min time"

        xv = root.createVariable("x", "f8", ("x",))
        yv = root.createVariable("y", "f8", ("y",))
        xv[:] = xgrid
        yv[:] = ygrid
        xv.units = "degrees_east"
        yv.units = "degrees_north"
        xv.long_name = "longitude"
        yv.long_name = "latitude"

        ny, nx = len(ygrid), len(xgrid)
        tchunk = max(1, min(int(args.time_chunk), len(times)))
        ychunk = min(ny, int(args.y_chunk)) if args.y_chunk else ny
        xchunk = min(nx, int(args.x_chunk)) if args.x_chunk else nx

        self.R = root.createVariable(
            "R",
            "f4",
            ("time", "y", "x"),
            zlib=not args.no_compression,
            complevel=int(args.complevel),
            shuffle=True,
            chunksizes=(tchunk, ychunk, xchunk),
            fill_value=np.float32(np.nan),
        )
        self.R.units = "mm h-1"
        self.R.long_name = "CML rain rate gridded from filtered publication R_mean_5min"

        root.title = "Filtered ESSD publication CML rainfall gridded with IDW"
        root.source_rainrate_file = str(args.input_nc)
        root.source_coordinate_file = str(args.coords_nc)
        root.rain_variable = args.rain_var
        root.baseline_method = args.baseline_method
        root.kr_model = args.kr_model
        root.waa_method = args.waa_method
        root.interpolation = f"IDW (nnear={args.nnear}, p={args.p}, max_dist={args.max_dist})"
        root.grid_resolution_deg = float(args.resolution)
        root.note = (
            "Input opened once and output streamed sequentially to avoid repeated "
            "netCDF4/HDF5 opens on shared HPC filesystems."
        )

    def write(self, start: int, data: np.ndarray) -> None:
        end = start + data.shape[0]
        self.R[start:end, :, :] = np.asarray(data, dtype="float32")
        self.root.sync()

    def close(self) -> None:
        self.root.close()


def run(args) -> None:
    t0 = time.perf_counter()

    args.input_nc = Path(args.input_nc).expanduser().resolve()
    args.coords_nc = Path(args.coords_nc).expanduser().resolve()
    args.output_nc = Path(args.output_nc).expanduser().resolve()

    for path, label in [
        (args.input_nc, "input rain-rate file"),
        (args.coords_nc, "coordinate file"),
    ]:
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")

    engine = choose_input_engine(args.input_engine)

    print("=" * 100)
    print("SAFE SINGLE-OPEN ESSD CML IDW GRIDDING")
    print("=" * 100)
    print(f"Input     : {args.input_nc}")
    print(f"Coordinates: {args.coords_nc}")
    print(f"Output    : {args.output_nc}")
    print(f"Engine    : {engine}")
    print(f"Combo     : {args.baseline_method} + {args.kr_model} + {args.waa_method}")
    print(f"IDW       : nnear={args.nnear}, p={args.p}, max_dist={args.max_dist}, resolution={args.resolution}")
    print("=" * 100)

    ds = open_input(args.input_nc, engine)
    writer = None

    try:
        if "cml_id" not in ds.coords:
            raise KeyError("Input file must contain cml_id coordinate.")

        cml_ids = np.asarray(ds["cml_id"].values).astype(str)
        rain = select_combo(
            ds,
            args.rain_var,
            args.baseline_method,
            args.kr_model,
            args.waa_method,
        )

        times_all = pd.DatetimeIndex(ds["time_5min"].values)
        file_i0, file_i1, times = resolve_time_window(
            times_all,
            args.tstart,
            args.tend,
        )

        lon_all, lat_all = load_midpoints(args.coords_nc, cml_ids, args)
        good = np.isfinite(lon_all) & np.isfinite(lat_all)
        good_idx = np.flatnonzero(good).astype(np.int64)
        if good_idx.size == 0:
            raise RuntimeError("No finite CML midpoint coordinates.")

        lon = lon_all[good]
        lat = lat_all[good]

        xg2d, yg2d = build_grid(lon, lat, args.resolution)
        xgrid = xg2d[0, :]
        ygrid = yg2d[:, 0]
        ny, nx = xg2d.shape

        print(f"Links for IDW : {good_idx.size}/{len(cml_ids)}")
        print(f"Grid          : {ny} x {nx} = {ny*nx:,} cells")
        print(f"Time          : {times[0]} -> {times[-1]}")
        print(f"Time steps    : {len(times):,}")
        print(f"Time chunk    : {args.time_chunk}")

        one_chunk_gb = min(args.time_chunk, len(times)) * ny * nx * 4 / 1024**3
        print(f"Output chunk  : ~{one_chunk_gb:.3f} GB float32")

        writer = GridWriter(
            args.output_nc,
            times,
            xgrid,
            ygrid,
            args,
        )

        starts = range(0, len(times), int(args.time_chunk))
        pbar = tqdm(starts, total=math.ceil(len(times) / int(args.time_chunk)), desc="gridding", unit="chunk")

        for out_start in pbar:
            out_end = min(out_start + int(args.time_chunk), len(times))
            file_start = file_i0 + out_start
            file_end = file_i0 + out_end

            read_t0 = time.perf_counter()

            # IMPORTANT: read from the SAME already-open dataset handle.
            Z = rain.isel(
                cml_id=good_idx,
                time_5min=slice(file_start, file_end),
            ).transpose(
                "time_5min",
                "cml_id",
            ).values

            Z = np.ascontiguousarray(Z, dtype="float32")

            G = idw_fast_series(
                lon,
                lat,
                Z,
                xg2d,
                yg2d,
                nnear=args.nnear,
                p=args.p,
                max_distance=args.max_dist,
                batch=args.batch,
                cell_chunk=args.cell_chunk,
                workers=args.query_workers,
            )

            writer.write(out_start, G)

            elapsed = time.perf_counter() - read_t0
            finite_fraction = float(np.isfinite(G).mean()) if G.size else np.nan
            pbar.set_postfix_str(
                f"{out_start}:{out_end} {elapsed:.1f}s finite={finite_fraction:.3f}"
            )

            del Z, G
            gc.collect()

        writer.close()
        writer = None

    finally:
        if writer is not None:
            writer.close()
        ds.close()

    print(f"\nSaved: {args.output_nc}")
    print(f"Elapsed: {time.perf_counter() - t0:.1f} s")


def parse_args():
    p = argparse.ArgumentParser(
        description="Safe single-open IDW gridding of filtered ESSD CML publication rainfall.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--input-nc",
        default="/home/jovyan/ESSD_publication_dataset/CZ_CML_RAINRATE_5MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--coords-nc",
        default="/home/jovyan/ESSD_publication_dataset/CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc",
    )
    p.add_argument(
        "--output-nc",
        default="/home/jovyan/ESSD_publication_dataset/CZ_CML_GRID_5MIN_PYCML_CONSTANT_ITU2003_PASTOREK2021_IDW.nc",
    )

    p.add_argument("--rain-var", default="R_mean_5min")
    p.add_argument("--baseline-method", default="pycml_constant")
    p.add_argument("--kr-model", default="ITU_2003")
    p.add_argument("--waa-method", default="pastorek_2021")

    p.add_argument("--tstart", default=None)
    p.add_argument("--tend", default=None)

    p.add_argument("--resolution", type=float, default=0.01)
    p.add_argument("--nnear", type=int, default=15)
    p.add_argument("--p", type=int, default=2)
    p.add_argument("--max-dist", type=float, default=0.3)

    p.add_argument("--time-chunk", type=int, default=288)
    p.add_argument("--batch", type=int, default=32)
    p.add_argument("--cell-chunk", type=int, default=20000)
    p.add_argument("--query-workers", type=int, default=1)

    p.add_argument("--input-engine", choices=["auto", "h5netcdf", "netcdf4"], default="auto")
    p.add_argument("--complevel", type=int, default=4)
    p.add_argument("--no-compression", action="store_true")
    p.add_argument("--y-chunk", type=int, default=None)
    p.add_argument("--x-chunk", type=int, default=None)
    p.add_argument("--overwrite", action="store_true")

    p.add_argument("--site-a-lat", default="site_a_latitude")
    p.add_argument("--site-a-lon", default="site_a_longitude")
    p.add_argument("--site-b-lat", default="site_b_latitude")
    p.add_argument("--site-b-lon", default="site_b_longitude")

    args = p.parse_args()

    if args.time_chunk < 1:
        p.error("--time-chunk must be >=1")
    if args.batch < 1:
        p.error("--batch must be >=1")
    if args.cell_chunk < 1:
        p.error("--cell-chunk must be >=1")
    if args.query_workers < 1:
        p.error("--query-workers must be >=1")
    if args.nnear < 1:
        p.error("--nnear must be >=1")
    if args.resolution <= 0:
        p.error("--resolution must be >0")

    return args


if __name__ == "__main__":
    run(parse_args())
