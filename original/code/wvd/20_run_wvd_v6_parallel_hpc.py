#!/usr/bin/env python3
"""
Parallel HPC launcher + merger for:
    19_wvd_eband_progressive_baseline_v6.py

SCIENCE IS UNCHANGED.

The V6 retrieval is independent by directional CML sublink. This launcher:
  1) reproduces the same 71-76 / 81-86 GHz, 2-10 km population ordering;
  2) splits the physical CMLs into contiguous chunks;
  3) runs the UNMODIFIED V6 script on those chunks in parallel;
  4) concatenates the resulting NetCDF along `series`;
  5) merges diagnostics/regimes/events/reinitialization CSVs and renumbers
     `series_index` to the original full-population order.

No V6 thresholds, P.676 physics, wet/dry masks, baseline rules, temporal
resolution, or retrieval equations are changed.

Designed for the Czech ESSD workflow: expected 676 physical CMLs /
1352 directional sublinks by default.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Run unmodified WVD V6 in parallel chunks and merge exactly by series."
    )
    p.add_argument("--v6-script", required=True)
    p.add_argument("--helper-script", required=True)
    p.add_argument("--cml-nc", required=True)
    p.add_argument("--gauge-nc", required=True)
    p.add_argument("--temperature-nc", required=True)

    p.add_argument("--out-nc", required=True)
    p.add_argument("--diagnostics-csv", required=True)
    p.add_argument("--states-csv", required=True)
    p.add_argument("--events-csv", required=True)
    p.add_argument("--reinit-csv", required=True)

    p.add_argument("--workers", type=int, default=40)
    p.add_argument("--expected-cmls", type=int, default=676)
    p.add_argument("--expected-series", type=int, default=1352)

    p.add_argument(
        "--tmp-dir",
        default=None,
        help="Worker-output directory. Default: hidden directory beside out-nc.",
    )
    p.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep worker NetCDF/CSV/log files after successful merge.",
    )
    p.add_argument("--overwrite", action="store_true")

    # These are ONLY population discovery settings. They match V6 defaults.
    p.add_argument("--length-min-km", type=float, default=2.0)
    p.add_argument("--length-max-km", type=float, default=10.0)
    p.add_argument("--band1-min-ghz", type=float, default=71.0)
    p.add_argument("--band1-max-ghz", type=float, default=76.0)
    p.add_argument("--band2-min-ghz", type=float, default=81.0)
    p.add_argument("--band2-max-ghz", type=float, default=86.0)

    return p.parse_args()


# ---------------------------------------------------------------------
# Population discovery: same order as V6
# ---------------------------------------------------------------------

def _decode_strings(a):
    a = np.asarray(a)

    # Classic netCDF fixed-width character arrays.
    if a.dtype.kind in ("S", "U") and a.ndim >= 2 and a.shape[-1] > 1:
        try:
            z = netCDF4.chartostring(a)
            return [str(x).strip() for x in np.asarray(z).reshape(-1)]
        except Exception:
            pass

    out = []
    for x in a.reshape(-1):
        if isinstance(x, bytes):
            out.append(x.decode("utf-8", errors="replace").strip())
        else:
            out.append(str(x).strip())
    return out


def _to_ghz(values, units=""):
    x = np.asarray(values, dtype=float)
    u = str(units or "").lower().replace(" ", "")
    if "ghz" in u:
        return x
    if "mhz" in u:
        return x / 1e3
    if "khz" in u:
        return x / 1e6
    if "hz" in u:
        return x / 1e9

    finite = np.abs(x[np.isfinite(x)])
    med = float(np.nanmedian(finite)) if finite.size else np.nan
    if np.isfinite(med) and med > 1e7:
        return x / 1e9
    if np.isfinite(med) and med > 1e3:
        return x / 1e3
    return x


def _to_km(values, units=""):
    x = np.asarray(values, dtype=float)
    u = str(units or "").lower().replace(" ", "")
    if "km" in u:
        return x
    if u in {"m", "meter", "metre", "meters", "metres"}:
        return x / 1e3

    finite = np.abs(x[np.isfinite(x)])
    med = float(np.nanmedian(finite)) if finite.size else np.nan
    if np.isfinite(med) and med > 100:
        return x / 1e3
    return x


@dataclass
class SeriesKey:
    cml_id: str
    sublink_id: str


@dataclass
class ChunkSpec:
    index: int
    cml_ids: list[str]
    series_keys: list[SeriesKey]
    series_offset: int


def discover_population(args):
    with netCDF4.Dataset(args.cml_nc, "r") as nc:
        if "cml_id" not in nc.variables:
            raise RuntimeError("Raw CML file has no `cml_id` variable")
        if "sublink_id" not in nc.variables:
            raise RuntimeError("Raw CML file has no `sublink_id` variable")
        if "frequency" not in nc.variables:
            raise RuntimeError("Raw CML file has no `frequency` variable")
        if "length" not in nc.variables:
            raise RuntimeError("Raw CML file has no `length` variable")

        cml_ids = _decode_strings(nc.variables["cml_id"][:])
        sub_ids = _decode_strings(nc.variables["sublink_id"][:])

        fv = nc.variables["frequency"]
        freq = _to_ghz(fv[:], getattr(fv, "units", ""))

        lv = nc.variables["length"]
        length = _to_km(lv[:], getattr(lv, "units", "")).reshape(-1)

    if freq.shape == (len(cml_ids), len(sub_ids)):
        freq = freq.T
    if freq.shape != (len(sub_ids), len(cml_ids)):
        raise RuntimeError(
            f"Unexpected frequency shape {freq.shape}; expected "
            f"({len(sub_ids)}, {len(cml_ids)}) after orientation."
        )

    selected_cmls = []
    series_keys = []
    keys_by_cml = {}

    for ci, cid in enumerate(cml_ids):
        L = float(length[ci])
        if not (
            np.isfinite(L)
            and args.length_min_km <= L <= args.length_max_km
        ):
            continue

        local = []
        for si, sid in enumerate(sub_ids):
            f = float(freq[si, ci])
            if not np.isfinite(f):
                continue
            f_ok = (
                args.band1_min_ghz <= f <= args.band1_max_ghz
                or args.band2_min_ghz <= f <= args.band2_max_ghz
            )
            if f_ok:
                key = SeriesKey(str(cid), str(sid))
                series_keys.append(key)
                local.append(key)

        if local:
            selected_cmls.append(str(cid))
            keys_by_cml[str(cid)] = local

    return selected_cmls, series_keys, keys_by_cml


def make_chunks(selected_cmls, keys_by_cml, n_workers):
    n_workers = max(1, min(int(n_workers), len(selected_cmls)))

    # Contiguous chunks preserve the exact serial V6 series ordering and make
    # the final NetCDF merge a simple series-axis concatenation.
    idx_chunks = np.array_split(np.arange(len(selected_cmls)), n_workers)

    chunks = []
    offset = 0
    for j, inds in enumerate(idx_chunks):
        ids = [selected_cmls[int(i)] for i in inds]
        keys = []
        for cid in ids:
            keys.extend(keys_by_cml[cid])
        chunks.append(
            ChunkSpec(
                index=j,
                cml_ids=ids,
                series_keys=keys,
                series_offset=offset,
            )
        )
        offset += len(keys)
    return chunks


# ---------------------------------------------------------------------
# Worker launch
# ---------------------------------------------------------------------

def worker_paths(tmp_dir: Path, idx: int):
    stem = tmp_dir / f"chunk_{idx:03d}"
    return {
        "nc": stem.with_suffix(".nc"),
        "diag": Path(str(stem) + "_DIAGNOSTICS.csv"),
        "states": Path(str(stem) + "_REGIMES.csv"),
        "events": Path(str(stem) + "_EVENTS.csv"),
        "reinit": Path(str(stem) + "_REINITIALIZATIONS.csv"),
        "log": Path(str(stem) + ".log"),
        "summary": stem.with_suffix(".summary.json"),
    }


def build_worker_command(args, chunk: ChunkSpec, paths):
    cmd = [
        sys.executable,
        "-u",
        str(Path(args.v6_script)),
        "--helper-script", str(Path(args.helper_script)),
        "--cml-nc", str(Path(args.cml_nc)),
        "--gauge-nc", str(Path(args.gauge_nc)),
        "--temperature-nc", str(Path(args.temperature_nc)),
        "--out-nc", str(paths["nc"]),
        "--diagnostics-csv", str(paths["diag"]),
        "--states-csv", str(paths["states"]),
        "--events-csv", str(paths["events"]),
        "--reinit-csv", str(paths["reinit"]),
        "--overwrite",
    ]
    for cid in chunk.cml_ids:
        cmd.extend(["--cml-id", cid])
    return cmd


def tail_text(path: Path, n=80):
    try:
        lines = path.read_text(errors="replace").splitlines()
        return "\n".join(lines[-n:])
    except Exception as e:
        return f"<could not read {path}: {e}>"


def run_workers(args, chunks, tmp_dir):
    env = os.environ.copy()

    # One Python process = one CPU. Prevent BLAS/OpenMP oversubscription.
    for k in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "BLIS_NUM_THREADS",
    ):
        env[k] = "1"

    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    env.setdefault("MALLOC_ARENA_MAX", "2")
    env["PYTHONUNBUFFERED"] = "1"

    running = []
    t0 = time.time()

    print("\n" + "=" * 110)
    print("PARALLEL V6 EXECUTION — SCIENCE UNCHANGED")
    print("=" * 110)
    print(f"Workers                    : {len(chunks)}")
    print(f"Physical CMLs              : {sum(len(c.cml_ids) for c in chunks)}")
    print(f"Directional sublinks       : {sum(len(c.series_keys) for c in chunks)}")
    print(f"Temporary directory        : {tmp_dir}")
    print("Threads per worker          : 1")
    print("V6 code                     : UNMODIFIED")
    print()

    manifest = []

    try:
        for chunk in chunks:
            paths = worker_paths(tmp_dir, chunk.index)
            cmd = build_worker_command(args, chunk, paths)
            log_fh = open(paths["log"], "w", buffering=1)
            p = subprocess.Popen(
                cmd,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                env=env,
                start_new_session=True,
            )
            running.append(
                {
                    "chunk": chunk,
                    "paths": paths,
                    "cmd": cmd,
                    "proc": p,
                    "log_fh": log_fh,
                    "reported": False,
                }
            )
            manifest.append(
                {
                    "chunk": chunk.index,
                    "n_cml": len(chunk.cml_ids),
                    "n_series": len(chunk.series_keys),
                    "series_offset": chunk.series_offset,
                    "cml_ids": chunk.cml_ids,
                    "command": cmd,
                }
            )

        (tmp_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2)
        )

        last_status = 0.0
        n_done = 0

        while n_done < len(running):
            time.sleep(2.0)

            for item in running:
                if item["reported"]:
                    continue
                rc = item["proc"].poll()
                if rc is None:
                    continue

                item["log_fh"].flush()
                item["log_fh"].close()
                item["reported"] = True
                n_done += 1

                if rc != 0:
                    # Terminate all other workers immediately.
                    for other in running:
                        if other["proc"].poll() is None:
                            try:
                                os.killpg(
                                    os.getpgid(other["proc"].pid),
                                    signal.SIGTERM,
                                )
                            except Exception:
                                try:
                                    other["proc"].terminate()
                                except Exception:
                                    pass

                    failed_log = item["paths"]["log"]
                    raise RuntimeError(
                        f"Worker chunk {item['chunk'].index:03d} failed "
                        f"with exit code {rc}.\n\n"
                        f"--- tail of {failed_log} ---\n"
                        f"{tail_text(failed_log)}"
                    )

                elapsed = time.time() - t0
                print(
                    f"[worker {item['chunk'].index+1:02d}/{len(chunks):02d}] "
                    f"DONE  CMLs={len(item['chunk'].cml_ids):2d} "
                    f"series={len(item['chunk'].series_keys):2d} "
                    f"elapsed={elapsed/60:.1f} min"
                )

            now = time.time()
            if now - last_status >= 60:
                alive = sum(x["proc"].poll() is None for x in running)
                done = len(running) - alive
                print(
                    f"[status] {done}/{len(running)} workers finished, "
                    f"{alive} running, elapsed={(now-t0)/60:.1f} min"
                )
                last_status = now

        return time.time() - t0

    finally:
        for item in running:
            try:
                if not item["log_fh"].closed:
                    item["log_fh"].close()
            except Exception:
                pass


# ---------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------

def verify_worker_order(chunks, tmp_dir):
    print("\nVerifying worker population/order against serial V6 ordering...")

    total = 0
    for chunk in chunks:
        paths = worker_paths(tmp_dir, chunk.index)
        df = pd.read_csv(paths["diag"])
        if "series_index" not in df:
            raise RuntimeError(f"{paths['diag']} has no series_index")

        df = df.sort_values("series_index", kind="stable")
        got = [
            (str(a), str(b))
            for a, b in zip(df["cml_id"].astype(str), df["sublink_id"].astype(str))
        ]
        exp = [(k.cml_id, k.sublink_id) for k in chunk.series_keys]

        if got != exp:
            raise RuntimeError(
                f"Worker {chunk.index:03d} series ordering differs from "
                "the serial V6 population ordering."
            )

        with netCDF4.Dataset(paths["nc"], "r") as nc:
            n = len(nc.dimensions["series"])
        if n != len(exp):
            raise RuntimeError(
                f"Worker {chunk.index:03d}: NetCDF series={n}, expected={len(exp)}"
            )

        total += n

    print(f"Verified directional sublinks: {total}")


# ---------------------------------------------------------------------
# Exact-by-series NetCDF merge
# ---------------------------------------------------------------------

def _creation_kwargs(src_var):
    kwargs = {}

    if "_FillValue" in src_var.ncattrs():
        kwargs["fill_value"] = src_var.getncattr("_FillValue")

    try:
        filters = src_var.filters() or {}
    except Exception:
        filters = {}

    # Preserve ordinary zlib compression used by V6.
    if filters.get("zlib", False):
        kwargs["zlib"] = True
        kwargs["complevel"] = int(filters.get("complevel", 4))
        kwargs["shuffle"] = bool(filters.get("shuffle", True))

    if filters.get("fletcher32", False):
        kwargs["fletcher32"] = True

    try:
        chunks = src_var.chunking()
        if isinstance(chunks, (list, tuple)):
            kwargs["chunksizes"] = tuple(int(x) for x in chunks)
    except Exception:
        pass

    try:
        endian = src_var.endian()
        if endian in ("little", "big"):
            kwargs["endian"] = endian
    except Exception:
        pass

    return kwargs


def merge_netcdf(chunks, tmp_dir, out_nc: Path):
    print("\nMerging worker NetCDF files along the unchanged `series` dimension...")

    src_paths = [worker_paths(tmp_dir, c.index)["nc"] for c in chunks]
    sources = [netCDF4.Dataset(p, "r") for p in src_paths]

    try:
        for s in sources:
            s.set_auto_mask(False)
            s.set_auto_scale(False)

        first = sources[0]
        total_series = sum(len(s.dimensions["series"]) for s in sources)

        tmp_out = Path(str(out_nc) + ".merging")
        if tmp_out.exists():
            tmp_out.unlink()

        dst = netCDF4.Dataset(tmp_out, "w", format=first.data_model)
        dst.set_auto_mask(False)
        dst.set_auto_scale(False)

        try:
            # Dimensions
            for name, dim in first.dimensions.items():
                if name == "series":
                    dst.createDimension(name, total_series)
                else:
                    dst.createDimension(
                        name,
                        None if dim.isunlimited() else len(dim)
                    )

            # Global attributes
            for a in first.ncattrs():
                dst.setncattr(a, first.getncattr(a))
            dst.setncattr(
                "parallel_execution",
                f"Unmodified V6 executed as {len(chunks)} independent chunks; "
                "merged only along the independent series dimension."
            )
            dst.setncattr(
                "parallel_science_changed",
                "false"
            )

            # Variables
            for name, sv0 in first.variables.items():
                kwargs = _creation_kwargs(sv0)

                # A source chunk size can never exceed its local dimension,
                # but is also valid in the larger merged series dimension.
                dv = dst.createVariable(
                    name,
                    sv0.datatype,
                    sv0.dimensions,
                    **kwargs,
                )
                dv.set_auto_mask(False)
                dv.set_auto_scale(False)

                for a in sv0.ncattrs():
                    if a == "_FillValue":
                        continue
                    dv.setncattr(a, sv0.getncattr(a))

                if "series" not in sv0.dimensions:
                    dv[...] = sv0[...]
                    continue

                axis = sv0.dimensions.index("series")
                offset = 0
                for src in sources:
                    sv = src.variables[name]
                    n_local = len(src.dimensions["series"])

                    dslice = [slice(None)] * dv.ndim
                    dslice[axis] = slice(offset, offset + n_local)
                    dv[tuple(dslice)] = sv[...]

                    offset += n_local

                if offset != total_series:
                    raise RuntimeError(
                        f"Variable {name}: merged {offset} series, "
                        f"expected {total_series}"
                    )

            dst.sync()

        finally:
            dst.close()

        os.replace(tmp_out, out_nc)

    finally:
        for s in sources:
            s.close()

    print(f"Merged NetCDF: {out_nc}")


# ---------------------------------------------------------------------
# CSV / summary merge
# ---------------------------------------------------------------------

def merge_csv(chunks, tmp_dir, key, out_path: Path):
    frames = []
    for chunk in chunks:
        p = worker_paths(tmp_dir, chunk.index)[key]
        try:
            df = pd.read_csv(p)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame()

        if len(df) and "series_index" in df.columns:
            df["series_index"] = (
                pd.to_numeric(df["series_index"], errors="raise").astype(int)
                + int(chunk.series_offset)
            )
        frames.append(df)

    nonempty_schema = [x for x in frames if len(x.columns)]
    if nonempty_schema:
        columns = list(nonempty_schema[0].columns)
        fixed = []
        for x in frames:
            if len(x.columns) == 0:
                fixed.append(pd.DataFrame(columns=columns))
            else:
                fixed.append(x)
        merged = pd.concat(fixed, ignore_index=True)
    else:
        merged = pd.DataFrame()

    merged.to_csv(out_path, index=False)
    return merged


def merge_summary(args, chunks, tmp_dir, out_nc, frames):
    summaries = []
    for c in chunks:
        p = worker_paths(tmp_dir, c.index)["summary"]
        if p.exists():
            summaries.append(json.loads(p.read_text()))

    base = dict(summaries[0]) if summaries else {}
    base["n_series"] = int(sum(len(c.series_keys) for c in chunks))
    base["n_physical_cml"] = int(sum(len(c.cml_ids) for c in chunks))
    base["n_hardware_events"] = int(len(frames["events"]))
    base["n_baseline_reinitializations"] = int(len(frames["reinit"]))
    base["n_state_rows_common"] = int(len(frames["states"]))
    base["parallel_workers"] = int(len(chunks))
    base["parallel_execution"] = True
    base["science_changed_by_parallelization"] = False

    out_summary = out_nc.with_suffix(".summary.json")
    out_summary.write_text(json.dumps(base, indent=2))
    return out_summary


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    args = parse_args()

    out_nc = Path(args.out_nc)
    diag_csv = Path(args.diagnostics_csv)
    states_csv = Path(args.states_csv)
    events_csv = Path(args.events_csv)
    reinit_csv = Path(args.reinit_csv)

    for p in (
        Path(args.v6_script),
        Path(args.helper_script),
        Path(args.cml_nc),
        Path(args.gauge_nc),
        Path(args.temperature_nc),
    ):
        if not p.exists():
            raise FileNotFoundError(p)

    finals = [
        out_nc,
        diag_csv,
        states_csv,
        events_csv,
        reinit_csv,
        out_nc.with_suffix(".summary.json"),
    ]

    if any(p.exists() for p in finals) and not args.overwrite:
        existing = [str(p) for p in finals if p.exists()]
        raise FileExistsError(
            "Final output(s) already exist. Use --overwrite:\n" + "\n".join(existing)
        )

    for p in finals:
        p.parent.mkdir(parents=True, exist_ok=True)

    tmp_dir = (
        Path(args.tmp_dir)
        if args.tmp_dir
        else out_nc.parent / f".{out_nc.stem}_parallel_tmp"
    )

    if tmp_dir.exists():
        if not args.overwrite:
            raise FileExistsError(
                f"Temporary directory exists: {tmp_dir}; use --overwrite"
            )
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 110)
    print("DISCOVERING THE EXACT V6 E-BAND POPULATION")
    print("=" * 110)

    selected_cmls, series_keys, keys_by_cml = discover_population(args)

    print(f"Physical CMLs         : {len(selected_cmls):,}")
    print(f"Directional sublinks  : {len(series_keys):,}")

    if args.expected_cmls is not None and len(selected_cmls) != args.expected_cmls:
        raise RuntimeError(
            f"Expected {args.expected_cmls} physical CMLs but discovered "
            f"{len(selected_cmls)}. Aborting rather than risk a population mismatch."
        )

    if args.expected_series is not None and len(series_keys) != args.expected_series:
        raise RuntimeError(
            f"Expected {args.expected_series} directional sublinks but discovered "
            f"{len(series_keys)}. Aborting rather than risk a population mismatch."
        )

    chunks = make_chunks(selected_cmls, keys_by_cml, args.workers)

    print(f"Parallel chunks       : {len(chunks)}")
    print(
        "Series/chunk          : "
        f"{min(len(c.series_keys) for c in chunks)}–"
        f"{max(len(c.series_keys) for c in chunks)}"
    )

    t_all = time.time()

    worker_seconds = run_workers(args, chunks, tmp_dir)

    verify_worker_order(chunks, tmp_dir)

    # Remove/replace old final outputs only after all workers succeeded.
    if args.overwrite:
        for p in finals:
            if p.exists():
                p.unlink()

    t_merge = time.time()

    merge_netcdf(chunks, tmp_dir, out_nc)

    frames = {
        "diag": merge_csv(chunks, tmp_dir, "diag", diag_csv),
        "states": merge_csv(chunks, tmp_dir, "states", states_csv),
        "events": merge_csv(chunks, tmp_dir, "events", events_csv),
        "reinit": merge_csv(chunks, tmp_dir, "reinit", reinit_csv),
    }

    summary_path = merge_summary(
        args, chunks, tmp_dir, out_nc, frames
    )

    # Final hard checks.
    with netCDF4.Dataset(out_nc, "r") as nc:
        merged_series = len(nc.dimensions["series"])

    if merged_series != len(series_keys):
        raise RuntimeError(
            f"Final NetCDF has {merged_series} series; expected {len(series_keys)}"
        )

    if len(frames["diag"]) != len(series_keys):
        raise RuntimeError(
            f"Final diagnostics has {len(frames['diag'])} rows; "
            f"expected {len(series_keys)}"
        )

    total_seconds = time.time() - t_all
    merge_seconds = time.time() - t_merge

    print("\n" + "=" * 110)
    print("DONE — PARALLEL V6")
    print("=" * 110)
    print(f"NetCDF          : {out_nc}")
    print(f"Diagnostics     : {diag_csv}")
    print(f"Regimes         : {states_csv}")
    print(f"Events          : {events_csv}")
    print(f"Reinitializations: {reinit_csv}")
    print(f"Summary         : {summary_path}")
    print(f"Workers         : {len(chunks)}")
    print(f"Series          : {merged_series}")
    print(f"Worker walltime : {worker_seconds/60:.1f} min")
    print(f"Merge walltime  : {merge_seconds/60:.1f} min")
    print(f"Total walltime  : {total_seconds/60:.1f} min")
    print("SCIENCE CHANGED : NO")
    print("=" * 110)

    if not args.keep_temp:
        shutil.rmtree(tmp_dir)
        print(f"Removed temporary worker directory: {tmp_dir}")
    else:
        print(f"Kept temporary worker directory: {tmp_dir}")


if __name__ == "__main__":
    main()
