#!/usr/bin/env python3
"""
Parallel launcher for 30c_v6_quantization_chunk.py.

Science is identical to the serial quantization characterization. It only
partitions the exact V6 diagnostics population across workers, then merges the
per-worker CSV/NetCDF summaries in original series_index order.

Each worker opens the raw NetCDF once and processes a contiguous block, which
avoids repeatedly reopening the large raw file for every individual link.
"""

from __future__ import annotations
import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    root = "/home/jovyan/ESSD_publication_dataset"
    p.add_argument("--chunk-script", default=f"{root}/30c_v6_quantization_chunk.py")
    p.add_argument("--v6-diagnostics-csv", default=f"{root}/CZ_CML_WVD_EBAND_PROGRESSIVE_V6_DIAGNOSTICS.csv")
    p.add_argument("--raw-cml-nc", default=f"{root}/CZ_CML_RAW_1MIN_20200101_20211231.nc")
    p.add_argument("--v6-script", default=f"{root}/19_wvd_eband_progressive_baseline_v6.py")
    p.add_argument("--physics-helper-script", default=f"{root}/13_wvd_raw_60_100_timevarying_T.py")
    p.add_argument("--out-nc", default=f"{root}/CZ_CML_WVD_V6_LINK_QUANTIZATION.nc")
    p.add_argument("--series-csv", default=f"{root}/CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES.csv")
    p.add_argument("--states-csv", default=f"{root}/CZ_CML_WVD_V6_LINK_QUANTIZATION_STATES.csv")
    p.add_argument("--workers", type=int, default=40)
    p.add_argument("--length-min-km", type=float, default=None)
    p.add_argument("--length-max-km", type=float, default=None)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def run_worker(k, start, stop, args, tmp):
    prefix = tmp / f"worker_{k:03d}"
    out_nc = prefix.with_suffix(".nc")
    series_csv = Path(str(prefix) + "_series.csv")
    states_csv = Path(str(prefix) + "_states.csv")
    log = Path(str(prefix) + ".log")

    cmd = [
        sys.executable, "-u", str(Path(args.chunk_script).resolve()),
        "--v6-diagnostics-csv", str(Path(args.v6_diagnostics_csv).resolve()),
        "--raw-cml-nc", str(Path(args.raw_cml_nc).resolve()),
        "--v6-script", str(Path(args.v6_script).resolve()),
        "--physics-helper-script", str(Path(args.physics_helper_script).resolve()),
        "--out-nc", str(out_nc),
        "--series-csv", str(series_csv),
        "--states-csv", str(states_csv),
        "--series-start", str(start),
        "--series-stop", str(stop),
        "--overwrite",
    ]
    if args.length_min_km is not None:
        cmd += ["--length-min-km", str(args.length_min_km)]
    if args.length_max_km is not None:
        cmd += ["--length-max-km", str(args.length_max_km)]

    t0 = time.time()
    with log.open("w", encoding="utf-8") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        tail = log.read_text(encoding="utf-8", errors="ignore").splitlines()[-30:]
        raise RuntimeError(
            f"worker {k} failed ({start}:{stop})\n" + "\n".join(tail)
        )
    return {
        "worker": k, "start": start, "stop": stop,
        "out_nc": out_nc, "series_csv": series_csv,
        "states_csv": states_csv, "log": log,
        "elapsed": time.time() - t0,
    }


def write_merged_netcdf(path, sdf, args):
    path = Path(path).resolve()
    if path.exists():
        if not args.overwrite:
            raise FileExistsError(f"{path} exists; use --overwrite")
        path.unlink()

    with netCDF4.Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("series", len(sdf))
        nc.createVariable("series_index", "i4", ("series",))[:] = sdf["series_index"].to_numpy(np.int32)
        nc.createVariable("cml_id", str, ("series",))[:] = sdf["cml_id"].astype(str).to_numpy(object)
        nc.createVariable("sublink_id", str, ("series",))[:] = sdf["sublink_id"].astype(str).to_numpy(object)

        for name, units in (
            ("frequency_ghz", "GHz"),
            ("length_km", "km"),
            ("tx_quantization_db", "dB"),
            ("rx_quantization_db", "dB"),
            ("tl_quantization_db", "dB"),
            ("specific_attenuation_quantization_db_per_km", "dB km-1"),
            ("wvd_quantization_resolution_ref20c_rho10_gm3", "g m-3"),
            ("wvd_quantization_resolution_standard_min_gm3", "g m-3"),
            ("wvd_quantization_resolution_standard_median_gm3", "g m-3"),
            ("wvd_quantization_resolution_standard_max_gm3", "g m-3"),
        ):
            v = nc.createVariable(name, "f4", ("series",), fill_value=np.nan)
            v[:] = pd.to_numeric(sdf[name], errors="coerce").to_numpy(np.float32)
            v.units = units

        nc.title = "Parallel link-derived quantization characterization for unchanged V6 WVD retrieval"
        nc.population_source = str(Path(args.v6_diagnostics_csv).resolve())
        nc.source_raw_cml = str(Path(args.raw_cml_nc).resolve())
        nc.v6_retrieval_modified = "no"
        nc.humidity_used = "no"
        nc.quantization_is_not_radio_noise = "true"
        nc.parallelization_changed_science = "no"


def main():
    args = parse_args()
    t0 = time.time()

    for p in (args.chunk_script, args.v6_diagnostics_csv, args.raw_cml_nc,
              args.v6_script, args.physics_helper_script):
        if not Path(p).expanduser().is_file():
            raise FileNotFoundError(Path(p).expanduser())

    # The chunk script applies optional length filters BEFORE series slicing.
    diag = pd.read_csv(args.v6_diagnostics_csv).sort_values("series_index").reset_index(drop=True)
    if args.length_min_km is not None:
        diag = diag[pd.to_numeric(diag["length_km"], errors="coerce") >= args.length_min_km]
    if args.length_max_km is not None:
        diag = diag[pd.to_numeric(diag["length_km"], errors="coerce") <= args.length_max_km]
    diag = diag.reset_index(drop=True)

    n = len(diag)
    if n == 0:
        raise RuntimeError("No V6 directional links selected")

    workers = max(1, min(int(args.workers), n))
    chunk = int(math.ceil(n / workers))
    slices = []
    for k in range(workers):
        a = k * chunk
        b = min(n, (k + 1) * chunk)
        if a < b:
            slices.append((k + 1, a, b))

    out_nc = Path(args.out_nc).expanduser().resolve()
    tmp = out_nc.parent / f".{out_nc.stem}_quant_parallel_tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)

    print("=" * 108)
    print("PARALLEL V6 LINK-QUANTIZATION CHARACTERIZATION")
    print("=" * 108)
    print(f"Directional links : {n:,}")
    print(f"Workers           : {len(slices)}")
    print("V6 retrieval      : UNCHANGED")
    print("Humidity          : NEVER READ")
    print(f"Temporary dir     : {tmp}")
    print("=" * 108)

    done = []
    last_status = 0.0
    with ThreadPoolExecutor(max_workers=len(slices)) as ex:
        futs = {
            ex.submit(run_worker, k, a, b, args, tmp): (k, a, b)
            for k, a, b in slices
        }
        while futs:
            completed = []
            for fut in list(futs):
                if fut.done():
                    completed.append(fut)
            for fut in completed:
                info = fut.result()
                done.append(info)
                futs.pop(fut)
                print(
                    f"[worker {info['worker']:02d}/{len(slices):02d}] DONE "
                    f"rows={info['stop']-info['start']} "
                    f"elapsed={info['elapsed']/60:.1f} min"
                )
            now = time.time()
            if futs and now - last_status >= 60:
                print(
                    f"[status] {len(done)}/{len(slices)} finished, "
                    f"{len(futs)} running, elapsed={(now-t0)/60:.1f} min"
                )
                last_status = now
            if futs:
                time.sleep(2)

    done = sorted(done, key=lambda x: x["start"])

    sdf = pd.concat(
        [pd.read_csv(x["series_csv"]) for x in done],
        ignore_index=True
    ).sort_values("series_index").reset_index(drop=True)

    stf = pd.concat(
        [pd.read_csv(x["states_csv"]) for x in done],
        ignore_index=True
    ).sort_values(["series_index", "state_id"]).reset_index(drop=True)

    # Exact population verification.
    expected = diag["series_index"].astype(int).to_numpy()
    got = sdf["series_index"].astype(int).to_numpy()
    if len(got) != len(expected) or not np.array_equal(got, expected):
        raise RuntimeError(
            "Merged worker series_index does not match selected V6 population"
        )

    series_out = Path(args.series_csv).expanduser().resolve()
    states_out = Path(args.states_csv).expanduser().resolve()
    series_out.parent.mkdir(parents=True, exist_ok=True)
    states_out.parent.mkdir(parents=True, exist_ok=True)
    sdf.to_csv(series_out, index=False)
    stf.to_csv(states_out, index=False)

    write_merged_netcdf(out_nc, sdf, args)

    summary = {
        "n_series": int(len(sdf)),
        "n_states": int(len(stf)),
        "n_series_with_tl_quantization": int(np.isfinite(pd.to_numeric(sdf["tl_quantization_db"], errors="coerce")).sum()),
        "fraction_with_tl_quantization": float(np.isfinite(pd.to_numeric(sdf["tl_quantization_db"], errors="coerce")).mean()),
        "workers": len(slices),
        "retrieval_modified": False,
        "humidity_used": False,
        "parallelization_changed_science": False,
        "elapsed_seconds": time.time() - t0,
    }
    out_nc.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    shutil.rmtree(tmp)

    print("\n" + "=" * 108)
    print("DONE — PARALLEL V6 QUANTIZATION")
    print("=" * 108)
    print("NetCDF     :", out_nc)
    print("Series CSV :", series_out)
    print("States CSV :", states_out)
    print(
        "Detected qTL:",
        f"{summary['n_series_with_tl_quantization']}/{summary['n_series']}"
    )
    print(f"Elapsed    : {summary['elapsed_seconds']/60:.1f} min")
    print("SCIENCE CHANGED BY PARALLELIZATION: NO")


if __name__ == "__main__":
    main()
