#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys
import pandas as pd

ROOT = Path("/home/jovyan/ESSD_publication_dataset")
APR = ROOT / "cmlonly_graf_nearby_comparison" / "apr_oct_2021"

PYTHON = Path("/home/jovyan/my-conda-envs/tf-env/bin/python")

S3 = ROOT / "03_blettner_event_maps_timeseries.py"

CACHE = (
    APR / "article_style" /
    "multiple_events_radar" /
    "selected_instances_cache.csv"
)

GEOM = ROOT / "CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc"

RADAR_ROOT = Path("/home/jovyan/dwd_radolan/2021/suplement")

GRIDS = {
    "GRAF": (
        APR / "grids" /
        "CZ_CML_GRID_GRAF_5MIN_APR_OCT_2021_IDW.nc"
    ),
    "NEARBY": (
        APR / "grids" /
        "CZ_CML_RAINRATE_NEARBY_5MIN_APR_OCT_2021_IDW.nc"
    ),
}

BASE_OUT = (
    APR / "article_style" /
    "multiple_events_radar"
)

if not CACHE.is_file():
    raise FileNotFoundError(
        f"Missing event cache: {CACHE}"
    )

for label, grid in GRIDS.items():
    if not grid.is_file():
        raise FileNotFoundError(
            f"Missing {label} grid: {grid}"
        )

events = pd.read_csv(
    CACHE,
    parse_dates=["time", "start_time"],
)

print("=" * 110)
print("MULTIPLE EVENT MAPS")
print("Every comparison is CML METHOD versus RADAR")
print("=" * 110)

for _, row in events.iterrows():

    rank = int(row["event_rank"])
    peak = pd.Timestamp(row["time"])
    start = pd.Timestamp(row["start_time"])

    tag = peak.strftime("%Y%m%d_%H%M")

    print()
    print("=" * 110)
    print(f"EVENT {rank}: {peak}")
    print("=" * 110)

    for method in ("GRAF", "NEARBY"):

        out = (
            BASE_OUT /
            f"event_{rank:02d}_{tag}" /
            f"{method}_vs_RADAR"
        )

        out.mkdir(
            parents=True,
            exist_ok=True,
        )

        cmd = [
            str(PYTHON),
            "-u",
            str(S3),

            "--gridded-nc",
            str(GRIDS[method]),

            "--cml-var",
            "R",

            "--cml-row-label",
            method,

            "--cml-geom-nc",
            str(GEOM),

            "--radar-root",
            str(RADAR_ROOT),

            "--selection-cache-csv",
            str(CACHE),

            "--event-rank",
            str(rank),

            "--instance-rank",
            "1",

            "--start-time",
            start.strftime("%Y-%m-%d %H:%M:%S"),

            "--ncols",
            "8",

            "--panel-step-minutes",
            "60",

            "--dpi",
            "600",

            "--out-dir",
            str(out),
        ]

        # Time series are produced separately from the current 2219-link
        # products. Do not invoke the obsolete Blettner V2 time-series routine.
        cmd += [
            "--skip-timeseries",
        ]

        print()
        print(f"{method} vs RADAR")
        print("Output:", out)

        logfile = out / "run.log"

        with logfile.open("w") as log:

            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            for line in process.stdout:
                sys.stdout.write(line)
                log.write(line)

            code = process.wait()

        if code != 0:
            raise RuntimeError(
                f"{method} vs RADAR failed "
                f"for event {rank}. "
                f"See {logfile}"
            )

print()
print("=" * 110)
print("ALL EVENTS FINISHED")
print("=" * 110)
