#!/usr/bin/env python3
"""
CHMI historical CSV -> CF/ACDD NetCDF downloader/converter
==================================================

Designed for:
https://opendata.chmi.cz/meteorology/climate/historical_csv/data/

Main use case:
    python chmi_historical_csv_to_netcdf_v5.py humidity \
        --resolution 10min \
        --years 2024 2025 \
        --root-dir opendata.chmi.cz \
        --stage all \
        --workers 16

The script:
1. discovers/downloads CHMI CSV files with one pipeline progress bar,
2. downloads/reads meta1..meta4,
3. scans local CSV files without concatenating all data into RAM,
4. writes values, FLAG and QUALITY directly into a compressed NetCDF,
5. attaches representative station metadata,
6. preserves the complete date-valid station metadata history (meta1),
7. preserves the complete date-valid element/sensor metadata history (meta2).

Important metadata design:
- STATION in data CSV is joined to WSI in meta1/meta2.
- ELEMENT in data CSV is joined to EG_EL_ABBREVIATION in meta2.
- meta1/meta2 BEGIN_DATE and END_DATE are respected.
- CHMI uses END_DATE values such as year 3999; these are NOT parsed with
  pandas Timestamp because pandas cannot safely represent year 3999.
- A representative static lat/lon/elevation is included for convenience,
  selected as the meta1 record with the greatest temporal overlap with
  the requested years.
- The complete metadata history is also stored in the NetCDF so station
  relocations are not lost.

Dependencies:
    pip install requests beautifulsoup4 tqdm pandas numpy netCDF4

Example humidity:
    python chmi_historical_csv_to_netcdf_v5.py humidity \
        --resolution 10min \
        --years 2024 2025 \
        --elements H \
        --root-dir opendata.chmi.cz \
        --stage all \
        --workers 16

List available categories:
    python chmi_historical_csv_to_netcdf_v5.py \
        --resolution 10min \
        --list-categories
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

try:
    from netCDF4 import Dataset
except ImportError as exc:
    raise SystemExit(
        "ERROR: netCDF4 is required.\n"
        "Install dependencies with:\n"
        "  pip install requests beautifulsoup4 tqdm pandas numpy netCDF4"
    ) from exc


BASE = "https://opendata.chmi.cz/meteorology/climate/historical_csv/"
DATA_BASE = urljoin(BASE, "data/")
METADATA_BASE = urljoin(BASE, "metadata/")

RESOLUTION_INFO = {
    "10min": {
        "obs_type": "10M",
        "step_seconds": 10 * 60,
        "label": "10-minute",
    },
    "1hour": {
        "obs_type": "1H",
        "step_seconds": 60 * 60,
        "label": "1-hour",
    },
}

METADATA_FILES = ("meta1.csv", "meta2.csv", "meta3.csv", "meta4.csv")

_thread_local = threading.local()


def make_progress_bar(total: int):
    """Create the single tqdm progress bar used by the pipeline."""
    return tqdm(
        total=total,
        desc="Progress",
        unit="step",
        dynamic_ncols=True,
        leave=True,
    )


# ---------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------

def make_session() -> requests.Session:
    """Create a retry-enabled requests session."""
    session = requests.Session()
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        status=5,
        backoff_factor=0.8,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET", "HEAD"]),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=32, pool_maxsize=32)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update(
        {
            "User-Agent": (
                "Mozilla/5.0 CHMI-historical-csv-netcdf/2.0 "
                "(research data downloader)"
            )
        }
    )
    return session


def get_thread_session() -> requests.Session:
    """One Session per worker thread."""
    if not hasattr(_thread_local, "session"):
        _thread_local.session = make_session()
    return _thread_local.session


def get_html(url: str, timeout: int) -> str:
    session = make_session()
    try:
        r = session.get(url, timeout=timeout)
        r.raise_for_status()
        return r.text
    finally:
        session.close()


def list_subdirectories(url: str, timeout: int) -> list[str]:
    """Return directory names from an Apache-style directory listing."""
    html = get_html(url, timeout)
    soup = BeautifulSoup(html, "html.parser")

    out = []
    for a in soup.find_all("a"):
        href = a.get("href")
        if not href or href in ("../", "./"):
            continue
        if href.endswith("/"):
            out.append(href.rstrip("/"))

    return sorted(set(out))


def list_csv_urls(url: str, timeout: int) -> list[tuple[str, str]]:
    """Return [(filename, absolute_url), ...] from a directory listing."""
    html = get_html(url, timeout)
    soup = BeautifulSoup(html, "html.parser")

    files: list[tuple[str, str]] = []
    for a in soup.find_all("a"):
        href = a.get("href")
        if href and href.lower().endswith(".csv"):
            filename = os.path.basename(href)
            files.append((filename, urljoin(url, href)))

    # Stable unique ordering.
    return sorted(dict(files).items())


def atomic_download(
    url: str,
    target: Path,
    timeout: int,
    overwrite: bool,
) -> str:
    """
    Download one file atomically.

    Returns:
        "downloaded" or "skipped"
    """
    target.parent.mkdir(parents=True, exist_ok=True)

    if target.exists() and target.stat().st_size > 0 and not overwrite:
        return "skipped"

    part = target.with_name(target.name + ".part")
    if part.exists():
        part.unlink()

    session = get_thread_session()

    with session.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with part.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

    if part.stat().st_size == 0:
        part.unlink(missing_ok=True)
        raise IOError(f"Downloaded zero-byte file: {url}")

    os.replace(part, target)
    return "downloaded"


# ---------------------------------------------------------------------
# Metadata readers
# ---------------------------------------------------------------------

def read_standard_metadata_csv(path: Path) -> pd.DataFrame:
    """
    Read CHMI metadata with the trailing empty CSV column removed.
    Suitable for meta1, meta2, meta3.
    """
    df = pd.read_csv(path, dtype=str)
    df = df.dropna(axis=1, how="all")
    return df


def read_meta4(path: Path) -> pd.DataFrame:
    """
    Read meta4 robustly.

    CHMI meta4 descriptions can themselves contain commas and some historical
    copies are not quoted in a way pandas' normal CSV parser accepts.
    The first comma separates FLAG2 from DESCRIPTION; a final trailing comma
    is ignored.
    """
    rows: list[tuple[str, str]] = []

    with path.open("r", encoding="utf-8-sig", errors="replace") as f:
        first_line = True
        for raw in f:
            line = raw.rstrip("\r\n")
            if first_line:
                first_line = False
                continue
            if not line:
                continue
            if "," not in line:
                continue

            code, description = line.split(",", 1)
            if description.endswith(","):
                description = description[:-1]

            rows.append((code.strip(), description.strip()))

    return pd.DataFrame(rows, columns=["FLAG2", "DESCRIPTION"])


def parse_chmi_datetime(value: object) -> Optional[datetime]:
    """
    Parse CHMI metadata timestamps with Python datetime.

    Python datetime supports year 3999, unlike pandas Timestamp.
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None

    s = str(value).strip()
    if not s or s.lower() == "nan":
        return None

    if s.endswith("Z"):
        s = s[:-1] + "+00:00"

    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)

    return dt.astimezone(timezone.utc)


def periods_for_years(years: list[int]) -> list[tuple[datetime, datetime]]:
    """UTC closed intervals for selected years."""
    periods = []
    for year in sorted(set(years)):
        start = datetime(year, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
        end = datetime(year, 12, 31, 23, 59, 59, tzinfo=timezone.utc)
        periods.append((start, end))
    return periods


def overlaps_any_period(
    begin: Optional[datetime],
    end: Optional[datetime],
    periods: list[tuple[datetime, datetime]],
) -> bool:
    if begin is None or end is None:
        return False
    return any(begin <= p_end and end >= p_start for p_start, p_end in periods)


def total_overlap_seconds(
    begin: Optional[datetime],
    end: Optional[datetime],
    periods: list[tuple[datetime, datetime]],
) -> float:
    if begin is None or end is None:
        return -1.0

    total = 0.0
    for p_start, p_end in periods:
        left = max(begin, p_start)
        right = min(end, p_end)
        if right >= left:
            total += (right - left).total_seconds() + 1.0
    return total


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Download CHMI historical CSV data and build a metadata-aware "
            "compressed NetCDF without concatenating all CSVs into RAM."
        ),
    )

    parser.add_argument(
        "category",
        nargs="?",
        help=(
            "CHMI category, e.g. humidity, precipitation, temperature, "
            "wind, air_pressure."
        ),
    )
    parser.add_argument(
        "--resolution",
        choices=sorted(RESOLUTION_INFO),
        default="10min",
        help="Temporal-resolution directory.",
    )
    parser.add_argument(
        "--years",
        nargs="+",
        type=int,
        help="Years to process, e.g. --years 2024 2025",
    )
    parser.add_argument(
        "--elements",
        nargs="+",
        default=None,
        help=(
            "Optional CHMI ELEMENT codes to retain, e.g. --elements H. "
            "If omitted, all elements found in the selected category are used."
        ),
    )
    parser.add_argument(
        "--root-dir",
        default="opendata.chmi.cz",
        help="Local root directory.",
    )
    parser.add_argument(
        "--stage",
        choices=("download", "build", "all"),
        default="all",
        help="Run only download, only NetCDF build, or both.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output NetCDF path. A descriptive default is generated if omitted.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=12,
        help="Parallel download workers.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=90,
        help="HTTP timeout in seconds.",
    )
    parser.add_argument(
        "--overwrite-downloads",
        action="store_true",
        help="Redownload existing non-empty CSV files.",
    )
    parser.add_argument(
        "--refresh-metadata",
        action="store_true",
        help="Redownload meta1..meta4 even if local copies exist.",
    )
    parser.add_argument(
        "--overwrite-nc",
        action="store_true",
        help="Overwrite an existing output NetCDF.",
    )
    parser.add_argument(
        "--compression",
        type=int,
        default=4,
        choices=range(0, 10),
        metavar="[0-9]",
        help="NetCDF zlib compression level.",
    )
    parser.add_argument(
        "--duplicate-policy",
        choices=("error", "last"),
        default="error",
        help=(
            "What to do if a local CSV contains duplicate timestamps for the "
            "same station/element."
        ),
    )
    parser.add_argument(
        "--list-categories",
        action="store_true",
        help="List categories available for --resolution and exit.",
    )

    args = parser.parse_args()

    if args.list_categories:
        return args

    if not args.category:
        parser.error("category is required unless --list-categories is used")

    if not args.years:
        parser.error("--years is required")

    if args.workers < 1:
        parser.error("--workers must be >= 1")

    return args


# ---------------------------------------------------------------------
# Directory layout
# ---------------------------------------------------------------------

def local_historical_root(root: Path) -> Path:
    return root / "meteorology" / "climate" / "historical_csv"


def local_data_dir(
    root: Path,
    resolution: str,
    category: str,
    year: int,
) -> Path:
    return (
        local_historical_root(root)
        / "data"
        / resolution
        / category
        / str(year)
    )


def local_metadata_dir(root: Path) -> Path:
    return local_historical_root(root) / "metadata"


def default_output_path(
    root: Path,
    resolution: str,
    category: str,
    years: list[int],
) -> Path:
    years_sorted = sorted(set(years))
    if len(years_sorted) == 1:
        year_tag = str(years_sorted[0])
    else:
        year_tag = f"{years_sorted[0]}_{years_sorted[-1]}"

    name = (
        f"CHMI_{resolution.upper()}_"
        f"{category.upper()}_{year_tag}.nc"
    )
    return root / "processed" / name


# ---------------------------------------------------------------------
# Download stage
# ---------------------------------------------------------------------

def ensure_metadata(
    root: Path,
    timeout: int,
    refresh: bool,
) -> dict[str, Path]:
    metadata_dir = local_metadata_dir(root)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}

    print("\nMetadata")
    print("-" * 70)

    for name in METADATA_FILES:
        target = metadata_dir / name
        paths[name] = target

        if target.exists() and target.stat().st_size > 0 and not refresh:
            print(f"[metadata] existing: {target}")
            continue

        print(f"[metadata] downloading: {name}")
        atomic_download(
            url=urljoin(METADATA_BASE, name),
            target=target,
            timeout=timeout,
            overwrite=True,
        )

    return paths


def filter_urls_for_elements(
    files: list[tuple[str, str]],
    elements: Optional[list[str]],
) -> list[tuple[str, str]]:
    """
    Filter from filename suffix where possible.

    Actual station and element identities are still verified from CSV content
    during the build stage.
    """
    if not elements:
        return files

    selected = []
    for filename, url in files:
        # Monthly files currently end as -ELEMENT-YYYYMM.csv.
        # Use re.escape so elements containing punctuation are safe.
        if any(
            re.search(
                rf"-{re.escape(element)}-\d{{6}}\.csv$",
                filename,
                flags=re.IGNORECASE,
            )
            for element in elements
        ):
            selected.append((filename, url))

    return selected


def element_codes_from_filenames(
    files: list[tuple[str, str]],
) -> set[str]:
    """Extract CHMI ELEMENT codes from monthly CSV filenames."""
    codes: set[str] = set()
    for filename, _ in files:
        match = re.search(r"-([^-]+)-\d{6}\.csv$", filename, flags=re.IGNORECASE)
        if match:
            codes.add(match.group(1))
    return codes


def download_data(
    root: Path,
    resolution: str,
    category: str,
    years: list[int],
    elements: Optional[list[str]],
    workers: int,
    timeout: int,
    overwrite: bool,
    progress=None,
    pipeline_total_multiplier: int = 1,
):
    resolution_url = urljoin(DATA_BASE, f"{resolution}/")
    categories = list_subdirectories(resolution_url, timeout)

    if category not in categories:
        raise ValueError(
            f"Category {category!r} is not available under {resolution!r}.\n"
            f"Available categories: {', '.join(categories)}"
        )

    category_url = urljoin(resolution_url, f"{category}/")
    available_years = set(list_subdirectories(category_url, timeout))

    # Discover first, validate ELEMENT codes second, then download.
    # This prevents confusing errors such as using humidity code H with
    # category=temperature.
    discovered_by_year: dict[int, list[tuple[str, str]]] = {}
    available_elements: set[str] = set()

    for year in sorted(set(years)):
        if str(year) not in available_years:
            print(
                f"[WARNING] {category}/{year} is not present on the server; "
                "skipping."
            )
            continue

        year_url = urljoin(category_url, f"{year}/")
        raw_files = list_csv_urls(year_url, timeout)
        discovered_by_year[year] = raw_files
        available_elements.update(element_codes_from_filenames(raw_files))

    if not discovered_by_year:
        raise RuntimeError(
            f"No requested years are available for category={category!r}."
        )

    if elements:
        requested_upper = {str(e).upper() for e in elements}
        available_upper = {str(e).upper() for e in available_elements}
        invalid = sorted(requested_upper - available_upper)

        if invalid:
            available_text = ", ".join(sorted(available_elements)) or "none detected"
            raise ValueError(
                "Requested CHMI ELEMENT code(s) are not available for the "
                f"selected category/year(s): {', '.join(invalid)}\n"
                f"Category: {category}\n"
                f"Years: {', '.join(map(str, sorted(discovered_by_year)))}\n"
                f"Available ELEMENT codes: {available_text}\n"
                "Either use one of those codes or omit --elements to download "
                "all elements in the category."
            )

    all_tasks: list[tuple[str, Path]] = []

    for year, raw_files in sorted(discovered_by_year.items()):
        files = filter_urls_for_elements(raw_files, elements)

        year_dir = local_data_dir(root, resolution, category, year)
        year_dir.mkdir(parents=True, exist_ok=True)

        print(
            f"[{year}] server files selected: {len(files):,} "
            f"-> {year_dir}"
        )

        for filename, url in files:
            all_tasks.append((url, year_dir / filename))

    if not all_tasks:
        available_text = ", ".join(sorted(available_elements)) or "none detected"
        raise RuntimeError(
            "No CSV files were selected for download.\n"
            f"Available ELEMENT codes: {available_text}"
        )

    skipped_initial = sum(
        1
        for _, path in all_tasks
        if path.exists() and path.stat().st_size > 0 and not overwrite
    )

    pending = [
        (url, path)
        for url, path in all_tasks
        if overwrite or not (path.exists() and path.stat().st_size > 0)
    ]

    print("\nData download")
    print("-" * 70)
    print(f"Selected files : {len(all_tasks):,}")
    print(f"Already present: {skipped_initial:,}")
    print(f"To download    : {len(pending):,}")
    print(f"Workers        : {workers}")

    failures: list[tuple[str, str]] = []

    # Exactly ONE tqdm instance for --stage all. The same object is passed
    # through download -> CSV scan -> NetCDF write. No nested tqdm bars.
    if progress is None:
        bar = make_progress_bar(len(all_tasks) * pipeline_total_multiplier)
    else:
        bar = progress
        if not bar.total:
            bar.reset(total=len(all_tasks) * pipeline_total_multiplier)

    if skipped_initial:
        bar.update(skipped_initial)

    if pending:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(
                    atomic_download,
                    url,
                    path,
                    timeout,
                    True,
                ): (url, path)
                for url, path in pending
            }

            for future in as_completed(futures):
                url, path = futures[future]
                try:
                    future.result()
                except Exception as exc:
                    failures.append((url, str(exc)))
                    tqdm.write(f"[FAILED] {path.name}: {exc}")
                finally:
                    bar.update(1)

    if failures:
        print("\nDownload failures:")
        for url, reason in failures[:20]:
            print(f"  {url}\n    {reason}")
        if len(failures) > 20:
            print(f"  ... plus {len(failures) - 20} more failures")

        raise RuntimeError(
            f"{len(failures)} file(s) failed. Re-run the command; completed "
            "files will be skipped."
        )

    return len(all_tasks), bar


# ---------------------------------------------------------------------
# Local CSV scan
# ---------------------------------------------------------------------

def local_csv_files(
    root: Path,
    resolution: str,
    category: str,
    years: list[int],
) -> list[Path]:
    paths: list[Path] = []
    for year in sorted(set(years)):
        d = local_data_dir(root, resolution, category, year)
        if not d.exists():
            continue
        paths.extend(sorted(d.glob("*.csv")))
    return sorted(paths)


def read_csv_identity(path: Path) -> tuple[str, str]:
    """
    Read only the header + first data row to identify STATION and ELEMENT.

    This avoids loading each CSV during the planning scan.
    """
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        row = next(reader, None)

    if not header:
        raise ValueError("empty CSV header")
    if not row:
        raise ValueError("CSV has no data rows")

    # CHMI CSVs often have a trailing empty column.
    header = [h.strip() for h in header]

    try:
        i_station = header.index("STATION")
        i_element = header.index("ELEMENT")
    except ValueError as exc:
        raise ValueError(
            f"required STATION/ELEMENT columns not found; header={header}"
        ) from exc

    station = row[i_station].strip()
    element = row[i_element].strip()

    if not station or not element:
        raise ValueError("blank STATION or ELEMENT in first data row")

    return station, element


def scan_local_files(
    files: list[Path],
    requested_elements: Optional[list[str]],
    progress=None,
) -> tuple[list[Path], list[str], list[str]]:
    selected: list[Path] = []
    stations: set[str] = set()
    elements: set[str] = set()

    requested = set(requested_elements or [])

    failures = []

    for path in files:
        try:
            station, element = read_csv_identity(path)
        except Exception as exc:
            failures.append((path, str(exc)))
        else:
            if not requested or element in requested:
                selected.append(path)
                stations.add(station)
                elements.add(element)
        finally:
            if progress is not None:
                progress.update(1)

    if failures:
        print("\n[WARNING] CSV identity scan failures:")
        for p, reason in failures[:20]:
            print(f"  {p}: {reason}")
        if len(failures) > 20:
            print(f"  ... plus {len(failures) - 20} more")

    if not selected:
        raise RuntimeError("No usable local CSV files remained after scanning.")

    return selected, sorted(stations), sorted(elements)


# ---------------------------------------------------------------------
# Time axis
# ---------------------------------------------------------------------

def build_time_seconds(years: list[int], step_seconds: int) -> np.ndarray:
    """
    Build a sorted UTC epoch-second time axis for exactly the requested years.

    Non-consecutive years are supported without filling intervening years.
    """
    chunks = []

    for year in sorted(set(years)):
        start = int(
            datetime(year, 1, 1, tzinfo=timezone.utc).timestamp()
        )
        stop = int(
            datetime(year + 1, 1, 1, tzinfo=timezone.utc).timestamp()
        )
        chunks.append(
            np.arange(start, stop, step_seconds, dtype=np.int64)
        )

    if not chunks:
        return np.empty(0, dtype=np.int64)

    return np.concatenate(chunks)


def expected_chunk_time(resolution: str) -> int:
    """Chunk approximately one month of time per station."""
    if resolution == "10min":
        return 31 * 24 * 6       # 4464
    if resolution == "1hour":
        return 31 * 24           # 744
    return 1024


# ---------------------------------------------------------------------
# Metadata matching
# ---------------------------------------------------------------------

def prepare_meta1(meta1: pd.DataFrame) -> pd.DataFrame:
    m = meta1.copy()

    for c in ("GEOGR1", "GEOGR2", "ELEVATION"):
        m[c] = pd.to_numeric(m[c], errors="coerce")

    m["_begin_dt"] = m["BEGIN_DATE"].map(parse_chmi_datetime)
    m["_end_dt"] = m["END_DATE"].map(parse_chmi_datetime)

    return m


def prepare_meta2(meta2: pd.DataFrame) -> pd.DataFrame:
    m = meta2.copy()

    m["HEIGHT"] = pd.to_numeric(m["HEIGHT"], errors="coerce")
    m["_begin_dt"] = m["BEGIN_DATE"].map(parse_chmi_datetime)
    m["_end_dt"] = m["END_DATE"].map(parse_chmi_datetime)

    return m


def filter_metadata_history(
    df: pd.DataFrame,
    stations: set[str],
    periods: list[tuple[datetime, datetime]],
    *,
    elements: Optional[set[str]] = None,
    obs_type: Optional[str] = None,
) -> pd.DataFrame:
    mask = df["WSI"].isin(stations)

    if elements is not None:
        mask &= df["EG_EL_ABBREVIATION"].isin(elements)

    if obs_type is not None:
        mask &= df["OBS_TYPE"].eq(obs_type)

    subset = df.loc[mask].copy()

    overlap = [
        overlaps_any_period(b, e, periods)
        for b, e in zip(subset["_begin_dt"], subset["_end_dt"])
    ]
    subset = subset.loc[overlap].copy()

    return subset.reset_index(drop=True)


def choose_representative_rows(
    df: pd.DataFrame,
    group_cols: list[str],
    periods: list[tuple[datetime, datetime]],
) -> pd.DataFrame:
    """
    Choose the record with maximum temporal overlap with selected years.

    This is only for convenient static metadata coordinates. Full metadata
    histories are retained separately in the output.
    """
    if df.empty:
        return df.copy()

    d = df.copy()
    d["_overlap"] = [
        total_overlap_seconds(b, e, periods)
        for b, e in zip(d["_begin_dt"], d["_end_dt"])
    ]

    # Latest begin wins a tie.
    d = d.sort_values(
        group_cols + ["_overlap", "_begin_dt"],
        ascending=[True] * len(group_cols) + [False, False],
        kind="stable",
    )

    return d.drop_duplicates(group_cols, keep="first").copy()


# ---------------------------------------------------------------------
# NetCDF helpers
# ---------------------------------------------------------------------

def write_vlen_string_variable(
    nc: Dataset,
    name: str,
    dim: str,
    values: Iterable[object],
    *,
    long_name: Optional[str] = None,
):
    var = nc.createVariable(name, str, (dim,))
    vals = [
        "" if x is None or (isinstance(x, float) and math.isnan(x)) else str(x)
        for x in values
    ]
    var[:] = np.asarray(vals, dtype=object)
    if long_name:
        var.long_name = long_name
    return var


def safe_element_var_name(element: str) -> str:
    """
    Create a NetCDF-safe variable name while retaining original ELEMENT in attrs.
    """
    name = re.sub(r"[^0-9A-Za-z_]", "_", element)
    if not name:
        name = "element"
    if name[0].isdigit():
        name = "el_" + name
    return name


def unique_var_names(elements: list[str]) -> dict[str, str]:
    mapping = {}
    used = set()

    for element in elements:
        base = safe_element_var_name(element)
        name = base
        i = 2
        while name in used:
            name = f"{base}_{i}"
            i += 1
        used.add(name)
        mapping[element] = name

    return mapping


def metadata_value_map(meta4: pd.DataFrame) -> dict[int, str]:
    out = {}
    for _, row in meta4.iterrows():
        try:
            code = int(float(str(row["FLAG2"])))
        except Exception:
            continue
        out[code] = str(row["DESCRIPTION"])
    return out


def element_flag_description_map(
    meta3: pd.DataFrame,
    element: str,
) -> dict[str, str]:
    subset = meta3[
        meta3["EL_ABBREVIATION"].astype(str).eq(element)
    ]
    return {
        str(row["FLAG"]): str(row["DESCRIPTION"])
        for _, row in subset.iterrows()
    }


def add_global_attributes(
    nc: Dataset,
    *,
    resolution: str,
    category: str,
    years: list[int],
    elements: list[str],
):
    nc.title = (
        f"CHMI historical {resolution} {category} station observations"
    )
    nc.institution = "Czech Hydrometeorological Institute (CHMI)"
    nc.source = (
        "https://opendata.chmi.cz/meteorology/climate/historical_csv/"
    )
    nc.history = (
        f"{datetime.now(timezone.utc).isoformat()} - "
        "downloaded/converted with chmi_historical_csv_to_netcdf_v5.py"
    )
    # CF + ACDD conventions: CF handles scientific semantics/coordinates;
    # ACDD adds discovery-oriented dataset metadata commonly used by Earth-science archives.
    nc.Conventions = "CF-1.13, ACDD-1.3"
    nc.featureType = "timeSeries"
    nc.cdm_data_type = "TimeSeries"
    nc.standard_name_vocabulary = "CF Standard Name Table v94"
    nc.summary = (
        f"Station-based {resolution} {category} observations from the Czech "
        "Hydrometeorological Institute (CHMI), converted from the CHMI historical "
        "CSV archive to a CF discrete-sampling-geometry time-series NetCDF. "
        "Observation values and CHMI quality/flag fields are preserved; no "
        "scientific quality filtering is applied during conversion."
    )
    nc.keywords = (
        f"meteorology, climatology, station observations, Czech Republic, CHMI, {category}"
    )
    nc.license = "Creative Commons Attribution 4.0 International (CC BY 4.0)"
    nc.references = (
        "https://opendata.chmi.cz/meteorology/climate/Klimatologicka_data_popis.pdf"
    )
    nc.metadata_link = "https://opendata.chmi.cz/meteorology/climate/historical_csv/metadata/"
    nc.processing_level = "Format conversion and metadata harmonization; observation values unchanged"
    nc.date_created = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    nc.chmi_resolution = resolution
    nc.chmi_category = category
    nc.chmi_years = ",".join(str(y) for y in sorted(set(years)))
    nc.chmi_elements = ",".join(elements)
    nc.metadata_join = (
        "data.STATION -> meta1/meta2.WSI; "
        "data.ELEMENT -> meta2.EG_EL_ABBREVIATION"
    )
    nc.metadata_validity_note = (
        "Representative static station metadata uses the record with greatest "
        "temporal overlap with the requested years. Complete date-valid meta1 "
        "and meta2 histories are stored in NetCDF groups "
        "metadata/station_history and metadata/element_history."
    )
    nc.time_zone = "UTC"


def create_netcdf(
    output: Path,
    *,
    stations: list[str],
    elements: list[str],
    time_seconds: np.ndarray,
    resolution: str,
    category: str,
    years: list[int],
    meta1_hist: pd.DataFrame,
    meta1_repr: pd.DataFrame,
    meta2_hist: pd.DataFrame,
    meta2_repr: pd.DataFrame,
    meta3: pd.DataFrame,
    meta4: pd.DataFrame,
    compression: int,
) -> tuple[Dataset, dict[str, dict[str, object]], dict[str, str]]:
    """
    Create output structure and return open Dataset plus data variable handles.
    """
    output.parent.mkdir(parents=True, exist_ok=True)

    nc = Dataset(output, "w", format="NETCDF4")
    add_global_attributes(
        nc,
        resolution=resolution,
        category=category,
        years=years,
        elements=elements,
    )

    n_station = len(stations)
    n_time = len(time_seconds)

    nc.createDimension("station", n_station)
    nc.createDimension("time", n_time)

    # Coordinates
    station_var = nc.createVariable("station", str, ("station",))
    station_var[:] = np.asarray(stations, dtype=object)
    station_var.long_name = "WIGOS Station Identifier (WSI)"
    station_var.standard_name = "platform_id"
    station_var.cf_role = "timeseries_id"

    time_var = nc.createVariable("time", "i8", ("time",))
    time_var[:] = time_seconds
    time_var.units = "seconds since 1970-01-01 00:00:00 UTC"
    time_var.calendar = "proleptic_gregorian"
    time_var.standard_name = "time"
    time_var.long_name = "Time of observation"
    time_var.axis = "T"

    if n_time:
        first_dt = datetime.fromtimestamp(int(time_seconds[0]), tz=timezone.utc)
        last_dt = datetime.fromtimestamp(int(time_seconds[-1]), tz=timezone.utc)
        nc.time_coverage_start = first_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
        nc.time_coverage_end = last_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
        if resolution == "10min":
            nc.time_coverage_resolution = "PT10M"
        elif resolution == "1hour":
            nc.time_coverage_resolution = "PT1H"

    station_to_i = {wsi: i for i, wsi in enumerate(stations)}

    # Representative station metadata
    rep1 = meta1_repr.set_index("WSI", drop=False) if not meta1_repr.empty else None

    lat = np.full(n_station, np.nan, dtype=np.float64)
    lon = np.full(n_station, np.nan, dtype=np.float64)
    elevation = np.full(n_station, np.nan, dtype=np.float32)
    names = [""] * n_station
    ghids = [""] * n_station
    meta_available = np.zeros(n_station, dtype=np.uint8)

    if rep1 is not None:
        for wsi, i in station_to_i.items():
            if wsi not in rep1.index:
                continue

            row = rep1.loc[wsi]
            # set_index may return DataFrame only if duplicate survived unexpectedly.
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]

            lon[i] = pd.to_numeric(row["GEOGR1"], errors="coerce")
            lat[i] = pd.to_numeric(row["GEOGR2"], errors="coerce")
            elevation[i] = pd.to_numeric(row["ELEVATION"], errors="coerce")
            names[i] = "" if pd.isna(row["FULL_NAME"]) else str(row["FULL_NAME"])
            ghids[i] = "" if pd.isna(row["GH_ID"]) else str(row["GH_ID"])
            meta_available[i] = 1

    v = nc.createVariable("lat", "f8", ("station",), fill_value=np.nan)
    v[:] = lat
    v.standard_name = "latitude"
    v.long_name = "Station latitude"
    v.units = "degrees_north"
    v.comment = "Representative station metadata; full date-valid history is in metadata/station_history."

    v = nc.createVariable("lon", "f8", ("station",), fill_value=np.nan)
    v[:] = lon
    v.standard_name = "longitude"
    v.long_name = "Station longitude"
    v.units = "degrees_east"
    v.comment = "Representative station metadata; full date-valid history is in metadata/station_history."

    v = nc.createVariable(
        "elevation",
        "f4",
        ("station",),
        fill_value=np.float32(np.nan),
    )
    v[:] = elevation
    v.standard_name = "surface_altitude"
    v.long_name = "Station surface altitude"
    v.units = "m"
    v.positive = "up"
    v.comment = "Representative CHMI ELEVATION metadata; full date-valid history is in metadata/station_history."

    station_name_var = write_vlen_string_variable(
        nc, "station_name", "station", names, long_name="Station name"
    )
    station_name_var.standard_name = "platform_name"

    gh_id_var = write_vlen_string_variable(
        nc, "GH_ID", "station", ghids, long_name="CHMI GH station identifier"
    )
    gh_id_var.comment = "CHMI station identifier supplied in meta1."

    # ACDD discovery extents derived from the representative station metadata.
    finite_lat = lat[np.isfinite(lat)]
    finite_lon = lon[np.isfinite(lon)]
    finite_elev = elevation[np.isfinite(elevation)]
    if finite_lat.size:
        nc.geospatial_lat_min = float(finite_lat.min())
        nc.geospatial_lat_max = float(finite_lat.max())
        nc.geospatial_lat_units = "degrees_north"
    if finite_lon.size:
        nc.geospatial_lon_min = float(finite_lon.min())
        nc.geospatial_lon_max = float(finite_lon.max())
        nc.geospatial_lon_units = "degrees_east"
    if finite_elev.size:
        nc.geospatial_vertical_min = float(finite_elev.min())
        nc.geospatial_vertical_max = float(finite_elev.max())
        nc.geospatial_vertical_units = "m"
        nc.geospatial_vertical_positive = "up"

    # -----------------------------------------------------------------
    # Full metadata histories are stored in NetCDF groups so the root
    # xarray Dataset stays simple. Each history record carries its WSI
    # directly; no integer pointer/index variable is needed.
    # -----------------------------------------------------------------
    metadata_group = nc.createGroup("metadata")

    station_history = metadata_group.createGroup("station_history")
    station_history.createDimension("record", len(meta1_hist))
    station_history.description = (
        "Date-valid CHMI meta1 station history for stations represented "
        "in the root dataset. WSI identifies the station directly."
    )

    if len(meta1_hist):
        write_vlen_string_variable(
            station_history, "WSI", "record", meta1_hist["WSI"],
            long_name="WIGOS Station Identifier",
        )
        write_vlen_string_variable(
            station_history, "begin_date", "record", meta1_hist["BEGIN_DATE"]
        )
        write_vlen_string_variable(
            station_history, "end_date", "record", meta1_hist["END_DATE"]
        )
        write_vlen_string_variable(
            station_history, "GH_ID", "record", meta1_hist["GH_ID"]
        )
        write_vlen_string_variable(
            station_history, "station_name", "record", meta1_hist["FULL_NAME"]
        )

        for name, source, dtype, units in (
            ("lon", "GEOGR1", "f8", "degrees_east"),
            ("lat", "GEOGR2", "f8", "degrees_north"),
            ("elevation", "ELEVATION", "f4", "m"),
        ):
            vals = pd.to_numeric(meta1_hist[source], errors="coerce").to_numpy()
            vv = station_history.createVariable(
                name, dtype, ("record",), fill_value=np.nan
            )
            vv[:] = vals
            vv.units = units

    element_history = metadata_group.createGroup("element_history")
    element_history.createDimension("record", len(meta2_hist))
    element_history.description = (
        "Date-valid CHMI meta2 element/sensor history for observations "
        "represented in the root dataset. WSI and element identify each record."
    )

    if len(meta2_hist):
        for name, source in (
            ("WSI", "WSI"),
            ("obs_type", "OBS_TYPE"),
            ("element", "EG_EL_ABBREVIATION"),
            ("begin_date", "BEGIN_DATE"),
            ("end_date", "END_DATE"),
            ("name", "NAME"),
            ("unit", "UN_DESCRIPTION"),
            ("schedule", "SCHEDULE"),
        ):
            write_vlen_string_variable(
                element_history, name, "record", meta2_hist[source]
            )

        vv = element_history.createVariable(
            "sensor_height",
            "f4",
            ("record",),
            fill_value=np.float32(np.nan),
        )
        vv[:] = pd.to_numeric(
            meta2_hist["HEIGHT"], errors="coerce"
        ).to_numpy(dtype=np.float32)
        vv.units = "m"

    # Convenience representative sensor height per element/station.
    rep2_lookup = {}
    if not meta2_repr.empty:
        for _, row in meta2_repr.iterrows():
            rep2_lookup[
                (str(row["WSI"]), str(row["EG_EL_ABBREVIATION"]))
            ] = row

    # Main observation variables
    chunk_time = min(expected_chunk_time(resolution), max(1, n_time))
    zlib = compression > 0
    varmap: dict[str, dict[str, object]] = {}
    element_to_varname = unique_var_names(elements)

    qmap = metadata_value_map(meta4)

    for element in elements:
        varname = element_to_varname[element]

        # Derive stable element metadata from meta2.
        em = meta2_hist[
            meta2_hist["EG_EL_ABBREVIATION"].astype(str).eq(element)
        ]

        names_unique = sorted(
            set(
                str(x)
                for x in em["NAME"].dropna().tolist()
                if str(x).strip()
            )
        )
        units_unique = sorted(
            set(
                str(x)
                for x in em["UN_DESCRIPTION"].dropna().tolist()
                if str(x).strip()
            )
        )

        value_var = nc.createVariable(
            varname,
            "f4",
            ("station", "time"),
            fill_value=np.float32(np.nan),
            zlib=zlib,
            complevel=compression,
            shuffle=True,
            chunksizes=(1, chunk_time),
        )
        value_var.chmi_element = element
        value_var.coverage_content_type = "physicalMeasurement"

        if len(names_unique) == 1:
            value_var.long_name = names_unique[0]
        elif names_unique:
            value_var.long_name = " / ".join(names_unique)

        if len(units_unique) == 1:
            value_var.units = units_unique[0]
        elif units_unique:
            value_var.chmi_units_seen = json.dumps(
                units_unique, ensure_ascii=False
            )

        # CF mapping for CHMI relative humidity.  Preserve the original CHMI
        # name as provenance while using an English scientific long_name.
        if element == "H":
            if hasattr(value_var, "long_name"):
                value_var.chmi_long_name = value_var.long_name
            value_var.standard_name = "relative_humidity"
            value_var.long_name = "Relative humidity"
            if not units_unique:
                value_var.units = "%"

        value_var.comment = (
            f"CHMI element {element}; values are preserved from the source CSV. "
            "No quality filtering is applied during conversion."
        )

        quality_var = nc.createVariable(
            f"{varname}_quality",
            "u1",
            ("station", "time"),
            fill_value=np.uint8(255),
            zlib=zlib,
            complevel=compression,
            shuffle=True,
            chunksizes=(1, chunk_time),
        )
        quality_var.long_name = f"Quality code for {getattr(value_var, 'long_name', element)}"
        quality_var.units = "1"
        quality_var.coverage_content_type = "qualityInformation"
        quality_var.chmi_quality_code_descriptions = json.dumps(
            qmap, ensure_ascii=False, sort_keys=True
        )

        if qmap:
            codes = sorted(qmap)
            quality_var.flag_values = np.asarray(codes, dtype=np.uint8)
            # NetCDF flag_meanings cannot contain spaces inside one meaning.
            quality_var.flag_meanings = " ".join(
                re.sub(r"\s+", "_", qmap[c].strip())
                for c in codes
            )

        flag_var = nc.createVariable(
            f"{varname}_flag",
            "S1",
            ("station", "time"),
            zlib=zlib,
            complevel=compression,
            shuffle=True,
            chunksizes=(1, chunk_time),
        )
        flag_var.long_name = f"CHMI element-specific flag for {getattr(value_var, 'long_name', element)}"
        flag_var.coverage_content_type = "qualityInformation"
        flag_var.chmi_flag_descriptions = json.dumps(
            element_flag_description_map(meta3, element),
            ensure_ascii=False,
            sort_keys=True,
        )

        height = np.full(n_station, np.nan, dtype=np.float32)
        for wsi, i in station_to_i.items():
            row = rep2_lookup.get((wsi, element))
            if row is not None:
                height[i] = pd.to_numeric(row["HEIGHT"], errors="coerce")

        height_var = nc.createVariable(
            f"{varname}_sensor_height",
            "f4",
            ("station",),
            fill_value=np.float32(np.nan),
        )
        height_var[:] = height
        height_var.units = "m"
        height_var.long_name = (
            f"CHMI sensor height for {getattr(value_var, 'long_name', element)}"
        )
        height_var.chmi_metadata_field = "HEIGHT"
        height_var.comment = (
            "Representative value from CHMI meta2. CHMI documentation names the "
            "field HEIGHT but does not define its vertical datum/reference in the "
            "published description; therefore no CF standard_name or positive "
            "attribute is asserted here. Full date-valid history is in "
            "metadata/element_history."
        )

        # CF auxiliary coordinates: these station-level variables describe where
        # each observation series is measured.  Because the sensor height is
        # listed here, xarray promotes it to Coordinates on open_dataset().
        coordinate_names = (
            f"lat lon elevation station_name GH_ID {varname}_sensor_height"
        )
        value_var.coordinates = coordinate_names
        quality_var.coordinates = coordinate_names
        flag_var.coordinates = coordinate_names

        # CF ancillary_variables explicitly links the data to its QC variables.
        value_var.ancillary_variables = f"{varname}_quality {varname}_flag"

        varmap[element] = {
            "value": value_var,
            "quality": quality_var,
            "flag": flag_var,
            "varname": varname,
        }

    return nc, varmap, element_to_varname




def iter_with_progress(items, progress=None):
    """Yield items and update a shared tqdm bar exactly once per item."""
    for item in items:
        try:
            yield item
        finally:
            if progress is not None:
                progress.update(1)


# ---------------------------------------------------------------------
# Streaming CSV -> NetCDF
# ---------------------------------------------------------------------

CSV_USECOLS = ["STATION", "ELEMENT", "DT", "VALUE", "FLAG", "QUALITY"]


def read_observation_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        usecols=CSV_USECOLS,
        dtype={
            "STATION": "string",
            "ELEMENT": "string",
            "DT": "string",
            "FLAG": "string",
        },
        low_memory=False,
    )

    df["VALUE"] = pd.to_numeric(df["VALUE"], errors="coerce")
    df["QUALITY"] = pd.to_numeric(df["QUALITY"], errors="coerce")
    return df


def write_indexed(
    var,
    station_i: int,
    indices: np.ndarray,
    values: np.ndarray,
) -> None:
    """
    Prefer fast slice assignment when indices are contiguous.
    """
    if len(indices) == 0:
        return

    if len(indices) == 1:
        var[station_i, int(indices[0])] = values[0]
        return

    if np.all(np.diff(indices) == 1):
        start = int(indices[0])
        stop = int(indices[-1]) + 1
        var[station_i, start:stop] = values
    else:
        var[station_i, indices] = values


def stream_csvs_into_netcdf(
    nc: Dataset,
    varmap: dict[str, dict[str, object]],
    files: list[Path],
    stations: list[str],
    time_seconds: np.ndarray,
    *,
    duplicate_policy: str,
    progress=None,
) -> dict[str, int]:
    station_to_i = {wsi: i for i, wsi in enumerate(stations)}

    stats = {
        "files": 0,
        "rows": 0,
        "rows_written": 0,
        "rows_invalid_time": 0,
        "rows_outside_selected_time": 0,
        "duplicate_rows": 0,
    }

    for path in iter_with_progress(files, progress):
        df = read_observation_csv(path)

        if df.empty:
            continue

        station_values = df["STATION"].dropna().astype(str).unique()
        element_values = df["ELEMENT"].dropna().astype(str).unique()

        if len(station_values) != 1 or len(element_values) != 1:
            raise ValueError(
                f"{path}: expected one STATION and one ELEMENT per file, "
                f"got stations={station_values}, elements={element_values}"
            )

        station = station_values[0]
        element = element_values[0]

        if station not in station_to_i:
            raise KeyError(f"{path}: station {station!r} not in scan plan")
        if element not in varmap:
            # File may have been excluded by --elements.
            continue

        stats["files"] += 1
        stats["rows"] += len(df)

        parsed = pd.to_datetime(df["DT"], utc=True, errors="coerce")
        valid_time = parsed.notna().to_numpy()
        stats["rows_invalid_time"] += int((~valid_time).sum())

        if not valid_time.any():
            continue

        df = df.loc[valid_time].copy()
        parsed = parsed.loc[valid_time]

        # Duplicate timestamp policy is explicit; no silent averaging.
        dup = parsed.duplicated(keep=False).to_numpy()
        if dup.any():
            ndup = int(dup.sum())
            stats["duplicate_rows"] += ndup

            if duplicate_policy == "error":
                example = df.loc[dup, "DT"].head(5).tolist()
                raise ValueError(
                    f"{path}: duplicate timestamps detected ({ndup} rows). "
                    f"Examples: {example}. "
                    "Use --duplicate-policy last only if that is intended."
                )

            keep = ~parsed.duplicated(keep="last").to_numpy()
            df = df.loc[keep].copy()
            parsed = parsed.loc[keep]

        # Convert UTC timestamps to integer Unix seconds explicitly.
        #
        # IMPORTANT (v5): do not use
        #     parsed.astype("int64") // 1_000_000_000
        # here.  The integer unit of pandas datetimes is version/backend
        # dependent (for example datetime64[us, UTC] in newer environments).
        # Dividing a microsecond integer representation by 1e9 makes the
        # timestamps wrong by a factor of 1000 and causes every observation
        # to miss the NetCDF time grid.  Converting through datetime64[s]
        # makes the requested unit explicit and is independent of pandas'
        # internal datetime resolution.
        seconds = (
            parsed.dt.tz_convert("UTC")
            .dt.tz_localize(None)
            .to_numpy(dtype="datetime64[s]")
            .astype(np.int64)
        )

        idx = np.searchsorted(time_seconds, seconds)
        in_bounds = idx < len(time_seconds)

        exact = np.zeros(len(idx), dtype=bool)
        if in_bounds.any():
            exact[in_bounds] = (
                time_seconds[idx[in_bounds]] == seconds[in_bounds]
            )

        stats["rows_outside_selected_time"] += int((~exact).sum())

        if not exact.any():
            continue

        df = df.iloc[np.flatnonzero(exact)].copy()
        idx = idx[exact].astype(np.int64)

        values = df["VALUE"].to_numpy(dtype=np.float32)

        quality_numeric = pd.to_numeric(
            df["QUALITY"], errors="coerce"
        ).to_numpy(dtype=np.float64)
        quality = np.full(len(df), 255, dtype=np.uint8)
        qvalid = np.isfinite(quality_numeric)
        if qvalid.any():
            qvals = quality_numeric[qvalid]
            if np.any((qvals < 0) | (qvals > 254)):
                raise ValueError(
                    f"{path}: QUALITY contains values outside uint8 range."
                )
            quality[qvalid] = qvals.astype(np.uint8)

        flag_series = df["FLAG"].fillna("").astype(str)
        too_long = flag_series.str.len() > 1
        if too_long.any():
            examples = flag_series.loc[too_long].unique()[:10].tolist()
            raise ValueError(
                f"{path}: FLAG values longer than one character were found: "
                f"{examples}. The output schema stores CHMI FLAG as S1."
            )

        flags = flag_series.to_numpy(dtype="S1")

        si = station_to_i[station]
        handles = varmap[element]

        write_indexed(handles["value"], si, idx, values)
        write_indexed(handles["quality"], si, idx, quality)
        write_indexed(handles["flag"], si, idx, flags)

        stats["rows_written"] += len(idx)

    if stats["rows"] > 0 and stats["rows_written"] == 0:
        raise RuntimeError(
            "CSV files contained observation rows, but zero rows matched the "
            "selected NetCDF time grid. Check timestamp conversion/alignment "
            "before using the output."
        )

    nc.sync()
    return stats


# ---------------------------------------------------------------------
# Build stage
# ---------------------------------------------------------------------

def build_netcdf(
    root: Path,
    resolution: str,
    category: str,
    years: list[int],
    requested_elements: Optional[list[str]],
    output: Path,
    compression: int,
    overwrite_nc: bool,
    duplicate_policy: str,
    progress=None,
) -> None:
    files = local_csv_files(root, resolution, category, years)
    if not files:
        raise FileNotFoundError(
            "No local CSV files found. Run with --stage all or --stage download first."
        )

    metadata_dir = local_metadata_dir(root)
    meta_paths = {name: metadata_dir / name for name in METADATA_FILES}

    missing_meta_files = [
        str(p) for p in meta_paths.values() if not p.exists()
    ]
    if missing_meta_files:
        raise FileNotFoundError(
            "Missing metadata file(s):\n  " + "\n  ".join(missing_meta_files)
        )

    print("\nBuild plan")
    print("-" * 70)
    print(f"Local CSV files found: {len(files):,}")

    own_progress = progress is None
    if progress is None:
        progress = make_progress_bar(2 * len(files))
    elif not progress.total:
        progress.reset(total=2 * len(files))

    selected_files, stations, elements = scan_local_files(
        files,
        requested_elements=requested_elements,
        progress=progress,
    )

    if progress is not None:
        # Correct the remaining total if --elements excluded local files or
        # the local directory contains files from an earlier broader run.
        progress.total = progress.n + len(selected_files)
        progress.refresh()

    print(f"Selected CSV files   : {len(selected_files):,}")
    print(f"Stations             : {len(stations):,}")
    print(f"Elements             : {', '.join(elements)}")

    # Read metadata.
    meta1 = prepare_meta1(
        read_standard_metadata_csv(meta_paths["meta1.csv"])
    )
    meta2 = prepare_meta2(
        read_standard_metadata_csv(meta_paths["meta2.csv"])
    )
    meta3 = read_standard_metadata_csv(meta_paths["meta3.csv"])
    meta4 = read_meta4(meta_paths["meta4.csv"])

    periods = periods_for_years(years)
    station_set = set(stations)
    element_set = set(elements)
    obs_type = RESOLUTION_INFO[resolution]["obs_type"]

    meta1_hist = filter_metadata_history(
        meta1,
        station_set,
        periods,
    )
    meta2_hist = filter_metadata_history(
        meta2,
        station_set,
        periods,
        elements=element_set,
        obs_type=obs_type,
    )

    meta1_repr = choose_representative_rows(
        meta1_hist,
        ["WSI"],
        periods,
    )
    meta2_repr = choose_representative_rows(
        meta2_hist,
        ["WSI", "EG_EL_ABBREVIATION"],
        periods,
    )

    stations_with_meta1 = set(meta1_hist["WSI"].astype(str))
    missing_station_meta = sorted(station_set - stations_with_meta1)

    expected_pairs = {(w, e) for w in stations for e in elements}
    observed_pairs = set(
        zip(
            meta2_hist["WSI"].astype(str),
            meta2_hist["EG_EL_ABBREVIATION"].astype(str),
        )
    )
    missing_element_meta = sorted(expected_pairs - observed_pairs)

    print(f"meta1 history records : {len(meta1_hist):,}")
    print(f"meta2 history records : {len(meta2_hist):,}")

    if missing_station_meta:
        print(
            f"[WARNING] {len(missing_station_meta)} downloaded station(s) "
            "have no overlapping meta1 record."
        )
        for wsi in missing_station_meta[:20]:
            print(f"  missing meta1: {wsi}")
        if len(missing_station_meta) > 20:
            print(f"  ... plus {len(missing_station_meta) - 20} more")

    # Only report meta2 gaps for station/element combinations that actually
    # occur in selected files, not the full cross-product.
    file_pairs = set()
    for p in selected_files:
        station, element = read_csv_identity(p)
        file_pairs.add((station, element))

    missing_file_pair_meta2 = sorted(file_pairs - observed_pairs)
    if missing_file_pair_meta2:
        print(
            f"[WARNING] {len(missing_file_pair_meta2)} downloaded "
            "station/element pair(s) have no overlapping meta2 record."
        )
        for pair in missing_file_pair_meta2[:20]:
            print(f"  missing meta2: {pair[0]} / {pair[1]}")
        if len(missing_file_pair_meta2) > 20:
            print(
                f"  ... plus {len(missing_file_pair_meta2) - 20} more"
            )

    step_seconds = RESOLUTION_INFO[resolution]["step_seconds"]
    time_seconds = build_time_seconds(years, step_seconds)

    print(f"Time steps           : {len(time_seconds):,}")

    if output.exists() and not overwrite_nc:
        raise FileExistsError(
            f"Output already exists: {output}\n"
            "Use --overwrite-nc to replace it."
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    temp_output = output.with_name(output.name + ".part")

    if temp_output.exists():
        temp_output.unlink()

    print(f"Output               : {output}")
    print(f"Compression           : zlib level {compression}")

    nc = None
    try:
        nc, varmap, element_names = create_netcdf(
            temp_output,
            stations=stations,
            elements=elements,
            time_seconds=time_seconds,
            resolution=resolution,
            category=category,
            years=years,
            meta1_hist=meta1_hist,
            meta1_repr=meta1_repr,
            meta2_hist=meta2_hist,
            meta2_repr=meta2_repr,
            meta3=meta3,
            meta4=meta4,
            compression=compression,
        )

        stats = stream_csvs_into_netcdf(
            nc,
            varmap,
            selected_files,
            stations,
            time_seconds,
            duplicate_policy=duplicate_policy,
            progress=progress,
        )

        nc.processing_statistics = json.dumps(stats, sort_keys=True)
        nc.element_variable_mapping = json.dumps(
            element_names,
            ensure_ascii=False,
            sort_keys=True,
        )
        nc.sync()
        nc.close()
        nc = None

        os.replace(temp_output, output)

    except Exception:
        if nc is not None:
            try:
                nc.close()
            except Exception:
                pass
        raise

    print("\nNetCDF complete")
    print("-" * 70)
    for k, v in stats.items():
        print(f"{k:28s}: {v:,}")
    print(f"Saved: {output}")

    if own_progress:
        progress.close()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    if args.list_categories:
        url = urljoin(DATA_BASE, f"{args.resolution}/")
        categories = list_subdirectories(url, args.timeout)

        print(
            f"Available CHMI historical_csv categories under "
            f"{args.resolution}:"
        )
        for item in categories:
            print(f"  {item}")
        return

    years = sorted(set(args.years))
    root = Path(args.root_dir).expanduser().resolve()

    print("=" * 78)
    print("CHMI HISTORICAL CSV -> CF/ACDD NETCDF")
    print("=" * 78)
    print(f"Resolution : {args.resolution}")
    print(f"Category   : {args.category}")
    print(f"Years      : {', '.join(map(str, years))}")
    print(f"Elements   : {args.elements or 'ALL in category'}")
    print(f"Root       : {root}")
    print(f"Stage      : {args.stage}")
    print("=" * 78)

    if args.stage in ("download", "all"):
        ensure_metadata(
            root,
            timeout=args.timeout,
            refresh=args.refresh_metadata,
        )

    progress = None

    try:
        if args.stage in ("download", "all"):
            _, progress = download_data(
                root=root,
                resolution=args.resolution,
                category=args.category,
                years=years,
                elements=args.elements,
                workers=args.workers,
                timeout=args.timeout,
                overwrite=args.overwrite_downloads,
                progress=None,
                pipeline_total_multiplier=(3 if args.stage == "all" else 1),
            )

        if args.stage == "build":
            # Metadata are required even when only building.
            required = [
                local_metadata_dir(root) / name
                for name in METADATA_FILES
            ]
            if any(not p.exists() for p in required):
                print(
                    "[INFO] Local metadata are incomplete; downloading metadata "
                    "before build."
                )
                ensure_metadata(
                    root,
                    timeout=args.timeout,
                    refresh=False,
                )

        if args.stage in ("build", "all"):
            output = (
                Path(args.output).expanduser().resolve()
                if args.output
                else default_output_path(
                    root,
                    args.resolution,
                    args.category,
                    years,
                )
            )

            build_netcdf(
                root=root,
                resolution=args.resolution,
                category=args.category,
                years=years,
                requested_elements=args.elements,
                output=output,
                compression=args.compression,
                overwrite_nc=args.overwrite_nc,
                duplicate_policy=args.duplicate_policy,
                progress=progress,
            )
    finally:
        if progress is not None:
            progress.close()


if __name__ == "__main__":
    main()
