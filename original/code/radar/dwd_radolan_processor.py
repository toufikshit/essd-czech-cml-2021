#!/usr/bin/env python3
"""
DWD RADOLAN Radar Data Processor for CML Intersections
=======================================================

Downloads RADOLAN radar data from DWD Open Data, extracts archives,
and computes rain rate time series along Commercial Microwave Link (CML) paths.

Data source:
    https://opendata.dwd.de/climate_environment/CDC/grids_germany/5_minutes/radolan/reproc/2017_002/netCDF/

Pipeline:
    1. Download RADOLAN .tar.gz files (with 2021 supplement fix)
    2. Extract all archives
    3. Calculate intersection weights between CML paths and radar grid
    4. Extract rain rate time series along each CML path
    5. Merge into single NetCDF output

Usage:
    python dwd_radolan_processor.py

Author: [Your Name]
License: MIT
"""

import os
import sys
import glob
import tarfile
import requests
import numpy as np
import xarray as xr
from bs4 import BeautifulSoup
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ============================================================================
# CONFIGURATION
# ============================================================================

# DWD Open Data base URL for RADOLAN 5-min reprocessed data
BASE_URL = "https://opendata.dwd.de/climate_environment/CDC/grids_germany/5_minutes/radolan/reproc/2017_002/netCDF/"

# Years to download
YEARS = range(2021, 2022)

# Input/Output paths
ROOT_DIR = "dwd_radolan"
CMLS_PATH = "txrx_CZ_June2021.nc"          # CML dataset with coordinates
OUTPUT_DIR = "dwd_radolan/RADKLIM_5min/"    # Processed intersection files
FINAL_OUTPUT = "dwd_radolan/cml_rr_RadlongS2021E2022_CZ.nc"  # Merged output

# Processing settings
MAX_WORKERS = 6  # Parallel workers for intersection processing


# ============================================================================
# STEP 1: DOWNLOAD RADOLAN FILES
# ============================================================================

def list_links(url, extension=None):
    """
    Parse DWD HTML directory listing and extract file links.

    Parameters
    ----------
    url : str
        URL of the DWD directory page.
    extension : str, optional
        Filter by file extension (e.g., '.tar.gz').

    Returns
    -------
    list of str
        List of filenames found.
    """
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"[ERROR] Cannot open: {url} — {e}")
        return []

    soup = BeautifulSoup(response.text, "html.parser")
    files = []
    for anchor in soup.find_all("a"):
        href = anchor.get("href")
        if href:
            if extension is None or href.endswith(extension):
                files.append(href)
    return files


def download_file(url, local_path, pbar=None):
    """
    Download file with progress and overwrite protection.

    Parameters
    ----------
    url : str
        URL to download.
    local_path : str
        Local path to save file.
    pbar : tqdm, optional
        Progress bar to update.
    """
    if os.path.exists(local_path):
        if pbar:
            pbar.set_postfix_str(f"Skip: {os.path.basename(local_path)}")
        return True

    try:
        response = requests.get(url, stream=True, timeout=120)
        response.raise_for_status()
        total_size = int(response.headers.get("content-length", 0))

        with open(local_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        if pbar:
            pbar.set_postfix_str(f"OK: {os.path.basename(local_path)}")
        return True

    except Exception as e:
        if pbar:
            pbar.set_postfix_str(f"FAILED: {e}")
        return False


def download_radolan_data(years=YEARS, base_url=BASE_URL, root_dir=ROOT_DIR):
    """
    Download RADOLAN .tar.gz files for specified years.

    Handles the 2021 supplement fix automatically.

    Parameters
    ----------
    years : iterable of int
        Years to download.
    base_url : str
        DWD base URL.
    root_dir : str
        Local root directory.

    Returns
    -------
    int
        Number of files downloaded.
    """
    os.makedirs(root_dir, exist_ok=True)
    total_downloaded = 0

    for year in years:
        print(f"\n{'='*60}")
        print(f"📅 YEAR: {year}")
        print("=" * 60)

        year_dir = os.path.join(root_dir, str(year))
        os.makedirs(year_dir, exist_ok=True)

        # Special handling for 2021 (faulty main file)
        if year == 2021:
            print("⚠️  Year 2021 has faulty main file — fetching from supplement/")
            supp_url = f"{base_url}supplement/"
            supp_files = list_links(supp_url, ".tar.gz")
            supp_corrected = [f for f in supp_files if "2021" in f]

            if not supp_corrected:
                print("[ERROR] No corrected supplement file found for 2021!")
                continue

            pbar = tqdm(supp_corrected, desc=f"Downloading {year}", unit="file")
            for filename in pbar:
                success = download_file(
                    url=supp_url + filename,
                    local_path=os.path.join(year_dir, filename),
                    pbar=pbar,
                )
                if success:
                    total_downloaded += 1
            continue

        # Normal years
        year_url = f"{base_url}{year}/"
        year_files = list_links(year_url, ".tar.gz")

        if not year_files:
            print(f"No .tar.gz files found for {year}")
            continue

        pbar = tqdm(year_files, desc=f"Downloading {year}", unit="file")
        for filename in pbar:
            success = download_file(
                url=year_url + filename,
                local_path=os.path.join(year_dir, filename),
                pbar=pbar,
            )
            if success:
                total_downloaded += 1

    print(f"\n✅ Download complete. Total files: {total_downloaded}")
    return total_downloaded


# ============================================================================
# STEP 2: EXTRACT TAR.GZ ARCHIVES
# ============================================================================

def extract_tar_gz(file_path, extract_to):
    """
    Extract a .tar.gz archive.

    Parameters
    ----------
    file_path : str
        Path to the archive.
    extract_to : str
        Directory to extract into.

    Returns
    -------
    bool
        True if successful.
    """
    try:
        with tarfile.open(file_path, "r:gz") as tar:
            tar.extractall(path=extract_to)
        return True
    except Exception as e:
        print(f"[ERROR] Failed extracting {file_path}: {e}")
        return False


def extract_all_archives(root_dir=ROOT_DIR):
    """
    Extract all .tar.gz files in the directory tree.

    Parameters
    ----------
    root_dir : str
        Root directory to search.

    Returns
    -------
    int
        Number of archives extracted.
    """
    # Find all tar.gz files
    archives = []
    for root, dirs, files in os.walk(root_dir):
        for fname in files:
            if fname.endswith(".tar.gz"):
                archives.append(os.path.join(root, fname))

    if not archives:
        print("No .tar.gz files found to extract.")
        return 0

    print(f"\n📦 Found {len(archives)} archives to extract")

    extracted = 0
    pbar = tqdm(archives, desc="Extracting archives", unit="file")
    for archive_path in pbar:
        pbar.set_postfix_str(os.path.basename(archive_path))
        extract_dir = os.path.dirname(archive_path)
        if extract_tar_gz(archive_path, extract_dir):
            extracted += 1

    print(f"✅ Extracted {extracted} archives")
    return extracted


# ============================================================================
# STEP 3: CALCULATE INTERSECTION WEIGHTS
# ============================================================================

def calculate_intersection_weights(cmls_path, sample_radar_path):
    """
    Calculate sparse intersection weights between CML paths and radar grid.

    Parameters
    ----------
    cmls_path : str
        Path to CML dataset with site coordinates.
    sample_radar_path : str
        Path to a sample RADOLAN file (for grid coordinates).

    Returns
    -------
    xr.Dataset
        Intersection weights dataset.
    """
    import poligrain as plg

    print("\n🔗 Calculating CML-radar intersection weights...")

    cmls = xr.open_dataset(cmls_path)
    ds_sample = xr.open_dataset(sample_radar_path)

    n_cmls = len(cmls.cml_id)
    print(f"   CMLs: {n_cmls}")
    print(f"   Radar grid: {ds_sample.lon.shape} × {ds_sample.lat.shape}")

    weights = plg.spatial.calc_sparse_intersect_weights_for_several_cmls(
        x1_line=cmls.site_a_longitude.values,
        y1_line=cmls.site_a_latitude.values,
        x2_line=cmls.site_b_longitude.values,
        y2_line=cmls.site_b_latitude.values,
        cml_id=cmls.cml_id.values,
        x_grid=ds_sample.lon.values,
        y_grid=ds_sample.lat.values,
        grid_point_location="lower_left",
    )

    print(f"✅ Weights calculated for {n_cmls} CMLs")
    return weights, cmls


# ============================================================================
# STEP 4: PROCESS RADAR FILES → CML INTERSECTIONS
# ============================================================================

def process_single_radar_file(args):
    """
    Process a single radar file (worker function for parallel execution).

    Parameters
    ----------
    args : tuple
        (file_path, weights, output_dir)

    Returns
    -------
    dict
        Result with status and message.
    """
    import poligrain as plg

    file_path, weights_path, output_dir = args
    basename = os.path.basename(file_path).replace(".nc", "_CML_RR.nc")
    outfile = os.path.join(output_dir, basename)

    result = {
        "file": file_path,
        "output": outfile,
        "status": "unknown",
        "message": "",
    }

    # Skip if already processed
    if os.path.exists(outfile):
        result["status"] = "skipped"
        result["message"] = "Already exists"
        return result

    try:
        ds = xr.open_dataset(file_path)

        # Check for rain rate variable
        if "RR" not in ds:
            result["status"] = "error"
            result["message"] = "Missing RR variable"
            return result

        # Load weights
        weights = xr.open_dataset(weights_path)

        # Extract time series along CML paths
        da_ts = plg.spatial.get_grid_time_series_at_intersections(
            grid_data=ds["RR"],
            intersect_weights=weights,
        )

        # Save
        da_ts.to_netcdf(outfile)

        result["status"] = "success"
        result["message"] = f"Saved {basename}"
        return result

    except Exception as e:
        result["status"] = "error"
        result["message"] = str(e)
        return result


def process_radar_files_parallel(
    radar_root=ROOT_DIR,
    weights=None,
    output_dir=OUTPUT_DIR,
    max_workers=MAX_WORKERS,
):
    """
    Process all radar files to extract CML intersection time series.

    Shows detailed progress bar with per-CML status.

    Parameters
    ----------
    radar_root : str
        Root directory with radar .nc files.
    weights : xr.Dataset
        Intersection weights.
    output_dir : str
        Output directory for processed files.
    max_workers : int
        Number of parallel workers.

    Returns
    -------
    dict
        Summary statistics.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Save weights temporarily for worker processes
    weights_path = os.path.join(output_dir, "_temp_weights.nc")
    weights.to_netcdf(weights_path)

    # Find all radar files
    radar_files = sorted(glob.glob(f"{radar_root}/**/*.nc", recursive=True))
    # Exclude our output files and temp files
    radar_files = [f for f in radar_files if "RADKLIM" not in f and "_temp_" not in f]

    if not radar_files:
        print("No radar .nc files found!")
        return {"total": 0, "success": 0, "skipped": 0, "error": 0}

    print(f"\n🛰️  Processing {len(radar_files)} radar files → CML intersections")
    print(f"   Workers: {max_workers}")
    print(f"   Output: {output_dir}")
    print()

    # Prepare arguments
    args_list = [(f, weights_path, output_dir) for f in radar_files]

    # Statistics
    stats = {"total": len(radar_files), "success": 0, "skipped": 0, "error": 0}

    # ══════════════════════════════════════════════════════════════════
    # MAIN PROGRESS BAR
    # ══════════════════════════════════════════════════════════════════
    pbar = tqdm(
        total=len(radar_files),
        desc="🌧️  Radar → CML",
        unit="file",
        ncols=100,
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}",
    )

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_radar_file, args): args[0] for args in args_list}

        for future in as_completed(futures):
            result = future.result()

            # Update statistics
            stats[result["status"]] += 1

            # Update progress bar
            status_icon = {"success": "✓", "skipped": "○", "error": "✗"}.get(result["status"], "?")
            pbar.set_postfix_str(
                f"{status_icon} {os.path.basename(result['file'])} | "
                f"✓{stats['success']} ○{stats['skipped']} ✗{stats['error']}"
            )
            pbar.update(1)

    pbar.close()

    # Cleanup temp weights file
    if os.path.exists(weights_path):
        os.remove(weights_path)

    # Print summary
    print(f"\n{'─'*60}")
    print(f"📊 Processing Summary:")
    print(f"   ✓ Success: {stats['success']}")
    print(f"   ○ Skipped: {stats['skipped']}")
    print(f"   ✗ Errors:  {stats['error']}")
    print(f"{'─'*60}")

    return stats


# ============================================================================
# STEP 5: MERGE OUTPUTS INTO SINGLE NETCDF
# ============================================================================

def merge_intersection_files(input_dir=OUTPUT_DIR, output_file=FINAL_OUTPUT):
    """
    Merge all processed intersection files into a single NetCDF.

    Parameters
    ----------
    input_dir : str
        Directory with processed _CML_RR.nc files.
    output_file : str
        Output merged NetCDF path.

    Returns
    -------
    xr.Dataset
        Merged dataset.
    """
    print(f"\n📎 Merging intersection files...")

    pattern = os.path.join(input_dir, "*_CML_RR.nc")
    files = sorted(glob.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No files found matching {pattern}")

    print(f"   Found {len(files)} files to merge")

    # Progress bar for loading
    chunks = []
    for f in tqdm(files, desc="Loading files", unit="file"):
        try:
            ds = xr.open_dataset(f)
            chunks.append(ds)
        except Exception as e:
            print(f"   [WARNING] Could not load {f}: {e}")

    if not chunks:
        raise ValueError("No valid files to merge!")

    print("   Concatenating along time dimension...")
    merged = xr.concat(chunks, dim="time")
    merged = merged.sortby("time")

    # Add metadata
    merged.attrs["title"] = "RADOLAN rain rate along CML paths"
    merged.attrs["source"] = "DWD RADOLAN reprocessed 2017.002"
    merged.attrs["processing"] = "poligrain intersection weights"

    print(f"   Saving to {output_file}...")
    merged.to_netcdf(
        output_file,
        encoding={"rr": {"dtype": "float32", "zlib": True, "complevel": 4}},
    )

    print(f"✅ Merged output: {output_file}")
    print(f"   Shape: {dict(merged.dims)}")
    print(f"   Time range: {merged.time.values[0]} → {merged.time.values[-1]}")

    return merged


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def find_sample_radar_file(root_dir=ROOT_DIR):
    """Find any .nc radar file to use as grid reference."""
    for root, dirs, files in os.walk(root_dir):
        for f in files:
            if f.endswith(".nc") and "RADKLIM" not in f:
                return os.path.join(root, f)
    return None


def main():
    """
    Main pipeline: download → extract → weights → process → merge.
    """
    print("=" * 70)
    print("  DWD RADOLAN Radar Data Processor for CML Intersections")
    print("=" * 70)

    # Step 1: Download
    print("\n" + "━" * 70)
    print("STEP 1/5: Downloading RADOLAN data")
    print("━" * 70)
    download_radolan_data()

    # Step 2: Extract
    print("\n" + "━" * 70)
    print("STEP 2/5: Extracting archives")
    print("━" * 70)
    extract_all_archives()

    # Step 3: Calculate weights
    print("\n" + "━" * 70)
    print("STEP 3/5: Calculating intersection weights")
    print("━" * 70)

    # Find a sample radar file for grid coordinates
    sample_file = find_sample_radar_file()
    if not sample_file:
        raise FileNotFoundError("No radar .nc files found! Run download/extract first.")

    print(f"   Using grid from: {sample_file}")

    if not os.path.exists(CMLS_PATH):
        raise FileNotFoundError(f"CML dataset not found: {CMLS_PATH}")

    weights, cmls = calculate_intersection_weights(CMLS_PATH, sample_file)

    # Step 4: Process intersections
    print("\n" + "━" * 70)
    print("STEP 4/5: Processing radar → CML intersections")
    print("━" * 70)
    process_radar_files_parallel(weights=weights)

    # Step 5: Merge
    print("\n" + "━" * 70)
    print("STEP 5/5: Merging output files")
    print("━" * 70)
    merge_intersection_files()

    print("\n" + "=" * 70)
    print("🎉 PIPELINE COMPLETE!")
    print("=" * 70)


if __name__ == "__main__":
    main()
