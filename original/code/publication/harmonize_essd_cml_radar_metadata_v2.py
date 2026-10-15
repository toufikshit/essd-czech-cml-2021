#!/usr/bin/env python3
"""
Harmonize metadata for the three ESSD CML/radar publication NetCDF files.

This script is designed for:

  ESSD_publication_dataset/CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc
  ESSD_publication_dataset/CZ_CML_RAINRATE_5MIN_20200101_20211231.nc
  ESSD_publication_dataset/CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc

What it does
------------
1. DOES NOT change the values or dimensions of the existing scientific arrays.
2. Adds/updates CF-1.13 + ACDD-1.3 style descriptive metadata.
3. Adds endpoint geometry coordinates to the CML rain-rate file by exact cml_id
   matching against the preprocessed CML file.
4. Adds CML path length + endpoint geometry coordinates to the radar file by
   exact cml_id matching against the preprocessed CML file.
5. Adds clear descriptions, units, coordinate relationships, method-axis
   descriptions, provenance, temporal coverage, and geographic coverage.
6. Uses exactly ONE tqdm progress bar.

No CF standard_name is invented where no clean standard name exists.

The script uses safe temporary working copies in the same directory, validates
them, and then atomically replaces the originals. This avoids HDF5 write-lock
errors when the original files are still open read-only in a Jupyter session.
The existing scientific arrays are copied byte-for-byte and are not recomputed.

Run
---
    python harmonize_essd_cml_radar_metadata_v2.py

Optional:
    python harmonize_essd_cml_radar_metadata_v2.py \
        --root ESSD_publication_dataset

Validation only:
    python harmonize_essd_cml_radar_metadata_v2.py --dry-run
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date
from tqdm import tqdm


CF_VERSION = "CF-1.13"
ACDD_VERSION = "ACDD-1.3"
STANDARD_NAME_VOCABULARY = "CF Standard Name Table v94"

DEFAULT_PREPROCESSED = "CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc"
DEFAULT_RAINRATE = "CZ_CML_RAINRATE_5MIN_20200101_20211231.nc"
DEFAULT_RADAR = "CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc"

ENDPOINT_COORDS = (
    "site_a_latitude",
    "site_a_longitude",
    "site_b_latitude",
    "site_b_longitude",
)

CML_GEOMETRY_COORDS = (
    "length",
    "site_a_latitude",
    "site_a_longitude",
    "site_b_latitude",
    "site_b_longitude",
)

PREPROCESSED_AUX = (
    "length frequency polarization "
    "site_a_latitude site_a_longitude "
    "site_b_latitude site_b_longitude detection_limit_mmhr"
)

RAINRATE_AUX_SUBLINK = (
    "length frequency polarization "
    "site_a_latitude site_a_longitude "
    "site_b_latitude site_b_longitude"
)

RAINRATE_AUX_LINK = (
    "length site_a_latitude site_a_longitude "
    "site_b_latitude site_b_longitude"
)

RADAR_AUX = (
    "length site_a_latitude site_a_longitude "
    "site_b_latitude site_b_longitude"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Make the three ESSD CML/radar publication NetCDF files "
            "self-describing without changing their scientific data arrays."
        ),
    )
    p.add_argument(
        "--root",
        default="ESSD_publication_dataset",
        help="Directory containing the three publication NetCDF files.",
    )
    p.add_argument("--preprocessed", default=DEFAULT_PREPROCESSED)
    p.add_argument("--rainrate", default=DEFAULT_RAINRATE)
    p.add_argument("--radar", default=DEFAULT_RADAR)
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate structures and cml_id matching only; do not modify files.",
    )

    # Optional discovery metadata: intentionally not invented.
    p.add_argument("--institution", default=None)
    p.add_argument("--creator-name", default=None)
    p.add_argument("--creator-email", default=None)
    p.add_argument("--creator-institution", default=None)
    p.add_argument("--project", default=None)
    p.add_argument("--license", default=None)
    p.add_argument("--references", default=None)
    p.add_argument("--doi", default=None)

    return p.parse_args()


def as_strings(values) -> np.ndarray:
    """Read a NetCDF string/char coordinate robustly as 1-D Python strings."""
    arr = np.asarray(values)

    if arr.dtype.kind in {"U", "O"}:
        return np.asarray([str(x) for x in arr], dtype=object)

    if arr.dtype.kind == "S":
        if arr.ndim == 1:
            return np.asarray(
                [
                    x.decode("utf-8", errors="replace")
                    if isinstance(x, (bytes, np.bytes_))
                    else str(x)
                    for x in arr
                ],
                dtype=object,
            )

        # Character array, e.g. (n, string_length)
        rows = []
        for row in arr:
            parts = []
            for x in np.ravel(row):
                if isinstance(x, (bytes, np.bytes_)):
                    parts.append(x.decode("utf-8", errors="replace"))
                else:
                    parts.append(str(x))
            rows.append("".join(parts).rstrip("\x00").strip())
        return np.asarray(rows, dtype=object)

    return arr.astype(str).astype(object)


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")


def require_vars(nc: Dataset, names: tuple[str, ...], label: str) -> None:
    missing = [name for name in names if name not in nc.variables]
    if missing:
        raise ValueError(
            f"{label}: missing required variable(s): {', '.join(missing)}"
        )


def append_history(nc: Dataset, text: str) -> None:
    old = str(getattr(nc, "history", "")).strip()
    nc.history = text if not old else f"{old}\n{text}"


def set_optional_global_metadata(nc: Dataset, args: argparse.Namespace) -> None:
    values = {
        "institution": args.institution,
        "creator_name": args.creator_name,
        "creator_email": args.creator_email,
        "creator_institution": args.creator_institution,
        "project": args.project,
        "license": args.license,
        "references": args.references,
        "id": args.doi,
    }

    for name, value in values.items():
        if value:
            setattr(nc, name, value)


def set_attr(var, name: str, value) -> None:
    """Explicit helper so all metadata assignment is easy to audit."""
    setattr(var, name, value)


def time_info(nc: Dataset, time_name: str) -> tuple[str, str, str | None]:
    var = nc.variables[time_name]

    if len(var) == 0:
        raise ValueError(f"{time_name} is empty")

    units = getattr(var, "units", None)
    if not units:
        raise ValueError(
            f"{time_name} has no encoded NetCDF units; cannot derive coverage."
        )

    calendar = getattr(var, "calendar", "standard")

    first_num = np.asarray(var[0]).item()
    last_num = np.asarray(var[-1]).item()

    first = num2date(
        first_num,
        units=units,
        calendar=calendar,
        only_use_cftime_datetimes=True,
    )
    last = num2date(
        last_num,
        units=units,
        calendar=calendar,
        only_use_cftime_datetimes=True,
    )

    def iso(dt) -> str:
        return (
            f"{int(dt.year):04d}-{int(dt.month):02d}-{int(dt.day):02d}"
            f"T{int(dt.hour):02d}:{int(dt.minute):02d}:"
            f"{int(dt.second):02d}Z"
        )

    resolution = None
    if len(var) > 1:
        a_num = np.asarray(var[0]).item()
        b_num = np.asarray(var[1]).item()

        a = num2date(
            a_num,
            units=units,
            calendar=calendar,
            only_use_cftime_datetimes=True,
        )
        b = num2date(
            b_num,
            units=units,
            calendar=calendar,
            only_use_cftime_datetimes=True,
        )

        seconds = int(round((b - a).total_seconds()))

        if seconds > 0:
            if seconds % 3600 == 0:
                resolution = f"PT{seconds // 3600}H"
            elif seconds % 60 == 0:
                resolution = f"PT{seconds // 60}M"
            else:
                resolution = f"PT{seconds}S"

    return iso(first), iso(last), resolution


def finite_minmax(var) -> tuple[float | None, float | None]:
    values = np.asarray(var[:], dtype=np.float64)
    good = np.isfinite(values)

    if not np.any(good):
        return None, None

    return float(np.nanmin(values[good])), float(np.nanmax(values[good]))


def set_common_global_metadata(
    nc: Dataset,
    *,
    title: str,
    summary: str,
    keywords: str,
    source: str,
    processing_level: str,
    time_name: str,
    args: argparse.Namespace,
) -> None:
    nc.Conventions = f"{CF_VERSION}, {ACDD_VERSION}"
    nc.standard_name_vocabulary = STANDARD_NAME_VOCABULARY
    nc.title = title
    nc.summary = summary
    nc.keywords = keywords
    nc.source = source
    nc.processing_level = processing_level

    t0, t1, tres = time_info(nc, time_name)
    nc.time_coverage_start = t0
    nc.time_coverage_end = t1
    if tres:
        nc.time_coverage_resolution = tres

    if all(name in nc.variables for name in ENDPOINT_COORDS):
        lat_min_a, lat_max_a = finite_minmax(nc.variables["site_a_latitude"])
        lat_min_b, lat_max_b = finite_minmax(nc.variables["site_b_latitude"])
        lon_min_a, lon_max_a = finite_minmax(nc.variables["site_a_longitude"])
        lon_min_b, lon_max_b = finite_minmax(nc.variables["site_b_longitude"])

        lat_values = [
            x for x in (lat_min_a, lat_max_a, lat_min_b, lat_max_b)
            if x is not None
        ]
        lon_values = [
            x for x in (lon_min_a, lon_max_a, lon_min_b, lon_max_b)
            if x is not None
        ]

        if lat_values:
            nc.geospatial_lat_min = min(lat_values)
            nc.geospatial_lat_max = max(lat_values)
            nc.geospatial_lat_units = "degrees_north"

        if lon_values:
            nc.geospatial_lon_min = min(lon_values)
            nc.geospatial_lon_max = max(lon_values)
            nc.geospatial_lon_units = "degrees_east"

    set_optional_global_metadata(nc, args)

    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    nc.date_metadata_modified = now
    append_history(
        nc,
        f"{now}: publication metadata harmonized by "
        "harmonize_essd_cml_radar_metadata_v2.py; existing scientific "
        "data values were not modified.",
    )


def set_time_attrs(var, long_name: str, resolution_text: str) -> None:
    set_attr(var, "standard_name", "time")
    set_attr(var, "long_name", long_name)
    set_attr(var, "axis", "T")
    set_attr(var, "comment", f"UTC observation time at {resolution_text} resolution.")


def set_cml_id_attrs(var) -> None:
    set_attr(var, "long_name", "Commercial microwave link identifier")
    set_attr(
        var,
        "comment",
        "Unique identifier used to match the same commercial microwave link "
        "across the publication NetCDF products.",
    )


def set_sublink_attrs(var) -> None:
    set_attr(var, "long_name", "Directional microwave sublink index")
    set_attr(var, "units", "1")
    set_attr(
        var,
        "comment",
        "Index distinguishing the two directional radio channels associated "
        "with a physical commercial microwave link. The numeric index itself "
        "should not be interpreted as a geographic direction.",
    )


def set_length_attrs(var) -> None:
    set_attr(var, "long_name", "Commercial microwave link path length")
    set_attr(var, "units", "km")
    set_attr(
        var,
        "comment",
        "Distance between endpoints A and B of the corresponding microwave link path.",
    )


def set_endpoint_attrs(nc: Dataset) -> None:
    specs = {
        "site_a_latitude": (
            "latitude",
            "Latitude of commercial microwave link endpoint A",
            "degrees_north",
        ),
        "site_a_longitude": (
            "longitude",
            "Longitude of commercial microwave link endpoint A",
            "degrees_east",
        ),
        "site_b_latitude": (
            "latitude",
            "Latitude of commercial microwave link endpoint B",
            "degrees_north",
        ),
        "site_b_longitude": (
            "longitude",
            "Longitude of commercial microwave link endpoint B",
            "degrees_east",
        ),
    }

    for name, (standard_name, long_name, units) in specs.items():
        var = nc.variables[name]
        set_attr(var, "standard_name", standard_name)
        set_attr(var, "long_name", long_name)
        set_attr(var, "units", units)


def set_frequency_attrs(var) -> None:
    set_attr(var, "standard_name", "radiation_frequency")
    set_attr(var, "long_name", "Microwave carrier frequency")
    # Preserve existing unit where it is already correct; otherwise file-specific
    # caller will set the unit before this helper.
    set_attr(
        var,
        "comment",
        "Carrier frequency of the directional commercial microwave sublink.",
    )


def set_polarization_attrs(var) -> None:
    set_attr(var, "long_name", "Microwave polarization")
    set_attr(
        var,
        "comment",
        "Polarization code associated with the directional commercial microwave "
        "sublink. Values are preserved from the source metadata.",
    )


def load_source_geometry(preprocessed: Path) -> dict[str, dict[str, object]]:
    """
    Load only the small coordinate arrays needed to enrich the other files.
    """
    with Dataset(preprocessed, "r") as nc:
        require_vars(
            nc,
            (
                "cml_id",
                "length",
                "site_a_latitude",
                "site_a_longitude",
                "site_b_latitude",
                "site_b_longitude",
            ),
            "preprocessed CML",
        )

        ids = as_strings(nc.variables["cml_id"][:])

        if len(np.unique(ids)) != len(ids):
            raise ValueError("Preprocessed CML cml_id values are not unique.")

        geometry = {}

        for i, cml_id in enumerate(ids):
            geometry[str(cml_id)] = {
                "length": float(nc.variables["length"][i]),
                "site_a_latitude": float(nc.variables["site_a_latitude"][i]),
                "site_a_longitude": float(nc.variables["site_a_longitude"][i]),
                "site_b_latitude": float(nc.variables["site_b_latitude"][i]),
                "site_b_longitude": float(nc.variables["site_b_longitude"][i]),
            }

        return geometry


def geometry_arrays_for_ids(
    target_ids: np.ndarray,
    geometry: dict[str, dict[str, object]],
    label: str,
) -> dict[str, np.ndarray]:
    missing = [str(x) for x in target_ids if str(x) not in geometry]

    if missing:
        preview = ", ".join(missing[:20])
        more = "" if len(missing) <= 20 else f" ... +{len(missing)-20} more"
        raise ValueError(
            f"{label}: {len(missing)} cml_id value(s) are not present in the "
            f"preprocessed CML file: {preview}{more}"
        )

    arrays = {
        "length": np.asarray(
            [geometry[str(x)]["length"] for x in target_ids],
            dtype=np.float32,
        ),
        "site_a_latitude": np.asarray(
            [geometry[str(x)]["site_a_latitude"] for x in target_ids],
            dtype=np.float32,
        ),
        "site_a_longitude": np.asarray(
            [geometry[str(x)]["site_a_longitude"] for x in target_ids],
            dtype=np.float32,
        ),
        "site_b_latitude": np.asarray(
            [geometry[str(x)]["site_b_latitude"] for x in target_ids],
            dtype=np.float32,
        ),
        "site_b_longitude": np.asarray(
            [geometry[str(x)]["site_b_longitude"] for x in target_ids],
            dtype=np.float32,
        ),
    }

    return arrays


def create_or_validate_1d_coordinate(
    nc: Dataset,
    name: str,
    values: np.ndarray,
    *,
    dimension: str = "cml_id",
) -> None:
    if name in nc.variables:
        var = nc.variables[name]

        if tuple(var.dimensions) != (dimension,):
            raise ValueError(
                f"{name}: existing dimensions {var.dimensions} do not match {(dimension,)}"
            )

        existing = np.asarray(var[:], dtype=np.float64)
        expected = np.asarray(values, dtype=np.float64)

        same = np.allclose(existing, expected, equal_nan=True, rtol=0, atol=1e-6)
        if not same:
            raise ValueError(
                f"{name}: existing values do not match exact cml_id-matched "
                "geometry from the preprocessed CML dataset."
            )
        return

    var = nc.createVariable(
        name,
        "f4",
        (dimension,),
        zlib=False,
        fill_value=np.float32(np.nan),
    )
    var[:] = values


def add_geometry_to_file(
    path: Path,
    geometry: dict[str, dict[str, object]],
    *,
    include_length: bool,
) -> int:
    """
    Add exact cml_id-matched geometry. Returns number of matched CMLs.
    """
    with Dataset(path, "r+") as nc:
        require_vars(nc, ("cml_id",), path.name)
        ids = as_strings(nc.variables["cml_id"][:])

        if len(np.unique(ids)) != len(ids):
            raise ValueError(f"{path.name}: cml_id values are not unique.")

        arrays = geometry_arrays_for_ids(ids, geometry, path.name)

        names = list(ENDPOINT_COORDS)
        if include_length:
            names.insert(0, "length")

        for name in names:
            create_or_validate_1d_coordinate(nc, name, arrays[name])

        nc.sync()

    return len(ids)


def harmonize_preprocessed(path: Path, args: argparse.Namespace) -> None:
    with Dataset(path, "r+") as nc:
        require_vars(
            nc,
            (
                "time",
                "cml_id",
                "sublink_id",
                "tx",
                "rx",
                "trsl",
                "detection_limit_mmhr",
                "length",
                "site_a_latitude",
                "site_a_longitude",
                "site_b_latitude",
                "site_b_longitude",
                "frequency",
                "polarization",
            ),
            "preprocessed CML",
        )

        set_time_attrs(nc.variables["time"], "Observation time", "1 min")
        set_cml_id_attrs(nc.variables["cml_id"])
        set_sublink_attrs(nc.variables["sublink_id"])
        set_length_attrs(nc.variables["length"])
        set_endpoint_attrs(nc)

        freq = nc.variables["frequency"]
        set_attr(freq, "units", "MHz")
        set_frequency_attrs(freq)

        set_polarization_attrs(nc.variables["polarization"])

        dlim = nc.variables["detection_limit_mmhr"]
        set_attr(dlim, "long_name", "CML rain-rate detection limit")
        set_attr(dlim, "units", "mm h-1")
        set_attr(dlim, "coverage_content_type", "auxiliaryInformation")
        set_attr(
            dlim,
            "comment",
            "Link-specific minimum detectable rain-rate threshold used by the "
            "source preprocessing workflow; derived using the ITU-R 2003 "
            "attenuation-rain-rate parameterization.",
        )

        tx = nc.variables["tx"]
        set_attr(tx, "long_name", "Transmitted signal level")
        set_attr(tx, "units", "dBm")
        set_attr(tx, "coverage_content_type", "physicalMeasurement")
        set_attr(tx, "coordinates", PREPROCESSED_AUX)
        set_attr(
            tx,
            "comment",
            "Transmitted microwave signal power level for each directional "
            "sublink. dBm denotes decibels relative to 1 mW.",
        )

        rx = nc.variables["rx"]
        set_attr(rx, "long_name", "Preprocessed received signal level")
        set_attr(rx, "units", "dBm")
        set_attr(rx, "coverage_content_type", "physicalMeasurement")
        set_attr(rx, "coordinates", PREPROCESSED_AUX)
        set_attr(
            rx,
            "comment",
            "Received microwave signal power level after the preprocessing "
            "steps documented by this variable and by the global preprocessing "
            "attributes. dBm denotes decibels relative to 1 mW.",
        )

        trsl = nc.variables["trsl"]
        set_attr(
            trsl,
            "long_name",
            "Transmitted-minus-received signal level",
        )
        set_attr(trsl, "units", "dB")
        set_attr(trsl, "coverage_content_type", "physicalMeasurement")
        set_attr(trsl, "coordinates", PREPROCESSED_AUX)
        set_attr(
            trsl,
            "comment",
            "Total link signal loss computed as transmitted signal level minus "
            "received signal level (TX - RX) after preprocessing. This quantity "
            "contains the baseline link loss together with weather-related "
            "attenuation and should not be interpreted as rain attenuation alone.",
        )

        set_common_global_metadata(
            nc,
            title=(
                "Preprocessed commercial microwave link signal-level "
                "observations over the Czech Republic, 2020-2021"
            ),
            summary=(
                "One-minute commercial microwave link observations prepared for "
                "precipitation-sensing analyses. Each physical link is represented "
                "by two directional sublinks. The file contains transmitted signal "
                "level, preprocessed received signal level, transmitted-minus-"
                "received signal loss, link geometry, carrier frequency, "
                "polarization, and a link-specific rain-rate detection limit. "
                "The preprocessing provenance and filtering steps are retained in "
                "the global and variable attributes."
            ),
            keywords=(
                "commercial microwave links, CML, microwave propagation, "
                "opportunistic sensing, signal level, attenuation, "
                "precipitation sensing, Czech Republic"
            ),
            source="Commercial microwave link signal-level observations",
            processing_level="Quality-controlled and preprocessed CML signal data",
            time_name="time",
            args=args,
        )

        nc.cdm_data_type = "Other"
        nc.dataset_representation = "fixed-path commercial microwave link network"
        nc.comment = (
            "The native representation uses dimensions time, sublink_id, and "
            "cml_id. Auxiliary coordinates describe the fixed path geometry and "
            "radio properties. Existing publication_* attributes preserve the "
            "original publication-preparation provenance."
        )
        nc.sync()


def harmonize_rainrate(path: Path, args: argparse.Namespace) -> None:
    with Dataset(path, "r+") as nc:
        require_vars(
            nc,
            (
                "cml_id",
                "sublink_id",
                "baseline_method",
                "kr_model",
                "waa_method",
                "time_5min",
                "length",
                "frequency",
                "polarization",
                "site_a_latitude",
                "site_a_longitude",
                "site_b_latitude",
                "site_b_longitude",
                "R_sublink_5min",
                "R_mean_5min",
                "R_radar",
                "waa_mean",
                "A_obs_mean",
                "A_rain_mean",
            ),
            "CML rain-rate",
        )

        set_cml_id_attrs(nc.variables["cml_id"])
        set_sublink_attrs(nc.variables["sublink_id"])
        set_time_attrs(nc.variables["time_5min"], "Five-minute observation time", "5 min")
        set_length_attrs(nc.variables["length"])
        set_endpoint_attrs(nc)

        freq = nc.variables["frequency"]
        set_attr(freq, "units", "GHz")
        set_frequency_attrs(freq)

        set_polarization_attrs(nc.variables["polarization"])

        baseline = nc.variables["baseline_method"]
        set_attr(baseline, "long_name", "Dry-reference baseline estimation method")
        set_attr(
            baseline,
            "comment",
            "Method identifier used to estimate the dry-reference signal-loss "
            "baseline before excess attenuation is derived. The publication file "
            "retains only the pycml_constant baseline method.",
        )

        kr = nc.variables["kr_model"]
        set_attr(kr, "long_name", "Attenuation-rain-rate relation model")
        set_attr(
            kr,
            "comment",
            "Identifier of the ITU attenuation-rain-rate coefficient set used "
            "to convert rain-induced specific attenuation to rain rate. "
            "Available publication values are ITU_2003 and ITU_2005.",
        )

        waa = nc.variables["waa_method"]
        set_attr(waa, "long_name", "Wet-antenna attenuation correction method")
        set_attr(
            waa,
            "comment",
            "Wet-antenna attenuation treatment used in the CML rain-rate "
            "retrieval. Publication values are none, leijnse_2008, "
            "pastorek_2021, and schleiss_2013; 'none' means no WAA correction.",
        )

        rs = nc.variables["R_sublink_5min"]
        set_attr(rs, "long_name", "CML-retrieved rain rate for each directional sublink")
        set_attr(rs, "units", "mm h-1")
        set_attr(rs, "coverage_content_type", "physicalMeasurement")
        set_attr(rs, "coordinates", RAINRATE_AUX_SUBLINK)
        set_attr(
            rs,
            "comment",
            "Five-minute rain rate retrieved independently for each directional "
            "CML sublink. Select baseline_method, kr_model, and waa_method to "
            "identify the retrieval configuration.",
        )

        rm = nc.variables["R_mean_5min"]
        set_attr(rm, "long_name", "Link-mean CML-retrieved rain rate")
        set_attr(rm, "units", "mm h-1")
        set_attr(rm, "coverage_content_type", "physicalMeasurement")
        set_attr(rm, "coordinates", RAINRATE_AUX_LINK)
        set_attr(
            rm,
            "comment",
            "Five-minute link-level rain rate obtained after the rain-rate "
            "retrieval is performed separately for the directional sublinks and "
            "the resulting sublink rain-rate estimates are combined at link level. "
            "Select baseline_method, kr_model, and waa_method to identify the "
            "retrieval configuration.",
        )

        rr = nc.variables["R_radar"]
        set_attr(rr, "long_name", "Path-averaged weather-radar rain rate")
        set_attr(rr, "units", "mm h-1")
        set_attr(rr, "coverage_content_type", "physicalMeasurement")
        set_attr(rr, "coordinates", RAINRATE_AUX_LINK)
        set_attr(
            rr,
            "comment",
            "Weather-radar precipitation rate spatially averaged along the "
            "corresponding CML path and sampled on the five-minute time grid. "
            "This is an independent radar reference, not a microwave-link "
            "rain-rate retrieval. The source radar product provides the actual "
            "availability of finite reference values.",
        )

        waa_mean = nc.variables["waa_mean"]
        set_attr(waa_mean, "long_name", "Mean wet-antenna attenuation diagnostic")
        set_attr(waa_mean, "units", "dB")
        set_attr(waa_mean, "coverage_content_type", "auxiliaryInformation")
        set_attr(waa_mean, "coordinates", RAINRATE_AUX_SUBLINK)
        set_attr(
            waa_mean,
            "comment",
            "Summary mean wet-antenna attenuation diagnostic produced by the "
            "source retrieval workflow for each CML sublink and WAA method. "
            "No time dimension is retained for this summary quantity.",
        )

        aobs = nc.variables["A_obs_mean"]
        set_attr(aobs, "long_name", "Mean non-negative observed excess attenuation")
        set_attr(aobs, "units", "dB")
        set_attr(aobs, "coverage_content_type", "auxiliaryInformation")
        set_attr(aobs, "coordinates", RAINRATE_AUX_SUBLINK)
        set_attr(
            aobs,
            "comment",
            "Summary mean of the non-negative observed excess attenuation "
            "diagnostic. The source workflow defines instantaneous observed "
            "excess attenuation as max(TRSL - baseline, 0). No time dimension "
            "is retained for this summary quantity.",
        )

        arain = nc.variables["A_rain_mean"]
        set_attr(arain, "long_name", "Mean rain-associated attenuation diagnostic")
        set_attr(arain, "units", "dB")
        set_attr(arain, "coverage_content_type", "auxiliaryInformation")
        set_attr(arain, "coordinates", RAINRATE_AUX_SUBLINK)
        set_attr(
            arain,
            "comment",
            "Summary mean rain-associated attenuation diagnostic produced by the "
            "source retrieval workflow after applying the selected wet-antenna "
            "attenuation treatment. No time dimension is retained for this "
            "summary quantity.",
        )

        set_common_global_metadata(
            nc,
            title=(
                "Commercial microwave link rain-rate retrievals and "
                "path-averaged radar reference over the Czech Republic, 2020-2021"
            ),
            summary=(
                "Five-minute rain-rate estimates derived from commercial microwave "
                "link signal attenuation. Retrievals are provided per directional "
                "sublink and as a link-level mean for combinations of baseline, "
                "attenuation-rain-rate, and wet-antenna attenuation methods. "
                "A path-averaged weather-radar rain-rate reference is included for "
                "the retained CML paths. Link geometry, carrier frequency, and "
                "polarization are included as auxiliary coordinates, while "
                "summary attenuation diagnostics document important intermediate "
                "retrieval quantities."
            ),
            keywords=(
                "commercial microwave links, CML, rainfall retrieval, rain rate, "
                "wet antenna attenuation, attenuation-rain-rate relation, "
                "weather radar, path-averaged radar, Czech Republic"
            ),
            source=(
                "Derived from preprocessed commercial microwave link observations "
                "and path-averaged weather-radar precipitation"
            ),
            processing_level="Derived five-minute CML rain-rate product",
            time_name="time_5min",
            args=args,
        )

        nc.cdm_data_type = "Other"
        nc.dataset_representation = "fixed-path commercial microwave link retrieval ensemble"
        nc.comment = (
            "Method coordinates are part of the scientific data model: "
            "baseline_method selects the dry-reference baseline, kr_model selects "
            "the attenuation-rain-rate relation, and waa_method selects the "
            "wet-antenna attenuation treatment. Existing processing and "
            "publication_* attributes retain detailed provenance from the source "
            "workflow."
        )
        nc.sync()


def harmonize_radar(path: Path, args: argparse.Namespace) -> None:
    with Dataset(path, "r+") as nc:
        require_vars(
            nc,
            (
                "time",
                "cml_id",
                "rr",
                "length",
                "site_a_latitude",
                "site_a_longitude",
                "site_b_latitude",
                "site_b_longitude",
            ),
            "radar CML-path",
        )

        set_time_attrs(nc.variables["time"], "Radar observation time", "5 min")
        set_cml_id_attrs(nc.variables["cml_id"])
        set_length_attrs(nc.variables["length"])
        set_endpoint_attrs(nc)

        rr = nc.variables["rr"]
        set_attr(rr, "long_name", "Path-averaged weather-radar precipitation rate")
        set_attr(rr, "units", "mm h-1")
        set_attr(rr, "coverage_content_type", "physicalMeasurement")
        set_attr(rr, "coordinates", RADAR_AUX)
        set_attr(
            rr,
            "comment",
            "Five-minute weather-radar precipitation rate spatially averaged "
            "along the straight path between endpoints A and B of each CML. "
            "The cml_id coordinate provides exact alignment with the corresponding "
            "commercial microwave links in the other publication products.",
        )

        set_common_global_metadata(
            nc,
            title=(
                "Path-averaged weather-radar precipitation rate along "
                "commercial microwave links over the Czech Republic, 2021"
            ),
            summary=(
                "Five-minute weather-radar precipitation rates spatially averaged "
                "along the paths of commercial microwave links retained in the "
                "publication CML rain-rate product. Each cml_id is accompanied by "
                "the microwave-link path length and endpoint coordinates so that "
                "the radar reference can be used without consulting another file "
                "for geometry. The actual time coordinate spans the available "
                "radar period in this publication file."
            ),
            keywords=(
                "weather radar, precipitation, rain rate, commercial microwave "
                "links, CML paths, path average, Czech Republic"
            ),
            source="Weather-radar precipitation spatially sampled along CML paths",
            processing_level="Five-minute path-averaged radar reference product",
            time_name="time",
            args=args,
        )

        nc.cdm_data_type = "Other"
        nc.dataset_representation = "weather-radar precipitation sampled along fixed CML paths"
        nc.comment = (
            "The filename and legacy publication_* crop-window attributes refer "
            "to the common publication window, whereas time_coverage_start and "
            "time_coverage_end are calculated from the actual radar time "
            "coordinate. For this file the observed time coordinate begins in "
            "2021. Existing publication_* attributes are retained as provenance."
        )
        nc.sync()


def validate_final(
    preprocessed: Path,
    rainrate: Path,
    radar: Path,
) -> None:
    with Dataset(preprocessed, "r") as nc:
        for name in (
            "tx", "rx", "trsl", "time", "cml_id", "length",
            "site_a_latitude", "site_a_longitude",
            "site_b_latitude", "site_b_longitude",
            "frequency", "polarization", "detection_limit_mmhr",
        ):
            if name not in nc.variables:
                raise RuntimeError(f"Final validation: missing {name} in preprocessed.")

        if "CF-1.13" not in str(getattr(nc, "Conventions", "")):
            raise RuntimeError("Final validation: preprocessed missing CF-1.13.")

    with Dataset(rainrate, "r") as nc:
        for name in (
            "R_sublink_5min", "R_mean_5min", "R_radar",
            "waa_mean", "A_obs_mean", "A_rain_mean",
            "site_a_latitude", "site_a_longitude",
            "site_b_latitude", "site_b_longitude",
        ):
            if name not in nc.variables:
                raise RuntimeError(f"Final validation: missing {name} in rain-rate.")

        for name in ("R_sublink_5min", "R_mean_5min", "R_radar"):
            if not hasattr(nc.variables[name], "long_name"):
                raise RuntimeError(
                    f"Final validation: {name} has no long_name."
                )

    with Dataset(radar, "r") as nc:
        for name in (
            "rr", "length",
            "site_a_latitude", "site_a_longitude",
            "site_b_latitude", "site_b_longitude",
        ):
            if name not in nc.variables:
                raise RuntimeError(f"Final validation: missing {name} in radar.")

        if not hasattr(nc.variables["rr"], "coordinates"):
            raise RuntimeError("Final validation: radar rr has no coordinates attribute.")



def make_safe_work_copy(src: Path) -> Path:
    """
    Create a same-directory working copy.

    GNU cp --reflink=auto is attempted first so filesystems supporting copy-on-
    write can avoid physically duplicating the multi-GB CML file. A normal
    shutil.copy2 fallback is used otherwise.
    """
    work = src.with_name(f".{src.name}.metadata_v2_work_{os.getpid()}")
    work.unlink(missing_ok=True)

    try:
        subprocess.run(
            [
                "cp",
                "--reflink=auto",
                "--preserve=mode,timestamps",
                "--",
                str(src),
                str(work),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        work.unlink(missing_ok=True)
        shutil.copy2(src, work)

    if not work.exists() or work.stat().st_size != src.stat().st_size:
        work.unlink(missing_ok=True)
        raise RuntimeError(f"Safe working copy failed for {src}")

    return work


def install_work_copy(work: Path, target: Path) -> None:
    """Atomically replace target with the already validated same-directory copy."""
    if work.parent != target.parent:
        raise RuntimeError("Working copy must be in the same directory as target.")
    os.replace(work, target)


def main() -> None:
    args = parse_args()

    root = Path(args.root).expanduser().resolve()
    preprocessed = root / args.preprocessed
    rainrate = root / args.rainrate
    radar = root / args.radar

    originals = (preprocessed, rainrate, radar)
    for path in originals:
        require_file(path)

    # Exactly one tqdm object in the whole script.
    total_steps = 10 if not args.dry_run else 4
    work_paths: dict[Path, Path] = {}

    try:
        with tqdm(
            total=total_steps,
            desc="ESSD metadata",
            unit="step",
            dynamic_ncols=True,
        ) as bar:

            # 1. Read the small source geometry arrays from the original file.
            bar.set_postfix_str("reading source CML geometry", refresh=True)
            geometry = load_source_geometry(preprocessed)
            bar.update(1)

            # 2. Validate rain-rate cml_id matching before any modification.
            bar.set_postfix_str("checking rain-rate cml_id matching", refresh=True)
            with Dataset(rainrate, "r") as nc:
                require_vars(nc, ("cml_id",), "rain-rate")
                rain_ids = as_strings(nc.variables["cml_id"][:])
            geometry_arrays_for_ids(rain_ids, geometry, "rain-rate")
            bar.update(1)

            # 3. Validate radar cml_id matching before any modification.
            bar.set_postfix_str("checking radar cml_id matching", refresh=True)
            with Dataset(radar, "r") as nc:
                require_vars(nc, ("cml_id",), "radar")
                radar_ids = as_strings(nc.variables["cml_id"][:])
            geometry_arrays_for_ids(radar_ids, geometry, "radar")
            bar.update(1)

            # 4. Validate the expected scientific variables.
            bar.set_postfix_str("checking existing variables", refresh=True)
            with Dataset(preprocessed, "r") as nc:
                require_vars(nc, ("tx", "rx", "trsl"), "preprocessed")
            with Dataset(rainrate, "r") as nc:
                require_vars(
                    nc,
                    ("R_sublink_5min", "R_mean_5min", "R_radar"),
                    "rain-rate",
                )
            with Dataset(radar, "r") as nc:
                require_vars(nc, ("rr",), "radar")
            bar.update(1)

            if args.dry_run:
                bar.set_postfix_str("validation complete", refresh=True)
                return

            # 5. Make safe same-directory copies. We never open the originals r+.
            bar.set_postfix_str("creating safe working copies", refresh=True)
            for original in originals:
                work_paths[original] = make_safe_work_copy(original)
            pre_work = work_paths[preprocessed]
            rain_work = work_paths[rainrate]
            radar_work = work_paths[radar]
            bar.update(1)

            # 6. Add/verify geometry in the CML rain-rate working copy.
            bar.set_postfix_str("adding geometry to CML rain-rate", refresh=True)
            add_geometry_to_file(
                rain_work,
                geometry,
                include_length=False,
            )
            with Dataset(rain_work, "r") as nc:
                rain_ids = as_strings(nc.variables["cml_id"][:])
                expected = geometry_arrays_for_ids(rain_ids, geometry, "rain-rate")
                existing = np.asarray(nc.variables["length"][:], dtype=np.float64)
                if not np.allclose(
                    existing,
                    expected["length"].astype(np.float64),
                    equal_nan=True,
                    rtol=0,
                    atol=1e-6,
                ):
                    raise ValueError(
                        "Rain-rate length values do not match preprocessed CML geometry."
                    )
            bar.update(1)

            # 7. Add/verify geometry in the radar working copy.
            bar.set_postfix_str("adding geometry to radar CML paths", refresh=True)
            add_geometry_to_file(
                radar_work,
                geometry,
                include_length=True,
            )
            bar.update(1)

            # 8. Write all self-describing metadata to the working copies.
            bar.set_postfix_str("writing self-describing metadata", refresh=True)
            harmonize_preprocessed(pre_work, args)
            harmonize_rainrate(rain_work, args)
            harmonize_radar(radar_work, args)
            bar.update(1)

            # 9. Validate all working copies before touching originals.
            bar.set_postfix_str("validating working copies", refresh=True)
            validate_final(pre_work, rain_work, radar_work)
            bar.update(1)

            # 10. Atomically install only after every working copy validates.
            bar.set_postfix_str("installing validated files", refresh=True)
            install_work_copy(pre_work, preprocessed)
            work_paths.pop(preprocessed, None)
            install_work_copy(rain_work, rainrate)
            work_paths.pop(rainrate, None)
            install_work_copy(radar_work, radar)
            work_paths.pop(radar, None)
            bar.update(1)
            bar.set_postfix_str("complete", refresh=True)

        print("\nMetadata harmonization complete.")
        print(f"Preprocessed CML : {preprocessed}")
        print(f"CML rain rate    : {rainrate}")
        print(f"Radar CML paths  : {radar}")
        print("Scientific array values and dimensions were not recomputed or changed.")
        print("Validated temporary copies were atomically installed over the originals.")

    finally:
        # On any failure, remove unfinished temporary copies and leave the
        # corresponding originals untouched.
        for work in list(work_paths.values()):
            try:
                work.unlink(missing_ok=True)
            except Exception:
                pass


if __name__ == "__main__":
    main()
