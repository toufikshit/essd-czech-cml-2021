#!/usr/bin/env python3
"""
Inspect the three ESSD publication NetCDF files without loading the large data arrays.

Run:
    python inspect_essd_three_netcdf.py

It prints:
- dimensions
- coordinates
- data variables
- shapes and dtypes
- variable attributes
- global attributes
- NetCDF groups
- compression/chunking information

It does NOT modify any file.
"""

from pathlib import Path
import xarray as xr
from netCDF4 import Dataset

FILES = [
    Path("ESSD_publication_dataset/CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc"),
    Path("ESSD_publication_dataset/CZ_CML_RAINRATE_5MIN_20200101_20211231.nc"),
    Path("ESSD_publication_dataset/CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc"),
]


def short_value(value, maxlen=220):
    text = repr(value)
    return text if len(text) <= maxlen else text[:maxlen] + " ..."


def print_attrs(attrs, indent="      "):
    if not attrs:
        print(f"{indent}(none)")
        return

    for key in sorted(attrs):
        print(f"{indent}{key}: {short_value(attrs[key])}")


def inspect_file(path: Path):
    print("\n" + "=" * 100)
    print(path)
    print("=" * 100)

    if not path.exists():
        print("ERROR: FILE NOT FOUND")
        return

    print(f"Size: {path.stat().st_size / 1024**2:,.2f} MiB")

    # ------------------------------------------------------------
    # Xarray-level structure
    # ------------------------------------------------------------
    ds = xr.open_dataset(path)

    try:
        print("\nDIMENSIONS")
        print("-" * 100)
        for name, size in ds.sizes.items():
            print(f"{name:35s} {size:,}")

        print("\nCOORDINATES")
        print("-" * 100)
        if not ds.coords:
            print("(none)")
        else:
            for name, var in ds.coords.items():
                print(
                    f"\n{name}"
                    f"\n    dims : {var.dims}"
                    f"\n    shape: {var.shape}"
                    f"\n    dtype: {var.dtype}"
                )
                print("    attrs:")
                print_attrs(var.attrs, indent="        ")

        print("\nDATA VARIABLES")
        print("-" * 100)
        if not ds.data_vars:
            print("(none)")
        else:
            for name, var in ds.data_vars.items():
                print(
                    f"\n{name}"
                    f"\n    dims : {var.dims}"
                    f"\n    shape: {var.shape}"
                    f"\n    dtype: {var.dtype}"
                )
                print("    attrs:")
                print_attrs(var.attrs, indent="        ")

        print("\nGLOBAL ATTRIBUTES")
        print("-" * 100)
        print_attrs(ds.attrs, indent="    ")

        # Only inspect time endpoints; don't load main arrays.
        if "time" in ds.coords and ds.sizes.get("time", 0) > 0:
            print("\nTIME COVERAGE")
            print("-" * 100)
            print("start:", ds.time.values[0])
            print("end  :", ds.time.values[-1])

    finally:
        ds.close()

    # ------------------------------------------------------------
    # NetCDF storage details
    # ------------------------------------------------------------
    print("\nNETCDF STORAGE / GROUPS")
    print("-" * 100)

    with Dataset(path, "r") as nc:
        print("data_model:", nc.data_model)
        print("groups    :", list(nc.groups.keys()) or "(none)")

        print("\nVariable storage:")
        for name, var in nc.variables.items():
            try:
                filters = var.filters()
            except Exception:
                filters = None

            try:
                chunks = var.chunking()
            except Exception:
                chunks = None

            print(
                f"  {name:35s}"
                f" dtype={str(var.dtype):10s}"
                f" dims={var.dimensions}"
                f" chunks={chunks}"
                f" filters={filters}"
            )


def main():
    print("Inspecting ESSD publication NetCDF files")
    print("No files will be modified.")

    for path in FILES:
        inspect_file(path)


if __name__ == "__main__":
    main()
