from pathlib import Path

import numpy as np
import xarray as xr


FILES = {
    "Raw CML dataset": Path(
        "/home/jovyan/ESSD_publication_dataset/"
        "CZ_CML_RAW_1MIN_20200101_20211231.nc"
    ),
    "Quality controlled CML dataset": Path(
        "/home/jovyan/ESSD_publication_dataset/"
        "CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc"
    ),
}

E_BAND_THRESHOLD_GHZ = 70.0


def calculate_e_band_fraction(path: Path) -> dict:
    with xr.open_dataset(path) as ds:
        if "frequency" not in ds:
            raise KeyError(f"{path}: variable 'frequency' was not found.")

        frequency = ds["frequency"].astype("float64")

        # Convert MHz to GHz when necessary.
        finite_values = frequency.values[np.isfinite(frequency.values)]
        if finite_values.size == 0:
            raise ValueError(f"{path}: no finite frequency values.")

        if np.nanmedian(finite_values) > 1000:
            frequency_ghz = frequency / 1000.0
        else:
            frequency_ghz = frequency

        if "sublink_id" in frequency_ghz.dims:
            # A physical CML is classified as E band when at least one
            # directional sublink operates above 70 GHz.
            finite_per_cml = frequency_ghz.notnull().any(dim="sublink_id")
            e_band_per_cml = (
                frequency_ghz > E_BAND_THRESHOLD_GHZ
            ).any(dim="sublink_id")
        else:
            finite_per_cml = frequency_ghz.notnull()
            e_band_per_cml = frequency_ghz > E_BAND_THRESHOLD_GHZ

        total_cmls = int(ds.sizes["cml_id"])
        cmls_with_frequency = int(finite_per_cml.sum().item())
        e_band_cmls = int((e_band_per_cml & finite_per_cml).sum().item())

        percentage_all = 100.0 * e_band_cmls / total_cmls
        percentage_valid = (
            100.0 * e_band_cmls / cmls_with_frequency
            if cmls_with_frequency
            else np.nan
        )

        return {
            "total_cmls": total_cmls,
            "cmls_with_frequency": cmls_with_frequency,
            "e_band_cmls": e_band_cmls,
            "percentage_all": percentage_all,
            "percentage_valid": percentage_valid,
        }


for label, path in FILES.items():
    if not path.exists():
        print(f"\n{label}: file not found\n  {path}")
        continue

    result = calculate_e_band_fraction(path)

    print(f"\n{label}")
    print("=" * len(label))
    print(f"Total CMLs                    : {result['total_cmls']:,}")
    print(f"CMLs with frequency metadata : {result['cmls_with_frequency']:,}")
    print(f"E band CMLs above 70 GHz     : {result['e_band_cmls']:,}")
    print(f"E band percentage, all CMLs  : {result['percentage_all']:.2f}%")
    print(
        "E band percentage, CMLs with frequency metadata: "
        f"{result['percentage_valid']:.2f}%"
    )
