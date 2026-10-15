#!/usr/bin/env python3
"""
cml_qc_filter_standalone_v4_3.py
================================
Create a minimal, publication-ready CML NetCDF from an existing QC dataset.

The scientific QC is NOT recomputed. The script applies the existing flags:

    qc_final_link_pass
    qc_final_sublink_pass
    qc_final_valid

Output layout
-------------
The original rectangular structure is preserved:

    time × sublink_id × cml_id

Only these data variables are written:

    tx(time, sublink_id, cml_id)
    rx(time, sublink_id, cml_id)
    trsl_qc(time, sublink_id, cml_id)

Only these auxiliary coordinates are written:

    frequency(sublink_id, cml_id)
    polarization(sublink_id, cml_id)
    length(cml_id)

The dimension coordinates time, sublink_id and cml_id are also retained.
All other source variables are dropped.

Selection rules
---------------
1. Physical links with qc_final_link_pass != 1 are physically removed.
2. A physical link is retained when at least one selected directional sublink
   remains, unless --require-all-sublinks is used.
3. Failed directional sublinks inside retained physical links remain as the
   shared sublink_id slot, but tx, rx, trsl_qc, frequency and polarization are
   missing for that direction. This preserves the original rectangular layout.
4. At sample level, tx, rx and trsl_qc are finite only where qc_final_valid == 1
   and the directional sublink is selected. Invalid samples are written as
   missing values.
5. Optional minimum valid-count/fraction criteria can further restrict already
   passing sublinks. Use --min-valid-fraction 0 --min-valid-observations 1 to
   apply the existing QC without a new global availability threshold.

The output signal variables carry a NetCDF ``coordinates`` attribute listing
frequency, polarization and length. Xarray therefore opens these variables as
auxiliary coordinates rather than ordinary data variables.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np

try:
    from netCDF4 import Dataset, chartostring
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "netCDF4 is required. Install it with `python -m pip install netCDF4`. "
        f"Import error: {exc}"
    ) from exc

SCRIPT_VERSION = "4.3"


def _safe_attrs(var) -> Dict[str, object]:
    return {
        name: var.getncattr(name)
        for name in var.ncattrs()
        if name not in {"_FillValue", "coordinates"}
    }


def _fill_value(var):
    return var.getncattr("_FillValue") if "_FillValue" in var.ncattrs() else None


def _masked_to_array(data) -> np.ndarray:
    if np.ma.isMaskedArray(data):
        kind = np.asarray(data.data).dtype.kind
        if kind == "S":
            return np.asarray(data.filled(b""))
        if kind == "U":
            return np.asarray(data.filled(""))
        return np.asarray(data.filled(np.nan))
    return np.asarray(data)


def _guess_name(
    available: Iterable[str],
    requested: Optional[str],
    candidates: Sequence[str],
    role: str,
) -> str:
    names = list(available)
    if requested:
        if requested not in names:
            raise KeyError(
                f"Requested {role} variable {requested!r} not found. Available: {names}"
            )
        return requested
    lower = {name.lower(): name for name in names}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    raise KeyError(f"Could not auto-detect {role}. Available: {names}")


def _guess_dim(
    dimensions: Sequence[str],
    requested: Optional[str],
    candidates: Sequence[str],
    role: str,
) -> str:
    dims = list(dimensions)
    if requested:
        if requested not in dims:
            raise KeyError(f"Requested {role} dimension {requested!r} not found: {dims}")
        return requested
    lower = {name.lower(): name for name in dims}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    raise KeyError(f"Could not auto-detect {role} dimension. Available: {dims}")


def _read_link_flag(var, cml_dim: str) -> np.ndarray:
    if tuple(var.dimensions) != (cml_dim,):
        raise ValueError(
            f"{var.name} must have dimensions ({cml_dim},); found {var.dimensions}"
        )
    return np.asarray(_masked_to_array(var[:]) == 1, dtype=bool)


def _read_sublink_flag(var, sublink_dim: str, cml_dim: str) -> np.ndarray:
    dims = list(var.dimensions)
    if len(dims) != 2 or sublink_dim not in dims or cml_dim not in dims:
        raise ValueError(
            f"{var.name} must contain exactly ({sublink_dim}, {cml_dim}); "
            f"found {var.dimensions}"
        )
    arr = _masked_to_array(var[:])
    arr = np.transpose(arr, (dims.index(sublink_dim), dims.index(cml_dim)))
    return np.asarray(arr == 1, dtype=bool)


def _read_valid_canonical(
    var,
    time_dim: str,
    sublink_dim: str,
    cml_dim: str,
    t0: int,
    t1: int,
    source_cml_indices: np.ndarray,
) -> np.ndarray:
    """Read final validity as (time, sublink, selected-cml)."""
    dims = list(var.dimensions)
    expected = {time_dim, sublink_dim, cml_dim}
    if len(dims) != 3 or set(dims) != expected:
        raise ValueError(
            f"{var.name} must contain exactly {sorted(expected)}; found {var.dimensions}"
        )
    sl = [slice(None)] * 3
    sl[dims.index(time_dim)] = slice(t0, t1)
    sl[dims.index(cml_dim)] = source_cml_indices
    arr = _masked_to_array(var[tuple(sl)])
    arr = np.transpose(
        arr,
        (dims.index(time_dim), dims.index(sublink_dim), dims.index(cml_dim)),
    )
    return np.asarray(arr == 1, dtype=bool)


def _read_signal_canonical(
    var,
    time_dim: str,
    sublink_dim: str,
    cml_dim: str,
    t0: int,
    t1: int,
    source_cml_indices: np.ndarray,
) -> np.ndarray:
    """Read a signal variable as (time, sublink, selected-cml)."""
    dims = list(var.dimensions)
    expected = {time_dim, sublink_dim, cml_dim}
    if len(dims) != 3 or set(dims) != expected:
        raise ValueError(
            f"{var.name} must contain exactly {sorted(expected)}; found {var.dimensions}"
        )
    sl = [slice(None)] * 3
    sl[dims.index(time_dim)] = slice(t0, t1)
    sl[dims.index(cml_dim)] = source_cml_indices
    arr = _masked_to_array(var[tuple(sl)])
    return np.transpose(
        arr,
        (dims.index(time_dim), dims.index(sublink_dim), dims.index(cml_dim)),
    )


def _decode_fixed_strings(
    arr: np.ndarray,
    dims: Tuple[str, ...],
    sublink_dim: str,
    cml_dim: str,
) -> Tuple[np.ndarray, Tuple[str, ...]]:
    """Collapse a trailing S1/U1 character dimension into full strings."""
    if (
        arr.ndim >= 2
        and len(dims) == arr.ndim
        and dims[-1] not in {sublink_dim, cml_dim}
        and arr.dtype.kind in "SU"
        and arr.dtype.itemsize in (1, 4)
    ):
        if arr.dtype.kind == "S":
            try:
                decoded = np.asarray(chartostring(arr), dtype=object)
            except Exception:
                rows = arr.reshape(-1, arr.shape[-1])
                decoded = np.asarray(
                    [
                        b"".join(np.asarray(row, dtype="S1").tolist())
                        .decode("utf-8", errors="ignore")
                        .rstrip("\x00")
                        for row in rows
                    ],
                    dtype=object,
                ).reshape(arr.shape[:-1])
        else:
            rows = arr.reshape(-1, arr.shape[-1])
            decoded = np.asarray(
                ["".join(str(x) for x in row).rstrip("\x00") for row in rows],
                dtype=object,
            ).reshape(arr.shape[:-1])
        return decoded, tuple(dims[:-1])
    return arr, dims


def _read_sublink_coordinate(
    var,
    sublink_dim: str,
    cml_dim: str,
    n_sublink: int,
    n_cml: int,
) -> np.ndarray:
    """Return frequency/polarization as canonical (sublink, cml)."""
    arr = _masked_to_array(var[:])
    dims = tuple(var.dimensions)
    arr, dims = _decode_fixed_strings(arr, dims, sublink_dim, cml_dim)

    if dims == (sublink_dim, cml_dim):
        out = arr
    elif dims == (cml_dim, sublink_dim):
        out = arr.T
    elif dims == (cml_dim,):
        out = np.broadcast_to(arr[None, :], (n_sublink, n_cml))
    elif dims == (sublink_dim,):
        out = np.broadcast_to(arr[:, None], (n_sublink, n_cml))
    elif arr.ndim == 0:
        out = np.full((n_sublink, n_cml), arr.item(), dtype=object)
    else:
        raise ValueError(
            f"Coordinate {var.name} has unsupported dimensions {dims}; expected "
            f"({sublink_dim},{cml_dim}), ({cml_dim},{sublink_dim}), "
            f"({cml_dim},) or ({sublink_dim},)."
        )
    return np.asarray(out)


def _read_link_coordinate(
    var,
    sublink_dim: str,
    cml_dim: str,
    n_sublink: int,
    n_cml: int,
) -> np.ndarray:
    """Return link length as canonical (cml,)."""
    arr = _masked_to_array(var[:])
    dims = tuple(var.dimensions)
    arr, dims = _decode_fixed_strings(arr, dims, sublink_dim, cml_dim)

    if dims == (cml_dim,):
        return np.asarray(arr)
    if dims == (sublink_dim, cml_dim):
        values = np.asarray(arr, dtype=np.float64)
    elif dims == (cml_dim, sublink_dim):
        values = np.asarray(arr.T, dtype=np.float64)
    elif arr.ndim == 0:
        return np.full(n_cml, arr.item())
    else:
        raise ValueError(
            f"Coordinate {var.name} has unsupported dimensions {dims}; expected "
            f"({cml_dim},) or a sublink-by-link array with identical path lengths."
        )

    # Link length should be identical for both directions. Use the first finite
    # value after verifying consistency where two finite values exist.
    out = np.full(n_cml, np.nan, dtype=np.float64)
    for c in range(n_cml):
        finite = values[:, c][np.isfinite(values[:, c])]
        if finite.size == 0:
            continue
        if finite.size > 1 and np.nanmax(finite) - np.nanmin(finite) > 1.0e-6:
            raise ValueError(
                f"{var.name} differs between sublinks for source cml index {c}: "
                f"{finite.tolist()}"
            )
        out[c] = finite[0]
    return out


def _numeric_fill(var):
    fill = _fill_value(var)
    if fill is not None:
        return fill
    dtype = np.dtype(var.dtype)
    if dtype.kind == "f":
        return np.nan
    if dtype.kind == "i":
        return np.iinfo(dtype).min
    if dtype.kind == "u":
        return np.iinfo(dtype).max
    raise TypeError(f"Variable {var.name} is not numeric")


def _create_numeric_var(
    dst: Dataset,
    name: str,
    dtype,
    dims: Tuple[str, ...],
    attrs: Dict[str, object],
    *,
    fill_value,
    compression_level: int,
    chunksizes: Optional[Tuple[int, ...]] = None,
):
    kwargs = {"fill_value": fill_value}
    if compression_level > 0:
        kwargs.update(zlib=True, complevel=int(compression_level), shuffle=True)
    if chunksizes is not None:
        kwargs["chunksizes"] = chunksizes
    var = dst.createVariable(name, dtype, dims, **kwargs)
    var.setncatts(attrs)
    return var


def _create_string_var(
    dst: Dataset,
    name: str,
    dims: Tuple[str, ...],
    attrs: Dict[str, object],
):
    var = dst.createVariable(name, str, dims)
    var.setncatts(attrs)
    return var


def _write_vlen_strings(var, data: np.ndarray) -> None:
    arr = np.asarray(data, dtype=object)
    try:
        var[:] = arr
        return
    except Exception:
        pass
    for idx in np.ndindex(arr.shape):
        value = arr[idx]
        if isinstance(value, (bytes, np.bytes_)):
            value = value.decode("utf-8", errors="ignore").rstrip("\x00")
        else:
            value = str(value)
        var[idx] = value


def _copy_dimension_coordinate(
    src: Dataset,
    dst: Dataset,
    dim_name: str,
    retained_indices: Optional[np.ndarray] = None,
) -> None:
    """Copy a one-dimensional dimension coordinate if it exists."""
    if dim_name not in src.variables:
        # NetCDF dimensions may exist without an explicit coordinate variable.
        return
    source = src.variables[dim_name]
    if tuple(source.dimensions) != (dim_name,):
        raise ValueError(
            f"Dimension coordinate {dim_name} must have dimensions ({dim_name},); "
            f"found {source.dimensions}"
        )
    data = _masked_to_array(source[:])
    if retained_indices is not None:
        data = data[retained_indices]

    attrs = _safe_attrs(source)
    if np.asarray(data).dtype.kind in "OSU" or source.datatype is str:
        out = _create_string_var(dst, dim_name, (dim_name,), attrs)
        _write_vlen_strings(out, np.asarray(data, dtype=object))
    else:
        fill = _fill_value(source)
        kwargs = {"fill_value": fill} if fill is not None else {}
        out = dst.createVariable(dim_name, source.datatype, (dim_name,), **kwargs)
        out.setncatts(attrs)
        out[:] = data


def _write_summary_csv(path: Path, summary: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value"])
        for key, value in summary.items():
            writer.writerow([key, value])


def process(args: argparse.Namespace) -> Dict[str, object]:
    started = time.perf_counter()
    src_path = Path(args.qc_nc).expanduser().resolve()
    out_path = Path(args.out_nc).expanduser().resolve()
    temp_path = out_path.with_name(out_path.name + ".partial")

    if not src_path.exists():
        raise FileNotFoundError(src_path)
    if src_path == out_path:
        raise ValueError("Input and output paths must differ")
    if out_path.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists: {out_path}. Add --overwrite to replace it.")
    if temp_path.exists():
        temp_path.unlink()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with Dataset(src_path, "r") as src:
        src.set_auto_mask(True)
        dims_all = list(src.dimensions)
        time_dim = _guess_dim(dims_all, args.time_dim, ["time"], "time")
        sublink_dim = _guess_dim(
            dims_all, args.sublink_dim, ["sublink_id", "sublink"], "sublink"
        )
        cml_dim = _guess_dim(dims_all, args.cml_dim, ["cml_id", "cml"], "CML")

        link_name = _guess_name(
            src.variables, args.link_pass_var, ["qc_final_link_pass"], "link-pass"
        )
        sublink_name = _guess_name(
            src.variables,
            args.sublink_pass_var,
            ["qc_final_sublink_pass"],
            "sublink-pass",
        )
        valid_name = _guess_name(
            src.variables, args.valid_var, ["qc_final_valid"], "sample-valid"
        )
        tx_name = _guess_name(src.variables, args.tx_var, ["tx"], "TX")
        rx_name = _guess_name(src.variables, args.rx_var, ["rx"], "RX")
        trsl_name = _guess_name(
            src.variables, args.trsl_var, ["trsl_qc"], "quality-controlled TRSL"
        )
        frequency_name = _guess_name(
            src.variables, args.frequency_var, ["frequency"], "frequency"
        )
        polarization_name = _guess_name(
            src.variables, args.polarization_var, ["polarization"], "polarization"
        )
        length_name = _guess_name(src.variables, args.length_var, ["length"], "length")

        n_time = len(src.dimensions[time_dim])
        n_sublink = len(src.dimensions[sublink_dim])
        n_cml = len(src.dimensions[cml_dim])

        link_pass = _read_link_flag(src.variables[link_name], cml_dim)
        sublink_pass = _read_sublink_flag(
            src.variables[sublink_name], sublink_dim, cml_dim
        )
        if sublink_pass.shape != (n_sublink, n_cml):
            raise ValueError(
                f"Unexpected sublink flag shape {sublink_pass.shape}; "
                f"expected {(n_sublink, n_cml)}"
            )

        valid_var = src.variables[valid_name]
        print("Scanning sample-level validity ...", flush=True)
        valid_counts = np.zeros((n_sublink, n_cml), dtype=np.int64)
        all_cml_indices = np.arange(n_cml, dtype=np.int64)
        for t0 in range(0, n_time, args.time_chunk):
            t1 = min(t0 + args.time_chunk, n_time)
            valid = _read_valid_canonical(
                valid_var,
                time_dim,
                sublink_dim,
                cml_dim,
                t0,
                t1,
                all_cml_indices,
            )
            valid_counts += valid.sum(axis=0, dtype=np.int64)
            print(f"  validity {t1:,}/{n_time:,}", flush=True)

        valid_fraction = valid_counts / float(n_time)
        extra_pass = (
            (valid_counts >= int(args.min_valid_observations))
            & (valid_fraction >= float(args.min_valid_fraction))
        )
        selected_sublink = sublink_pass & extra_pass & link_pass[None, :]

        if args.require_all_sublinks:
            retained_link_mask = link_pass & selected_sublink.all(axis=0)
        else:
            retained_link_mask = link_pass & selected_sublink.any(axis=0)

        retained_indices = np.flatnonzero(retained_link_mask).astype(np.int64)
        if retained_indices.size == 0:
            raise RuntimeError(
                "No links remain. Check QC flags and optional valid-count/fraction thresholds."
            )

        selected_out = selected_sublink[:, retained_indices]
        valid_counts_out = valid_counts[:, retained_indices]
        retained_links = int(retained_indices.size)
        retained_sublinks = int(selected_out.sum())
        failed_slots = int((~selected_out).sum())
        valid_observations = int((valid_counts_out * selected_out).sum())
        possible_cells = int(n_time * n_sublink * retained_links)
        missing_cells = int(possible_cells - valid_observations)

        # Read compact coordinates before creating the destination.
        frequency = _read_sublink_coordinate(
            src.variables[frequency_name], sublink_dim, cml_dim, n_sublink, n_cml
        )[:, retained_indices]
        polarization = _read_sublink_coordinate(
            src.variables[polarization_name], sublink_dim, cml_dim, n_sublink, n_cml
        )[:, retained_indices]
        length = _read_link_coordinate(
            src.variables[length_name], sublink_dim, cml_dim, n_sublink, n_cml
        )[retained_indices]

        # Failed directional slots carry no scientific metadata.
        if frequency.dtype.kind in "fc":
            frequency = np.where(selected_out, frequency, np.nan)
        elif frequency.dtype.kind in "iu":
            frequency = frequency.astype(np.float64)
            frequency = np.where(selected_out, frequency, np.nan)
        else:
            frequency = np.asarray(frequency, dtype=object)
            frequency[~selected_out] = ""

        polarization = np.asarray(polarization, dtype=object)
        polarization[~selected_out] = ""

        summary: Dict[str, object] = {
            "script": Path(__file__).name,
            "script_version": SCRIPT_VERSION,
            "source_qc_file": str(src_path),
            "output_file": str(out_path),
            "output_dimensions": f"{time_dim},{sublink_dim},{cml_dim}",
            "output_data_variables": f"{tx_name},{rx_name},{trsl_name}",
            "output_auxiliary_coordinates": (
                f"{frequency_name},{polarization_name},{length_name}"
            ),
            "original_links": int(n_cml),
            "retained_links": retained_links,
            "dropped_links": int(n_cml - retained_links),
            "original_directional_sublinks": int(n_sublink * n_cml),
            "retained_directional_sublinks": retained_sublinks,
            "dropped_directional_sublinks_total": int(
                n_sublink * n_cml - retained_sublinks
            ),
            "failed_directional_slots_inside_retained_links": failed_slots,
            "valid_observations_retained": valid_observations,
            "invalid_or_failed_output_cells_written_missing": missing_cells,
            "minimum_valid_fraction": float(args.min_valid_fraction),
            "minimum_valid_observations": int(args.min_valid_observations),
            "require_all_sublinks": bool(args.require_all_sublinks),
        }

        print("=" * 92)
        print("CML QC minimal same-structure output")
        print(f"Source QC file       : {src_path}")
        print(f"Output               : {out_path}")
        print(f"Output structure      : {time_dim} × {sublink_dim} × {cml_dim}")
        print(f"Data variables        : {tx_name}, {rx_name}, {trsl_name}")
        print(
            f"Auxiliary coordinates : {frequency_name}, {polarization_name}, {length_name}"
        )
        print(f"Original links        : {n_cml:,}")
        print(f"Retained links        : {retained_links:,}")
        print(f"Dropped links         : {n_cml - retained_links:,}")
        print(f"Retained sublinks     : {retained_sublinks:,}")
        print(f"Failed slots retained : {failed_slots:,} (all scientific values missing)")
        print(f"Valid observations    : {valid_observations:,}")
        print(f"Missing output cells  : {missing_cells:,}")
        print("=" * 92)

        with Dataset(temp_path, "w", format="NETCDF4") as dst:
            dst.createDimension(time_dim, len(src.dimensions[time_dim]))
            dst.createDimension(sublink_dim, n_sublink)
            dst.createDimension(cml_dim, retained_links)

            # Copy source global metadata, then document the compact product.
            dst.setncatts({name: src.getncattr(name) for name in src.ncattrs()})
            dst.setncattr("strict_qc_compaction_program", Path(__file__).name)
            dst.setncattr("strict_qc_compaction_version", SCRIPT_VERSION)
            dst.setncattr("strict_qc_source_file", str(src_path))
            dst.setncattr(
                "strict_qc_content",
                "Minimal QC product containing only tx, rx and trsl_qc, plus time, sublink_id, cml_id, frequency, polarization and length coordinates.",
            )
            dst.setncattr(
                "strict_qc_link_rule",
                "Physical links failing qc_final_link_pass are removed. A retained link has at least one selected directional sublink."
                if not args.require_all_sublinks
                else "Physical links are retained only when every directional sublink is selected.",
            )
            dst.setncattr(
                "strict_qc_sublink_rule",
                "A failed directional sublink inside a retained link remains in the shared sublink_id slot, but tx, rx, trsl_qc, frequency and polarization are missing for that direction.",
            )
            dst.setncattr(
                "strict_qc_observation_rule",
                "tx, rx and trsl_qc are retained only where qc_final_valid equals 1 and the directional sublink is selected; all other samples are missing.",
            )
            dst.setncattr("strict_qc_summary_json", json.dumps(summary, sort_keys=True))
            dst.setncattr(
                "strict_qc_created_utc",
                time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            )

            _copy_dimension_coordinate(src, dst, time_dim)
            _copy_dimension_coordinate(src, dst, sublink_dim)
            _copy_dimension_coordinate(src, dst, cml_dim, retained_indices)

            # Auxiliary coordinate: frequency(sublink_id, cml_id)
            freq_src = src.variables[frequency_name]
            freq_attrs = _safe_attrs(freq_src)
            if frequency.dtype.kind in "fc":
                freq_dtype = "f4" if np.dtype(freq_src.dtype).itemsize <= 4 else "f8"
                freq_var = _create_numeric_var(
                    dst,
                    frequency_name,
                    freq_dtype,
                    (sublink_dim, cml_dim),
                    freq_attrs,
                    fill_value=np.nan,
                    compression_level=args.compression_level,
                    chunksizes=(n_sublink, min(retained_links, args.link_chunk)),
                )
                freq_var[:, :] = frequency
            else:
                freq_var = _create_string_var(
                    dst, frequency_name, (sublink_dim, cml_dim), freq_attrs
                )
                _write_vlen_strings(freq_var, frequency)

            # Auxiliary coordinate: polarization(sublink_id, cml_id)
            pol_var = _create_string_var(
                dst,
                polarization_name,
                (sublink_dim, cml_dim),
                _safe_attrs(src.variables[polarization_name]),
            )
            _write_vlen_strings(pol_var, polarization)

            # Auxiliary coordinate: length(cml_id)
            length_src = src.variables[length_name]
            length_dtype = "f4" if np.dtype(length_src.dtype).itemsize <= 4 else "f8"
            length_var = _create_numeric_var(
                dst,
                length_name,
                length_dtype,
                (cml_dim,),
                _safe_attrs(length_src),
                fill_value=np.nan,
                compression_level=args.compression_level,
                chunksizes=(min(retained_links, args.link_chunk),),
            )
            length_var[:] = length

            coordinates_attr = f"{frequency_name} {polarization_name} {length_name}"
            signal_names = [tx_name, rx_name, trsl_name]
            signal_outputs = {}
            for name in signal_names:
                source = src.variables[name]
                attrs = _safe_attrs(source)
                attrs["coordinates"] = coordinates_attr
                attrs["strict_qc_mask"] = (
                    f"Finite only where {valid_name}=1 and {sublink_name}=1 for retained links."
                )
                chunks = (
                    min(n_time, args.time_chunk),
                    n_sublink,
                    min(retained_links, args.link_chunk),
                )
                signal_outputs[name] = _create_numeric_var(
                    dst,
                    name,
                    source.datatype,
                    (time_dim, sublink_dim, cml_dim),
                    attrs,
                    fill_value=_numeric_fill(source),
                    compression_level=args.compression_level,
                    chunksizes=chunks,
                )

            # Stream only the three signal variables. The QC flags are read from
            # the source but are intentionally not copied to the minimal output.
            for oc0 in range(0, retained_links, args.link_batch):
                oc1 = min(oc0 + args.link_batch, retained_links)
                source_c = retained_indices[oc0:oc1]
                selected_batch = selected_out[:, oc0:oc1]
                for t0 in range(0, n_time, args.time_chunk):
                    t1 = min(t0 + args.time_chunk, n_time)
                    valid = _read_valid_canonical(
                        valid_var,
                        time_dim,
                        sublink_dim,
                        cml_dim,
                        t0,
                        t1,
                        source_c,
                    )
                    keep = valid & selected_batch[None, :, :]

                    output_slice = (
                        slice(t0, t1),
                        slice(None),
                        slice(oc0, oc1),
                    )
                    for name in signal_names:
                        source = src.variables[name]
                        data = _read_signal_canonical(
                            source,
                            time_dim,
                            sublink_dim,
                            cml_dim,
                            t0,
                            t1,
                            source_c,
                        )
                        missing = _numeric_fill(source)
                        data = np.where(keep, data, missing)
                        signal_outputs[name][output_slice] = data

                print(
                    f"Written retained links {oc0:5d}:{oc1:5d}/{retained_links}",
                    flush=True,
                )

            dst.sync()

    if out_path.exists():
        out_path.unlink()
    os.replace(temp_path, out_path)

    elapsed = float(time.perf_counter() - started)
    summary["elapsed_seconds"] = elapsed
    summary["output_size_gib"] = float(out_path.stat().st_size / (1024**3))

    json_path = (
        Path(args.summary_json).expanduser().resolve()
        if args.summary_json
        else out_path.with_suffix(".summary.json")
    )
    csv_path = (
        Path(args.summary_csv).expanduser().resolve()
        if args.summary_csv
        else out_path.with_suffix(".summary.csv")
    )
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    _write_summary_csv(csv_path, summary)

    print("=" * 92)
    print(f"Saved minimal QC dataset: {out_path}")
    print(f"Saved summary JSON      : {json_path}")
    print(f"Saved summary CSV       : {csv_path}")
    print(f"Elapsed                 : {elapsed:.1f} s")
    print(f"Output size             : {summary['output_size_gib']:.3f} GiB")
    print("=" * 92)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create a minimal same-structure QC NetCDF containing only tx, rx, "
            "trsl_qc and the frequency, polarization and length coordinates."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    parser.add_argument("--qc-nc", required=True)
    parser.add_argument("--out-nc", required=True)
    parser.add_argument("--summary-json", default=None)
    parser.add_argument("--summary-csv", default=None)

    parser.add_argument("--link-pass-var", default="qc_final_link_pass")
    parser.add_argument("--sublink-pass-var", default="qc_final_sublink_pass")
    parser.add_argument("--valid-var", default="qc_final_valid")
    parser.add_argument("--tx-var", default="tx")
    parser.add_argument("--rx-var", default="rx")
    parser.add_argument("--trsl-var", default="trsl_qc")
    parser.add_argument("--frequency-var", default="frequency")
    parser.add_argument("--polarization-var", default="polarization")
    parser.add_argument("--length-var", default="length")

    parser.add_argument("--time-dim", default="time")
    parser.add_argument("--sublink-dim", default="sublink_id")
    parser.add_argument("--cml-dim", default="cml_id")

    parser.add_argument("--min-valid-fraction", type=float, default=0.0)
    parser.add_argument("--min-valid-observations", type=int, default=1)
    parser.add_argument(
        "--require-all-sublinks",
        action="store_true",
        help="Retain only physical links for which every directional sublink passes.",
    )

    parser.add_argument("--time-chunk", type=int, default=10080)
    parser.add_argument("--link-batch", type=int, default=16)
    parser.add_argument("--link-chunk", type=int, default=16)
    parser.add_argument("--compression-level", type=int, choices=range(0, 10), default=1)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0.0 <= args.min_valid_fraction <= 1.0:
        raise SystemExit("--min-valid-fraction must be between 0 and 1")
    if args.min_valid_observations < 0:
        raise SystemExit("--min-valid-observations must be >= 0")
    if args.time_chunk < 1 or args.link_batch < 1 or args.link_chunk < 1:
        raise SystemExit("Chunk sizes must be >= 1")
    print(f"Running {Path(__file__).name} version {SCRIPT_VERSION}")
    try:
        process(args)
    except KeyboardInterrupt:
        print("Interrupted. Any .partial output can be removed safely.", file=sys.stderr)
        raise SystemExit(130)


if __name__ == "__main__":
    main()
