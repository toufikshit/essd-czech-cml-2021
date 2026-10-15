"""Generate the Zenodo staging inventory: SHA256SUMS, file+variable inventory JSON,
and a Markdown data dictionary extracted from the actual cropped files.

Usage: python -I zenodo_inventory.py <staging_root> <out_prefix>
  staging_root: directory containing primary/ and supplementary_wvd/
Writes <staging_root>/SHA256SUMS.txt, <out_prefix>_inventory.json, and
<staging_root>/DATA_DICTIONARY.md.
"""
import hashlib
import json
import os
import sys
import numpy as np
import netCDF4

LICENCE = {
    "CZ_CML_RAW_1MIN": "CC BY-NC 4.0 (intended; pending rights sign-off)",
    "CZ_CML_PREPROCESSED_1MIN": "CC BY-NC 4.0 (intended; pending rights sign-off)",
    "CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN": "CC BY-NC 4.0 (intended; pending rights sign-off)",
    "CZ_RAIN_GAUGE_10MIN": "CC BY 4.0 (CHMI; pending attribution sign-off)",
    "CHMI_10MIN_TEMPERATURE": "CC BY 4.0 (CHMI; pending attribution sign-off)",
    "CHMI_10MIN_HUMIDITY": "CC BY 4.0 (CHMI; pending attribution sign-off)",
    "CZ_RADAR_CMLPATH_5MIN": "CC BY 4.0 (derived from DWD RADKLIM; pending attribution sign-off)",
    "CZ_CML_WVD_EBAND_PROGRESSIVE_V6": "CC BY-NC 4.0 (intended; pending rights sign-off)",
    "CZ_CML_WVD_V6_LINK_QUANTIZATION_SERIES": "CC BY-NC 4.0 (intended; pending rights sign-off)",
}


def licence_for(fn):
    for k, v in LICENCE.items():
        if fn.startswith(k):
            return v
    return "TBD"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(16 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def nc_info(path):
    info = {"dimensions": {}, "variables": [], "global_attributes": {}}
    with netCDF4.Dataset(path) as ds:
        info["dimensions"] = {k: len(v) for k, v in ds.dimensions.items()}
        info["global_attributes"] = {k: str(ds.getncattr(k)) for k in ds.ncattrs()}
        for name, v in ds.variables.items():
            e = {"name": name, "dtype": str(v.dtype), "dims": list(v.dimensions),
                 "units": getattr(v, "units", None), "long_name": getattr(v, "long_name", None),
                 "fill_value": None if not hasattr(v, "_FillValue") else str(v._FillValue)}
            if "time" in name and v.ndim == 1:
                x = v[:]
                try:
                    e["coverage"] = [str(netCDF4.num2date(float(np.asarray(x)[i]), v.units)) for i in (0, -1)]
                    e["n"] = int(np.asarray(x).size)
                except Exception:
                    pass
            info["variables"].append(e)
    return info


def main():
    root, prefix = sys.argv[1], sys.argv[2]
    inventory = []
    sums = []
    for sub in ("primary", "supplementary_wvd"):
        base = os.path.join(root, sub)
        if not os.path.isdir(base):
            continue
        for fn in sorted(os.listdir(base)):
            p = os.path.join(base, fn)
            if not os.path.isfile(p):
                continue
            entry = {"group": sub, "file": fn, "bytes": os.path.getsize(p),
                     "sha256": sha256(p), "licence": licence_for(fn)}
            if fn.endswith(".nc"):
                entry["netcdf"] = nc_info(p)
            inventory.append(entry)
            sums.append(f"{entry['sha256']}  {sub}/{fn}")
    with open(os.path.join(root, "SHA256SUMS.txt"), "w", newline="\n") as f:
        f.write("\n".join(sums) + "\n")
    with open(prefix + "_inventory.json", "w", encoding="utf-8") as f:
        json.dump(inventory, f, indent=2)

    md = ["# Data dictionary — ESSD Czech CML 2021 Zenodo staging",
          "",
          "Generated directly from the staged files. Sizes, checksums and licences per",
          "file are in `SHA256SUMS.txt` and the staging README; licences are *intended*",
          "terms pending written sign-off.", ""]
    for e in inventory:
        md.append(f"## `{e['group']}/{e['file']}`")
        md.append(f"*Size:* {e['bytes']:,} bytes · *SHA256:* `{e['sha256']}` · *Licence:* {e['licence']}")
        if "netcdf" in e:
            n = e["netcdf"]
            md.append(f"*Dimensions:* " + ", ".join(f"`{k}`={v:,}" for k, v in n["dimensions"].items()))
            md.append("")
            md.append("| variable | dtype | dims | units | long_name / note |")
            md.append("|---|---|---|---|---|")
            for v in n["variables"]:
                note = v.get("long_name") or ""
                if v.get("coverage"):
                    note = (note + f" covers {v['coverage'][0]} .. {v['coverage'][1]} ({v['n']:,} steps)").strip()
                md.append(f"| `{v['name']}` | {v['dtype']} | {','.join(v['dims'])} | "
                          f"{v.get('units') or ''} | {note} |")
        md.append("")
    with open(os.path.join(root, "DATA_DICTIONARY.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(md) + "\n")
    print(f"inventory: {len(inventory)} files")


if __name__ == "__main__":
    main()
