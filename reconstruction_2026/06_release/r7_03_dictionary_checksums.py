"""Stage 7: write SHA256SUMS.txt and DATA_DICTIONARY.md for the corrected package, from the files themselves.

Same layout as tools/zenodo_inventory.py (used for the earlier staging); the licence texts are the decided ones
(CML products CC BY-NC 4.0; CHMI and DWD derived products CC BY 4.0) and a short note is added for the CML files.
Usage: python -I r7_03_dictionary_checksums.py <package dir>
"""
import hashlib, os, sys
import numpy as np
import netCDF4

LICENCE = [("CZ_CML_", "CC BY-NC 4.0"), ("CZ_RAIN_GAUGE_10MIN", "CC BY 4.0 (CHMI)"),
           ("CHMI_10MIN_", "CC BY 4.0 (CHMI)"), ("CZ_RADAR_CMLPATH_5MIN", "CC BY 4.0 (derived from DWD RADKLIM)")]
NOTE = {
    "CZ_CML_RAW_1MIN": "Raw transmitted (`tx`) and received (`rx`) signal levels of the 2,574 CMLs of the raw dataset.",
    "CZ_CML_PREPROCESSED_1MIN": ("Quality controlled signals of the 1,987 CMLs retained after quality control "
                                 "(3,731 sublinks with `qc_sublink_status` = 0). `qc_flag` bits: 1 invalid value, 2 mixed "
                                 "placeholder minute, 4 low RSL plateau, 8 total signal loss reconstructed during a "
                                 "blackout, 16 total signal loss interpolated. Parameters in `qc_parameters_json`."),
    "CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN": ("Native 1-minute CML rainfall of the 1,987 quality controlled CMLs "
                                                   "(Graf RSTD, constant baseline, Pastorek wet antenna correction, "
                                                   "ITU-R P.838-3 2005). `trsl_source_flag_1min` bit 8/16 marks rates "
                                                   "from reconstructed/interpolated total signal loss."),
}


def licence_for(fn):
    return next((v for k, v in LICENCE if fn.startswith(k)), "CC BY-NC 4.0")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(16 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def nc_vars(path):
    with netCDF4.Dataset(path) as ds:
        dims = {k: len(v) for k, v in ds.dimensions.items()}
        rows = []
        for name, v in ds.variables.items():
            note = getattr(v, "long_name", "") or ""
            if "time" in name and v.ndim == 1:
                x = np.asarray(v[:])
                try:
                    c = [str(netCDF4.num2date(float(x[i]), v.units)) for i in (0, -1)]
                    note = (note + f" covers {c[0]} .. {c[1]} ({x.size:,} steps)").strip()
                except Exception:
                    pass
            rows.append(f"| `{name}` | {v.dtype} | {','.join(v.dimensions)} | {getattr(v, 'units', '') or ''} | {note} |")
    return dims, rows


root = sys.argv[1]
sums, md = [], ["# Data dictionary — ESSD Czech CML 2021 dataset", "",
                "Generated directly from the files in the archive. Checksums are also in `SHA256SUMS.txt`.", ""]
for sub in ("primary", "supplementary_wvd"):
    for fn in sorted(os.listdir(os.path.join(root, sub))):
        p = os.path.join(root, sub, fn)
        h = sha256(p); sums.append(f"{h}  {sub}/{fn}")
        md += [f"## `{sub}/{fn}`", f"*Size:* {os.path.getsize(p):,} bytes · *SHA256:* `{h}` · *Licence:* {licence_for(fn)}"]
        for k, v in NOTE.items():
            if fn.startswith(k + "_"):
                md += ["", v]
        if fn.endswith(".nc"):
            dims, rows = nc_vars(p)
            md += ["", "*Dimensions:* " + ", ".join(f"`{k}`={v:,}" for k, v in dims.items()), "",
                   "| variable | dtype | dims | units | long_name / note |", "|---|---|---|---|---|"] + rows
        md.append("")
with open(os.path.join(root, "DATA_DICTIONARY.md"), "w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(md) + "\n")
for fn in ("README_DATA.md", "DATA_DICTIONARY.md"):
    sums.append(f"{sha256(os.path.join(root, fn))}  {fn}")
with open(os.path.join(root, "SHA256SUMS.txt"), "w", newline="\n") as f:
    f.write("\n".join(sums) + "\n")
print("\n".join(sums))
