"""Build the GitHub release-candidate repository from the verified handoff tree.

Usage: python -I build_github_repo.py <handoff_repository_root> <dest_repo_root> <report_json>

Every file is byte-copied from the handoff `repository/` tree except:
  * config/active_profile.json           -> regenerated portable default (profile "user")
  * notebooks/06_WVD_Sensitivity_Length_Frequency.ipynb
                                         -> bytes of the newest HPC variant (_REFRESHED_V6)
  * notebooks/02_Rainfall_Validation.ipynb
                                         -> one fallback literal removed (documented patch)
Historical/non-release items are preserved under attic/ instead of being deleted.
Large executed outputs (NetCDF, npz, per-sensor CSV dumps) stay out of git; selected
small executed evidence is copied under docs/evidence/hpc_run_2026-10/ with a manifest.
The script verifies, after building, that every byte-copied file hashes equal to its
source, and writes a machine-readable build report.
"""
import hashlib
import json
import os
import shutil
import sys

BS = chr(92)


def wl(p):
    p = os.path.abspath(p)
    if os.name == "nt" and not p.startswith(BS * 2 + "?" + BS):
        p = BS * 2 + "?" + BS + p
    return p


def sha256(path):
    h = hashlib.sha256()
    with open(wl(path), "rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def copy_file(src, dst):
    os.makedirs(wl(os.path.dirname(dst)), exist_ok=True)
    shutil.copyfile(wl(src), wl(dst))


# ---------------------------------------------------------------- copy plan
BYTE_COPY_DIRS = ["src", "original", "tests", "docs"]
# tests/test_connected_chain.py is byte-copied by the walk, then patched in place.

ROOT_BYTE_COPY = [
    "run.py",
    "tools_verify_originals.py",
    "requirements.txt",
    "requirements-notebooks.txt",
    "requirements-test.txt",
    "README_CONNECTED_WORKFLOW.md",
    ".gitignore",  # replaced later by generated one, but keep original copy in attic
]

CONFIG_BYTE_COPY = [
    "config/hpc_source.example.json",
    "config/user_paths.example.json",
    "config/downloaded_release.example.json",
]

NOTEBOOKS_BYTE_COPY = [
    "notebooks/00_Explore_Data.ipynb",
    "notebooks/01_Graf_Rainfall.ipynb",
    "notebooks/02_Rainfall_Validation.ipynb",
    "notebooks/03_Grid_CML_Rainfall.ipynb",
    "notebooks/04_Spatial_Radar_Maps.ipynb",
    "notebooks/05_Eband_Water_Vapour.ipynb",
    # 06 canonicalised separately
    "notebooks/07_WVD_Paper_Population_Hexbin.ipynb",
    "notebooks/08_CHMI_CSV_to_NetCDF4.ipynb",
    "notebooks/outputs/.gitkeep",
]

# moved, byte-identical, into docs/status/
STATUS_MOVES = {
    "NOTEBOOK_STATUS.md": "docs/status/NOTEBOOK_STATUS.md",
    "TEST_STATUS.md": "docs/status/TEST_STATUS.md",
    "TEST_LOG_NOTEBOOKS.md": "docs/status/TEST_LOG_NOTEBOOKS.md",
}

# moved, byte-identical, into attic/ (historical repo-management artifacts)
ATTIC_MOVES = {
    "README.md": "attic/README_handoff_root.md",
    "README_NOTEBOOKS.md": "attic/README_NOTEBOOKS_handoff.md",
    "RUN_FIRST.txt": "attic/RUN_FIRST_handoff.txt",
    ".gitignore": "attic/gitignore_handoff.txt",
    "NOTEBOOK_ADDED_SHA256.json": "attic/NOTEBOOK_ADDED_SHA256.json",
    "NOTEBOOK_RELEASE_MANIFEST.txt": "attic/NOTEBOOK_RELEASE_MANIFEST.txt",
    "upgrade_ESSD_wrapper_python_ONLY.py": "attic/upgrade_ESSD_wrapper_python_ONLY.py",
    "notebooks/06_WVD_Sensitivity_Length_Frequency.ipynb":
        "attic/notebook_versions/06_WVD_Sensitivity_Length_Frequency__V2_original.ipynb",
    "notebooks/06_WVD_Sensitivity_Length_Frequency_FIXED_NETCDF4.ipynb":
        "attic/notebook_versions/06_WVD_Sensitivity_Length_Frequency__FIXED_NETCDF4.ipynb",
    "config/active_profile.json": "attic/active_profile_hpc_handoff.json",
}

CANONICAL_06_SOURCE = "notebooks/06_WVD_Sensitivity_Length_Frequency_REFRESHED_V6.ipynb"
CANONICAL_06_DEST = "notebooks/06_WVD_Sensitivity_Length_Frequency.ipynb"

# small executed HPC evidence kept in git (paths relative to results/connected_workflow)
EVIDENCE_FILES = [
    "pipeline_manifest.json",
    "00_exploration/network/reviewer_f01_network_geometry.png",
    "00_exploration/cml/reviewer_f03_cml_characteristics.png",
    "00_exploration/monthly_availability/reviewer_f02_monthly_availability.csv",
    "00_exploration/monthly_availability/reviewer_f02_monthly_availability.png",
    "00_exploration/native_availability/availability_native_timestep_percentage.png",
    "00_exploration/native_availability/availability_native_timestep_percentage.pdf",
    "01_rainfall_rebuild/rebuild_summary.json",
    "01_rainfall_rebuild/rainfall_retrieval_per_sublink.csv",
    "02_rainfall_figures/reviewer_cml_rain_timeseries.png",
    "02_gauge_validation/validation_workspace/figures/FINAL_GRAF_validation_summary.json",
    "02_gauge_validation/validation_workspace/figures/FINAL_GRAF_confusion_GAUGE_same_population_202104_202110.pdf",
    "02_gauge_validation/validation_workspace/figures/FINAL_GRAF_hexbin_gauge_202104_202110.pdf",
    "02_gauge_validation/validation_workspace/figures/FINAL_GRAF_hexbin_gauge_WETONLY_202104_202110.pdf",
    "02_gauge_validation/validation_workspace/figures/FINAL_GRAF_gauge_shared_population_audit.csv",
    "04_maps/event04_instance01_multitime_summary.csv",
    "04_maps/hourly_1h/event04_instance01_multitime_hourly1h.png",
    "04_maps/native_5min/event04_instance01_multitime_native5min.png",
    "05_wvd/reviewer_wvd_v6_small.summary.json",
    "05_wvd/reviewer_wvd_v6_small_diagnostics.csv",
    "05_wvd/reviewer_wvd_v6_small_events.csv",
    "05_wvd/reviewer_wvd_v6_small_reinit.csv",
    "05_wvd/reviewer_wvd_v6_small_states.csv",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/run_summary.json",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/population_counts.csv",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/binned_validation_metrics.csv",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/validation_population_links.csv",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/station_assignments.csv",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/indexing_compat_provenance.json",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/33_v6_binned_validation_hexbins_NETCDF4_COMPAT.py",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/figures/01_hexbin_by_sensitivity_graf.pdf",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/figures/02_hexbin_by_frequency_graf.pdf",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/figures/03_hexbin_by_length_graf.pdf",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/figures/04_hexbin_frequency_x_length_graf.pdf",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/figures/05_hexbin_frequency_x_sensitivity_graf.pdf",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/figures/06_hexbin_joint_71_76GHz_graf.pdf",
    "WVD_V6_QUANTIZATION_BINNED_VALIDATION/figures/07_hexbin_joint_81_86GHz_graf.pdf",
    "07_wvd_paper_population/wvp_support_metrics.csv",
    "07_wvd_paper_population/wvp_support_70GHz.pdf",
    "07_wvd_paper_population/wvp_support_80GHz.pdf",
    "07_wvd_paper_population/wvp_support_all_E_band.pdf",
    "07_wvd_paper_population/wvp_support_three_panel.pdf",
]
EVIDENCE_DEST = "docs/evidence/hpc_run_2026-10"

PORTABLE_ACTIVE_PROFILE = {
    "profile": "user",
    "workspace": "results/connected_workflow_user",
    "science_python_by_profile": {
        "hpc": "/home/jovyan/my-conda-envs/tf-env/bin/python",
        "user": None,
    },
}

# exact, minimal compatibility patch for the reordered-IDs regression test:
# the synthetic rain file stores uint8 flags when netCDF4 is present, and current
# xarray refuses to write uint8 through the scipy NETCDF3 backend. Mirror the
# fixture's own engine selection instead of forcing scipy.
TEST_CHAIN = "tests/test_connected_chain.py"
TEST_CHAIN_OLD = "    swap=tmp_path/'rain_swapped.nc'; reordered.to_netcdf(swap,engine='scipy')\n"
TEST_CHAIN_NEW = (
    "    try:import netCDF4;engine='netcdf4'  # match mini_fixture: uint8 flags need NETCDF4\n"
    "    except ImportError:engine='scipy'\n"
    "    swap=tmp_path/'rain_swapped.nc'; reordered.to_netcdf(swap,engine=engine)\n"
)

# Exact, minimal in-cell notebook patches. Each entry must match exactly once.
# They fall into two classes, both required by the repository's own
# tests/test_notebook_adapters.py safe-default and portability assertions:
#   * portability: no hardcoded HPC interpreter path, no '!pip' shell magic;
#   * safe defaults: the author's HPC walkthrough left expensive/unverified
#     toggles enabled in the executed notebooks; public copies must default off.
NOTEBOOK_PATCHES = {
    "notebooks/01_Graf_Rainfall.ipynb": [(
        "RECOMPUTE = True       # False: register released rain without recomputing; safest first run",
        "RECOMPUTE = False      # False: register released rain without recomputing; safest first run",
        "safe default: full 2,219-link recompute off by default (HPC run left it True; "
        "executed evidence preserved under docs/evidence/hpc_run_2026-10/)",
    )],
    "notebooks/02_Rainfall_Validation.ipynb": [
        (
            "SCIENCE_PYTHON = (\n"
            "        \"/home/jovyan/my-conda-envs/tf-env/bin/python\"\n"
            "        if PROFILE == \"hpc\" else sys.executable\n"
            "    )",
            "# Set science_python_by_profile in config/active_profile.json instead of\n"
            "    # editing this notebook; the notebook kernel's interpreter is the fallback.\n"
            "    SCIENCE_PYTHON = sys.executable",
            "portability: hardcoded /home/jovyan interpreter fallback removed; the science "
            "interpreter comes from config/active_profile.json",
        ),
        (
            "RUN_GAUGE_VALIDATION = True  # computationally heavy; original publication detector unresolved",
            "RUN_GAUGE_VALIDATION = False  # computationally heavy; original publication detector unresolved",
            "safe default: heavy gauge validation off by default (HPC run left it True)",
        ),
        (
            "ALLOW_UNVERIFIED_V3 = True  # explicit opt-in only; this v3 is NOT proven historical source",
            "ALLOW_UNVERIFIED_V3 = False  # explicit opt-in only; this v3 is NOT proven historical source",
            "safe default: unverified v3 gauge detector requires explicit opt-in "
            "(HPC run left it True)",
        ),
    ],
    "notebooks/04_Spatial_Radar_Maps.ipynb": [
        (
            "! pip install cartopy\n",
            "# Requires cartopy; if missing, install it in this notebook's environment first (pip install cartopy).\n",
            "portability: '!pip' IPython shell magic is not valid Python and broke the "
            "repository's own notebook-syntax test and non-IPython execution",
        ),
        (
            "RUN_MAPS=True     # True to execute original map renderer; requires daily native RADKLIM and Cartopy",
            "RUN_MAPS=False    # True to execute original map renderer; requires daily native RADKLIM and Cartopy",
            "safe default: external native-RADKLIM map rendering off by default "
            "(HPC run left it True)",
        ),
    ],
    "notebooks/06_WVD_Sensitivity_Length_Frequency.ipynb": [(
        "RUN_SENSITIVITY = True  # safe first pass; set True AFTER 1352-series check passes",
        "RUN_SENSITIVITY = False  # set True AFTER the 1,352-series input check above passes",
        "safe default: expensive full-population sensitivity pass off by default, matching "
        "the notebook's own markdown instructions (HPC run left it True with a stale comment)",
    )],
}


def main():
    src_root, dst_root, report_out = sys.argv[1], sys.argv[2], sys.argv[3]
    src_root = os.path.abspath(src_root)
    dst_root = os.path.abspath(dst_root)
    if os.path.exists(wl(dst_root)):
        raise SystemExit(f"Destination exists, refusing to overwrite: {dst_root}")

    copied = []    # (src_rel, dst_rel) byte-identical pairs
    generated = []  # dst_rel
    patched = []   # dicts

    def byte_copy(rel_src, rel_dst=None):
        rel_dst = rel_dst or rel_src
        copy_file(os.path.join(src_root, rel_src), os.path.join(dst_root, rel_dst))
        copied.append((rel_src, rel_dst))

    # 1) whole directories
    for d in BYTE_COPY_DIRS:
        base = os.path.join(src_root, d)
        for root, _dirs, files in os.walk(wl(base)):
            for fn in files:
                full = os.path.join(root, fn).replace(BS * 2 + "?" + BS, "")
                rel = os.path.relpath(full, src_root).replace(os.sep, "/")
                byte_copy(rel)

    # 2) root files, configs, notebooks
    for rel in ROOT_BYTE_COPY:
        if rel == ".gitignore":
            continue  # handled via attic + generated replacement
        byte_copy(rel)
    for rel in CONFIG_BYTE_COPY:
        byte_copy(rel)
    for rel in NOTEBOOKS_BYTE_COPY:
        byte_copy(rel)

    # 3) status moves and attic moves (byte-identical, new location)
    for rel_src, rel_dst in {**STATUS_MOVES, **ATTIC_MOVES}.items():
        byte_copy(rel_src, rel_dst)

    # 4) canonical notebook 06 (byte-identical content of REFRESHED_V6)
    byte_copy(CANONICAL_06_SOURCE, CANONICAL_06_DEST)

    # 4b) patch the chain test copied by the directory walk (engine selection)
    tc_src = os.path.join(src_root, TEST_CHAIN)
    tc_dst = os.path.join(dst_root, TEST_CHAIN)
    text = open(wl(tc_dst), encoding="utf-8").read()
    if text.count(TEST_CHAIN_OLD) != 1:
        raise SystemExit("Expected exactly 1 occurrence of the scipy-engine line in test_connected_chain.py")
    with open(wl(tc_dst), "w", encoding="utf-8", newline="") as f:
        f.write(text.replace(TEST_CHAIN_OLD, TEST_CHAIN_NEW))
    copied[:] = [(s, d) for s, d in copied if s != TEST_CHAIN]
    patched.append({
        "file": TEST_CHAIN,
        "change": "reordered-IDs regression test now writes its temporary swapped file with "
                  "the same NetCDF engine selection as tests/mini_fixture.py (netCDF4 when "
                  "available) instead of always forcing the scipy NETCDF3 backend, which "
                  "cannot store the uint8 wet/dry flags under current xarray",
        "source_sha256": sha256(tc_src),
        "dest_sha256": sha256(tc_dst),
    })

    # 5) table-driven notebook patches, applied to the DEST copies (06 has already
    #    been canonicalised, so patch targets are the shipped notebook names)
    for nb_rel, rules in NOTEBOOK_PATCHES.items():
        nb_dst = os.path.join(dst_root, nb_rel)
        before_sha = sha256(nb_dst)
        with open(wl(nb_dst), encoding="utf-8") as f:
            nb = json.load(f)
        changes = []
        for old, new, why in rules:
            hits = 0
            for cell in nb["cells"]:
                if cell["cell_type"] != "code":
                    continue
                src_text = "".join(cell["source"])
                if old in src_text:
                    cell["source"] = src_text.replace(old, new).splitlines(keepends=True)
                    hits += 1
            if hits != 1:
                raise SystemExit(f"{nb_rel}: expected exactly 1 occurrence of patch "
                                 f"({why.split(':')[0]}), found {hits}")
            changes.append(why)
        with open(wl(nb_dst), "w", encoding="utf-8", newline="\n") as f:
            json.dump(nb, f, indent=1, ensure_ascii=False)
            f.write("\n")
        copied[:] = [(s, d) for s, d in copied if d != nb_rel]
        patched.append({
            "file": nb_rel,
            "changes": changes,
            "pre_patch_sha256": before_sha,
            "dest_sha256": sha256(nb_dst),
        })

    # 6) portable active profile
    ap = os.path.join(dst_root, "config/active_profile.json")
    os.makedirs(wl(os.path.dirname(ap)), exist_ok=True)
    with open(wl(ap), "w", encoding="utf-8", newline="\n") as f:
        json.dump(PORTABLE_ACTIVE_PROFILE, f, indent=2)
        f.write("\n")
    generated.append("config/active_profile.json")

    # 7) executed evidence subset
    ev_manifest = []
    for rel in EVIDENCE_FILES:
        s = os.path.join(src_root, "results/connected_workflow", rel)
        d = os.path.join(dst_root, EVIDENCE_DEST, rel)
        copy_file(s, d)
        ev_manifest.append({
            "evidence_path": f"{EVIDENCE_DEST}/{rel}",
            "original_result_path": f"results/connected_workflow/{rel}",
            "sha256": sha256(s),
            "bytes": os.path.getsize(wl(s)),
        })
        copied.append((f"results/connected_workflow/{rel}", f"{EVIDENCE_DEST}/{rel}"))
    with open(wl(os.path.join(dst_root, EVIDENCE_DEST, "EVIDENCE_MANIFEST.json")), "w",
              encoding="utf-8", newline="\n") as f:
        json.dump({
            "description": "Byte-identical copies of selected small outputs of the author's "
                           "real HPC execution (October 2026) of the connected workflow. "
                           "Large NetCDF/npz outputs are not in git; see docs/EVIDENCE.md.",
            "entries": ev_manifest,
        }, f, indent=2)
        f.write("\n")
    generated.append(f"{EVIDENCE_DEST}/EVIDENCE_MANIFEST.json")

    # 8) placeholder dirs
    for d in ("results", "notebooks/outputs", "examples", "figures"):
        os.makedirs(wl(os.path.join(dst_root, d)), exist_ok=True)
    for keep in ("results/.gitkeep", "figures/.gitkeep"):
        open(wl(os.path.join(dst_root, keep)), "w").close()
        generated.append(keep)

    # ---------------------------------------------------------- verification
    mismatched = []
    for rel_src, rel_dst in copied:
        if sha256(os.path.join(src_root, rel_src)) != sha256(os.path.join(dst_root, rel_dst)):
            mismatched.append({"src": rel_src, "dst": rel_dst})

    report = {
        "source_root": src_root,
        "dest_root": dst_root,
        "n_byte_copied": len(copied),
        "n_generated_by_script": len(generated),
        "patched": patched,
        "generated": generated,
        "byte_copy_mismatches": mismatched,
        "copied_pairs": [{"src": s, "dst": d} for s, d in copied],
    }
    with open(wl(report_out), "w", encoding="utf-8", newline="\n") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: report[k] for k in ("n_byte_copied", "n_generated_by_script")}, indent=2))
    print("patched:", [p["file"] for p in patched])
    print("BYTE-COPY VERIFY:", "OK" if not mismatched else f"{len(mismatched)} MISMATCHES")
    if mismatched:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
