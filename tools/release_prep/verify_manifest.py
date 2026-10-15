"""Verify extracted handoff trees against handoff/FILE_MANIFEST_SHA256.json.

Usage: python -I verify_manifest.py <manifest.json> <extract_root> <report_out.json> [--skip-hash-above BYTES]

extract_root is the directory containing repository/, data/, handoff/ (i.e. the
ESSD_CZ_2021_HANDOFF folder). Writes a JSON report with matched / mismatched /
missing / extra files and a size inventory by group.
"""
import hashlib
import json
import os
import sys
import time


def win_long(p):
    p = os.path.abspath(p)
    if os.name == "nt" and not p.startswith("\\\\?\\"):
        p = "\\\\?\\" + p
    return p


def sha256_file(path, bufsize=16 * 1024 * 1024):
    h = hashlib.sha256()
    with open(win_long(path), "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main():
    manifest_path, extract_root, report_out = sys.argv[1], sys.argv[2], sys.argv[3]
    with open(win_long(manifest_path), "r", encoding="utf-8") as f:
        manifest = json.load(f)
    entries = manifest["entries"]

    t0 = time.time()
    matched, mismatched, missing = [], [], []
    size_by_group = {}
    total_bytes = 0
    for i, e in enumerate(entries):
        rel = e["archive_path"]
        p = os.path.join(extract_root, rel.replace("/", os.sep))
        g = e.get("group", "?")
        if not os.path.exists(win_long(p)):
            missing.append(rel)
            continue
        actual_bytes = os.path.getsize(win_long(p))
        actual_hash = sha256_file(p)
        total_bytes += actual_bytes
        size_by_group[g] = size_by_group.get(g, 0) + actual_bytes
        if actual_hash == e["sha256"] and actual_bytes == e["bytes"]:
            matched.append(rel)
        else:
            mismatched.append({
                "archive_path": rel,
                "expected_sha256": e["sha256"], "actual_sha256": actual_hash,
                "expected_bytes": e["bytes"], "actual_bytes": actual_bytes,
            })
        if (i + 1) % 25 == 0:
            print(f"  hashed {i+1}/{len(entries)} ({total_bytes/1e9:.2f} GB, {time.time()-t0:.0f}s)", flush=True)

    # extra files present on disk but not in manifest
    manifest_set = {e["archive_path"] for e in entries}
    extra = []
    for root, _dirs, files in os.walk(win_long(extract_root)):
        for fn in files:
            full = os.path.join(root, fn)
            rel = os.path.relpath(full.replace("\\\\?\\", ""), os.path.abspath(extract_root)).replace(os.sep, "/")
            if rel not in manifest_set:
                extra.append(rel)

    report = {
        "manifest": os.path.abspath(manifest_path),
        "extract_root": os.path.abspath(extract_root),
        "n_entries": len(entries),
        "n_matched": len(matched),
        "n_mismatched": len(mismatched),
        "n_missing": len(missing),
        "n_extra_on_disk": len(extra),
        "mismatched": mismatched,
        "missing": missing,
        "extra_on_disk": sorted(extra),
        "bytes_hashed": total_bytes,
        "size_by_group_bytes": size_by_group,
        "declared_missing_or_optional": manifest.get("missing_or_optional"),
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(win_long(report_out), "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: report[k] for k in (
        "n_entries", "n_matched", "n_mismatched", "n_missing", "n_extra_on_disk", "bytes_hashed", "elapsed_s")}, indent=2))
    print("OK" if not mismatched and not missing else "PROBLEMS FOUND")


if __name__ == "__main__":
    main()
