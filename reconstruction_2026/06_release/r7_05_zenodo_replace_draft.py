"""Stage 7: replace the outdated files and description of the unpublished Zenodo draft 23398309.

Replaces the four files of the draft (ESSD_CZ_CML_2021_dataset.zip, README_DATA.md, DATA_DICTIONARY.md,
SHA256SUMS.txt) with the corrected package, using the same steps as tools/zenodo_replace_small_file.py
(delete entry, start, curl PUT, commit, MD5 check), and updates only the `description` metadata field.
Draft only: there is no publish call and no new DOI. Refuses to run if the record is not an unsubmitted draft.
Without --apply it prints what it would do.
Usage: python -I r7_05_zenodo_replace_draft.py <release dir> [--apply]
"""
import hashlib, json, os, subprocess, sys
import requests

API = "https://zenodo.org/api"
REC = 23398309
rel = sys.argv[1]; apply = '--apply' in sys.argv
LOCAL = {"ESSD_CZ_CML_2021_dataset.zip": os.path.join(rel, "ESSD_CZ_CML_2021_dataset.zip"),
         "README_DATA.md": os.path.join(rel, "package", "README_DATA.md"),
         "DATA_DICTIONARY.md": os.path.join(rel, "package", "DATA_DICTIONARY.md"),
         "SHA256SUMS.txt": os.path.join(rel, "package", "SHA256SUMS.txt")}
DESCRIPTION = (
    "Seven data products for the Czech Republic, April to October 2021, at native resolution: raw commercial microwave "
    "link (CML) signals at 1 min (2,574 links), quality controlled CML signals at 1 min (1,987 links), final CML "
    "rainfall at its native 1 min resolution for the 1,987 quality controlled links (Graf wet and dry detection, "
    "constant baseline, Pastorek wet antenna correction, ITU R P.838 2005 coefficients), CHMI rain gauges at 10 min "
    "(424 stations), RADKLIM radar rainfall averaged along CML paths at 5 min (1,132 paths), CHMI temperature at 10 min "
    "(251 stations) and CHMI relative humidity at 10 min (247 stations). Optional supplement: E band water vapour "
    "retrieval products (1,352 directional series). The CML rainfall is published only at 1 min; the 5 min series used "
    "in the companion article are derived from this file by averaging five complete valid 1-minute directional samples "
    "per block (wet when any valid minute is wet, physical-link mean over the finite directional 5-minute rates), and "
    "the accompanying software performs this derivation. The deposit includes a README, a data dictionary and SHA256 "
    "checksums. This record uses per file licences: the CML derived files are CC BY-NC 4.0 (provided under an agreement "
    "with T-Mobile Czech Republic, available on request; noncommercial use only) and the CHMI and DWD derived files are "
    "CC BY 4.0. The record licence field is set to other (open) for this reason; the deposit README and SHA256SUMS list "
    "the licence per file. The wet and dry thresholds and the dry baselines are computed on the April to October 2021 "
    "record and are stored inside the rainfall file. There was no calibration of the CML rainfall against gauges or "
    "radar.")


def md5_of(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


tok = open(os.path.join(os.environ["USERPROFILE"], ".zenodo_token")).read().strip()
H = {"Authorization": f"Bearer {tok}"}
hdr = os.path.join(os.environ["USERPROFILE"], ".zenodo_curl_header")
dep_url = f"{API}/deposit/depositions/{REC}"
dep = requests.get(dep_url, headers=H, timeout=60).json()
if dep.get("submitted") or dep.get("state") != "unsubmitted":
    raise SystemExit(f"unexpected record state: state={dep.get('state')} submitted={dep.get('submitted')}")
print("draft state:", dep.get("state"), "| reserved DOI:", dep["metadata"].get("prereserve_doi", {}).get("doi"))
files = f"{API}/records/{REC}/draft/files"
entries = {e["key"]: e for e in requests.get(files, headers=H, timeout=60).json()["entries"]}
plan = {}
for name, p in LOCAL.items():
    if name not in entries:
        raise SystemExit(f"{name} is not on the draft; refusing to add a new file")
    plan[name] = (p, os.path.getsize(p), md5_of(p))
    print(f"replace {name}: remote {entries[name].get('size')} {entries[name].get('checksum')} -> local {plan[name][1]} md5:{plan[name][2]}")
print("description: 2,580/2,219 -> 2,574/1,987; full-history thresholds sentence -> April to October 2021")
if not apply:
    raise SystemExit("dry run only (add --apply)")

for name, (p, size, md5) in plan.items():
    r = requests.delete(f"{files}/{name}", headers=H, timeout=120)
    if r.status_code >= 400:
        raise SystemExit(f"delete {name} failed ({r.status_code}): {r.text[:300]}")
    r = requests.post(files, headers=H, json=[{"key": name}], timeout=120)
    if r.status_code >= 400:
        raise SystemExit(f"start {name} failed ({r.status_code}): {r.text[:300]}")
    rr = subprocess.run(["curl.exe", "--fail", "-sS", "-o", "NUL", "-H", f"@{hdr}",
                         "-H", "Content-Type: application/octet-stream", "-X", "PUT", "-T", p,
                         f"{files}/{name}/content"], capture_output=True, text=True)
    if rr.returncode != 0:
        raise SystemExit(f"upload {name} failed: {rr.stderr[:300]}")
    r = requests.post(f"{files}/{name}/commit", headers=H, timeout=600)
    if r.status_code >= 400:
        raise SystemExit(f"commit {name} failed: {r.text[:300]}")
    remote = r.json().get("checksum", "").replace("md5:", "")
    if remote != md5:
        raise SystemExit(f"MD5 mismatch for {name}: remote {remote} local {md5}")
    print(f"replaced and verified: {name} size {size} md5 {md5}", flush=True)

meta = dict(dep["metadata"]); meta["description"] = DESCRIPTION
r = requests.put(dep_url, headers={**H, "Content-Type": "application/json"}, data=json.dumps({"metadata": meta}), timeout=120)
if r.status_code >= 400:
    raise SystemExit(f"metadata update failed ({r.status_code}): {r.text[:500]}")
after = requests.get(dep_url, headers=H, timeout=60).json()
changed = [k for k in set(after["metadata"]) | set(dep["metadata"]) if after["metadata"].get(k) != dep["metadata"].get(k)]
print("metadata fields changed:", changed)
print("state:", after.get("state"), "| submitted:", after.get("submitted"), "| DOI:", after["metadata"].get("prereserve_doi", {}).get("doi"))
for e in requests.get(files, headers=H, timeout=60).json()["entries"]:
    print(f"  {e['key']}  {e.get('size')}  {e.get('checksum')}  {e.get('status')}")
