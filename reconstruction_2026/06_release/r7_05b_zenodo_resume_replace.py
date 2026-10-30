"""Stage 7: resume the replacement of r7_05 after an interrupted upload (curl error 55 on the ZIP).

Per file of the corrected package: an entry already completed with the local MD5 is kept; a pending entry (started but
without content) gets the content again; a completed entry with another MD5 is deleted, started and uploaded. Each
upload is retried up to four times and every commit is checked against the local MD5. Afterwards only the
`description` metadata field is set (the text of r7_05). Draft only: no publish call, no new DOI, no new file keys.
Usage: python -I r7_05b_zenodo_resume_replace.py <release dir>
"""
import ast, hashlib, json, os, subprocess, sys, time
import requests

API = "https://zenodo.org/api"
REC = 23398309
rel = sys.argv[1]
LOCAL = {"ESSD_CZ_CML_2021_dataset.zip": os.path.join(rel, "ESSD_CZ_CML_2021_dataset.zip"),
         "README_DATA.md": os.path.join(rel, "package", "README_DATA.md"),
         "DATA_DICTIONARY.md": os.path.join(rel, "package", "DATA_DICTIONARY.md"),
         "SHA256SUMS.txt": os.path.join(rel, "package", "SHA256SUMS.txt")}
# the description text is taken unchanged from r7_05 (not executed: its state check fails while an entry is pending)
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "r7_05_zenodo_replace_draft.py"), encoding="utf-8").read()
DESCRIPTION = next(ast.literal_eval(n.value) for n in ast.parse(src).body
                   if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "DESCRIPTION")


def md5_of(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


tok = open(os.path.join(os.environ["USERPROFILE"], ".zenodo_token")).read().strip()
H = {"Authorization": f"Bearer {tok}"}
hdr = os.path.join(os.environ["USERPROFILE"], ".zenodo_curl_header")
files = f"{API}/records/{REC}/draft/files"


def entries():
    r = requests.get(files, headers=H, timeout=60)
    r.raise_for_status()
    return {e["key"]: e for e in r.json()["entries"]}


def upload(name, p):
    for attempt in range(1, 5):
        rr = subprocess.run(["curl.exe", "--fail", "-sS", "-o", "NUL", "-H", f"@{hdr}",
                             "-H", "Content-Type: application/octet-stream", "-X", "PUT", "-T", p,
                             f"{files}/{name}/content"], capture_output=True, text=True)
        if rr.returncode == 0:
            return
        print(f"  upload {name} attempt {attempt} failed: {rr.stderr.strip()[:200]}", flush=True)
        time.sleep(60)
    raise SystemExit(f"upload {name} failed after 4 attempts")


ent = entries()
# the interrupted ZIP entry was removed by Zenodo; only that key of the earlier draft may be missing
if not set(ent) <= set(LOCAL) or set(LOCAL) - set(ent) - {"ESSD_CZ_CML_2021_dataset.zip"}:
    raise SystemExit(f"unexpected draft files: {sorted(ent)}")
for name, p in LOCAL.items():
    md5 = md5_of(p); e = ent.get(name, {})
    remote = (e.get("checksum") or "").replace("md5:", "")
    if e.get("status") == "completed" and remote == md5:
        print(f"already replaced: {name} md5 {md5}", flush=True); continue
    if not e:
        r = requests.post(files, headers=H, json=[{"key": name}], timeout=120)
        if r.status_code >= 400:
            raise SystemExit(f"start {name} failed ({r.status_code}): {r.text[:300]}")
    elif e.get("status") == "completed":
        r = requests.delete(f"{files}/{name}", headers=H, timeout=120)
        if r.status_code >= 400:
            raise SystemExit(f"delete {name} failed ({r.status_code}): {r.text[:300]}")
        r = requests.post(files, headers=H, json=[{"key": name}], timeout=120)
        if r.status_code >= 400:
            raise SystemExit(f"start {name} failed ({r.status_code}): {r.text[:300]}")
    upload(name, p)
    r = requests.post(f"{files}/{name}/commit", headers=H, timeout=600)
    if r.status_code >= 400:
        raise SystemExit(f"commit {name} failed: {r.text[:300]}")
    remote = r.json().get("checksum", "").replace("md5:", "")
    if remote != md5:
        raise SystemExit(f"MD5 mismatch for {name}: remote {remote} local {md5}")
    print(f"replaced and verified: {name} size {os.path.getsize(p)} md5 {md5}", flush=True)

dep_url = f"{API}/deposit/depositions/{REC}"
dep = requests.get(dep_url, headers=H, timeout=60).json()
if dep.get("submitted") or dep.get("state") != "unsubmitted":
    raise SystemExit(f"unexpected record state: state={dep.get('state')} submitted={dep.get('submitted')}")
meta = dict(dep["metadata"]); meta["description"] = DESCRIPTION
r = requests.put(dep_url, headers={**H, "Content-Type": "application/json"}, data=json.dumps({"metadata": meta}), timeout=120)
if r.status_code >= 400:
    raise SystemExit(f"metadata update failed ({r.status_code}): {r.text[:500]}")
after = requests.get(dep_url, headers=H, timeout=60).json()
changed = [k for k in set(after["metadata"]) | set(dep["metadata"]) if after["metadata"].get(k) != dep["metadata"].get(k)]
print("metadata fields changed:", changed)
print("state:", after.get("state"), "| submitted:", after.get("submitted"), "| DOI:", after["metadata"].get("prereserve_doi", {}).get("doi"))
for e in entries().values():
    print(f"  {e['key']}  {e.get('size')}  {e.get('checksum')}  {e.get('status')}")
