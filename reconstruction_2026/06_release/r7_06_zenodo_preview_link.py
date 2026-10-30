"""Stage 7: create a "Can preview" share link for the unpublished Zenodo draft 23398309 and test it without login.

Uses the access-links endpoint of the record (parent) (permission "preview"); no publish call, no metadata change.
The test opens the draft page and the draft file list with the link token only (no Authorization header) and
compares the files seen with the corrected package (names and MD5).
The link (it contains the access token) is written to a file in the user profile, not printed in full.
Usage: python -I r7_06_zenodo_preview_link.py [--create]
"""
import json, os, sys
import requests

API = "https://zenodo.org/api"
REC = 23398309
EXPECTED = {"ESSD_CZ_CML_2021_dataset.zip": "b807b0fb672784e63117312a0d37dddd"}
tok = open(os.path.join(os.environ["USERPROFILE"], ".zenodo_token")).read().strip()
H = {"Authorization": f"Bearer {tok}"}
out = os.path.join(os.environ["USERPROFILE"], ".zenodo_preview_link.txt")

dep = requests.get(f"{API}/deposit/depositions/{REC}", headers=H, timeout=60).json()
if dep.get("submitted") or dep.get("state") != "unsubmitted":
    raise SystemExit(f"unexpected record state: {dep.get('state')}")
links_url = f"{API}/records/{REC}/access/links"  # share links belong to the parent record
r = requests.get(links_url, headers=H, timeout=60)
print("existing links:", r.status_code, [(h.get("id"), h.get("permission"), h.get("expires_at")) for h in r.json().get("hits", {}).get("hits", [])] if r.ok else r.text[:300])
if "--create" in sys.argv:
    r = requests.post(links_url, headers={**H, "Content-Type": "application/json"}, data=json.dumps({"permission": "preview"}), timeout=60)
    if r.status_code >= 400:
        raise SystemExit(f"create failed ({r.status_code}): {r.text[:500]}")
    link = r.json()
    print("created link id", link.get("id"), "permission", link.get("permission"), "expires", link.get("expires_at"))
    token = link["token"]
    url = f"https://zenodo.org/records/{REC}?preview=1&token={token}"
    with open(out, "w") as f:
        f.write(url + "\n")
    print("link written to", out)
if not os.path.exists(out):
    raise SystemExit("no link yet (run with --create)")
url = open(out).read().strip(); token = url.split("token=")[1]

# test without login
s = requests.Session()
p = s.get(url, timeout=120, allow_redirects=True)
print("draft page without login:", p.status_code, "| title found:", "Operational microwave links" in p.text, "| 2,574 in page:", "2,574" in p.text)
f = s.get(f"{API}/records/{REC}/draft/files", params={"token": token}, timeout=120)
print("draft file list with token only:", f.status_code)
if f.ok:
    for e in f.json().get("entries", []):
        md5 = (e.get("checksum") or "").replace("md5:", "")
        tag = "matches corrected package" if EXPECTED.get(e["key"]) == md5 else ""
        print(f"  {e['key']}  {e.get('size')}  {md5}  {tag}")
    dl = s.get(f"{API}/records/{REC}/draft/files/README_DATA.md/content", params={"token": token}, timeout=120)
    print("download README_DATA.md with token only:", dl.status_code, "| corrected text (2,574 CMLs):", "2,574 CMLs" in dl.text)
    head = s.get(f"{API}/records/{REC}/draft/files/ESSD_CZ_CML_2021_dataset.zip/content", params={"token": token},
                 headers={"Range": "bytes=0-1023"}, timeout=120, stream=True)
    print("ZIP readable with token only (first 1 kB):", head.status_code, "| starts with PK:", head.raw.read(2) == b"PK")
nolink = requests.get(f"https://zenodo.org/records/{REC}", timeout=60)
print("same page without the token:", nolink.status_code, "(expected 404 while unpublished)")
