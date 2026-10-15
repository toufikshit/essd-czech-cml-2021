#!/usr/bin/env python3
"""Confirm original scripts/notebooks have not been changed by restructuring."""
import hashlib
import json
from pathlib import Path
root=Path(__file__).resolve().parent/'original'
manifest=json.loads((root/'CONSOLIDATED_SOURCE_MANIFEST.json').read_text())
bad=[]
for item in manifest['originals']:
    path=root/item['archive_path']
    h=hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else 'NOT FOUND'
    if h!=item['sha256']:bad.append(str(path))
print(f"Verified {len(manifest['originals'])-len(bad)}/{len(manifest['originals'])} original source files/notebooks")
if bad:
    print('FAILED:','\n'.join(bad))
    raise SystemExit(1)
print('PASS: all originals byte-identical to recovered files')
