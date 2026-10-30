#!/usr/bin/env python3
"""Check that the scripts in original/code/ are unchanged (SHA-256 against original/SHA256SUMS.txt)."""
import hashlib
from pathlib import Path

root = Path(__file__).resolve().parent / 'original'
bad = []
entries = [line.split('  ', 1) for line in (root / 'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines() if line]
for sha, rel in entries:
    path = root / rel
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != sha:
        bad.append(rel)
print(f'Verified {len(entries) - len(bad)}/{len(entries)} files in original/')
if bad:
    raise SystemExit('CHANGED OR MISSING: ' + ', '.join(bad))
print('PASS: all files in original/ are unchanged')
