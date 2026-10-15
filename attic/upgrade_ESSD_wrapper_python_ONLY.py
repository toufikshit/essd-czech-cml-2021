#!/usr/bin/env python3
"""Safely migrate a V2 connected ESSD repository to V2.1 interpreter selection.

Updates ONLY two wrapper .py files, one notebook code cell, and shared active
profile configuration. Preserves original/ code and all results/ and manifests.
Reads updated source bytes from ESSD_CZ_2021_CONNECTED_WALKTHROUGH_V2_1.zip.
"""
from __future__ import annotations
import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile

BASE = 'ESSD_CZ_2021_CONNECTED_WALKTHROUGH_V2/'
UPDATES = ('src/essd_cz/workflows.py','src/essd_cz/connected.py')
NB = 'notebooks/01_Graf_Rainfall.ipynb'
CFG = 'config/active_profile.json'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path.cwd())
    p.add_argument('--archive', type=Path, default=None)
    args = p.parse_args()
    repo = args.repo.expanduser().resolve()
    archive = (args.archive or repo.parent/'ESSD_CZ_2021_CONNECTED_WALKTHROUGH_V2_1.zip').expanduser().resolve()
    if not (repo/'tools_verify_originals.py').is_file():
        raise RuntimeError(f'Not an ESSD reviewer repository: {repo}')
    if not archive.is_file(): raise FileNotFoundError(f'Upload the V2.1 archive: {archive}')
    new = {}
    with ZipFile(archive) as z:
        if z.testzip() is not None: raise RuntimeError('Updated archive failed ZIP integrity')
        for name in (*UPDATES,NB,CFG):
            new[name] = z.read(BASE+name)
    # Compute new contents in memory; no writes on a failed preflight.
    updated={}
    for name in UPDATES:
        path=repo/name
        old=path.read_text(encoding='utf-8')
        if name.endswith('workflows.py'):
            expected = 'def rebuild(paths, output, workers=4):'
        else:
            expected = 'def recompute_rainfall(self, workers=4):'
        replacement = new[name].decode('utf-8')
        if replacement == old:
            print('ALREADY PATCHED',name)
            continue
        if expected not in old:
            raise RuntimeError(f'Refusing to replace custom or unknown wrapper: {path}')
        updated[name] = replacement.encode('utf-8')
    old_nb=json.loads((repo/NB).read_text(encoding='utf-8'))
    new_nb=json.loads(new[NB].decode('utf-8'))
    assert len(old_nb['cells'])>2 and len(new_nb['cells'])>2
    oldsrc=''.join(old_nb['cells'][2]['source'])
    newsrc=new_nb['cells'][2]['source']
    if 'SCIENCE_PYTHON' not in oldsrc:
        if 'project.recompute_rainfall(workers=WORKERS)' not in oldsrc:
            raise RuntimeError('Notebook 01 cell 2 has unexpected/custom code; refusing to alter')
        old_nb['cells'][2]['source']=newsrc
        # Discard stale execution outputs only for the edited cell; keep all other cells.
        old_nb['cells'][2]['execution_count']=None
        old_nb['cells'][2]['outputs']=[]
        updated[NB] = (json.dumps(old_nb,indent=1,ensure_ascii=False)+'\n').encode('utf-8')
    else:
        print('ALREADY PATCHED',NB)
    current=json.loads((repo/CFG).read_text(encoding='utf-8'))
    if 'science_python_by_profile' not in current:
        current['science_python_by_profile']={
            'hpc':'/home/jovyan/my-conda-envs/tf-env/bin/python',
            'user':None,
        }
        updated[CFG]=(json.dumps(current,indent=2)+'\n').encode('utf-8')
    else:
        print('ALREADY CONFIGURED',CFG)
    if not updated:
        print('NO CHANGES NEEDED. Existing results and manifest preserved.')
        return
    timestamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    backup=repo/'results'/f'V2_before_scientific_python_fix_{timestamp}'
    if backup.exists():raise FileExistsError(backup)
    for name in updated:
        dest=backup/name
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(repo/name,dest)
    for name,contents in updated.items():
        (repo/name).write_bytes(contents)
        print('PATCHED',name)
    print('BACKUPS',backup)
    print('The original scientific scripts, data, and pipeline_manifest.json were not touched.')
    print('Restart the notebook kernel, then rerun Notebook 01 from the start.')

if __name__=='__main__':
    try: main()
    except Exception as exc:
        print(type(exc).__name__+': '+str(exc),file=sys.stderr)
        sys.exit(1)
