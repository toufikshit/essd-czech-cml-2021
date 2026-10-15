#!/usr/bin/env python3
"""GATE 2b: unchanged five-sensor original figure on tiny artificial NetCDFs."""
from pathlib import Path
import argparse
import subprocess
import sys
import json
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tests'))
from mini_fixture import create,sha


def run(out):
    out=Path(out).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'Refusing to overwrite previous test fixture: {out}')
    inp=out/'input'
    paths=create(inp)
    hashes={r:sha(p) for r,p in paths.items()}
    exe=ROOT/'original/code/figures/availability_native_timestep_percentage.py'
    fig=out/'availability'
    cmd=[sys.executable,str(exe),
         '--cml-file',paths['qc'],'--radar-file',paths['radar'],
         '--gauge-file',paths['gauge'],'--temperature-file',paths['temperature'],
         '--humidity-file',paths['humidity'],
         '--start','2021-04-01','--end','2021-04-02',
         '--out-dir',str(fig),'--dpi','100']
    subprocess.run(cmd,check=True)
    for name in ['cml','radar','gauge','temperature','humidity']:
        assert (fig/f'availability_native_timestep_percentage_{name}.csv').is_file(),name
    for ext in ['png','pdf']:
        assert (fig/f'availability_native_timestep_percentage.{ext}').is_file(),ext
    if {r:sha(p) for r,p in paths.items()}!=hashes:
        raise AssertionError('Input data changed during availability plotting')
    print('PASS: 5 CSVs + PNG + PDF; all 7 input checksums unchanged')
    (out/'AVAILABILITY_MINI_PASS.json').write_text(json.dumps({'status':'PASS','input_sha256':hashes},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out-dir',type=Path,default=ROOT/'results/mini_availability')
    a=p.parse_args();run(a.out_dir)
