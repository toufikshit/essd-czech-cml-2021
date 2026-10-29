"""Stage 4: prove that the scientific code of s4_01_rain_1min.py is a verbatim copy of the final exporter.

Compares the block from 'from __future__ import annotations' to the end of meta2() in
essd-czech-cml-2021/src/essd_cz/export_1min.py with the same block in s4_01_rain_1min.py, and also checks the
exporter against the preserved verified MetaCentrum original (only the docstring and argument defaults differ).
Usage: python -I s4_00_check_verbatim.py [--repo DIR]
"""
import argparse, hashlib, os, difflib
ap = argparse.ArgumentParser()
ap.add_argument('--repo', default=r'C:\Users\taouf\ESSD_RELEASE\essd-czech-cml-2021')
a = ap.parse_args()
here = os.path.dirname(os.path.abspath(__file__))
def block(path):
    t = open(path, encoding='utf-8').read()
    i = t.index('from __future__ import annotations')
    end_marker = 'return decode(a) if strings else np.asarray(a,float)'
    j = t.index(end_marker) + len(end_marker)
    return t[i:j]
src = os.path.join(a.repo, 'src', 'essd_cz', 'export_1min.py')
ori = os.path.join(a.repo, 'docs', 'evidence', 'hpc_run_2026-10', 'cml_graf_export_1min_verified.py')
new = os.path.join(here, 's4_01_rain_1min.py')
for p in (src, ori, new):
    print(hashlib.sha256(open(p, 'rb').read()).hexdigest(), p)
bs, bo, bn = block(src), block(ori), block(new)
print('scientific block lines:', len(bs.splitlines()))
print('s4_01 block == maintained exporter block:', bn == bs, hashlib.sha256(bn.encode()).hexdigest())
print('maintained exporter block == verified MetaCentrum original block:', bs == bo)
for name in ('def estimate', 'def graf', 'def pol_hv', 'def call_kr'):
    print(f'  {name}: present', name in bn)
if bn != bs:
    print(''.join(difflib.unified_diff(bs.splitlines(True), bn.splitlines(True))))
