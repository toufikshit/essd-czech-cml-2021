"""Stage 5B: pixel comparison of a reproduced figure with the published image, and old|new side-by-side panels.
Usage: python -I s5b_02_pixel_compare.py"""
import os, json
import numpy as np
from PIL import Image
H = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
O = os.path.join(H, 'original_figures'); R = os.path.join(H, 'fig568', 'control_reproduction'); N = os.path.join(H, 'fig568', 'new')
C = os.path.join(H, 'fig568', 'comparison'); os.makedirs(C, exist_ok=True)
res = {}
for f in ('FINAL_GRAF_confusion_GAUGE_same_population_202104_202110.png', 'FINAL_GRAF_hexbin_gauge_WETONLY_202104_202110.png'):
    a = np.asarray(Image.open(os.path.join(O, f)).convert('RGB'), np.int16)
    b = np.asarray(Image.open(os.path.join(R, f)).convert('RGB'), np.int16)
    n = np.asarray(Image.open(os.path.join(N, f)).convert('RGB'), np.int16)
    same = a.shape == b.shape
    d = np.abs(a - b) if same else None
    res[f] = dict(published_px=list(a.shape[:2][::-1]), reproduced_px=list(b.shape[:2][::-1]), new_px=list(n.shape[:2][::-1]),
                  identical=bool(same and d.max() == 0),
                  frac_pixels_differing=float((d.max(2) > 0).mean()) if same else None,
                  max_abs_diff=int(d.max()) if same else None,
                  frac_pixels_diff_gt_30=float((d.max(2) > 30).mean()) if same else None,
                  new_vs_published_frac_differing=float((np.abs(a - n).max(2) > 0).mean()) if a.shape == n.shape else None)
    # side by side: published | new
    A = Image.open(os.path.join(O, f)).convert('RGB'); B = Image.open(os.path.join(N, f)).convert('RGB')
    h = max(A.height, B.height); canvas = Image.new('RGB', (A.width + B.width + 40, h), 'white')
    canvas.paste(A, (0, 0)); canvas.paste(B, (A.width + 40, 0))
    canvas.thumbnail((2400, 2400)); canvas.save(os.path.join(C, f.replace('.png', '_old_vs_new.png')))
print(json.dumps(res, indent=2))
