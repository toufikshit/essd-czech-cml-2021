#!/usr/bin/env python3
"""Stage 5B / Agent 3b - side-by-side image (original | new) for a figure, with a header line.
Usage: python -I f4_03_side_by_side.py <old.png> <new.png> <out.png> "<left title>" "<right title>" """
import sys
from PIL import Image, ImageDraw, ImageFont
old, new, out, lt, rt = sys.argv[1:6]
a = Image.open(old).convert('RGB'); b = Image.open(new).convert('RGB')
H = max(a.height, b.height); pad = 40; head = 110
c = Image.new('RGB', (a.width + b.width + pad, H + head), 'white')
c.paste(a, (0, head)); c.paste(b, (a.width + pad, head))
d = ImageDraw.Draw(c)
try:
    f = ImageFont.truetype('DejaVuSans.ttf', 60)
except Exception:
    try:
        f = ImageFont.truetype('arial.ttf', 60)
    except Exception:
        f = ImageFont.load_default()
d.text((20, 20), lt, fill='black', font=f); d.text((a.width + pad + 20, 20), rt, fill='black', font=f)
d.line([(a.width + pad // 2, head), (a.width + pad // 2, H + head)], fill=(150, 150, 150), width=4)
c.save(out)
print('WROTE', out, c.size)
