"""compare_images.py - side-by-side (or stacked) comparison image of an original and a new figure.

Usage: python -I -B compare_images.py ORIGINAL NEW OUT [--stack v|h] [--labels "Original (published)" "New"] [--width 2400]
Both images are scaled to the same width (stack v) or height (stack h) for display only; a title strip with the
label and the native pixel size is added above each image. Read-only on inputs.
"""
import argparse
from PIL import Image, ImageDraw, ImageFont
p = argparse.ArgumentParser()
p.add_argument("orig"); p.add_argument("new"); p.add_argument("out")
p.add_argument("--stack", default="v", choices=["v", "h"])
p.add_argument("--labels", nargs=2, default=["ORIGINAL (published)", "NEW (candidate)"])
p.add_argument("--width", type=int, default=2400)
a = p.parse_args()
ims = [Image.open(a.orig).convert("RGB"), Image.open(a.new).convert("RGB")]
try:
    from matplotlib import font_manager
    font = ImageFont.truetype(font_manager.findfont("DejaVu Sans"), 34)
except Exception:
    font = ImageFont.load_default()
tiles = []
for im, lab in zip(ims, a.labels):
    w0, h0 = im.size
    if a.stack == "v":
        s = a.width / w0
    else:
        s = (a.width * 0.45) / h0
    im2 = im.resize((max(1, int(w0 * s)), max(1, int(h0 * s))), Image.LANCZOS)
    t = Image.new("RGB", (im2.width, im2.height + 50), "white")
    d = ImageDraw.Draw(t); d.text((10, 8), f"{lab}   [{w0} x {h0} px]", fill="black", font=font)
    t.paste(im2, (0, 50)); tiles.append(t)
if a.stack == "v":
    W = max(t.width for t in tiles); H = sum(t.height for t in tiles) + 20
    out = Image.new("RGB", (W, H), (200, 200, 200)); y = 0
    for t in tiles: out.paste(t, (0, y)); y += t.height + 20
else:
    W = sum(t.width for t in tiles) + 20; H = max(t.height for t in tiles)
    out = Image.new("RGB", (W, H), (200, 200, 200)); x = 0
    for t in tiles: out.paste(t, (x, 0)); x += t.width + 20
out.save(a.out)
print("saved", a.out, out.size)
