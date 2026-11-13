"""Gallery images of README.md: single panels cut from the figures stored in the executed notebooks.

Output: docs/images/gallery_chNN.png.  Run from the repository root:
    python scripts/make_readme_gallery.py        (needs Pillow)
"""
import base64, glob, io, json
import numpy as np
from PIL import Image, ImageOps

# chapter -> (row, column) of the panel in the chapter's 2-by-2 figure
PANELS = {4: (1, 0), 7: (1, 0), 8: (1, 0), 9: (1, 0), 10: (0, 0), 11: (1, 0)}

def gutter(img, axis):
    """A blank strip between the panels: the widest one near the centre along x (axis=0), the lowest one
    along y (axis=1), which lies between the x labels of the upper panels and the lower panels."""
    blank = (np.asarray(img.convert("L")) > 250).all(axis=axis)
    lo, hi = int(0.40 * blank.size), int(0.60 * blank.size)
    runs, start = [], None
    for i in range(lo, hi + 1):
        if blank[i] and start is None:
            start = i
        if start is not None and (not blank[i] or i == hi):
            runs.append((i - start, (start + i) // 2))
            start = None
    return max(runs)[1] if axis == 0 else runs[-1][1]


for chapter, (row, col) in PANELS.items():
    nb = json.load(open(glob.glob(f"notebooks/*/{chapter:02d}_*.ipynb")[0]))
    pngs = [o["data"]["image/png"] for c in nb["cells"] if c["cell_type"] == "code"
            for o in c.get("outputs", []) if "image/png" in o.get("data", {})]
    img = Image.open(io.BytesIO(base64.b64decode(pngs[-1]))).convert("RGB")
    xc, yc = gutter(img, axis=0), gutter(img, axis=1)        # blank strips between the panels
    box = (0 if col == 0 else xc, 0 if row == 0 else yc, xc if col == 0 else img.width, yc if row == 0 else img.height)
    panel = img.crop(box)
    panel = panel.crop(ImageOps.invert(panel).getbbox())      # trim the white margin ...
    panel = ImageOps.expand(panel, border=24, fill="white")   # ... and put back an even one
    panel = panel.resize((800, round(panel.height * 800 / panel.width)), Image.LANCZOS)
    panel.save(f"docs/images/gallery_ch{chapter:02d}.png", optimize=True)
print("wrote", len(PANELS), "gallery images")
