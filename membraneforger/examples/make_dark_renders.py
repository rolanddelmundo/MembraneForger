#!/usr/bin/env python
"""Make dark-mode copies of the docs renders by replacing their white background.

    python membraneforger/examples/make_dark_renders.py IMAGE [IMAGE ...] [--bg 1e2129] [--min-area 25] [--outside-only]

Each IMAGE.jpg is written as IMAGE_dark.jpg. White areas larger than --min-area pixels (the background and the
gaps between helices) become the background colour; small white areas are specular highlights on the molecules
and stay. A 1-pixel band around the molecules is un-mixed from white (colour to alpha), so edges keep their
anti-aliasing without a white fringe. --outside-only replaces only the white that touches the image border; use it
for a render with a see-through water box, where the white gaps between lipids belong to the box. The default colour
is the docs site's dark-mode page background.
"""
import argparse
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage


def darken(path, bg, min_area=25, white=246, band_px=1, outside_only=False):
    rgb = np.asarray(Image.open(path).convert("RGB")).astype(float)
    near_white = rgb.min(axis=2) >= white
    labels, n = ndimage.label(near_white)
    sizes = ndimage.sum(near_white, labels, index=np.arange(1, n + 1))
    large = np.zeros(n + 1, bool)
    large[1:] = sizes >= min_area
    if outside_only:
        large[:] = False
        large[np.unique(np.r_[labels[0], labels[-1], labels[:, 0], labels[:, -1]])] = True
        large[0] = False
    background = large[labels]
    band = ndimage.binary_dilation(background, iterations=band_px) & ~background
    alpha = np.ones(rgb.shape[:2])
    alpha[background] = 0
    alpha[band] = ((255 - rgb[band]) / 255).max(axis=1).clip(0, 1)
    colour = rgb.copy()
    a = alpha[band][:, None]
    colour[band] = np.where(a > 0, (rgb[band] - 255 * (1 - a)) / np.where(a > 0, a, 1), 0).clip(0, 255)
    out = colour * alpha[..., None] + np.asarray(bg, float) * (1 - alpha[..., None])
    return Image.fromarray(out.clip(0, 255).astype(np.uint8))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("images", nargs="+", type=Path)
    parser.add_argument("--bg", default="1e2129", help="background colour as hex (default: the site's dark page colour)")
    parser.add_argument("--min-area", type=int, default=25, help="white areas at least this many pixels are background")
    parser.add_argument("--outside-only", action="store_true", help="replace only the white that touches the image border")
    parser.add_argument("--quality", type=int, default=85)
    args = parser.parse_args()
    bg = tuple(int(args.bg.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4))
    for path in args.images:
        target = path.with_name(f"{path.stem}_dark.jpg")
        darken(path, bg, args.min_area, outside_only=args.outside_only).save(target, quality=args.quality, optimize=True)
        print("wrote", target)


if __name__ == "__main__":
    main()
