"""Renderer logic that needs neither VMD nor Tachyon: palette rules, scene-file lights, 16-bit image I/O, chain colours."""
import struct
import tempfile
import unittest
import zlib
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

from membraneforger import visualization as viz

SCENE_FILE = """Begin_Scene
Resolution 100 80
Camera
  Projection Orthographic
  Zoom 0.5
  Center  0 0 -2
End_Camera
Directional_Light Direction 0.1 -0.1 1 Color 1 1 1
Directional_Light Direction -1 -2 0.5 Color 1 1 1

Background 1 1 1
Sphere
"""


class Palette(unittest.TestCase):
    def test_fifteen_distinct_colours_with_their_own_ids(self):
        viz.check_palette()
        self.assertEqual(len(viz.PALETTE), 15)
        self.assertEqual(len({viz.COLOUR_IDS[name] for name in viz.PALETTE}), 15)
        used = [v[2] for v in viz.LIPID_TYPES.values()] + viz.CHAIN_COLOURS
        self.assertEqual(sorted(used), sorted(viz.PALETTE))

    def test_a_repeated_colour_is_refused(self):
        with mock.patch.dict(viz.PALETTE, {"Aqua Jelly": viz.PALETTE["Cyber Celadon"]}):
            with self.assertRaises(SystemExit):
                viz.check_palette()

    def test_a_colour_id_shared_with_an_element_is_refused(self):
        with mock.patch.dict(viz.COLOUR_IDS, {"water": viz.ELEMENT_COLOURS["O"][0]}):
            with self.assertRaises(SystemExit):
                viz.check_palette()


class SceneFile(unittest.TestCase):
    def test_lights_are_replaced_and_key_intensity_is_conserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            dat = Path(tmp) / "scene.dat"
            dat.write_text(SCENE_FILE)
            lights = viz.tune_scene_file(dat, viz.MODES["publication"], 100, 80)
            written = [line for line in dat.read_text().splitlines() if line.startswith("Directional_Light")]
        self.assertEqual(len(written), viz.MODES["publication"].key_lights + len(viz.LIGHTS) - 1)
        key = sum(light["intensity"] for light in lights if light["light"] == "key")
        self.assertAlmostEqual(key, viz.LIGHTS["key"][1], places=3)
        self.assertNotIn("Color 1 1 1", "\n".join(written[1:]))  # VMD's own lights are gone

    def test_an_unknown_scene_layout_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            dat = Path(tmp) / "scene.dat"
            dat.write_text("Begin_Scene\nSphere\n")
            with self.assertRaises(SystemExit):
                viz.tune_scene_file(dat, viz.MODES["preview"], 100, 80)


class DeepImages(unittest.TestCase):
    def test_ppm48_to_png16_keeps_every_sixteen_bit_value(self):
        pixels = (np.arange(4 * 5 * 3, dtype=np.uint32) * 1093 % 65536).astype(np.uint16).reshape(4, 5, 3)
        with tempfile.TemporaryDirectory() as tmp:
            ppm, png = Path(tmp) / "a.ppm", Path(tmp) / "a.png"
            ppm.write_bytes(b"P6\n5 4\n65535\n" + pixels.astype(">u2").tobytes())
            read = viz.read_ppm48(ppm)
            viz.write_png16(png, read, 600)
            data = png.read_bytes()
        self.assertTrue(np.array_equal(read, pixels))
        self.assertEqual(struct.unpack(">IIBB", data[16:26]), (5, 4, 16, 2))  # width, height, 16 bit, RGB
        start = data.index(b"IDAT") + 4
        rows = np.frombuffer(zlib.decompress(data[start:start + struct.unpack(">I", data[start - 8:start - 4])[0]]),
                             dtype=np.uint8).reshape(4, 31)
        undone = np.cumsum(rows[:, 1:].reshape(4, 5, 6).astype(np.uint16), axis=1).astype(np.uint8)  # undo the Sub filter
        self.assertTrue(np.array_equal(undone.reshape(4, 30).view(">u2").reshape(4, 5, 3), pixels))

    def test_an_eight_bit_ppm_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            ppm = Path(tmp) / "a.ppm"
            ppm.write_bytes(b"P6\n1 1\n255\n" + bytes(3))
            with self.assertRaises(SystemExit):
                viz.read_ppm48(ppm)


class ChainColours(unittest.TestCase):
    reference = [["ALA"] * 30 + ["GLY"] * 10, ["SER"] * 12, ["LEU"] * 25 + ["TRP"] * 5]

    def test_a_chain_keeps_its_colour_when_the_file_order_changes(self):
        scene = SimpleNamespace(reference_chains=self.reference)
        chains = [{"sequence": self.reference[k]} for k in (2, 0, 1)]
        self.assertEqual(viz.chain_colours(scene, chains), [viz.CHAIN_COLOURS[k] for k in (2, 0, 1)])

    def test_a_chain_unlike_any_reference_chain_is_refused(self):
        scene = SimpleNamespace(reference_chains=self.reference)
        with self.assertRaises(SystemExit):
            viz.chain_colours(scene, [{"sequence": ["PRO"] * 20}])

    def test_more_chains_than_colours_is_refused(self):
        scene = SimpleNamespace(reference_chains=None)
        with self.assertRaises(SystemExit):
            viz.chain_colours(scene, [{"sequence": ["ALA"]}] * (len(viz.CHAIN_COLOURS) + 1))


if __name__ == "__main__":
    unittest.main()
