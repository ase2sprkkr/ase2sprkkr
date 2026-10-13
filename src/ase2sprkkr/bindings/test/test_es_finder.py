import unittest
from ase.build import bulk

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from .. import es_finder  # NOQA: E402


class TestEsFinder(TestCase):
    def test_radii_ratios_map_is_copied(self):
        original = dict(es_finder.winger_seitz_radii)

        radii = es_finder._prepare_radii_ratios_map(None)
        radii["DefinitelyNotAnElement"] = 1.0

        self.assertEqual(original, es_finder.winger_seitz_radii)
        self.assertIsNot(radii, es_finder.winger_seitz_radii)

    def test_radii_ratios_map_overrides_defaults(self):
        radii = es_finder._prepare_radii_ratios_map({"Co": 123.0})

        self.assertEqual(123.0, radii["Co"])
        self.assertEqual(es_finder.winger_seitz_radii.keys(), radii.keys())

    def test_es_finder(self):
        if not es_finder.is_enabled:
            raise unittest.SkipTest("ES finder cannot be imported, skipping the test")

        atoms = bulk("Co")
        self.assertEqual(0, len(es_finder.empty_spheres(atoms).positions))

        atoms = bulk("NaCl", "rocksalt", a=5.64)
        out = es_finder.empty_spheres(atoms)
        self.assertEqual(12, len(out.positions))
