if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import unittest
from .. import es_finder
from ase.build import bulk
from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
import numpy as np

class TestEsFinder(TestCase):

  def test_es_finder(self):

       if not es_finder.is_enabled:
         raise unittest.SkipTest("ES finder cannot be imported, skipping the test")

       atoms = bulk('Co')
       self.assertEqual(0, len(es_finder.empty_spheres(atoms).positions))

       atoms = bulk('NaCl',  "rocksalt", a=5.64)
       out=es_finder.empty_spheres(atoms)
       self.assertEqual(12, len(out.positions))
