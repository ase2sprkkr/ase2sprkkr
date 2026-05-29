# Test the SpecResult class without loading a file
from ase2sprkkr.outputs.readers.spec import SpecResult
import numpy as np

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)


class TestSpecResult(TestCase):
    def test_result(self):
        # Create a test class that inherits from SpecResult
        class TestSpecResult(SpecResult):
            def __init__(self, content):
                self._file_content = content
                self.output_file = None  # Needed by parent class

        # Sample content that matches the format of spec.out
        test_content = """
 parameter for potential barrier:
# ibar:  1
# epsx: 0.5000000E-02
# zparup: 0.0000000E+00 -0.1000000E+01 0.1000000E+01
# zpardn: 0.0000000E+00 -0.1000000E+01 0.1000000E+01
# bparp: 0.1340000E+01 -0.2000000E+01 0.1018000E+01
=======================
** start new input **
===================
lattice constant      :     4.6298300
real basis 1          :   4.62983   0.00000
real basis 2          :  -2.31492   4.00955
reciprocal basis 1    :   1.35711   0.78353
reciprocal basis 2    :   0.00000   1.56706

vacuum-wavevector of the photonfield:   0.00000000   0.00000000   0.00611428
vector potential  of the photonfield:
aa(i=1,3)=  (  1.00000000  0.00000000) (  0.00000000  0.00000000) (  0.00000000  0.00000000)
stookes-vector    of the photonfield  0.10000E+01  100.000%    0.000%    0.000%
1.8937E+00  -9.6825E-02  6.1567E-08 -3.9380E-04  3.0771E-08  3.0795E-08
1.9000E+00  -9.6825E-02  5.7837E-08 -4.2368E-04  2.8906E-08  2.8931E-08
"""

        # Create a test instance
        out = TestSpecResult(test_content)

        # Test the properties
        # print("Testing potential_barrier:")
        self.assertEqual(
            out.potential_barrier,
            {
                "ibar": 1,
                "epsx": 0.005,
                "zparup": [0.0, -1.0, 1.0],
                "zpardn": [0.0, -1.0, 1.0],
                "bparp": [1.34, -2.0, 1.018],
            },
        )
        self.assertEqual(out.lattice_constants, {"a": 4.62983})
        self.assertEqual(out.basis_vectors["real"], np.asarray([[4.62983, 0.0], [-2.31492, 4.00955]]))
        self.assertEqual(out.basis_vectors["reciprocal"], np.asarray([[1.35711, 0.78353], [0.0, 1.56706]]))
        self.assertEqual(
            out.spectral_data,
            np.asarray(
                [
                    [1.8937e00, -9.6825e-02, 6.1567e-08, -3.9380e-04, 3.0771e-08, 3.0795e-08],
                    [1.9000e00, -9.6825e-02, 5.7837e-08, -4.2368e-04, 2.8906e-08, 2.8931e-08],
                ]
            ),
        )
