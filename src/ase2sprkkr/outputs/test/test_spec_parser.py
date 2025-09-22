# Test the SpecResult class without loading a file
from ase2sprkkr.outputs.readers.spec import SpecResult

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
test_result = TestSpecResult(test_content)

# Test the properties
print("Testing potential_barrier:")
print(test_result.potential_barrier)
print("\nTesting lattice_constants:")
print(test_result.lattice_constants)
print("\nTesting basis_vectors:")
for name, vector in test_result.basis_vectors.items():
    print(f"{name}: {vector}")
print("\nTesting spectral_data:")
print(test_result.spectral_data)