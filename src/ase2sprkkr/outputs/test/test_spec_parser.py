import numpy as np

from ase2sprkkr.outputs.readers.spec import SpecResult


def test_spec_result_properties():
    class InMemorySpecResult(SpecResult):
        def __init__(self, content):
            self._file_content = content
            self.output_file = None

    content = """
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

    result = InMemorySpecResult(content)

    assert result.potential_barrier == {
        "ibar": 1,
        "epsx": 0.005,
        "zparup": [0.0, -1.0, 1.0],
        "zpardn": [0.0, -1.0, 1.0],
        "bparp": [1.34, -2.0, 1.018],
    }
    assert result.lattice_constants == {"a": 4.62983}
    np.testing.assert_allclose(
        result.basis_vectors["real"],
        [[4.62983, 0.0], [-2.31492, 4.00955]],
    )
    np.testing.assert_allclose(
        result.basis_vectors["reciprocal"],
        [[1.35711, 0.78353], [0.0, 1.56706]],
    )
    np.testing.assert_allclose(
        result.spectral_data,
        [
            [1.8937e00, -9.6825e-02, 6.1567e-08, -3.9380e-04, 3.0771e-08, 3.0795e-08],
            [1.9000e00, -9.6825e-02, 5.7837e-08, -4.2368e-04, 2.8906e-08, 2.8931e-08],
        ],
    )


def test_arpes_result_uses_native_spec_file_manifest(tmp_path):
    from ase2sprkkr.input_parameters.input_parameters import InputParameters
    from ase2sprkkr.outputs.readers.arpes import ArpesResult

    for task in ("ARPES", "CLPES"):
        for polarization, expected in (
            ("P", "sample_data.spc"),
            ("S", "sample_data_pol1.spc"),
        ):
            parameters = InputParameters.create(task)
            parameters.CONTROL.DATASET = "sample"
            parameters.SPEC_PH.POL_P = polarization
            result = ArpesResult(parameters, None, tmp_path)

            assert result.spc_filename == str(tmp_path / expected)
            assert result.output_values["SPC"].file_type == "spc"
            assert result.path_to("SPC") == str(tmp_path / expected)
