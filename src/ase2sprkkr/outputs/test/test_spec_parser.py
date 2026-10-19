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
    # Unstructured diagnostic rows inside the text log are not a native SPC ResultFile.


def test_arpes_result_uses_native_spec_file_manifest(tmp_path):
    from ase2sprkkr.input_parameters.input_parameters import InputParameters
    from ase2sprkkr.outputs.readers.arpes import ArpesResult

    for task in ("ARPES", "CLPES"):
        # POL_P selects photon polarization, not the number of output files.
        # The latter is controlled by SPEC_PH.NPOL.
        for polarization in ("P", "S"):
            expected = "sample_data.spc"
            parameters = InputParameters.create(task)
            parameters.CONTROL.DATASET = "sample"
            parameters.SPEC_PH.POL_P = polarization
            result = ArpesResult(parameters, None, tmp_path)

            assert result.spc_filename == str(tmp_path / expected)
            assert result.output_values["SPC"].file_type == "spc"
            assert result.path_to("SPC") == str(tmp_path / expected)

            # CLPES only supports NPOL=1/2; dichroic outputs belong to
            # photoemission tasks that allow NPOL=0/3.
            if task == "ARPES":
                parameters.SPEC_PH.NPOL = 3
                dichroic = ArpesResult(parameters, None, tmp_path)
                assert dichroic.spc_filenames == (
                    str(tmp_path / "sample_data_pol1.spc"),
                    str(tmp_path / "sample_data_pol2.spc"),
                )

def test_spec_reader_requires_a_native_structured_band_file(tmp_path):
    import pytest
    import pyparsing as pp
    from ase2sprkkr.outputs.readers.spec import SpecResult

    (tmp_path / "Fe_data.spc").write_text(
        "# Legacy three-column diagnostic: not a native BAND SPC file\\n"
        "NT 1\\nNP 1\\n1.0D+00 2.0D+00 3.0D+00\\n"
    )
    result = SpecResult(None, None, tmp_path, output_file="Fe_SPEC.out", task_name="BAND")
    with pytest.raises((pp.ParseBaseException, ValueError)):
        _ = result.spectral_data


def test_arpes_registers_both_polarizations_from_native_files(tmp_path):
    from ase2sprkkr.outputs.readers.arpes import ArpesResult

    for suffix in ("_data_pol1.spc", "_data_pol2.spc"):
        (tmp_path / ("Fe" + suffix)).write_text("# SPC\n")
    result = ArpesResult(None, None, tmp_path, output_file="Fe_SPEC.out", task_name="ARPES")
    assert result.spc_filename == str(tmp_path / "Fe_data_pol1.spc")
    assert "SPC_2" in result.output_values
    assert result.path_to("SPC_2") == str(tmp_path / "Fe_data_pol2.spc")

def test_spectra_unpolarized_is_a_tuple_and_spc_is_deprecated(tmp_path, monkeypatch):
    import pytest
    from ase2sprkkr.outputs.readers.arpes import ArpesResult

    (tmp_path / "Fe_data.spc").write_text("# example\\n")
    result = ArpesResult(None, None, tmp_path, output_file="Fe_SPEC.out", task_name="ARPES")
    spectrum = object()
    monkeypatch.setattr(result, "parse_result_file", lambda file, potential=None: spectrum)
    assert result.spectra == (spectrum,)
    assert result.spectra is result.spectra
    assert "SPC" in result.output_values
    with pytest.warns(DeprecationWarning, match="result.spectra"):
        assert result.spc is spectrum


def test_spectra_includes_both_polarizations_in_result_values(tmp_path, monkeypatch):
    from ase2sprkkr.outputs.readers.arpes import ArpesResult
    from ase2sprkkr.outputs.readers.spec import SpecResult

    for suffix in ("_data_pol1.spc", "_data_pol2.spc"):
        (tmp_path / ("Fe" + suffix)).write_text("# example\\n")
    for result_class in (ArpesResult, SpecResult):
        result = result_class(None, None, tmp_path, output_file="Fe_SPEC.out", task_name="ARPES")
        monkeypatch.setattr(result, "parse_result_file", lambda file, potential=None: file.path())
        assert result.spectra == (
            str(tmp_path / "Fe_data_pol1.spc"),
            str(tmp_path / "Fe_data_pol2.spc"),
        )
        assert tuple(k for k in result.output_values if k.startswith("SPC")) == ("SPC", "SPC_2")
        assert result.output_values["SPC_2"].path() == str(tmp_path / "Fe_data_pol2.spc")
