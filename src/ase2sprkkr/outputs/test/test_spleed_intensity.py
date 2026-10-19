"""SPLEED intensity is a structured ResultFile read from the SPEC.out log."""

import numpy as np
import pytest
from ase2sprkkr.result_files import ResultFile
from ase2sprkkr.result_files.definitions.spleed_intensity import SpleedIntensityNotFound, SpleedIntensityResultFile
from ase2sprkkr.outputs.readers.spec import SpecResult


def test_two_column_log_is_result_file(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        "diagnostic\n i(e)-diagram\n 0.11000D+02   0.10000D+00\n"
        " status = 1 2 3\n 0.12000D+02   0.20000D+00\n"
        " rotation diagram\n 0.13000D+02   0.90000D+00\n"
    )
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    spectrum = result.spleed_intensity
    assert isinstance(spectrum, SpleedIntensityResultFile)
    np.testing.assert_allclose(spectrum.ENERGY(), [11, 12])
    np.testing.assert_allclose(spectrum.INTENSITY(), [0.1, 0.2])
    assert spectrum.RAW_DATA().shape == (2, 2)
    assert result.output_values["SPLEED_INTENSITY"].file_type == "spleed_intensity"


def test_seven_column_output_is_not_two_column_format(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        " rotation diagram\n"
        " 0.11000E+02  0.45000E+02  0.50000E+00  "
        "0.10000E+00  0.20000E+00  0.30000E+00  0.40000E+00\n"
    )
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    assert result.spleed_intensity is None
    with pytest.raises(SpleedIntensityNotFound):
        ResultFile.from_file(filename, try_only="spleed_intensity")


def test_no_diagram_header_cannot_match(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text("unrelated\n 0.11000D+02   0.10000D+00\n")
    with pytest.raises(SpleedIntensityNotFound):
        ResultFile.from_file(filename, try_only="spleed_intensity")

def test_intensity_stops_before_unrelated_numeric_diagnostics(tmp_path):
    from ase2sprkkr.result_files.definitions.spleed_intensity import _read_diagram
    import numpy as np

    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        " i(e)-diagram\n"
        " 0.11000D+02   0.10000D+00\n"
        " finished SPLOUT\n"
        " 0.14000D+02   0.99900D+00\n"
    )
    np.testing.assert_allclose(_read_diagram(filename), [[11., .1]])


def test_intensity_can_be_refreshed_after_output_grows(tmp_path):
    from ase2sprkkr.outputs.readers.spec import SpecResult
    filename = tmp_path / "Fe_SPEC.out"
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    assert result.spleed_intensity is None
    assert "SPLEED_INTENSITY" not in result.output_values
    filename.write_text(" i(e)-diagram\n 0.11000D+02   0.10000D+00\n")
    assert result.spleed_intensity is not None
    assert result.output_values["SPLEED_INTENSITY"].file_type == "spleed_intensity"
    filename.write_text(" i(e)-diagram\n 0.11000D+02   0.20000D+00\n")
    np.testing.assert_allclose(result.spleed_intensity.INTENSITY(), [.2])


def test_embedded_result_is_read_only(tmp_path):
    import pytest
    from ase2sprkkr.result_files import ResultFile
    filename = tmp_path / "Fe_SPEC.out"
    content = " i(e)-diagram\n 0.11000D+02   0.10000D+00\n"
    filename.write_text(content)
    result = ResultFile.from_file(filename, try_only="spleed_intensity")
    assert result._definition.extension == "SPEC.out"
    with pytest.raises(TypeError, match="view of"):
        result.save_to_file(filename, atoms=False)
    assert filename.read_text() == content
