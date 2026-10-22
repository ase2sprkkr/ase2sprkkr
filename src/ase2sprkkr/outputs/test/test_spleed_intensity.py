"""SPLEED intensity is parsed inside SPEC.out, not a separate ResultFile."""

from io import StringIO

import numpy as np

from ase2sprkkr.outputs.readers.spec import (
    SpecOutputParser,
    SpecResult,
    SpleedIntensity,
)
from ase2sprkkr.result_files import ResultFile


def test_two_column_log_is_parsed_from_main_output(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        "diagnostic\n i(e)-diagram\n"
        " 0.11000D+02   0.10000D+00\n"
        " 0.12000D+02   0.20000D+00\n"
        " finished SPLOUT\n"
    )
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    spectrum = result.spleed_intensity
    assert isinstance(spectrum, SpleedIntensity)
    np.testing.assert_allclose(spectrum.ENERGY(), [11, 12])
    np.testing.assert_allclose(spectrum.INTENSITY(), [.1, .2])
    assert spectrum.RAW_DATA().shape == (2, 2)
    assert "SPLEED_INTENSITY" not in result.files
    assert "SPLEED_INTENSITY" not in result.output_values
    assert "spleed_intensity" not in ResultFile.definitions


def test_seven_column_log_is_not_two_column_diagram(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        " rotation diagram\n"
        " 0.11000E+02  0.45000E+02  0.50000E+00  "
        "0.10000E+00  0.20000E+00  0.30000E+00  0.40000E+00\n"
    )
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    assert result.spleed_intensity is None


def test_no_diagram_header_is_not_matched():
    content = "unrelated\n 0.11000D+02   0.10000D+00\n"
    assert SpecOutputParser.parse_spleed_intensity(StringIO(content)) is None


def test_parser_stops_at_first_unrelated_nonempty_line():
    content = (
        " i(e)-diagram\n"
        " 0.11000D+02   0.10000D+00\n"
        " finished SPLOUT\n"
        " 0.14000D+02   0.99900D+00\n"
    )
    data = SpecOutputParser.parse_spleed_intensity(StringIO(content))
    np.testing.assert_allclose(data.RAW_DATA(), [[11, .1]])


def test_parser_handles_native_blank_line_after_heading():
    content = " i(e)-diagram\n\n 0.11000D+02   0.10000D+00\n"
    data = SpecOutputParser.parse_spleed_intensity(StringIO(content))
    np.testing.assert_allclose(data.ENERGY(), [11])


def test_live_output_can_be_read_after_new_rows_appear(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    assert result.spleed_intensity is None
    assert "SPLEED_INTENSITY" not in result.output_values
    filename.write_text(" i(e)-diagram\n 0.11000D+02   0.10000D+00\n")
    assert result.spleed_intensity is not None
    np.testing.assert_allclose(result.spleed_intensity.INTENSITY(), [.1])
    filename.write_text(" i(e)-diagram\n 0.11000D+02   0.20000D+00\n")
    np.testing.assert_allclose(result.spleed_intensity.INTENSITY(), [.2])
    assert "SPLEED_INTENSITY" not in result.output_values


def test_spleed_header_only_spc_is_not_a_spectrum(tmp_path):
    from ase2sprkkr.input_parameters.input_parameters import InputParameters
    p = InputParameters.create("SPLEED")
    p.CONTROL.DATASET = "Fe"
    p.SPEC_EL.TYP = 0
    result = SpecResult(p, None, tmp_path, output_file="Fe_SPEC.out")
    assert result.spc_filenames == ()
    assert "SPC" not in result.output_values


def test_standalone_log_uses_parser_to_reject_header_only_spc(tmp_path):
    (tmp_path / "Fe_data.spc").write_text("header only")
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(" i(e)-diagram\n 0.11000D+02   0.10000D+00\n")
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    assert result.spc_filenames == ()
