"""The SPLOUT intensity diagram and SPLEED SPC are distinct formats."""

import numpy as np
import pytest

from ase2sprkkr.outputs.readers.spleed_intensity import SpleedIntensity
from ase2sprkkr.outputs.readers.spec import SpecResult


def test_two_column_intensity_is_extracted_from_spec_output(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        "Other KKRSPEC diagnostics\n"
        " i(e)-diagram\n"
        " 0.11000D+02  0.10000D+00\n"
        " status = 1 2 3\n"
        " 0.12000D+02  0.20000D+00\n"
    )
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    data = result.spleed_intensity
    assert isinstance(data, SpleedIntensity)
    np.testing.assert_allclose(data.ENERGY(), [11, 12])
    np.testing.assert_allclose(data.INTENSITY(), [0.1, 0.2])
    assert data.RAW_DATA().shape == (2, 2)
    assert not hasattr(data, "ASYMMETRY")


def test_seven_column_output_is_not_intensity_format(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        " rotation diagram\n"
        " 0.11000E+02  0.45000E+02  0.50000E+00  "
        "0.10000E+00  0.20000E+00  0.30000E+00  0.40000E+00\n"
    )
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    assert result.spleed_intensity is None


def test_no_intensity_without_native_diagram_header(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text("unrelated\n 0.11000D+02  0.10000D+00\n")
    assert SpleedIntensity.from_output(filename) is None


def test_non_spleed_result_does_not_expose_spleed_data(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(" i(e)-diagram\n 0.11000D+02  0.10000D+00\n")
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="BAND")
    assert result.spleed_intensity is None


def test_intensity_rejects_seven_column_arrays():
    with pytest.raises(ValueError, match="exactly two columns"):
        SpleedIntensity(np.ones((2, 7)))
