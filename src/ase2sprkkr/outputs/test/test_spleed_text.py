"""The NOUT1 SPLEED records are not standalone .spc files."""

import numpy as np
import pytest

from ase2sprkkr.outputs.readers.spleed_text import SpleedTextSpectrum
from ase2sprkkr.outputs.readers.spec import SpecResult


def test_two_column_nout1_spectrum(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        "Unrelated startup info\n"
        " i(e)-diagram\n"
        " 0.11000D+02  0.10000D+00\n"
        " status = 1 2 3\n"
        " 0.12000D+02  0.20000D+00\n"
    )
    spectrum = SpleedTextSpectrum.from_file(filename)
    np.testing.assert_allclose(spectrum.ENERGY(), [11, 12])
    np.testing.assert_allclose(spectrum.INTENSITY(), [0.1, 0.2])
    with pytest.raises(AttributeError, match="not present"):
        spectrum.ASYMMETRY()


def test_seven_column_nout1_spectrum(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        " rotation diagram\n"
        " 0.11000E+02  0.45000E+02  0.50000E+00  "
        "0.10000E+00  0.20000E+00  0.30000E+00  0.40000E+00\n"
    )
    result = SpecResult(None, None, tmp_path, output_file=filename, task_name="SPLEED")
    spectrum = result.spleed_text_spectrum
    assert spectrum is not None
    np.testing.assert_allclose(spectrum.ENERGY(), [11])
    np.testing.assert_allclose(spectrum.THETA(), [45])
    np.testing.assert_allclose(spectrum.INTENSITY(), [0.5])
    np.testing.assert_allclose(spectrum.ASYMMETRY(), [0.1])
    np.testing.assert_allclose(spectrum.P_X(), [0.2])
    np.testing.assert_allclose(spectrum.P_Y(), [0.3])
    np.testing.assert_allclose(spectrum.P_Z(), [0.4])


def test_nout1_without_splout_diagram_is_not_spectrum(tmp_path):
    filename = tmp_path / "Fe_SPEC.out"
    filename.write_text(
        "unrelated numbers\n"
        " 0.11000D+02  0.10000D+00\n"
    )
    assert SpleedTextSpectrum.from_file(filename) is None
