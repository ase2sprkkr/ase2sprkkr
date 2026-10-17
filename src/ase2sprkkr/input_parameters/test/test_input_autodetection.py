import io

import numpy as np

from ..input_parameters import InputParameters


def _roundtrip(parameters):
    stream = io.StringIO(parameters.to_string(validate=True))
    return InputParameters.from_file(stream)


def test_from_file_detects_bsfek_from_task_selector():
    parameters = InputParameters.create("BSFEK")
    parameters.CONTROL.POTFIL = "Fe.pot"
    parameters.TASK.KPATH = 6

    parsed = _roundtrip(parameters)

    assert parsed.definition_name == "BSF"
    assert parsed.task_name == "BSFEK"
    assert parsed.CONTROL.ADSI() == "BLOCHSF"
    assert parsed.TASK.KPATH() == "6"


def test_from_file_detects_bsfkk_from_task_selector():
    parameters = InputParameters.create("BSFKK")
    parameters.CONTROL.POTFIL = "Fe.pot"
    parameters.TASK.NK1 = 45
    parameters.TASK.NK2 = 55
    parameters.TASK.K1 = [1.0, 0.0, 0.0]
    parameters.TASK.K2 = [0.0, 1.0, 0.0]
    parameters.TASK.KA = [[0.0, 0.0, 0.0]]

    parsed = _roundtrip(parameters)

    assert parsed.definition_name == "BSF"
    assert parsed.task_name == "BSFKK"
    assert parsed.CONTROL.ADSI() == "BLOCHSF"
    assert parsed.TASK.NK1() == 45
    assert parsed.TASK.NK2() == 55
    np.testing.assert_array_equal(parsed.TASK.KA(), [[0.0, 0.0, 0.0]])


def test_from_file_without_task_section_uses_scf_schema():
    text = """CONTROL DATASET=Fe
        ADSI=CUSTOM
        POTFIL=Fe.pot
TAU BZINT=POINTS
    NKTAB=250
ENERGY GRID={5}
    NE={32}
SCF NITER=80
    MIX=0.2
SITES NL={3}
"""
    parsed = InputParameters.from_file(io.StringIO(text))
    assert parsed.definition_name == "SCF"
    assert parsed.task_name == "SCF"
    assert parsed.CONTROL.ADSI() == "CUSTOM"
    assert parsed.SCF.NITER() == 80


def test_from_file_jxc_example_with_preceding_dmi_flag():
    from pathlib import Path

    path = Path(__file__).parent.parent / "examples" / "jxc.in"
    parsed = InputParameters.from_file(path)
    assert parsed.definition_name == "JXC"
    assert parsed.TASK.TASK() == "JXC"
    assert parsed.TASK.DMI()


def test_jxc_combined_native_selectors_choose_jxc_input_definition():
    from pathlib import Path

    path = Path(__file__).parent.parent / "examples" / "jxc.in"
    params = InputParameters.from_file(path)
    assert params.definition_name == "JXC"
    assert params.TASK.TASK() == "JXC"
    assert params.TASK.DMI()
