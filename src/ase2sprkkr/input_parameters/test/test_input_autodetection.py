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



def test_jxc_combined_selectors_separated_by_other_task_options():
    # The unordered TASK section must recognize both native keywords, even
    # when a non-selector parameter appears between them.
    from ..input_parameters_definitions import _task_preparse_definition

    preview = _task_preparse_definition()
    parsed = preview.read_from_string("TASK DMI CLURAD=4.0 JXC")
    assert list(parsed.TASK.TASK()) == ["DMI", "JXC"]


def test_alias_is_exactly_the_same_task_and_defaults():
    jxc = InputParameters.create("JXC")
    dmi = InputParameters.create("DMI")
    assert jxc.definition_name == dmi.definition_name == "JXC"
    assert jxc.task_definition is dmi.task_definition
    assert jxc.task_name == dmi.task_name == "JXC"
    assert jxc.requested_task_name == dmi.requested_task_name == "JXC"
    assert jxc.TASK.JXC() and dmi.TASK.JXC()
    assert not jxc.TASK.DMI() and not dmi.TASK.DMI()
    assert jxc.TASK.TASK() == dmi.TASK.TASK() == "JXC"


def test_native_dmi_only_and_combined_flags_are_not_task_aliases():
    for dmi_only in (True, False):
        params = InputParameters.create("JXC")
        params.CONTROL.POTFIL = "Fe.pot"
        params.TASK.DMI = True
        if dmi_only:
            params.TASK.JXC = False
        restored = _roundtrip(params)
        assert restored.definition_name == "JXC"
        assert restored.TASK.DMI()
        assert restored.TASK.JXC() is not dmi_only
        assert restored.task_name == "JXC"
        assert restored.TASK.TASK() == "JXC"
        assert restored.task_definition is InputParameters.task_registry["DMI"]
