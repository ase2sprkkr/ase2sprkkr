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


def test_jxc_and_dmi_have_distinct_defaults_on_shared_schema():
    jxc = InputParameters.create("JXC")
    dmi = InputParameters.create("DMI")
    assert jxc.definition_name == dmi.definition_name == "JXC"
    assert jxc.task_definition is not dmi.task_definition
    assert jxc.requested_task_name == "JXC"
    assert dmi.requested_task_name == "DMI"
    assert jxc.task_name == "JXC"
    assert dmi.task_name == "DMI"
    assert jxc.TASK.JXC() and not jxc.TASK.DMI()
    assert not dmi.TASK.JXC() and dmi.TASK.DMI()
    assert jxc.TASK.TASK() == "JXC"
    assert dmi.TASK.TASK() == "DMI"


def test_adding_secondary_native_switch_keeps_primary_task_enabled():
    jxc = InputParameters.create("JXC")
    dmi = InputParameters.create("DMI")
    jxc.TASK.DMI = True
    dmi.TASK.JXC = True
    assert jxc.TASK.JXC() and jxc.TASK.DMI()
    assert dmi.TASK.JXC() and dmi.TASK.DMI()
    assert jxc.task_name == dmi.task_name == "JXC"
    assert jxc.requested_task_name == "JXC"
    assert dmi.requested_task_name == "DMI"


def test_native_dmi_only_and_combined_flags_roundtrip():
    for requested in ("DMI", "JXC"):
        params = InputParameters.create(requested)
        params.CONTROL.POTFIL = "Fe.pot"
        if requested == "JXC":
            params.TASK.DMI = True
        restored = _roundtrip(params)
        assert restored.definition_name == "JXC"
        assert restored.TASK.DMI()
        assert restored.TASK.JXC() == (requested == "JXC")
        assert restored.task_name == requested
        assert restored.task_definition is InputParameters.task_registry[requested]
        assert restored.TASK.TASK() == requested


def test_combined_native_flags_have_jxc_identity_regardless_of_origin():
    for requested in ("DMI", "JXC"):
        params = InputParameters.create(requested)
        params.CONTROL.POTFIL = "Fe.pot"
        params.TASK.DMI = True
        params.TASK.JXC = True
        restored = _roundtrip(params)
        assert restored.task_name == "JXC"
        assert restored.TASK.TASK() == "JXC"
        assert restored.task_definition is InputParameters.task_registry["JXC"]
        assert restored.TASK.JXC() and restored.TASK.DMI()


def test_changing_between_jxc_and_dmi_switches_defaults():
    params = InputParameters.create("JXC")
    params.change_task("DMI")
    assert params.requested_task_name == "DMI"
    assert params.task_name == "DMI"
    assert params.TASK.DMI() and not params.TASK.JXC()
    params.change_task("JXC")
    assert params.requested_task_name == "JXC"
    assert params.task_name == "JXC"
    assert params.TASK.JXC() and not params.TASK.DMI()
