import io

from ..input_parameters import InputParameters


def _roundtrip(parameters):
    stream = io.StringIO(parameters.to_string(validate=True))
    return InputParameters.from_file(stream)


def test_from_file_detects_bsfek_shared_schema_by_adsi():
    parameters = InputParameters.create("BSFEK")
    parameters.CONTROL.POTFIL = "Fe.pot"
    parameters.TASK.KPATH = 6

    parsed = _roundtrip(parameters)

    assert parsed.definition_name == "BSF"
    assert parsed.task_name == "BSFEK"
    assert parsed.CONTROL.ADSI() == "BLOCHSF"
    assert parsed.TASK.KPATH() == "6"


def test_from_file_detects_bsfkk_shared_schema_by_adsi():
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
    assert parsed.TASK.KA() == [0.0, 0.0, 0.0]
