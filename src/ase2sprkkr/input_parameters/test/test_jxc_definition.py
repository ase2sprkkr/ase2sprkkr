from ..input_parameters import InputParameters


def test_jxc_dmi_is_a_supported_flag():
    parameters = InputParameters.create("JXC")

    assert parameters.TASK.DMI() is False

    parameters.TASK.DMI = True
    output = parameters.to_string(validate=False)

    assert "DMI" in output
    parsed = InputParameters.definition("JXC").read_from_string(output)
    assert parsed.TASK.DMI() is True
