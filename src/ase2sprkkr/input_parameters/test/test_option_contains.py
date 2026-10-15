from ase2sprkkr.input_parameters.input_parameters import InputParameters


def test_repeated_option_contains_only_explicit_occurrences():
    parameters = InputParameters.create("scf")
    option = parameters.MODE.SOC

    assert option() == option.default_value
    assert "def" not in option
    assert 2 not in option

    option.set({"def": option.default_value, 2: 0.4})

    assert "def" in option
    assert 2 in option
    assert 3 not in option

    option[2] = None

    assert "def" in option
    assert 2 not in option
