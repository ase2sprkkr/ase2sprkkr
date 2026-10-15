from ase2sprkkr.sprkkr.calculator import SPRKKR


def test_persistent_is_calculator_parameter(capsys):
    calculator = SPRKKR(persistent=True)

    assert calculator.persistent is True
    assert "unknown arguments" not in capsys.readouterr().out


def test_persistent_defaults_to_false():
    assert SPRKKR().persistent is False
