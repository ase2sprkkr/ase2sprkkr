from ..input_parameters import InputParameters


def test_xas_defines_upper_energy_bound_with_expected_default():
    parameters = InputParameters.create("XAS")

    assert parameters.ENERGY.EMAX() == 4.0
    assert parameters.ENERGY.EMAXEV() is None
