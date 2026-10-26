"""Shared energy defaults must keep task values, precedence, and copy semantics."""
import pickle

import pytest

from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.input_parameters.definitions.sections import ENERGY, TASK
from ase2sprkkr.input_parameters.input_parameters_definitions import InputParametersDefinition
from ase2sprkkr.common.warnings import DataValidityError


def test_custom_energy_defaults_use_the_same_precedence_as_arpes():
    definition = InputParametersDefinition('energy-test', [
        TASK('DOS'), ENERGY(emin=(.1, 'Minimum', -6.), emax=(.7, 'Maximum', 4.))])
    p = definition.create_object()
    assert p.ENERGY.EMINEV() == -6.
    assert p.ENERGY.EMAXEV() == 4.
    assert p.ENERGY.EMIN() is None
    copied = definition.copy().create_object()
    assert copied.ENERGY.EMINEV() == -6.
    assert copied.ENERGY.EMAXEV() == 4.
    p.ENERGY.set({'EMIN': .2, 'EMAX': .8})
    for candidate in (p, p.copy(copy_values=True), pickle.loads(pickle.dumps(p))):
        assert candidate.ENERGY.EMINEV() is None
        assert candidate.ENERGY.EMAXEV() is None
        assert candidate.ENERGY.EMIN() == .2
        candidate.ENERGY.set({'EMINEV': -3., 'EMAXEV': 2.})
        assert candidate.ENERGY.EMINEV() == -3.
        assert candidate.ENERGY.EMAXEV() == 2.
        candidate.ENERGY.set({'EMIN': None, 'EMAX': None, 'EMINEV': None, 'EMAXEV': None})
        assert candidate.ENERGY.EMINEV() == -6.
        assert candidate.ENERGY.EMAXEV() == 4.


@pytest.mark.parametrize('task,emin,emax', [('scf', -.2, None), ('dos', -.2, 1.),
                                        ('bsfek', -.2, -1.), ('bsfkk', .7, .7)])
def test_absolute_task_defaults_are_preserved(task, emin, emax):
    p = InputParameters.create(task)
    assert p.ENERGY.EMIN() == emin
    if emax is not None:
        assert p.ENERGY.EMAX() == emax
    p.ENERGY.EMINEV = -2.
    assert p.ENERGY.EMIN() is None
    p.ENERGY.EMINEV.clear()
    assert p.ENERGY.EMIN() == emin


@pytest.mark.parametrize('task,theta,phi,spol', [('ARPES', 0., 0., 2), ('SPLEED', 45., 270., 1)])
def test_spectroscopy_defaults_follow_reader(task, theta, phi, spol):
    p = InputParameters.create('arpes')
    p.TASK.TASK = task
    assert tuple(p.SPEC_EL.THETA()) == (theta, theta)
    assert tuple(p.SPEC_EL.PHI()) == (phi, phi)
    assert p.SPEC_EL.NT() == p.SPEC_EL.NP() == 1
    assert p.SPEC_EL.SPOL() == spol


@pytest.mark.parametrize("task", ("DOS", "XAS"))
@pytest.mark.parametrize("omitted", (("EMIN",), ("EMAX",), ("EMIN", "EMAX")))
def test_absolute_energy_bounds_are_optional_in_parsed_input(task, omitted):
    p = InputParameters.create(task)
    p.CONTROL.POTFIL = "Fe.pot"
    original = p.to_string()
    lines = original.splitlines()
    for name in omitted:
        assert any(line.strip().startswith(name + "=") for line in lines)

    missing_bounds = "\n".join(
        line for line in lines
        if not any(line.strip().startswith(name + "=") for name in omitted)
    ) + "\n"
    restored = p._definition.read_from_string(missing_bounds)

    # Omitting a keyword must not change the Python task defaults or their
    # serialization policy; SPRKKR's native omitted-value behavior is separate.
    for name in ("EMIN", "EMAX"):
        assert restored.ENERGY[name]() == p.ENERGY[name]()
    assert restored.to_string() == original


def test_derived_upper_bound_is_optional_too():
    energy = ENERGY(emax="emin")
    assert energy["EMIN"].is_optional is True
    assert energy["EMAX"].is_optional is True

    p = InputParameters.create("WFPLOT")
    p.CONTROL.POTFIL = "Fe.pot"
    without_upper = "\n".join(
        line for line in p.to_string().splitlines()
        if not line.strip().startswith("EMAX=")
    ) + "\n"
    restored = p._definition.read_from_string(without_upper)
    assert restored.ENERGY.EMAX() == restored.ENERGY.EMIN()


def test_optional_absolute_bounds_do_not_disable_relative_pair_validation():
    p = InputParameters.create("DOS")
    p.ENERGY.EMINEV = -3.0
    with pytest.raises(DataValidityError, match="must be supplied together"):
        p.validate("save")
