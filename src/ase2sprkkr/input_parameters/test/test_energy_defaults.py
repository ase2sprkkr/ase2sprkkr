"""Shared energy defaults must keep task values, precedence, and copy semantics."""
import pickle

import pytest

from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.input_parameters.definitions.sections import ENERGY, TASK
from ase2sprkkr.input_parameters.input_parameters_definitions import InputParametersDefinition


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
