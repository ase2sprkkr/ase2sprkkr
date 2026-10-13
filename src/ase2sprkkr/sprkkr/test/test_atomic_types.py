from types import SimpleNamespace

import pytest

from ase2sprkkr.sprkkr.atomic_types import AtomicType


def make_atom():
    return AtomicType(
        "Fe",
        n_electrons=26,
        n_core=18,
        n_valence=8,
        n_semicore=0,
    )


def test_explicit_n_core_is_returned():
    atom = make_atom()
    assert atom.n_core == 18


@pytest.mark.parametrize(
    ("attribute", "value"),
    [
        ("n_electrons", 26),
        ("n_core", 18),
        ("n_valence", 8),
        ("n_semicore", 0),
    ],
)
def test_electron_count_setters_accept_consistent_values(attribute, value):
    atom = make_atom()

    setattr(atom, attribute, value)

    assert getattr(atom, attribute) == value


@pytest.mark.parametrize(
    ("attribute", "invalid_value", "expected_value"),
    [
        ("n_electrons", 25, 26),
        ("n_core", 17, 18),
        ("n_valence", 7, 8),
        ("n_semicore", 1, 0),
    ],
)
def test_electron_count_setters_reject_without_mutating_state(
    attribute,
    invalid_value,
    expected_value,
):
    atom = make_atom()

    with pytest.raises(ValueError):
        setattr(atom, attribute, invalid_value)

    assert getattr(atom, attribute) == expected_value
    assert (atom.n_electrons, atom.n_core, atom.n_valence, atom.n_semicore) == (26, 18, 8, 0)


def test_invalid_n_electrons_does_not_update_potential_z():
    atom = make_atom()
    atom._potential = SimpleNamespace(z=26)

    with pytest.raises(ValueError):
        atom.n_electrons = 25

    assert atom.n_electrons == 26
    assert atom._potential.z == 26
