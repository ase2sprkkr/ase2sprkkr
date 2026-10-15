import io
from collections.abc import Mapping, Sequence

import numpy as np
import pytest
from ase.build import bulk

from ..definitions.sections.lattice import LatticeSectionDefinition
from ..definitions.sections.sites import SitesSectionDefinition
from ..potential_definitions import PotentialDefinition
from ..potentials import Potential


def _assert_nested_almost_equal(actual, expected):
    if isinstance(expected, Mapping):
        assert actual.keys() == expected.keys()
        for key in expected:
            _assert_nested_almost_equal(actual[key], expected[key])
        return

    if isinstance(expected, np.ndarray):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
        return

    if isinstance(expected, Sequence) and not isinstance(expected, (str, bytes)):
        assert len(actual) == len(expected)
        for actual_item, expected_item in zip(actual, expected):
            _assert_nested_almost_equal(actual_item, expected_item)
        return

    if isinstance(expected, (float, np.floating)):
        assert actual == pytest.approx(expected, abs=1e-12)
        return

    assert actual == expected


@pytest.mark.parametrize("element", ["Co", "Ni", "U", "Ge", "As", "Na", "In", "Pr"])
def test_lattice_round_trip(element):
    lattice_definition = LatticeSectionDefinition("LATTICE")
    sites_definition = SitesSectionDefinition("SITES")
    definition = PotentialDefinition("pot", [lattice_definition, sites_definition])
    potential = Potential(definition=definition)
    lattice = potential.LATTICE

    atoms = bulk(element)
    potential.atoms = atoms

    output = io.StringIO()
    potential.save_to_file(output)
    output.seek(0)
    restored = definition.read_from_file(output)

    _assert_nested_almost_equal(restored.LATTICE.as_dict(), lattice.as_dict())
    assert restored.atoms is not None
    assert restored.atoms.cell is not atoms.cell
    np.testing.assert_allclose(restored.atoms.cell, atoms.cell, rtol=0, atol=1e-12)
    np.testing.assert_allclose(restored.atoms.positions, atoms.positions, rtol=0, atol=1e-12)
