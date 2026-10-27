"""Behaviour of the show-structure CLI without opening GUI windows."""

import argparse

import pytest

from ..commands import show_structure
from ...ase import visualise
from ...potentials.potentials import Potential
from ...sprkkr import structure


@pytest.fixture
def structures(monkeypatch):
    writes = []
    views = []
    conversions = []

    class FakeAtoms:
        def __init__(self, name):
            self.name = name

        def write(self, filename, format):
            writes.append((self.name, str(filename), format))

    potential_atoms = FakeAtoms("potential")
    input_atoms = FakeAtoms("input")

    class FakePotential:
        atoms = potential_atoms

    monkeypatch.setattr(
        Potential, "from_file", staticmethod(lambda filename: FakePotential())
    )

    def convert(filename, potential, n_bulk, vacuum_height):
        conversions.append((filename, potential, n_bulk, vacuum_height))
        return input_atoms

    monkeypatch.setattr(structure, "structure_file_to_atoms", convert)
    monkeypatch.setattr(
        visualise, "view",
        lambda atoms, scale_radii, blocking=False, close_after=None: views.append((atoms.name, scale_radii)),
    )
    return writes, views, conversions


def run(*arguments):
    parser = argparse.ArgumentParser()
    show_structure.parser(parser)
    show_structure.run(parser.parse_args(["test.pot", *arguments]), {})
    

@pytest.mark.parametrize(
    "arguments, expected_writes, expected_views",
    [
        ([], [], ["potential"]),
        (["-i", "struct.inp"], [], ["input"]),
        (["-a"], [], ["potential"]),
        (["-i", "struct.inp", "-a"], [], ["input"]),
        (["-i", "struct.inp", "-A"], [], ["potential"]),
        (["-i", "struct.inp", "-a", "-A"], [], ["potential", "input"]),
        (["-a", "-A"], [], ["potential"]),
        (["-o", "result.cif"], [("potential", "result.cif", "cif")], []),
        (["-i", "struct.inp", "-o", "result.cif"], [("input", "result.cif", "cif")], []),
        (["-O", "potential.xyz"], [("potential", "potential.xyz", "cif")], []),
        (["-i", "struct.inp", "-O", "potential.xyz"],
         [("potential", "potential.xyz", "cif")], ["input"]),
        (["-i", "struct.inp", "-o", "out.cif", "-O"],
         [("input", "out.cif", "cif"), ("potential", "out_pot.cif", "cif")], []),
        (["-o", "out.xyz", "-O", "-f", "xyz"],
         [("potential", "out.xyz", "xyz"), ("potential", "out_pot.xyz", "xyz")], []),
        (["-a", "-o", "out.cif"], [("potential", "out.cif", "cif")], ["potential"]),
        (["-i", "struct.inp", "-a", "-O", "pot.cif"],
         [("potential", "pot.cif", "cif")], ["input"]),
        (["--no-view"], [], []),
        (["-i", "struct.inp", "-a", "-A", "--no-view"], [], []),
    ],
)
def test_show_structure_decisions(structures, arguments, expected_writes, expected_views):
    writes, views, conversions = structures
    run(*arguments)
    assert writes == expected_writes
    assert [name for name, scale in views] == expected_views
    assert len(conversions) == int("-i" in arguments)


def test_show_structure_passes_builder_and_viewer_options(structures):
    _, views, conversions = structures
    run("-i", "struct.inp", "-b", "4", "-v", "12.5", "-s", "0.8")
    assert conversions[0][0] == "struct.inp"
    assert conversions[0][2:] == (4, 12.5)
    assert views == [("input", 0.8)]


def test_unnamed_potential_output_requires_primary_output(structures):
    with pytest.raises(ValueError, match="requires --output"):
        run("-O")
    writes, views, conversions = structures
    assert writes == views == conversions == []


def test_show_structure_does_not_generate_temp_files(monkeypatch, structures):
    def unexpected_temp(*args, **kwargs):
        pytest.fail("show-structure must not allocate temporary files")
    monkeypatch.setattr("tempfile.NamedTemporaryFile", unexpected_temp)
    run("-i", "struct.inp", "-a")
