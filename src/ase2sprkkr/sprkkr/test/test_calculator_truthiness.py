import numpy as np
from ase import Atoms
from ase.build import bulk
import pytest

from ...potentials.definitions.sections.sites import SitesSection
from ...potentials.potentials import Potential
from .. import calculation
from ..atoms_region import AtomsRegion
from ..calculator import SPRKKR
from ..configuration import ConfigurationSectionTrait
from ..io_data import ReadIoData, WriteIoData
from ..sprkkr_atoms import SPRKKRAtoms


class _PotentialLike:
    def __init__(self, atoms):
        self.atoms = atoms


class _Result:
    task_name = "SCF"


class _Runner:
    def __init__(self, result):
        self.result = result

    def run(self):
        return self.result


class _ConfigurationTraitRecorder(ConfigurationSectionTrait):
    def _set_from_atoms(self, atoms, io_data):
        self.seen = atoms, io_data


class _Option:
    def __init__(self, value):
        self.value = value

    def __call__(self):
        return self.value


class _SitesStub:
    _container = object()

    def __init__(self):
        self.options = {
            "SCALED_ATOMIC_POSITIONS": _Option(np.empty((0, 3))),
            "BASSCALE": _Option(np.ones(3)),
        }

    def __getitem__(self, name):
        return self.options[name]


class _ReadData(dict):
    def __init__(self):
        super().__init__({"lattice.alat": 1.0})
        self.updated = []

    def update_atoms(self, atoms):
        self.updated.append(atoms)


def test_empty_atoms_are_not_treated_as_missing_by_calculator(monkeypatch):
    atoms = Atoms()
    calculator = SPRKKR(atoms=atoms)
    seen = []
    generated = object()

    def from_atoms(value):
        seen.append(value)
        return generated

    monkeypatch.setattr(Potential, "from_atoms", staticmethod(from_atoms))
    calculator._potential = None

    assert calculator.atoms is not None
    assert len(calculator.atoms) == 0
    assert calculator.potential is generated
    assert seen == [calculator.atoms]


def test_empty_atoms_are_explicit_when_potential_is_also_supplied():
    with pytest.raises(ValueError, match="both atoms and potential"):
        SPRKKR(atoms=Atoms(), potential=_PotentialLike(bulk("Fe")))


def test_potential_with_empty_atoms_replaces_calculator_atoms():
    calculator = SPRKKR(atoms=bulk("Fe"))
    atoms = Atoms()
    potential = _PotentialLike(atoms)

    calculator.potential = potential

    assert calculator.atoms is atoms


def test_explicit_empty_atoms_override_default_atoms(monkeypatch, tmp_path):
    explicit_atoms = Atoms()
    default_atoms = bulk("Fe")
    generated = object()
    seen = []

    def from_atoms(value):
        seen.append(value)
        return generated

    monkeypatch.setattr(Potential, "from_atoms", staticmethod(from_atoms))

    potential, atoms, potential_file, save_input = (
        calculation._resolve_potential_and_atoms(
            potential=None,
            atoms=explicit_atoms,
            potential_file=None,
            default_potential=None,
            default_atoms=default_atoms,
            restart_scf=False,
            directory=tmp_path,
        )
    )

    assert potential is generated
    assert atoms is explicit_atoms
    assert potential_file is None
    assert save_input is None
    assert seen == [explicit_atoms]


def test_empty_atoms_are_not_skipped_by_empty_sphere_policy(monkeypatch):
    atoms = Atoms()
    calls = []

    monkeypatch.setattr(
        calculation,
        "add_empty_spheres",
        lambda value, **kwargs: calls.append((value, kwargs)),
    )

    calculation._resolve_empty_spheres(
        input_parameters=None,
        atoms=atoms,
        empty_spheres=True,
        default_empty_spheres=False,
    )

    assert calls == [(atoms, {})]


def test_empty_atoms_are_used_for_filename_symbols(monkeypatch):
    atoms = Atoms()
    monkeypatch.setattr(calculation, "filename_from_symbols", lambda symbols: "empty")

    assert (
        calculation._symbols_for_filenames(
            atoms,
            input_parameters=None,
            input_file=None,
            potential_file=None,
            default_input_file=None,
            default_potential_file=None,
        )
        == "empty"
    )


def test_configuration_trait_preserves_empty_atoms_and_explicit_io_data():
    atoms = Atoms()
    io_data = WriteIoData(atoms)
    recorder = _ConfigurationTraitRecorder()

    returned = recorder.set_from_atoms(atoms, io_data)

    promoted_atoms, seen_io_data = recorder.seen
    assert promoted_atoms is atoms
    assert isinstance(promoted_atoms, SPRKKRAtoms)
    assert returned is io_data
    assert seen_io_data is io_data


def test_read_io_data_applies_handler_to_empty_atoms():
    atoms = Atoms()
    seen = []
    io_data = ReadIoData()

    def handler(value):
        seen.append(value)

    io_data.apply_on_atoms(handler, atoms)

    assert seen == [atoms]
    assert io_data._postponed == [handler]


def test_atoms_region_treats_empty_atoms_as_supplied(monkeypatch):
    atoms = Atoms()
    seen = []

    def set_atoms(region, value, add=True):
        seen.append((region, value, add))

    monkeypatch.setattr(AtomsRegion, "set_atoms", set_atoms)

    region = AtomsRegion("empty", slice(None), np.eye(3), atoms=atoms)

    assert seen == [(region, atoms, True)]


def test_sites_section_updates_explicit_empty_atoms():
    atoms = Atoms()
    read_data = _ReadData()

    result = SitesSection._update_atoms(_SitesStub(), atoms, read_data)

    assert result is None
    assert read_data.updated == []
    assert len(atoms) == 0


def test_result_potential_is_adopted_only_when_physical_state_was_omitted(
    monkeypatch,
):
    calculator = SPRKKR()
    result = _Result()
    adopted = []

    # The mocked runner and result update leave the calculator unchanged.
    # Exercise all three argument policies without reinitializing SPRKKR.
    monkeypatch.setattr(calculator, "create_run", lambda **kwargs: _Runner(result))
    monkeypatch.setattr(
        calculator,
        "_update_from_result",
        lambda value, *, adopt_potential: adopted.append((value, adopt_potential)),
    )

    for kwargs, expected in (
        ({}, True),
        ({"atoms": Atoms()}, False),
        ({"potential": False}, False),
    ):
        adopted.clear()
        assert calculator.calculate(**kwargs) is result
        assert adopted == [(result, expected)], kwargs
