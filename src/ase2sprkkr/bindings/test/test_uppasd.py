if __package__:
    from .init_tests import patch_package
else:
    from ase2sprkkr.bindings.test.init_tests import patch_package
__package__, __name__ = patch_package(__package__, __name__)

from pathlib import Path
import tempfile

import numpy as np
from ase import Atoms

from .. import uppasd as uppasd_bindings
from ..uppasd import (
    DMFILE_FORMAT,
    JFILE_FORMAT,
    MOMFILE_FORMAT,
    POSFILE_FORMAT,
    write_dmfile,
    write_jfile,
    write_mom_file,
    write_pos_file,
)
from ...output_files.definitions.jxc import Coordinates, JXCOutputFile, Selector


class FakeMoments:
    def __init__(self, spin_moment):
        self.spin_moment = spin_moment


class FakeAtomicType:
    def __init__(self, symbol="Fe", spin_moment=None):
        self.symbol = symbol
        self.moments = FakeMoments(spin_moment) if spin_moment is not None else None


class FakeSite:
    def __init__(self, site_type, occupation=None):
        self.site_type = site_type
        self.occupation = occupation or {site_type: 1.0}


class FakeJxcOutput:
    def __init__(self, rows=None, selected_iqs=None, selector=None):
        self.rows = rows
        self.selected_iqs = selected_iqs
        self.selector = selector
        self.filtered_calls = []
        self.iq_selector_calls = []
        self.create_selector_calls = []

    def filtered_data(self, **kwargs):
        self.filtered_calls.append(kwargs)
        return self.rows

    def iq_selector(self, **kwargs):
        self.iq_selector_calls.append(kwargs)
        return self.selected_iqs

    def create_selector(self, **kwargs):
        self.create_selector_calls.append(kwargs)
        return self.selector


J_DATA_DTYPE = np.dtype(
    [
        ("IT", int),
        ("IQ", int),
        ("JT", int),
        ("JQ", int),
        ("N1", int),
        ("N2", int),
        ("N3", int),
        ("DRX", float),
        ("DRY", float),
        ("DRZ", float),
        ("DR", float),
        ("JXX", float),
        ("JYY", float),
        ("JXY", float),
        ("JYX", float),
    ]
)

DM_DATA_DTYPE = np.dtype(
    [
        ("IT", int),
        ("IQ", int),
        ("JT", int),
        ("JQ", int),
        ("N1", int),
        ("N2", int),
        ("N3", int),
        ("DRX", float),
        ("DRY", float),
        ("DRZ", float),
        ("DR", float),
        ("DX", float),
        ("DY", float),
        ("DZ", float),
    ]
)


def make_atoms():
    atoms = Atoms(
        symbols=["Fe", "Ni", "Fe"],
        positions=[(0.0, 0.1, 0.2), (1.0, 1.1, 1.2), (2.0, 2.1, 2.2)],
        cell=np.eye(3),
        pbc=True,
    )
    shared_type = FakeAtomicType("Fe", 2.5)
    second_type = FakeAtomicType("Ni", 1.5)
    atoms.sites = [FakeSite(shared_type), FakeSite(second_type), FakeSite(shared_type)]
    return atoms


def test_write_jfile_uses_selector_and_lattice_coordinates():
    rows = np.array([(1, 1, 1, 2, 3, 4, 5, 0.1, 0.2, 0.3, 0.75, 6.25, 0.0, 0.0, 0.0)], dtype=J_DATA_DTYPE)
    selector = Selector(iq=[1], it={1})
    output = FakeJxcOutput(rows=rows)

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "jfile.dat"
        written = write_jfile(output, path, selector=selector, coordinates=Coordinates.lattice)

        assert written is True
        assert output.filtered_calls == [
            {
                "selector": selector,
                "iq": None,
                "it": None,
                "exclude_it": None,
                "exclude_vc": True,
                "exchange_radius": None,
            }
        ]
        assert path.read_text() == JFILE_FORMAT.format(1, 2, 3, 4, 5, 6.25, 0.75)


def test_write_dmfile_uses_cartesian_coordinates():
    rows = np.array([(1, 2, 1, 3, 7, 8, 9, 0.25, 0.5, 0.75, 1.5, 2.5, 3.5, 4.5)], dtype=DM_DATA_DTYPE)
    selector = Selector(iq=[2], it={1})
    output = FakeJxcOutput(rows=rows)

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "dmfile.dat"
        written = write_dmfile(output, path, selector=selector, coordinates=Coordinates.cartesian)

        assert written is True
        assert path.read_text() == DMFILE_FORMAT.format(2, 3, 0.25, 0.5, 0.75, 2.5, 3.5, 4.5, 1.5)


def test_write_pos_file_uses_selected_iqs(monkeypatch):
    atoms = make_atoms()
    selector = Selector(iq=[1, 3], it={1})
    output = FakeJxcOutput(selected_iqs=[1, 3])
    monkeypatch.setattr(uppasd_bindings.SPRKKRAtoms, "promote_ase_atoms", staticmethod(lambda atoms: atoms))

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "posfile.dat"
        written = write_pos_file(atoms, path, selector=selector, jxc_ouput_file=output)

        assert written is True
        assert output.iq_selector_calls == [
            {"selector": selector, "iq": None, "it": None, "exclude_it": None, "exclude_vc": True}
        ]
        assert path.read_text() == (
            POSFILE_FORMAT.format(1, 1, 0.0, 0.1, 0.2) + POSFILE_FORMAT.format(3, 1, 2.0, 2.1, 2.2)
        )


def test_write_mom_file_uses_selector_and_atomic_type_order(monkeypatch):
    atoms = make_atoms()
    selector = Selector(iq=..., it=...)
    output = FakeJxcOutput(selector=selector)
    monkeypatch.setattr(uppasd_bindings.SPRKKRAtoms, "promote_ase_atoms", staticmethod(lambda atoms: atoms))

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "momfile.dat"
        written = write_mom_file(atoms, path, selector=selector, jxc_ouput_file=output)

        assert written is True
        assert output.create_selector_calls == [
            {"selector": selector, "iq": None, "it": None, "exclude_it": None, "exclude_vc": True}
        ]
        assert path.read_text() == (
            MOMFILE_FORMAT.format(1, 1, 2.5, 0.0, 0.0, 1.0)
            + MOMFILE_FORMAT.format(3, 1, 2.5, 0.0, 0.0, 1.0)
            + MOMFILE_FORMAT.format(2, 1, 1.5, 0.0, 0.0, 1.0)
        )


def test_jxc_output_file_selector_excludes_vacuum_sites(monkeypatch):
    atoms = Atoms(
        symbols=["Fe", "X", "Ni"],
        positions=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
        cell=np.eye(3),
        pbc=True,
    )
    fe = FakeAtomicType("Fe", 2.0)
    vc = FakeAtomicType("Vc")
    ni = FakeAtomicType("Ni", 1.0)
    atoms.sites = [FakeSite(fe), FakeSite(vc), FakeSite(ni)]
    monkeypatch.setattr(uppasd_bindings.SPRKKRAtoms, "promote_ase_atoms", staticmethod(lambda atoms: atoms))

    jxc_output = JXCOutputFile.from_atoms(atoms)
    selector = jxc_output.create_selector(it=[1, 2, 3], exclude_it=set(), exclude_vc=True)

    assert selector.it == {1, 3}
    assert jxc_output.iq_selector(selector=selector) == {1, 3}
