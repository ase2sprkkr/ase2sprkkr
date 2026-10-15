import pickle

from ..calculator import SPRKKR
from ..sprkkr_atoms import SPRKKRAtoms


class TestSPRKKRPickle:
    def test_calculator_pickles(self):
        pickle.dumps(SPRKKR())

    def test_atoms_roundtrip(self):
        atoms = SPRKKRAtoms("H")
        atoms.calc = SPRKKR()

        restored = pickle.loads(pickle.dumps(atoms))

        assert restored.symbols[0] == "H"
        assert isinstance(restored.calc, SPRKKR)

    def test_atoms_with_potential_roundtrip(self):
        atoms = SPRKKRAtoms(
            "Fe",
            positions=[[0.0, 0.0, 0.0]],
            cell=[2.87, 2.87, 2.87],
            pbc=True,
        )
        original_potential = atoms.potential

        restored = pickle.loads(pickle.dumps(atoms))

        assert restored.symbols[0] == "Fe"
        assert restored.has_potential()
        assert restored.potential is not original_potential
