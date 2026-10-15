import os
import re
from datetime import datetime
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk

from ...potentials.potentials import Potential
from ..build import semiinfinite_system
from ..calculator import FilenameTemplator, SPRKKR


def _fast_atoms(atoms, jrws=20, r1=2e-6):
    for site in atoms.sites:
        site.mesh.jrws = jrws
        site.mesh.r1 = r1


def _sprkkr_enabled():
    return not os.environ.get("DO_NOT_RUN_SPRKKR", "")


class TestCalculator:
    def test_filename_template_date(self):
        calculator = SPRKKR()
        templator = FilenameTemplator(calculator, lambda: "SCF", lambda: "Fe")

        rendered = templator("%d")
        datetime.strptime(rendered, "%Y-%m-%d_%H:%M")
        assert templator("%d") == rendered

    def test_atoms_lazy_loads_string_potential(self):
        atoms = bulk("Fe")
        loaded_potential = type("LoadedPotential", (), {"atoms": atoms})()
        calculator = SPRKKR(potential="x.pot")

        assert calculator._potential == "x.pot"
        with patch.object(Potential, "from_file", return_value=loaded_potential) as from_file:
            assert calculator.atoms is atoms
            from_file.assert_called_once_with("x.pot")
            assert calculator._potential is loaded_potential
            assert calculator.atoms is atoms
            from_file.assert_called_once_with("x.pot")

    def test_save_input_rejects_atoms_with_explicit_potential(self):
        atoms = bulk("Fe")
        potential = Potential.from_atoms(bulk("Cu"))
        calculator = SPRKKR()

        for explicit_potential in (potential, "explicit.pot"):
            with pytest.raises(ValueError, match="both a potential and atoms"):
                calculator.save_input(atoms=atoms, potential=explicit_potential)

    @pytest.mark.slow
    def test_2d(self, calc_args):
        atoms = Atoms(
            symbols="C",
            positions=[[0, 0, 0]],
            cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            pbc=[1, 1, 1],
        )
        slab = semiinfinite_system(atoms, repeat=2)
        assert np.max(slab.arrays["spacegroup_kinds"]) == 5

        calculator = SPRKKR(atoms=slab, **calc_args())
        _fast_atoms(slab)
        calculator.input_parameters.set_from_atoms(slab)
        calculator.input_parameters.CONTROL.POTFIL = "xxx"
        text = calculator.input_parameters.to_string()
        assert re.search("NKTAB3D=5", text)
        assert not re.match("NKTAB=", text)

        if not _sprkkr_enabled():
            pytest.skip("SPR-KKR execution disabled by DO_NOT_RUN_SPRKKR")

        output = SPRKKR().calculate(
            slab,
            **calc_args(options={"EMIN": 8.0, "NITER": 2, "NKTAB3D": 4}),
        )
        output.input_parameters.set_from_atoms(slab)
        output.input_parameters.CONTROL.POTFIL = "xxx"
        text = output.input_parameters.to_string()
        assert re.search("NKTAB3D=", text)
        assert not re.match("NKTAB=", text)

    @pytest.mark.slow
    def test_calculator(self, calc_args, tmp_path):
        def here(name):
            return str(tmp_path / name)

        atoms = bulk("Li")
        calculator = SPRKKR(atoms=atoms, **calc_args())
        _fast_atoms(atoms)

        calculator.input_parameters.set(NE=21111)
        assert calculator.input_parameters.get("NE") == 21111
        calculator.set(NE=31111)
        assert calculator.input_parameters.get("NE") == 31111
        assert calculator.input_parameters.get("ENERGY.NE") == 31111
        calculator.input_parameters.set({"ENERGY.NE": 41111})
        assert calculator.input_parameters.get("NE") == 41111
        calculator.set({"ENERGY.NE": 51111})
        assert calculator.input_parameters.get("NE") == 51111
        calculator.set("ENERGY.NE", 451111)
        assert calculator.input_parameters.get("NE") == 451111

        calculator.input_parameters.set(NE=11111)
        assert calculator.input_parameters.get("NE") == 11111
        assert calculator.get("NE") == 11111
        assert calculator.input_parameters.ENERGY.NE() == 11111
        assert calculator.input_parameters.TASK.TASK() == "SCF"
        calculator.change_task("DOS", retain_values=True)
        assert calculator.input_parameters.TASK.TASK() == "DOS"
        assert calculator.input_parameters.get("NE") == 11111
        calculator.input_parameters.change_task("JXC", retain_values=True)
        assert calculator.input_parameters.TASK.TASK() == "JXC"
        assert calculator.input_parameters.get("NE") == 11111
        calculator.input_parameters = "PHAGEN"
        assert calculator.input_parameters.TASK.TASK() == "PHAGEN"
        assert calculator.input_parameters.get("NE") != 11111

        calculator = SPRKKR(atoms=atoms, **calc_args())
        inp_file = here("output_test_calc.inp")
        pot_file = here("output_test_calc.pot")
        Path(inp_file).touch()
        Path(pot_file).touch()

        times = {"inp": 0, "pot": 0}

        def assert_change(*expected):
            for extension, should_change in zip(("inp", "pot"), expected):
                filename = here(f"output_test_calc.{extension}")
                mtime = os.path.getmtime(filename)
                changed = abs(mtime - times[extension]) > 2.0
                assert changed == should_change
                if changed:
                    mtime -= 10
                    times[extension] = mtime
                    os.utime(filename, (mtime, mtime))

        assert_change(True, True)
        assert_change(False, False)
        calculator.save_input()
        assert_change(True, True)
        assert_change(False, False)
        calculator.save_input(input_parameters=inp_file)
        assert_change(True, True)
        calculator.save_input(input_parameters=inp_file, potential=False)
        assert_change(False, False)
        calculator.save_input(input_parameters=inp_file, options={"NITER": 2}, potential=False)
        assert_change(True, False)
        calculator.save_input(input_parameters=inp_file, potential=atoms.potential)
        assert_change(True, True)
        calculator.save_input(input_parameters=inp_file, potential=pot_file)
        assert_change(True, False)

    @pytest.mark.slow
    def test_run(self, calc_args, tmp_path):
        if not _sprkkr_enabled():
            pytest.skip("SPR-KKR execution disabled by DO_NOT_RUN_SPRKKR")

        atoms = bulk("Li")
        calculator = SPRKKR(atoms=atoms, **calc_args(mpi=False))
        _fast_atoms(atoms)
        atoms.calc = calculator
        calculator.input_parameters.find("NITER").set(2)
        calculator.input_parameters.find("NKTAB").set(50)
        assert isinstance(atoms.get_potential_energy(), float)
        calculator.input_parameters.find("NITER").set(100)

        output = calculator.calculate(options={"NITER": 2}, mpi=4)
        assert len(output.iterations) == 2
        assert str(atoms.symbols) == str(output.atoms.symbols)
        assert not output.converged
        assert isinstance(output.potential, Potential)
        assert isinstance(output.calculator.potential_object, Potential)

        reread = SPRKKR.InputParameters.create("scf").read_output_from_file(output.output_file)
        assert output.output_file == reread.output_file
        assert len(reread.iterations) == 2
        output.plot(filename=str(tmp_path / "output_test_calc.png"))

        atoms.calc = calculator
        assert isinstance(atoms.get_potential_energy(), float)

        atoms = bulk("LiCl", "rocksalt", a=5.64) * (2, 1, 1)
        calculator = SPRKKR(atoms=atoms, **calc_args())
        output = calculator.calculate(options={"EMIN": -1.0, "NITER": 1, "NKTAB": 40, "NE": 20})
        assert len(output.iterations) == 1
        assert len(output.iterations[-1]["atomic_types"]) == 3
        for atomic_type in output.iterations[-1]["atomic_types"].values():
            assert len(atomic_type["orbitals"]) == 4
        assert str(atoms.symbols) == str(output.atoms.symbols)

        reread = SPRKKR.InputParameters.from_file(output.input_file).read_output_from_file(output.output_file)
        assert str(atoms.symbols) == str(reread.atoms.symbols)
        assert len(reread.iterations) == 1

    @pytest.mark.slow
    def test_phagen(self, calc_args):
        if not _sprkkr_enabled():
            pytest.skip("SPR-KKR execution disabled by DO_NOT_RUN_SPRKKR")

        atoms = bulk("Li")
        calculator = SPRKKR(atoms=atoms, **calc_args())
        _fast_atoms(atoms)
        input_parameters = SPRKKR.InputParameters.create("scf")
        input_parameters.SCF.NITER = 1
        input_parameters.TAU.NKTAB = 20

        output = calculator.calculate(input_parameters=input_parameters)
        assert output.atoms == output.potential.atoms
        assert len(output.iterations) == 1
        assert str(atoms.symbols) == str(output.atoms.symbols)
        assert not output.converged

        calculator = SPRKKR(**calc_args())
        output = calculator.calculate(potential=output.potential, input_parameters=input_parameters)
        assert str(atoms.symbols) == str(output.atoms.symbols)
        assert len(output.iterations) == 1

        calculator.calculate(
            input_parameters="PHAGEN",
            potential=output.potential_filename,
            options={"NKTAB": 30},
        )
        output.calculator.calculate(input_parameters="PHAGEN", options={"NKTAB": 10})
