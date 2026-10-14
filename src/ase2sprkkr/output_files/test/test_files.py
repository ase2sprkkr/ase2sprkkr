from os import listdir
from os.path import dirname, join as path_join, isfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import tempfile

from ase.units import Rydberg as Ry
import numpy as np
import pytest

from ..output_files import OutputFile
from ..definitions.jxc import JXCOutputFile, Coordinates, Selector
from ...common.configuration_containers import DisabledAttributeError
from ...common.configuration_containers import RootConfigurationContainer


class TestOutput:

    @pytest.mark.slow
    def test_output(self):
        directory = path_join(dirname(dirname(__file__)), "examples")
        for filename in listdir(directory):
            path = path_join(directory, filename)
            if not isfile(path):
                continue

            extension = path.rsplit(".", 1)[1]
            output = OutputFile.from_file(path, unknown=False)
            if extension == "spc":
                assert output.KEYWORD() == "ARPES"
                if output.NE() > 1:
                    if filename == "WSe2_ARPES.spc":
                        assert output.ENERGY().shape == (200, 160)
                    doubled = output + output
                    for option in "ENERGY", "THETA", "K", "DETERMINANT":
                        np.testing.assert_equal(output[option](), doubled[option]())
                    for option in "TOTAL", "POLARIZATION", "UP", "DOWN":
                        np.testing.assert_equal(2 * output[option](), doubled[option]())
            elif extension == "dos":
                if filename == "TaAs_DOS.dos":
                    assert output.n_orbitals(1) == 3
                    assert output.n_spins() == 2
                    assert output.dos_for_site_type("Ta").shape == (2, 3, 1200)
                    np.testing.assert_equal(
                        output.DOS["Ta"][5] / Ry,
                        output.dos_for_site_type("Ta", 1, 2)[:],
                    )
                    np.testing.assert_equal(output["Ta"].dos, output[0].dos)
                    np.testing.assert_equal(
                        output.total_dos().dos,
                        (output[0] * 2 + output[1] * 2 + output[2] * 4).dos,
                    )
            elif extension == "bsf":
                if output.MODE() == "EK-REL":
                    first = output.NE()
                    second = output.NK()
                else:
                    first = output.NK1()
                    second = output.NK2()
                assert output.I().shape == (output.NQ_EFF(), first, second)
                if output.KEYWORD() == "BSF":
                    assert output.I_UP().shape == (output.NQ_EFF(), first, second)
                    with pytest.raises(DisabledAttributeError):
                        output.I_X
                else:
                    assert output.I_X().shape == (output.NQ_EFF(), first, second)
                    with pytest.raises(DisabledAttributeError):
                        output.I_UP

                if output.MODE() == "EK-REL":
                    assert len(output.K()) == output.NK()
                    assert len(output.E()) == output.NE()
                else:
                    assert len(output.K1()) == output.NK1()
                    assert len(output.K2()) == output.NK2()

            if hasattr(output, "plot"):
                if isinstance(output, JXCOutputFile):
                    continue
                with tempfile.NamedTemporaryFile() as target:
                    kwargs = {"filename": target}
                    if "exclude_vc" in getattr(output, "plot_parameters", {}):
                        kwargs["exclude_vc"] = False
                    output.plot(**kwargs)

    def test_guess_potential_filename_and_manual_override(self, tmp_path):
        output_path = tmp_path / "Fe_JXC_XCPLTEN_Jij.dat"
        output_path.write_text("unused")
        guessed = tmp_path / "Fe.pot"
        guessed.write_text("dummy")

        class DummyDefinition:
            @staticmethod
            def parse_file(file, allow_dangerous=False):
                return {"parsed": True}

        class DummyOutput(OutputFile):
            def stage_clear(self, transaction, *, check_required=True):
                object.__setattr__(self, "cleared", not check_required)
                return True

            def stage(self, transaction, values, *, unknown="add"):
                object.__setattr__(
                    self, "set_values", (values, unknown, transaction.policy.why)
                )
                return True

            def _validate(self, why):
                pass

        output = object.__new__(DummyOutput)
        object.__setattr__(output, "_definition", DummyDefinition())
        object.__setattr__(output, "_container", None)
        object.__setattr__(output, "_filename", None)
        object.__setattr__(output, "_potential", None)
        object.__setattr__(output, "_potential_filename", None)
        object.__setattr__(output, "cleared", False)
        object.__setattr__(output, "set_values", None)
        RootConfigurationContainer.read_from_file(output, str(output_path))

        assert Path(output._filename).name == "Fe_JXC_XCPLTEN_Jij.dat"
        assert Path(output.potential_filename).name == "Fe.pot"
        assert output.cleared
        assert output.set_values == ({"parsed": True}, "add", "parse")

        output.set_potential_filename(tmp_path / "manual.pot_new")
        assert Path(output.potential_filename).name == "manual.pot_new"

    def test_potential_property_accepts_filename_or_potential(self, tmp_path):
        output_path = tmp_path / "Fe_JXC_XCPLTEN_Jij.dat"
        output_path.write_text("unused")
        guessed = tmp_path / "Fe.pot"
        guessed.write_text("dummy")

        class DummyOutput(OutputFile):
            pass

        output = object.__new__(DummyOutput)
        object.__setattr__(output, "_definition", None)
        object.__setattr__(output, "_filename", str(output_path))
        object.__setattr__(output, "_potential", None)
        object.__setattr__(output, "_potential_filename", None)

        with patch(
            "ase2sprkkr.output_files.output_files.Potential.from_file",
            return_value="loaded",
        ) as loader:
            assert output.potential == "loaded"
            loader.assert_called_once_with(str(guessed))

        class FakePotential:
            def __init__(self, filename):
                self._filename = filename
                self.atoms = "atoms"

        with patch("ase2sprkkr.output_files.output_files.Potential", FakePotential):
            fake_potential = FakePotential(str(tmp_path / "explicit.pot"))
            output.potential = fake_potential
            assert output.potential is fake_potential
            assert output.atoms == "atoms"
            assert Path(output.potential_filename).name == "explicit.pot"

    def test_jxc_labels_are_taken_from_atoms(self):
        class FakeAtomicType:
            def __init__(self, symbol):
                self.symbol = symbol

        type_fe_1 = FakeAtomicType("Fe")
        type_ni = FakeAtomicType("Ni")
        type_fe_2 = FakeAtomicType("Fe")
        atoms = SimpleNamespace(
            sites=[
                SimpleNamespace(occupation={type_fe_1: 1.0}),
                SimpleNamespace(occupation={type_ni: 1.0}),
                SimpleNamespace(occupation={type_fe_2: 1.0}),
            ]
        )

        output = JXCOutputFile.from_atoms(atoms)
        assert output.it_labels == {1: "Fe_1", 2: "Ni", 3: "Fe_2"}

    def test_dij_plot_accepts_axis_argument(self):
        filename = (
            Path(dirname(dirname(__file__)))
            / "examples"
            / "SrTiO3_JXC_XCPLTEN_Dij.dat"
        )
        output = OutputFile.from_file(str(filename), unknown=False)

        with tempfile.NamedTemporaryFile() as target:
            assert "axis" in output.plot_parameters
            output.plot(filename=target, exclude_vc=False, axis="y")

    def test_write_uppasd_file_dispatches_to_current_binding_writer(self):
        class DummyJxc(JXCOutputFile):
            def is_Jij(self):
                return True

        output = object.__new__(DummyJxc)
        selector = Selector(iq=[1], it={1})

        with patch("ase2sprkkr.bindings.uppasd.write_jfile") as write_jfile:
            output.write_uppasd_file(
                "jfile.dat",
                directory="exports",
                selector=selector,
                exchange_radius=6.5,
                coordinates=Coordinates.cartesian,
            )

        write_jfile.assert_called_once_with(
            output,
            "jfile.dat",
            directory="exports",
            selector=selector,
            iq=None,
            it=None,
            exclude_it=None,
            exclude_vc=True,
            exchange_radius=6.5,
            coordinates=Coordinates.cartesian,
        )

    def test_write_uppasd_file_uses_dm_writer_for_dij_output(self):
        class DummyDij(JXCOutputFile):
            def is_Jij(self):
                return False

        output = object.__new__(DummyDij)

        with patch("ase2sprkkr.bindings.uppasd.write_dmfile") as write_dmfile:
            output.write_uppasd_file("dmfile.dat", directory="exports")

        write_dmfile.assert_called_once_with(
            output,
            "dmfile.dat",
            directory="exports",
            selector=None,
            iq=None,
            it=None,
            exclude_it=None,
            exclude_vc=True,
            exchange_radius=None,
            coordinates=Coordinates.lattice,
        )
