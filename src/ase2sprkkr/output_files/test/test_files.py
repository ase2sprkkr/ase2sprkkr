from os import listdir
from os.path import dirname, join as path_join, isfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import tempfile
import numpy as np
from ase.units import Rydberg as Ry

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from ase2sprkkr.output_files.test.init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..output_files import OutputFile
    from ..definitions.jxc import JXCOutputFile, Coordinates, Selector
    from ...common.configuration_containers import DisabledAttributeError
    from ...common.configuration_containers import RootConfigurationContainer


class TestOutput(TestCase):
    def test_output(self):
        dire = path_join(dirname(dirname(__file__)), "examples")
        for i in listdir(dire):
            fname = path_join(dire, i)
            if isfile(fname):
                ext = fname.rsplit(".", 1)[1]
                out = OutputFile.from_file(fname, unknown=False)
                if ext == "spc":
                    self.assertEqual("ARPES", out.KEYWORD())
                    if out.NE() > 1:
                        self.assertEqual((200, 160), out.ENERGY().shape)
                        o2 = out + out
                        for i in "ENERGY", "THETA", "K", "DETERMINANT":
                            self.assertEqual(out[i](), o2[i]())
                        for i in "TOTAL", "POLARIZATION", "UP", "DOWN":
                            self.assertEqual(2 * out[i](), o2[i]())
                elif ext == "dos":
                    self.assertEqual(out.n_orbitals(1), 3)
                    self.assertEqual(out.n_spins(), 2)
                    self.assertEqual((2, 3, 1200), out.dos_for_site_type("Ta").shape)
                    self.assertEqual(
                        out.DOS["Ta"][5] / Ry, out.dos_for_site_type("Ta", 1, 2)[:]
                    )
                    self.assertEqual(out["Ta"].dos, out[0].dos)
                    self.assertEqual(out["Ta"].dos, out[0].dos)
                    self.assertEqual(
                        out.total_dos().dos, (out[0] * 2 + out[1] * 2 + out[2] * 4).dos
                    )

                elif ext == "bsf":
                    if out.MODE() == "EK-REL":
                        first = out.NE()
                        second = out.NK()
                    else:
                        first = out.NK1()
                        second = out.NK2()
                    self.assertEqual(out.I().shape, (out.NQ_EFF(), first, second))
                    if out.KEYWORD() in ("BSF"):
                        self.assertEqual(
                            out.I_UP().shape, (out.NQ_EFF(), first, second)
                        )
                        self.assertRaises(DisabledAttributeError, lambda: out.I_X)
                    else:
                        self.assertEqual(out.I_X().shape, (out.NQ_EFF(), first, second))
                        self.assertRaises(DisabledAttributeError, lambda: out.I_UP)

                    if out.MODE() == "EK-REL":
                        self.assertEqual(len(out.K()), out.NK())
                        self.assertEqual(len(out.E()), out.NE())
                    else:
                        self.assertEqual(len(out.K1()), out.NK1())
                        self.assertEqual(len(out.K2()), out.NK2())

                if hasattr(out, "plot"):
                    if isinstance(out, JXCOutputFile):
                        continue
                    with tempfile.NamedTemporaryFile() as name:
                        kwargs = {"filename": name}
                        if "exclude_vc" in getattr(out, "plot_parameters", {}):
                            kwargs["exclude_vc"] = False
                        out.plot(**kwargs)

    def test_guess_potential_filename_and_manual_override(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            output = directory / "Fe_JXC_XCPLTEN_Jij.dat"
            output.write_text("unused")
            guessed = directory / "Fe.pot"
            guessed.write_text("dummy")

            class DummyDefinition:
                @staticmethod
                def parse_file(file, allow_dangerous=False):
                    return {"parsed": True}

            class DummyOutput(OutputFile):
                def clear(self, deep):
                    object.__setattr__(self, "cleared", deep)

                def set(self, values, unknown="add"):
                    object.__setattr__(self, "set_values", (values, unknown))

            out = object.__new__(DummyOutput)
            object.__setattr__(out, "_definition", DummyDefinition())
            object.__setattr__(out, "_filename", None)
            object.__setattr__(out, "_potential", None)
            object.__setattr__(out, "_potential_filename", None)
            object.__setattr__(out, "cleared", False)
            object.__setattr__(out, "set_values", None)
            RootConfigurationContainer.read_from_file(out, str(output))

            self.assertEqual(Path(out._filename).name, "Fe_JXC_XCPLTEN_Jij.dat")
            self.assertEqual(Path(out.potential_filename).name, "Fe.pot")
            self.assertTrue(out.cleared)
            self.assertEqual(out.set_values, ({"parsed": True}, "add"))

            out.set_potential_filename(directory / "manual.pot_new")
            self.assertEqual(Path(out.potential_filename).name, "manual.pot_new")

    def test_potential_property_accepts_filename_or_potential(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            output = directory / "Fe_JXC_XCPLTEN_Jij.dat"
            output.write_text("unused")
            guessed = directory / "Fe.pot"
            guessed.write_text("dummy")

            class DummyOutput(OutputFile):
                pass

            out = object.__new__(DummyOutput)
            object.__setattr__(out, "_definition", None)
            object.__setattr__(out, "_filename", str(output))
            object.__setattr__(out, "_potential", None)
            object.__setattr__(out, "_potential_filename", None)

            with patch(
                "ase2sprkkr.output_files.output_files.Potential.from_file",
                return_value="loaded",
            ) as loader:
                self.assertEqual(out.potential, "loaded")
                loader.assert_called_once_with(str(guessed))

            class FakePotential:
                def __init__(self, filename):
                    self._filename = filename
                    self.atoms = "atoms"

            with patch("ase2sprkkr.output_files.output_files.Potential", FakePotential):
                fake_potential = FakePotential(str(directory / "explicit.pot"))
                out.potential = fake_potential
                self.assertEqual(out.potential, fake_potential)
                self.assertEqual(out.atoms, "atoms")
                self.assertEqual(Path(out.potential_filename).name, "explicit.pot")

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

        out = JXCOutputFile.from_atoms(atoms)
        self.assertEqual(out.it_labels, {1: "Fe_1", 2: "Ni", 3: "Fe_2"})

    def test_dij_plot_accepts_axis_argument(self):
        filename = (
            Path(dirname(dirname(__file__))) / "examples" / "SrTiO3_JXC_XCPLTEN_Dij.dat"
        )
        out = OutputFile.from_file(str(filename), unknown=False)

        with tempfile.NamedTemporaryFile() as name:
            self.assertTrue("axis" in out.plot_parameters)
            out.plot(filename=name, exclude_vc=False, axis="y")

    def test_write_uppasd_file_dispatches_to_current_binding_writer(self):
        class DummyJxc(JXCOutputFile):
            def is_Jij(self):
                return True

        out = object.__new__(DummyJxc)
        selector = Selector(iq=[1], it={1})

        with patch("ase2sprkkr.bindings.uppasd.write_jfile") as write_jfile:
            out.write_uppasd_file(
                "jfile.dat",
                directory="exports",
                selector=selector,
                exchange_radius=6.5,
                coordinates=Coordinates.cartesian,
            )

        write_jfile.assert_called_once_with(
            out,
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

        out = object.__new__(DummyDij)

        with patch("ase2sprkkr.bindings.uppasd.write_dmfile") as write_dmfile:
            out.write_uppasd_file("dmfile.dat", directory="exports")

        write_dmfile.assert_called_once_with(
            out,
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
