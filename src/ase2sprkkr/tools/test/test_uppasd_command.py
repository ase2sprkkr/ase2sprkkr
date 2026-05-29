if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from pathlib import Path
from types import SimpleNamespace
from ase import Atoms

from ..commands import uppasd
from ...bindings import uppasd as uppasd_bindings
from ...output_files.output_files import OutputFile
from ...output_files.definitions.jxc import Selector
from ...potentials.potentials import Potential


def make_atoms():
    atoms = Atoms(
        symbols=["Fe", "Ni", "Fe"],
        positions=[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (2.0, 0.0, 0.0),
        ],
    )
    atoms.sites = [object(), object(), object()]
    return atoms


class FakeOutput:
    def __init__(self, selector):
        self.selector = selector
        self.create_selector_calls = []
        self.write_calls = []
        self.potential = None
        self.it_labels = {1: "Fe", 2: "Ni"}

    def create_selector(self, **kwargs):
        self.create_selector_calls.append(kwargs)
        return self.selector

    def write_uppasd_file(self, path, **kwargs):
        self.write_calls.append((Path(path).name, kwargs))
        return True

    def is_Jij(self):
        return True

    def plot(self, **kwargs):
        return True


class TestUppasdCommand(TestCase):
    def test_run_passes_selector_to_current_exports(self, monkeypatch, tmp_path):
        atoms = make_atoms()
        selector = Selector(iq=[1, 3], it={1, 2})
        jxc_output = FakeOutput(selector)
        calls = {"pos": [], "mom": []}

        monkeypatch.setattr(
            Potential,
            "from_file",
            staticmethod(lambda filename: SimpleNamespace(atoms=atoms)),
        )
        monkeypatch.setattr(
            OutputFile, "from_file", staticmethod(lambda filename, **kwargs: jxc_output)
        )
        monkeypatch.setattr(
            uppasd_bindings,
            "write_pos_file",
            lambda atoms, path, **kwargs: (
                calls["pos"].append((Path(path).name, kwargs)) or True
            ),
        )
        monkeypatch.setattr(
            uppasd_bindings,
            "write_mom_file",
            lambda atoms, path, **kwargs: (
                calls["mom"].append((Path(path).name, kwargs)) or True
            ),
        )

        args = SimpleNamespace(
            pot_file="calc.pot",
            jxc_file="calc_XCPLTEN_Jij.dat",
            dmi_file=None,
            output_dir=str(tmp_path),
            no_write=False,
            plot=False,
            no_plot=False,
            separate_plots=False,
            axis="all",
            exchange_radius=7.5,
            coordinates="cartesian",
            exclude=["Ni"],
            include=["Fe"],
            include_vacuum=True,
            font_size=14,
        )

        uppasd.run(args, {})

        assert jxc_output.create_selector_calls == [
            {
                "iq": ["Fe"],
                "it": ["Fe"],
                "exclude_it": ["Ni"],
                "exclude_vc": False,
            }
        ]
        assert jxc_output.write_calls == [
            (
                "jfile.dat",
                {
                    "selector": selector,
                    "exchange_radius": 7.5,
                    "coordinates": uppasd_bindings.Coordinates.cartesian,
                },
            )
        ]
        assert calls["pos"] == [("posfile.dat", {"selector": selector})]
        assert calls["mom"] == [("momfile.dat", {"selector": selector})]

    def test_run_no_write_skips_export_helpers(self, monkeypatch, tmp_path):
        atoms = make_atoms()
        selector = Selector(iq=[1], it={1})
        jxc_output = FakeOutput(selector)
        calls = {"pos": 0, "mom": 0}

        monkeypatch.setattr(
            Potential,
            "from_file",
            staticmethod(lambda filename: SimpleNamespace(atoms=atoms)),
        )
        monkeypatch.setattr(
            OutputFile, "from_file", staticmethod(lambda filename, **kwargs: jxc_output)
        )
        monkeypatch.setattr(
            uppasd_bindings,
            "write_pos_file",
            lambda *args, **kwargs: calls.__setitem__("pos", calls["pos"] + 1) or True,
        )
        monkeypatch.setattr(
            uppasd_bindings,
            "write_mom_file",
            lambda *args, **kwargs: calls.__setitem__("mom", calls["mom"] + 1) or True,
        )

        args = SimpleNamespace(
            pot_file="calc.pot",
            jxc_file="calc_XCPLTEN_Jij.dat",
            dmi_file=None,
            output_dir=str(tmp_path),
            no_write=True,
            plot=False,
            no_plot=False,
            separate_plots=False,
            axis="all",
            exchange_radius=4.0,
            coordinates="lattice",
            exclude=None,
            include=None,
            include_vacuum=False,
            font_size=14,
        )

        uppasd.run(args)

        assert jxc_output.create_selector_calls == [
            {
                "iq": None,
                "it": None,
                "exclude_it": None,
                "exclude_vc": True,
            }
        ]
        assert jxc_output.write_calls == []
        assert calls == {"pos": 0, "mom": 0}
