from pathlib import Path
from types import SimpleNamespace

from ase import Atoms

from ..commands import uppasd
from ...bindings import uppasd as uppasd_bindings
from ...result_files import ResultFile
from ...result_files.definitions.jxc import Selector
from ...potentials.potentials import Potential


def make_atoms():
    atoms = Atoms(symbols=["Fe", "Ni", "Fe"], positions=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
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


class TestUppasdCommand:
    def test_run_skips_exports_then_passes_selector_to_writers(self, monkeypatch, tmp_path):
        # Exercise both paths with the same mocked dependencies. The no-write
        # path runs first, so it cannot benefit from an earlier export.
        atoms = make_atoms()
        calls = {"pos": [], "mom": []}
        current_output = [None]

        monkeypatch.setattr(
            Potential,
            "from_file",
            staticmethod(lambda filename: SimpleNamespace(atoms=atoms)),
        )
        monkeypatch.setattr(
            ResultFile,
            "from_file",
            staticmethod(lambda filename, **kwargs: current_output[0]),
        )
        monkeypatch.setattr(
            uppasd_bindings,
            "write_pos_file",
            lambda atoms, path, **kwargs: calls["pos"].append((Path(path).name, kwargs)) or True,
        )
        monkeypatch.setattr(
            uppasd_bindings,
            "write_mom_file",
            lambda atoms, path, **kwargs: calls["mom"].append((Path(path).name, kwargs)) or True,
        )

        args = SimpleNamespace(
            outfile=None,
            inpsd=None,
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
            quiet=True,
        )

        no_write_output = FakeOutput(Selector(iq=[1], it={1}))
        current_output[0] = no_write_output
        uppasd.run(args, {"debug": True})

        assert no_write_output.create_selector_calls == [
            {"iq": None, "it": None, "exclude_it": None, "exclude_vc": True}
        ]
        assert no_write_output.write_calls == []
        assert calls == {"pos": [], "mom": []}

        # Fresh output object keeps the two behavioural scenarios independent,
        # while the common dependency setup is reused.
        selector = Selector(iq=[1, 3], it={1, 2})
        write_output = FakeOutput(selector)
        current_output[0] = write_output
        args.no_write = False
        args.exchange_radius = 7.5
        args.coordinates = "cartesian"
        args.exclude = ["Ni"]
        args.include = ["Fe"]
        args.include_vacuum = True

        uppasd.run(args, {"debug": True})

        assert write_output.create_selector_calls == [
            {"iq": ["Fe"], "it": ["Fe"], "exclude_it": ["Ni"], "exclude_vc": False}
        ]
        assert write_output.write_calls == [
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
