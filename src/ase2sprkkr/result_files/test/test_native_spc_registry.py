"""Several SPC parsers can coexist and results select one by task."""
from ase2sprkkr.result_files import ResultFile
from ase2sprkkr.result_files.definition_registry import result_file_definition_registry

def test_native_spc_definitions_are_distinct():
    names = ResultFile.definitions
    assert names["band_spc"].extension == "spc"
    assert names["spleed_spc"].extension == "spc"
    assert names["spc"].extension == "spc"
    assert len(result_file_definition_registry.by_extension()["spc"]) == 3

def test_band_and_spleed_results_register_typed_spectra(tmp_path):
    from ase2sprkkr.outputs.readers.spec import SpecResult
    (tmp_path / "Fe_data.spc").write_text("dummy")
    for task, kind in (("BAND", "band_spc"), ("SPLEED", "spleed_spc")):
        result = SpecResult(None, None, tmp_path, output_file="Fe_SPEC.out", task_name=task)
        assert result.output_values["SPC"].file_type == kind


def test_band_plot_uses_one_scatter_plot_and_correct_k_axis(monkeypatch):
    """A BAND plot is not an ARPES intensity colormesh."""
    from contextlib import contextmanager
    from types import SimpleNamespace
    import numpy as np
    import pytest
    from ase2sprkkr.result_files.definitions import band_spc

    class FakeAxis:
        def __init__(self):
            self.points = []
            self.lines = []
            self.labels = []

        def scatter(self, x, y, **kwargs):
            self.points.append((np.asarray(x), np.asarray(y), kwargs))

        def axhline(self, *args, **kwargs):
            self.lines.append((args, kwargs))

        def set_xlabel(self, label):
            self.labels.append(("x", label))

        def set_ylabel(self, label):
            self.labels.append(("y", label))

    for nt, axis_symbol in ((1, "k_z"), (2, "k_{\\\\parallel}")):
        axes = FakeAxis()
        options = {}
        @contextmanager
        def fake_single_plot(**kwargs):
            options.update(kwargs)
            yield axes, {}

        monkeypatch.setattr(band_spc, "single_plot", fake_single_plot)
        fake = SimpleNamespace(
            NT=lambda: nt,
            RAW_DATA=lambda: np.array([[0.0, -1.0, 8.0], [0.5, 1.5, 12.0]]),
        )
        actual = band_spc.BandSPCResultFile.plot(
            fake, filename="band.png", show=False
        )
        assert actual is axes
        assert len(axes.points) == 1
        np.testing.assert_array_equal(axes.points[0][0], [0.0, 0.5])
        np.testing.assert_array_equal(axes.points[0][1], [-1.0, 1.5])
        assert axes.lines[0][0] == (0,)
        assert axis_symbol in axes.labels[0][1]
        assert "E" in axes.labels[1][1]
        assert options["filename"] == "band.png"
        assert options["show"] is False


def test_band_plot_can_skip_fermi_line(monkeypatch):
    from contextlib import contextmanager
    from types import SimpleNamespace
    import numpy as np
    from ase2sprkkr.result_files.definitions import band_spc

    class FakeAxis:
        def __init__(self):
            self.lines = []
            self.points = []
        def scatter(self, x, y, **kwargs):
            self.points.append((x, y, kwargs))
        def axhline(self, *args, **kwargs):
            self.lines.append(args)
        def set_xlabel(self, label):
            pass
        def set_ylabel(self, label):
            pass

    axis = FakeAxis()
    @contextmanager
    def fake_single_plot(**kwargs):
        yield axis, {}
    monkeypatch.setattr(band_spc, "single_plot", fake_single_plot)
    fake = SimpleNamespace(
        NT=lambda: 1,
        RAW_DATA=lambda: np.array([[0.0, 0.0, 0.0]]),
    )
    band_spc.BandSPCResultFile.plot(fake, axis=axis, show_fermi=False, s=3)
    assert len(axis.points) == 1
    assert axis.points[0][2]["s"] == 3
    assert axis.lines == []
