"""Command-line SPC -> NeXus conversion and optional-dependency regressions."""

import argparse
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

from ase2sprkkr_cli.commands import nexus
from ase2sprkkr_cli.errors import CLIUsageError
from ase2sprkkr_cli.main import run as cli_run


DATA = Path(__file__).resolve().parents[2] / "ase2sprkkr" / "result_files" / "test" / "data"


def parse_args(*argv):
    parser = argparse.ArgumentParser()
    nexus.parser(parser)
    return parser.parse_args(argv)


def test_export_command_and_default_output():
    args = parse_args("calculation.spc")
    assert args.spc == Path("calculation.spc")
    assert args.output is None
    assert args.input_file is None


def test_nexus_requires_spc_file():
    with pytest.raises(SystemExit):
        parse_args()


def test_nexus_help_shows_export_options():
    parser = argparse.ArgumentParser(prog="ase2sprkkr nexus")
    nexus.parser(parser)
    help_text = parser.format_help()
    assert "usage: ase2sprkkr nexus" in help_text
    assert "spc" in help_text
    assert "--output FILE" in help_text
    assert "--input-file FILE" in help_text


def test_real_cli_nexus_help_shows_export_options(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli_run(["nexus", "--help"])
    assert exit_info.value.code == 0
    output = capsys.readouterr().out
    assert "spc" in output
    assert "--output" in output
    assert "--input-file" in output


def test_missing_optional_dependency_is_cli_usage_error(tmp_path, monkeypatch, capsys):
    source = tmp_path / "sample.spc"
    source.write_text("placeholder", encoding="utf-8")
    monkeypatch.setitem(sys.modules, "h5py", None)

    # Exercise the actual dispatcher: a normal usage error, no traceback.
    with pytest.raises(SystemExit) as exc:
        cli_run(["nexus", str(source)])
    assert exc.value.code == -2
    printed = capsys.readouterr()
    assert "ase2sprkkr[nexus]" in printed.out
    assert "pip install" in printed.out
    assert "Traceback" not in printed.err
    assert not source.with_suffix(".nxs").exists()


def test_source_cannot_be_overwritten(tmp_path):
    source = tmp_path / "calculation.spc"
    source.write_text("native data", encoding="utf-8")
    args = parse_args(str(source), "-o", str(source))
    with pytest.raises(CLIUsageError, match="must not overwrite"):
        nexus.run(args)
    assert source.read_text(encoding="utf-8") == "native data"


def test_missing_input_metadata_file_is_reported(tmp_path):
    source = tmp_path / "sample.spc"
    source.touch()
    args = parse_args(str(source), "-i", str(tmp_path / "missing.inp"))
    with pytest.raises(CLIUsageError, match="input file not found"):
        nexus.run(args)


def test_missing_spc_file_is_reported(tmp_path):
    args = parse_args(str(tmp_path / "missing.spc"))
    with pytest.raises(CLIUsageError, match="SPC file not found"):
        nexus.run(args)


@pytest.mark.parametrize(
    ("fixture", "expected_axes", "use_output"),
    [
        ("Mo2Se4_ARPES_small.spc", ["energy", "theta"], False),
        ("Mo2Se4_ARPES_kx_ky_small.spc", ["k_y", "k_x"], True),
    ],
)
def test_export_native_spc_preserves_axis_order_and_values(
    tmp_path, fixture, expected_axes, use_output, capsys
):
    h5py = pytest.importorskip("h5py")
    from ase2sprkkr.result_files import ResultFile

    source = tmp_path / "spectrum.spc"
    shutil.copyfile(DATA / fixture, source)
    target = tmp_path / ("custom.nxs" if use_output else "spectrum.nxs")
    arguments = ["nexus", str(source)]
    if use_output:
        arguments.extend(["--output", str(target)])

    with pytest.warns(UserWarning, match="No SPR-KKR input"):
        cli_run(arguments)

    assert target.is_file()
    assert str(target) in capsys.readouterr().out
    original = ResultFile.from_file(source, try_only="spc")
    with h5py.File(target, "r") as file:
        data = file["entry/data"]
        assert list(data.attrs["axes"]) == expected_axes
        np.testing.assert_array_equal(data["data"], original.TOTAL())
        np.testing.assert_array_equal(
            file["entry/sprkkr_spc/raw_data"], original.RAW_DATA()
        )
        if expected_axes == ["k_y", "k_x"]:
            assert data["data"].shape == (original.NP(), original.NT())
            assert np.atleast_1d(data.attrs["k_y_indices"]).tolist() in ([0], [0, 1])
            assert np.atleast_1d(data.attrs["k_x_indices"]).tolist() in ([1], [0, 1])
            assert data["data"][1, 2] == original.TOTAL()[1, 2]
        else:
            assert np.atleast_1d(data.attrs["energy_indices"]).tolist() in ([0], [0, 1])
            assert np.atleast_1d(data.attrs["theta_indices"]).tolist() in ([1], [0, 1])


def test_explicit_input_file_is_forwarded_to_writer(tmp_path, monkeypatch, capsys):
    pytest.importorskip("h5py")
    from ase2sprkkr.result_files import ResultFile

    source = tmp_path / "spectrum.spc"
    source.touch()
    params = tmp_path / "experiment.inp"
    params.write_text("TASK ARPES\n", encoding="utf-8")
    target = tmp_path / "output.nxs"
    received = []

    class Spectrum:
        def to_nexus(self, path, *, input_file):
            received.append((path, input_file))
            return path

    monkeypatch.setattr(
        ResultFile, "from_file", staticmethod(
            lambda path, *, try_only: Spectrum()
        ),
    )
    cli_run(["nexus", str(source), "-o", str(target), "-i", str(params)])
    assert received == [(target, params)]
    assert str(target) in capsys.readouterr().out
