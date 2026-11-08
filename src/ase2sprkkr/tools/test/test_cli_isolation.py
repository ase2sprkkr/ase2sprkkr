"""Regression tests for the separate, lightweight CLI import boundary."""

import argparse
import subprocess
import sys

import pytest

from ase2sprkkr_cli.main import Parser
from ase2sprkkr_cli.commands import run as run_command


@pytest.mark.parametrize(
    "argv",
    [
        ["--help"],
        ["run", "--help"],
        ["config", "--help"],
        ["nomad", "--help"],
        ["info", "--help"],
        ["--version"],
    ],
)
def test_simple_cli_does_not_import_calculation_package(argv):
    # A subprocess is essential: pytest has already imported ase2sprkkr here.
    program = """
import importlib.abc
import sys

class ForbidASE2SPRKKR(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "ase2sprkkr" or fullname.startswith("ase2sprkkr.") or fullname == "kkrunner" or fullname.startswith("kkrunner."):
            raise AssertionError("Heavy package imported: " + fullname)
        return None

sys.meta_path.insert(0, ForbidASE2SPRKKR())
if sys.argv[1:] == ["--version"]:
    import importlib.metadata
    original_version = importlib.metadata.version
    importlib.metadata.version = lambda dist: "0.0.0" if dist == "ase2sprkkr" else original_version(dist)
from ase2sprkkr_cli.main import run
run(sys.argv[1:])
"""
    proc = subprocess.run(
        [sys.executable, "-c", program, *argv],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert proc.stdout.strip()


def test_extra_help_only_evaluated_when_formatted():
    calls = []
    parser = Parser(prog="ase2sprkkr example")
    parser.extra_help = lambda: calls.append("called") or "Generated details"
    parser.add_argument("--example", help="Ordinary help")

    assert calls == []
    assert parser.parse_args(["--example", "value"]).example == "value"
    assert calls == []
    text = parser.format_help()
    assert "Ordinary help" in text
    assert text.endswith("Generated details\n")
    assert calls == ["called"]


def test_run_arguments_remain_unresolved_until_execution():
    parser = argparse.ArgumentParser()
    run_command.parser(parser)
    args = parser.parse_args(["potential.pot", "--task", "SCF", "NKTAB=7"])
    assert args.task == "SCF"
    assert args.options == ["NKTAB=7"]


def test_real_subparser_appends_dynamic_help(monkeypatch, capsys):
    from ase2sprkkr_cli.main import run
    from ase2sprkkr_cli.commands import plot

    calls = []
    monkeypatch.setattr(plot, "extra_help", lambda: calls.append(True) or "DYNAMIC FORMAT LIST")
    with pytest.raises(SystemExit) as exc:
        run(["plot", "--help"])

    assert exc.value.code == 0
    assert calls == [True]
    assert "DYNAMIC FORMAT LIST" in capsys.readouterr().out


@pytest.mark.parametrize(
    "argv",
    [
        ["run", "missing.pot", "--task", "NONEXISTENT_CLI_TASK"],
        ["run", "missing.pot", "INVALID_PARAMETER_SYNTAX"],
        ["plot", "missing.bsf", "--set", "INVALID_PARAMETER_SYNTAX"],
        ["plot", "missing.bsf", "--sites", "not-a-number"],
        ["plot", "missing.bsf", "--plot_size", "not-a-size"],
        ["plot", "missing.bsf", "--plot_size", "10 cmmmm,4 cm"],
    ],
)
def test_invalid_arguments_report_cli_usage_error(argv, capsys):
    from ase2sprkkr_cli.main import run

    with pytest.raises(SystemExit) as exc:
        run(argv)
    assert exc.value.code == -2
    captured = capsys.readouterr()
    assert "Invalid" in captured.out or "Unknown task" in captured.out
    assert "Traceback" not in captured.err


def test_remote_help_delegates_only_when_selected(capsys):
    from ase2sprkkr_cli.main import run
    import pytest

    with pytest.raises(SystemExit) as exited:
        run(["remote", "--help"])
    assert exited.value.code == 0
    assert "remote" in capsys.readouterr().out.lower()
