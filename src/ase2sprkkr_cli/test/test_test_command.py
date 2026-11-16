import argparse
from types import SimpleNamespace

import pytest

from ase2sprkkr_cli.commands import test as test_command


@pytest.mark.smoke
def test_pytest_arguments_preserve_unknown_option_order():
    parser = argparse.ArgumentParser()
    test_command.parser(parser)

    args, remainder = parser.parse_known_args(
        [
            "-k",
            "remote_executor_over_ssh or remote_executor_runs_psij_slurm_over_ssh",
        ]
    )
    getattr(args, test_command.unknowns).extend(remainder)

    assert args.pytest_arguments == [
        "-k",
        "remote_executor_over_ssh or remote_executor_runs_psij_slurm_over_ssh",
    ]


def test_pytest_arguments_preserve_mixed_positionals_and_options():
    parser = argparse.ArgumentParser()
    test_command.parser(parser)

    args, remainder = parser.parse_known_args(
        ["test_remote.py", "-k", "remote", "--maxfail=1"]
    )
    getattr(args, test_command.unknowns).extend(remainder)

    assert args.pytest_arguments == [
        "test_remote.py",
        "-k",
        "remote",
        "--maxfail=1",
    ]


@pytest.mark.parametrize(
    ("options", "expected"),
    [
        ([], []),
        (["--smoke"], ["-m", "(smoke)"]),
        (["--no-smoke"], ["-m", "(not smoke)"]),
        (["--fast"], ["-m", "(not slow)"]),
        (["--smoke", "--fast"], ["-m", "(smoke) and (not slow)"]),
        (["--no-smoke", "--fast"], ["-m", "(not smoke) and (not slow)"]),
    ],
)
def test_test_command_smoke_and_fast_flags(options, expected):
    parser = argparse.ArgumentParser()
    test_command.parser(parser)
    args = parser.parse_args(options)
    assert test_command._pytest_selection(
        [], fast=args.fast, smoke=args.smoke,
    ) == expected


def test_smoke_and_no_smoke_are_mutually_exclusive():
    parser = argparse.ArgumentParser()
    test_command.parser(parser)
    with pytest.raises(SystemExit) as error:
        parser.parse_args(["--smoke", "--no-smoke"])
    assert error.value.code == 2


def test_smoke_selection_intersects_forwarded_pytest_marker():
    assert test_command._pytest_selection(
        ["-q", "-m", "not remote", "-k", "parsing", "--durations=10"],
        fast=True, smoke=True,
    ) == ["-q", "-k", "parsing", "--durations=10",
          "-m", "(smoke) and (not slow) and ((not remote))"]


def test_smoke_flags_reach_pytest_without_disabling_full_default(monkeypatch):
    calls = []

    def fake_run(argv):
        calls.append(argv)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("subprocess.run", fake_run)
    parser = argparse.ArgumentParser()
    test_command.parser(parser)
    for flags in ([], ["--smoke"], ["--no-smoke"], ["--fast", "--smoke"]):
        args = parser.parse_args(flags)
        with pytest.raises(SystemExit) as error:
            test_command.run(args, {})
        assert error.value.code == 0
    assert calls[0][1:3] == ["-m", "pytest"]
    assert "-m" not in calls[0][3:]
    assert calls[1][-2:] == ["-m", "(smoke)"]
    assert calls[2][-2:] == ["-m", "(not smoke)"]
    assert calls[3][-2:] == ["-m", "(smoke) and (not slow)"]
