import argparse

from ase2sprkkr.tools.commands import test as test_command


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
