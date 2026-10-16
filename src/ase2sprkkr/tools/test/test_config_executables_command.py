import argparse

from ase2sprkkr.tools.commands import config as config_command


def test_config_parser_accepts_executables_option():
    parser = argparse.ArgumentParser()
    config_command.parser(parser)

    args = parser.parse_args(["--executables"])

    assert args.executables is True


def test_config_executables_lists_reader_defaults(capsys):
    config_command._print_executables()

    output = capsys.readouterr().out
    assert "Executable" in output
    assert "Default reader" in output
    assert "kkrscf" in output
    assert "kkrspec" in output
    assert "SPEC" in output
