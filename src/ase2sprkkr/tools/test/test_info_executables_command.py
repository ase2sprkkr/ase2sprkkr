import argparse

from ase2sprkkr.tools.commands import info as info_command


def _parser():
    parser = argparse.ArgumentParser()
    info_command.parser(parser)
    return parser


def test_info_parser_accepts_executables_option():
    args = _parser().parse_args(["--executables"])

    assert args.executables is True


def test_info_parser_leaves_executables_unselected_by_default():
    args = _parser().parse_args([])

    assert args.executables is None


def test_info_executables_lists_only_executable_capabilities(capsys):
    args = _parser().parse_args(["--executables"])

    info_command.run(args, None)

    output = capsys.readouterr().out
    assert "Known SPR-KKR executables" in output
    assert "Executable" in output
    assert "Default reader" in output
    assert "kkrscf" in output
    assert "kkrspec" in output
    assert "SPEC" in output
    assert "The known SPRKKR tasks" not in output
    assert "Known SPR-KKR structured result files" not in output
