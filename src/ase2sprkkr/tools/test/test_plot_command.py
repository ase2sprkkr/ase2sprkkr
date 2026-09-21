import argparse

import pytest

from ..commands import plot
from ...output_files.definitions.sfn import SFNOutputFile
from ...output_files.output_files import OutputFile


def parse_args(*arguments):
    parser = argparse.ArgumentParser()
    plot.parser(parser)
    return parser, parser.parse_args(arguments)


def test_sfn_mode_flags_are_documented():
    parser, args = parse_args("file.sfn", "--shape")

    assert args.what == "shape"
    help_text = parser.format_help()
    for mode in plot.SFN_PLOT_MODES:
        assert f"--{mode}" in help_text


def test_sfn_plot_modes_are_mutually_exclusive():
    with pytest.raises(SystemExit):
        parse_args("file.sfn", "--shape", "--radial")


def test_sfn_mode_flag_is_rejected_for_non_sfn(monkeypatch):
    class OtherOutput:
        pass

    monkeypatch.setattr(
        OutputFile,
        "from_file",
        staticmethod(lambda *args, **kwargs: OtherOutput()),
    )
    _, args = parse_args("file.dos", "--shape")

    with pytest.raises(ValueError, match="only valid for SFN files"):
        plot.run(args, {"debug": True})


def test_sfn_mode_flag_raises_for_non_sfn(monkeypatch):
    class OtherOutput:
        pass

    monkeypatch.setattr(
        OutputFile,
        "from_file",
        staticmethod(lambda *args, **kwargs: OtherOutput()),
    )
    _, args = parse_args("file.dos", "--shape")

    with pytest.raises(ValueError, match="only valid for SFN files"):
        plot.run(args, {"debug": False})


def test_sfn_mode_cannot_be_set_through_generic_option():
    _, args = parse_args("file.sfn", "--mesh", "-S", "what=shape")

    with pytest.raises(ValueError, match="cannot be set using"):
        plot.run(args, {"debug": False})


def test_sfn_mode_flag_is_passed_as_what(monkeypatch):
    class FakeSFN(SFNOutputFile):
        def __init__(self):
            object.__setattr__(self, "plot_kwargs", None)

        def plot(self, **kwargs):
            object.__setattr__(self, "plot_kwargs", kwargs)

    output = FakeSFN()
    monkeypatch.setattr(
        OutputFile,
        "from_file",
        staticmethod(lambda *args, **kwargs: output),
    )
    _, args = parse_args("file.sfn", "--shape")

    assert plot.run(args, {"debug": True}) is None
    assert output.plot_kwargs["what"] == "shape"
