import argparse
from functools import cache

import pytest

from ase2sprkkr_cli.commands import plot
from ase2sprkkr.result_files import ResultFile
from ase2sprkkr.result_files.definitions.sfn import SFNResultFile


@cache
def _parser():
    parser = argparse.ArgumentParser()
    plot.parser(parser)
    return parser


def parse_args(*arguments):
    parser = _parser()
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
        ResultFile,
        "from_file",
        staticmethod(lambda *args, **kwargs: OtherOutput()),
    )
    _, args = parse_args("file.dos", "--shape")

    for debug in (True, False):
        with pytest.raises(ValueError, match="only valid for SFN files"):
            plot.run(args, {"debug": debug})


def test_sfn_mode_cannot_be_set_through_generic_option():
    _, args = parse_args("file.sfn", "--mesh", "-S", "what=shape")

    with pytest.raises(ValueError, match="cannot be set using"):
        plot.run(args, {"debug": False})


@pytest.mark.parametrize(
    "sites, expected",
    [
        ("2", 1),
        ("2,4", slice(1, 4)),
    ],
)
def test_sfn_mode_flag_and_deferred_sites_are_passed_to_plot(monkeypatch, sites, expected):
    class FakeSFN(SFNResultFile):
        def __init__(self):
            object.__setattr__(self, "plot_kwargs", None)

        def plot(self, **kwargs):
            object.__setattr__(self, "plot_kwargs", kwargs)

    output = FakeSFN()
    monkeypatch.setattr(
        ResultFile,
        "from_file",
        staticmethod(lambda *args, **kwargs: output),
    )
    _, args = parse_args("file.sfn", "--shape", "--sites", sites)

    assert args.sites == sites  # The parser preserves text for deferred conversion.
    assert plot.run(args, {"debug": True}) is None
    assert output.plot_kwargs["what"] == "shape"
    assert output.plot_kwargs["sites"] == expected
