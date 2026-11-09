"""Live SCF progress uses the full parser's iteration formatting."""

from pathlib import Path

from ase2sprkkr.outputs.readers.scf import ScfOutputParser, ScfProgressParser
from ase2sprkkr.outputs.task_result import TaskResult


EXAMPLE = Path(__file__).parent.parent / "examples" / "scf.out"


def test_latest_scf_iteration_matches_complete_parser():
    result = TaskResult.from_file(EXAMPLE)
    latest = ScfProgressParser.latest_from_file(EXAMPLE)
    last = result.last_iteration
    assert latest["iteration"] == last.iteration()
    assert latest["error"] == last.error()
    assert latest["moment"]["spin"] == last.moment.spin()
    assert latest["converged"] == last.converged()


def test_scf_progress_ignores_incomplete_iteration():
    output = []
    parser = ScfProgressParser(on_info=output.append)
    header = "  8 ERR 0.123E-04 0.456E-06 EF 0.6 0.01 D 12.00 M 1.2 0.3\n"
    footer = "    ETOT    -42.5    SCF - cycle not converged - continue\n"
    parser.feed(header[:20])
    assert not output
    parser.feed(header[20:] + footer[:20])
    assert not output
    parser.feed(footer[20:])
    assert len(output) == 1
    assert output[0] == ScfOutputParser.iteration_info(parser.last_iteration) + "\n"
    parser.feed("  9 ERR 0.111E-05 0.000E+00 EF 0.6 0.0 D 11 M 1.2 0.3\n")
    assert parser.last_iteration["iteration"] == 8


def test_scf_progress_tail_ignores_partial_cycle(tmp_path):
    content = ("  8 ERR 0.123E-04 0.456E-06 EF 0.6 0.01 D 12 M 1.2 0.3\n"
               " ETOT -42.5 SCF - cycle not converged - continue\n"
               "  9 ERR 0.222E-04 0.333E-06 EF 0.6 0.01 D 12 M 1.2 0.3\n")
    path = tmp_path / "unfinished.out"
    path.write_text(content)
    assert ScfProgressParser.latest_from_file(path)["iteration"] == 8
