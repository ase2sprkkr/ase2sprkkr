import os
from pathlib import Path

from ..readers import scf as scf_reader
from ..readers.scf import ScfOutputParser, ScfResult


def test_scf_result_calculator_without_original_calculator(monkeypatch, tmp_path):
    result = ScfResult(None, None, tmp_path)
    result.files.add_file("potential", "Fe.pot")
    expected_potential = os.path.join(os.fspath(tmp_path), "Fe.pot_new")
    marker = object()

    def make_calculator(*, potential):
        assert potential == expected_potential
        return marker

    monkeypatch.setattr(scf_reader, "SPRKKR", make_calculator)

    assert result.calculator is marker


def test_scf_info_output_includes_final_status(capsys):
    output = Path(__file__).parent.parent / "examples" / "scf.out"
    parser = ScfOutputParser("info")
    result = ScfResult(None, None, None)

    parser.read_from_file(output, read_args=[result])

    assert parser.print_info is True
    assert parser.print_output is False
    assert "OK: The computation converged." in capsys.readouterr().out
