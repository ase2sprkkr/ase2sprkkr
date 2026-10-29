import os
from pathlib import Path

import pytest

from ..readers import scf as scf_reader
from ..readers.scf import ScfOutputParser, ScfResult
from ..task_result import TaskResult


def test_scf_output_without_input_file_preserves_result_semantics(tmp_path):
    source = Path(__file__).parent.parent / "examples" / "scf.out"
    output = tmp_path / "scf.out"
    output.write_bytes(source.read_bytes())

    result = TaskResult.from_file(output)

    assert isinstance(result, ScfResult)
    assert result.input_parameters is None
    assert result.task_name == "scf"
    assert result.produces_potential is True


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


@pytest.mark.parametrize("has_original_calculator", (False, True))
def test_new_task_continuation_is_independent_of_original_calculator(
    tmp_path, has_original_calculator,
):
    class OriginalCalculator:
        def copy_with_potential(self, potential):
            raise AssertionError("continuation must not copy the original calculator")

    original = OriginalCalculator() if has_original_calculator else None
    result = TaskResult(None, original, tmp_path)
    result.files.add_file("potential", "Fe.pot")

    calculator = result.new_task("DOS")

    assert calculator._potential == os.path.join(os.fspath(tmp_path), "Fe.pot")
    assert calculator.input_parameters.task_definition.name == "DOS"


def test_scf_info_output_includes_final_status(capsys):
    output = Path(__file__).parent.parent / "examples" / "scf.out"
    parser = ScfOutputParser(print_info=True)
    result = ScfResult(None, None, None)

    parser.read_from_file(output, read_args=[result])

    assert parser.print_info is True
    assert "OK: The computation converged." in capsys.readouterr().out
