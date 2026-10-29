from types import SimpleNamespace

from ase2sprkkr import InputParameters
from ase2sprkkr.outputs.task_result import KkrOutputReader, TaskResult


class _Parser:
    def __init__(self, print_info):
        self.print_info = print_info


class _Reader(KkrOutputReader):
    parser_class = _Parser
    result_class = TaskResult


def test_explicit_input_file_is_registered(tmp_path):
    input_file = tmp_path / "run.inp"
    input_file.write_text("input")

    result = TaskResult(
        None,
        None,
        str(tmp_path),
        output_file="run.out",
        input_file=str(input_file),
    )

    assert result.path_to("input") == str(input_file)
    assert result.input_parameters_file == str(input_file)


def test_result_uses_effective_task_identity(tmp_path):
    parameters = InputParameters.create("ARPES")
    parameters.TASK.TASK = "SPLEED"

    assert parameters.requested_task_name == "ARPES"
    assert parameters.task_name == "SPLEED"
    assert parameters.task_definition.name == "SPLEED"

    result = TaskResult(parameters, None, str(tmp_path))
    assert result.task_name == "spleed"

    reader = _Reader(parameters, None, str(tmp_path))
    assert reader.task_name == "SPLEED"
    assert reader.create_result("run.out").task_name == "spleed"


def test_jxc_continuation_uses_source_result_energy(tmp_path):
    source_parameters = SimpleNamespace(
        CONTROL=SimpleNamespace(POTFIL=lambda: "Fe.pot")
    )
    result = TaskResult(
        source_parameters,
        None,
        str(tmp_path),
        task_name="SCF",
    )
    result.last_iteration = SimpleNamespace(
        energy=SimpleNamespace(EMIN=lambda: -0.35)
    )

    continuation = result.new_task("JXC")

    assert continuation.input_parameters.task_definition.name == "JXC"
    assert continuation.input_parameters.ENERGY.EMIN() == -0.35

def test_autodetect_single_matching_native_input_file(tmp_path):
    parameters = InputParameters.create("ARPES")
    parameters.CONTROL.DATASET = "sample"
    parameters.CONTROL.POTFIL = "sample.pot"
    source = tmp_path / "sample_SPEC.inp"
    source.write_text(parameters.to_string(validate=False))
    result = TaskResult(None, None, str(tmp_path), output_file="sample_SPEC.out", task_name="ARPES")
    assert result.input_parameters_file == str(source)
    assert result.input_parameters.task_name == "ARPES"


def test_autodetect_is_not_ambiguous_or_cross_task(tmp_path):
    for name in ("first.inp", "second.inp"):
        (tmp_path / name).write_text("invalid")
    result = TaskResult(None, None, str(tmp_path), task_name="ARPES")
    assert result.input_parameters is None

def test_no_result_task_identity_raises_clear_error(tmp_path):
    import pytest

    result = TaskResult(None, None, str(tmp_path))
    assert result.input_parameters is None
    with pytest.raises(ValueError, match="Cannot determine the task"):
        _ = result.task_name


def test_historical_input_autodetect_accepts_in_suffix(tmp_path):
    params = InputParameters.create("ARPES")
    params.CONTROL.DATASET = "sample"
    params.CONTROL.POTFIL = "sample.pot"
    source = tmp_path / "sample_SPEC.in"
    source.write_text(params.to_string(validate=False))
    result = TaskResult(None, None, str(tmp_path), output_file="sample_SPEC.out", task_name="ARPES")
    assert result.input_parameters_file == str(source)
    assert result.input_parameters.task_name == "ARPES"


def test_invalid_guessed_input_is_nonfatal_but_explicit_input_error_propagates(
    tmp_path, monkeypatch
):
    import pytest

    path = tmp_path / "legacy.inp"
    path.write_text("invalid input")

    def invalid_input(_cls, filename):
        raise ValueError("invalid historical input")

    monkeypatch.setattr(InputParameters, "from_file", classmethod(invalid_input))
    guessed = TaskResult(None, None, str(tmp_path), task_name="SCF")
    assert guessed.input_parameters is None
    assert guessed.input_parameters_file is None

    explicit = TaskResult(None, None, str(tmp_path), input_file=str(path))
    with pytest.raises(ValueError, match="invalid historical input"):
        _ = explicit.input_parameters
