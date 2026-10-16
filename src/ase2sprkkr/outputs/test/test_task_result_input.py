from types import SimpleNamespace

from ase2sprkkr import InputParameters
from ase2sprkkr.outputs.task_result import KkrOutputReader, TaskResult


class _Parser:
    def __init__(self, print_output):
        self.print_output = print_output


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
