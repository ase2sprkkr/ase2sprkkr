from ase2sprkkr.input_parameters import InputParameters
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

    assert parameters.task_name == "ARPES"
    assert parameters.task_definition.name == "SPLEED"

    result = TaskResult(parameters, None, str(tmp_path))
    assert result.task_name == "spleed"

    reader = _Reader(parameters, None, str(tmp_path))
    assert reader.task_name == "SPLEED"
    assert reader.create_result("run.out").task_name == "spleed"
