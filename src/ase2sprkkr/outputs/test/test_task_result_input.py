from ase2sprkkr.outputs.task_result import TaskResult


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
