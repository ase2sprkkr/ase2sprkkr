import pytest

from ..input_parameters import InputParameters


def _reader_setting(values):
    def read_from_file(self, _file, *, allow_dangerous=False):
        self.set(values, retain_invalid="all", report_invalid="ignore")

    return read_from_file


@pytest.mark.parametrize(
    ("requested", "values", "effective"),
    (
        ("ARPES", {"TASK": {"TASK": "SPLEED"}}, "SPLEED"),
        ("SPLEED", {"TASK": {"TASK": "ARPES"}}, "ARPES"),
        ("BSFEK", {"ENERGY": {"NE": [1]}}, "BSFKK"),
        ("BSFKK", {"ENERGY": {"NE": [200]}}, "BSFEK"),
    ),
)
def test_read_for_task_rejects_mismatched_effective_task(
    monkeypatch, requested, values, effective
):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry[requested]
    monkeypatch.setattr(InputParameters, "read_from_file", _reader_setting(values))

    with pytest.raises(
        ValueError,
        match=rf"read for task '{requested}'.*task '{effective}'",
    ):
        registry.read_for_task(task, "ignored.inp")


@pytest.mark.parametrize(
    ("requested", "values"),
    (
        ("SCF", {}),
        ("ARPES", {"TASK": {"TASK": "ARPES"}}),
        ("SPLEED", {"TASK": {"TASK": "SPLEED"}}),
        ("BSFEK", {"ENERGY": {"NE": [200]}}),
        ("BSFKK", {"ENERGY": {"NE": [1]}}),
    ),
)
def test_read_for_task_accepts_matching_effective_task(
    monkeypatch, requested, values
):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry[requested]
    monkeypatch.setattr(InputParameters, "read_from_file", _reader_setting(values))

    parameters = registry.read_for_task(task, "ignored.inp")

    assert parameters.task_name == requested
    assert parameters.task_definition.name == requested


def _write_parameters(tmp_path, task_name):
    parameters = InputParameters.create(task_name)
    path = tmp_path / f"{task_name.lower()}.inp"
    path.write_text(parameters.to_string(validate=False))
    return path


@pytest.mark.parametrize(
    ("requested", "file_task", "effective"),
    (
        ("SPLEED", "ARPES", "ARPES"),
        ("BSFKK", "BSFEK", "BSFEK"),
    ),
)
def test_read_for_task_rejects_real_shared_schema_mismatch(
    tmp_path, requested, file_task, effective
):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry[requested]
    filename = _write_parameters(tmp_path, file_task)

    with pytest.raises(
        ValueError,
        match=rf"read for task '{requested}'.*task '{effective}'",
    ):
        registry.read_for_task(task, filename)


@pytest.mark.parametrize("requested", ("ARPES", "BSFEK"))
def test_read_for_task_accepts_real_matching_shared_schema(tmp_path, requested):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry[requested]
    filename = _write_parameters(tmp_path, requested)

    parameters = registry.read_for_task(task, filename)

    assert parameters.task_name == requested
    assert parameters.task_definition.name == requested
