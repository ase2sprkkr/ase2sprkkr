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
        (
            "BSFEK",
            {"TASK": {"TASK": "BSF"}, "ENERGY": {"NE": [1]}},
            "BSFKK",
        ),
        (
            "BSFKK",
            {"TASK": {"TASK": "BSF"}, "ENERGY": {"NE": [200]}},
            "BSFEK",
        ),
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
        ("SCF", {"TASK": {"TASK": "SCF"}}),
        ("ARPES", {"TASK": {"TASK": "ARPES"}}),
        ("SPLEED", {"TASK": {"TASK": "SPLEED"}}),
        ("BSFEK", {"TASK": {"TASK": "BSF"}, "ENERGY": {"NE": [200]}}),
        ("BSFKK", {"TASK": {"TASK": "BSF"}, "ENERGY": {"NE": [1]}}),
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


def test_read_for_task_rejects_defaulted_task_discriminator(monkeypatch):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["SCF"]
    monkeypatch.setattr(InputParameters, "read_from_file", _reader_setting({}))

    with pytest.raises(ValueError, match="does not explicitly identify a task"):
        registry.read_for_task(task, "ignored.inp")
