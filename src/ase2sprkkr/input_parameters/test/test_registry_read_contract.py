from io import StringIO
from pathlib import Path
import re

import pyparsing as pp
import pytest

from ...common.warnings import DataValidityError
from ..input_parameters import InputParameters
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskMismatchError, TaskRegistry


EXAMPLES = Path(__file__).parent.parent / "examples"


def _reader_setting(values):
    def read_from_file(self, _file, *, allow_dangerous=False):
        self.set(values, retain_invalid="all", report_invalid="ignore")

    return read_from_file


@pytest.mark.parametrize(
    ("requested", "values", "task_name"),
    (
        ("ARPES", {"TASK": {"TASK": "SPLEED"}}, "SPLEED"),
        ("SPLEED", {"TASK": {"TASK": "ARPES"}}, "ARPES"),
        ("BSFEK", {"ENERGY": {"NE": [1]}}, "BSFKK"),
        ("BSFKK", {"ENERGY": {"NE": [200]}}, "BSFEK"),
    ),
)
def test_read_for_task_rejects_mismatched_task(
    monkeypatch, requested, values, task_name
):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry[requested]
    monkeypatch.setattr(InputParameters, "read_from_file", _reader_setting(values))

    with pytest.raises(
        TaskMismatchError,
        match=rf"read for task '{requested}'.*task '{task_name}'",
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
def test_read_for_task_accepts_matching_task(monkeypatch, requested, values):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry[requested]
    monkeypatch.setattr(InputParameters, "read_from_file", _reader_setting(values))

    parameters = registry.read_for_task(task, "ignored.inp")

    assert parameters.task_name == requested
    assert parameters.task_definition.name == requested


def test_task_name_is_schema_derived_without_task_registry_lookup(monkeypatch):
    scf = InputParameters.create("SCF")
    bsf = InputParameters.create("BSFEK")

    def must_not_resolve(_self, _name):
        raise AssertionError("task_name must not resolve TaskDefinition")

    monkeypatch.setattr(TaskRegistry, "__getitem__", must_not_resolve)

    assert scf._definition.task_definition_name(scf) == "SCF"
    assert scf.task_name == "SCF"
    assert bsf.task_name == "BSFEK"


def test_read_for_task_uses_definition_reader(monkeypatch):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["BSFEK"]
    definition = registry["BSF"]
    parameters = InputParameters.create("BSFEK")
    calls = []

    def read_from_file(file, *, allow_dangerous=False, **kwargs):
        calls.append((file, allow_dangerous, kwargs))
        return parameters

    monkeypatch.setattr(definition, "read_from_file", read_from_file)

    out = registry.read_for_task(
        task,
        "historical.inp",
        input_definition=definition,
        allow_dangerous=True,
    )

    assert out is parameters
    assert calls == [
        (
            "historical.inp",
            True,
            {
                "requested_task_name": "BSFEK",
                "inputfile": "historical.inp",
            },
        )
    ]


def test_read_for_task_accepts_explicit_unregistered_schema(monkeypatch):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["SCF"]
    definition = InputParametersDefinition("experimental", [])
    parameters = InputParameters.create("SCF")

    monkeypatch.setattr(
        definition,
        "read_from_file",
        lambda _file, **_kwargs: parameters,
    )

    assert (
        registry.read_for_task(
            task,
            "experimental.inp",
            input_definition=definition,
        )
        is parameters
    )


def _write_parameters(tmp_path, task_name):
    parameters = InputParameters.create(task_name)
    path = tmp_path / f"{task_name.lower()}.inp"
    path.write_text(parameters.to_string(validate=False))
    return path


@pytest.mark.parametrize(
    ("requested", "file_task", "task_name"),
    (
        ("SPLEED", "ARPES", "ARPES"),
        ("BSFKK", "BSFEK", "BSFEK"),
    ),
)
def test_read_for_task_rejects_real_shared_schema_mismatch(
    tmp_path, requested, file_task, task_name
):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry[requested]
    filename = _write_parameters(tmp_path, file_task)

    with pytest.raises(
        TaskMismatchError,
        match=rf"read for task '{requested}'.*task '{task_name}'",
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


def test_arpes_task_selector_is_required_in_parsed_input():
    text = (EXAMPLES / "arpes.in").read_text().replace(
        "TASK     ARPES\n", "TASK\n", 1
    )
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["ARPES"]

    with pytest.raises((pp.ParseBaseException, DataValidityError)):
        registry.read_for_task(task, StringIO(text))


def test_bsf_energy_count_is_required_in_parsed_input():
    text = (EXAMPLES / "bsf.inp").read_text().replace("\tNE={200}\n", "", 1)
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["BSFEK"]

    with pytest.raises((pp.ParseBaseException, DataValidityError)):
        registry.read_for_task(task, StringIO(text))


@pytest.mark.parametrize("requested", ("ARPES", "AIPES", "SPLEED", "BAND"))
def test_shared_task_default_is_materialized_when_writing(requested):
    parameters = InputParameters.create(requested)
    parameters.CONTROL.POTFIL = "case.pot"

    text = parameters.to_string(validate=False)

    assert re.search(rf"(?m)^TASK\s+{requested}\s*$", text)

    task = InputParameters.task_registry[requested]
    parsed = InputParameters.input_parameters_registry.read_for_task(
        task, StringIO(text)
    )
    assert parsed.task_name == requested


@pytest.mark.parametrize(
    ("requested", "ne"), (("BSFEK", 200), ("BSFKK", 1))
)
def test_bsf_task_default_is_materialized_when_writing(requested, ne):
    parameters = InputParameters.create(requested)
    parameters.CONTROL.POTFIL = "case.pot"

    text = parameters.to_string(validate=False)

    assert f"\tNE={{{ne}}}" in text

    task = InputParameters.task_registry[requested]
    parsed = InputParameters.input_parameters_registry.read_for_task(
        task, StringIO(text)
    )
    assert parsed.task_name == requested