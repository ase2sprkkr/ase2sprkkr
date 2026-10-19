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


def test_read_for_task_uses_definition_reader_without_expected_task_context(monkeypatch):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["BSFEK"]
    definition = registry["BSF"]
    parameters = InputParameters.create("BSFEK")
    parameters.ENERGY.NE = [200]
    parameters._requested_task_name = None
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
            {"inputfile": "historical.inp"},
        )
    ]
    assert out.requested_task_name == "BSFEK"


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


@pytest.mark.parametrize(
    ("requested", "filename", "task_name"),
    (
        ("SPLEED", "arpes.in", "ARPES"),
        ("BSFKK", "bsf.inp", "BSFEK"),
    ),
)
def test_read_for_task_rejects_real_shared_schema_mismatch(
    requested, filename, task_name
):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry[requested]

    with pytest.raises(
        TaskMismatchError,
        match=rf"read for task '{requested}'.*task '{task_name}'",
    ):
        registry.read_for_task(task, EXAMPLES / filename)


@pytest.mark.parametrize(
    ("requested", "filename"),
    (("ARPES", "arpes.in"), ("BSFEK", "bsf.inp")),
)
def test_read_for_task_accepts_real_matching_shared_schema(requested, filename):
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry[requested]

    parameters = registry.read_for_task(task, EXAMPLES / filename)

    assert parameters.task_name == requested
    assert parameters.task_definition.name == requested
    assert parameters.requested_task_name == requested


def test_expected_task_does_not_supply_missing_shared_task_selector():
    text = (EXAMPLES / "arpes.in").read_text().replace(
        "TASK     ARPES\n", "TASK\n", 1
    )
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["SPLEED"]

    with pytest.raises(
        TaskMismatchError,
        match=r"read for task 'SPLEED'.*task None",
    ):
        registry.read_for_task(task, StringIO(text))


def test_expected_task_does_not_supply_missing_bsf_mode():
    text = (EXAMPLES / "bsf.inp").read_text().replace("\tNE={200}\n", "", 1)
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["BSFEK"]

    with pytest.raises((pp.ParseBaseException, DataValidityError, TaskMismatchError)):
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
    if requested == "BSFEK":
        parameters.TASK.KPATH = 3
    else:
        parameters.TASK.NK1 = 10
        parameters.TASK.NK2 = 10
        parameters.TASK.K1 = [0.0, 0.0, 0.0]
        parameters.TASK.K2 = [1.0, 1.0, 0.0]

    text = parameters.to_string(validate=False)

    assert f"\tNE={{{ne}}}" in text

    task = InputParameters.task_registry[requested]
    parsed = InputParameters.input_parameters_registry.read_for_task(
        task, StringIO(text)
    )
    assert parsed.task_name == requested
