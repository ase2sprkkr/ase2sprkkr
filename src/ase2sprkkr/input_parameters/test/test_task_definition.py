import pickle
from types import SimpleNamespace

import pytest

from ..input_parameters import InputParameters
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TASKS, TaskDefinition


def test_runtime_input_parameters_resolve_effective_task_through_schema():
    parameters = InputParameters.create("SCF")
    task = TASKS.definition("SCF")

    assert isinstance(parameters.task_definition, TaskDefinition)
    assert parameters.task_definition is task
    assert parameters.definition is task.input_parameters_definition
    assert parameters.definition.task_definition(parameters) is task
    assert parameters.requested_task_definition is task
    assert parameters.task_name == "SCF"
    assert parameters.requested_task_name == "SCF"


def test_schema_can_create_runtime_object_without_requested_task():
    task = TASKS.definition("DOS")
    parameters = task.input_parameters_definition.create_object()

    assert isinstance(parameters, InputParameters)
    assert parameters.definition is task.input_parameters_definition
    assert parameters.task_definition is task
    assert parameters.requested_task_definition is None


def test_schema_can_read_runtime_object_without_requested_task():
    task = TASKS.definition("DOS")
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.POTFIL = "test.pot"
    text = parameters.to_string(validate=False)

    parsed = task.input_parameters_definition.read_from_string(
        text,
        allow_dangerous=True,
    )

    assert parsed.definition is task.input_parameters_definition
    assert parsed.task_definition is task
    assert parsed.requested_task_definition is None


def test_task_create_records_requested_task_definition():
    task = TASKS.definition("DOS")
    parameters = task.create_input_parameters()

    assert parameters.requested_task_definition is task
    assert parameters.task_definition is task


def test_copy_preserves_requested_task_definition():
    parameters = InputParameters.create("DOS")
    requested = parameters.requested_task_definition
    parameters.CONTROL.DATASET = "copy-test"

    copied = parameters.copy(copy_values=True)

    assert copied.requested_task_definition is requested
    assert copied.task_definition is parameters.task_definition
    assert copied.definition is parameters.definition
    assert copied.CONTROL.DATASET() == "copy-test"


def test_standard_schema_keeps_one_direct_task_reference():
    schema = InputParametersDefinition("single", [])
    first = TaskDefinition(schema, name="FIRST")

    assert schema._task_definition is first
    assert schema.task_definition(SimpleNamespace()) is first

    with pytest.raises(ValueError, match="already associated"):
        TaskDefinition(schema, name="SECOND")


def test_dynamic_schema_resolves_through_external_registry(monkeypatch):
    class SharedDefinition(InputParametersDefinition):
        dynamic_task = True

        def task_definition_name(self, parameters):
            return parameters._selected_task_name

    schema = SharedDefinition("shared", [])
    first = TaskDefinition(schema, name="FIRST")
    second = TaskDefinition(schema, name="SECOND")
    registry = {"FIRST": first, "SECOND": second}
    monkeypatch.setattr(InputParameters, "task_registry", registry)

    parameters = first.create_input_parameters()
    parameters._selected_task_name = "FIRST"

    assert schema._task_definition is None
    assert parameters.requested_task_definition is first
    assert parameters.task_definition is first

    parameters._selected_task_name = "SECOND"

    assert parameters.requested_task_definition is first
    assert parameters.task_definition is second
    assert parameters.requested_task_name == "FIRST"
    assert parameters.task_name == "SECOND"


def test_dynamic_schema_requires_external_registry(monkeypatch):
    class SharedDefinition(InputParametersDefinition):
        dynamic_task = True

        def task_definition_name(self, parameters):
            return "FIRST"

    schema = SharedDefinition("shared", [])
    TaskDefinition(schema, name="FIRST")
    monkeypatch.setattr(InputParameters, "task_registry", None)
    parameters = schema.create_object()

    with pytest.raises(RuntimeError, match="requires an external task registry"):
        _ = parameters.task_definition


def test_schema_copy_drops_direct_task_reference():
    task = TASKS.definition("SCF")
    schema = task.input_parameters_definition

    copied = schema.copy()

    assert copied is not schema
    assert copied._task_definition is None
    with pytest.raises(RuntimeError, match="must override task_definition_name"):
        copied.task_definition(SimpleNamespace(task_registry=None))

    copied_task = TaskDefinition(copied, name="SCF-COPY")
    assert copied.task_definition(SimpleNamespace()) is copied_task


def test_spec_schema_resolves_spleed_through_external_registry(monkeypatch):
    schema = TASKS.definition("ARPES").input_parameters_definition.copy()
    arpes = TaskDefinition(schema, name="ARPES", executable="kkrspec")
    spleed = TaskDefinition(schema, name="SPLEED", executable="kkrspec")
    monkeypatch.setattr(
        InputParameters,
        "task_registry",
        {"ARPES": arpes, "SPLEED": spleed},
    )
    parameters = arpes.create_input_parameters()

    assert parameters.task_definition is arpes
    assert parameters.requested_task_definition is arpes

    parameters.TASK.TASK = "SPLEED"

    assert parameters.task_definition is spleed
    assert parameters.requested_task_definition is arpes


def test_bsf_schema_resolves_ek_and_kk_through_external_registry(monkeypatch):
    schema = TASKS.definition("BSF").input_parameters_definition.copy()
    bsfek = TaskDefinition(schema, name="BSFEK", executable="kkrgen")
    bsfkk = TaskDefinition(schema, name="BSFKK", executable="kkrgen")
    monkeypatch.setattr(
        InputParameters,
        "task_registry",
        {"BSFEK": bsfek, "BSFKK": bsfkk},
    )
    parameters = bsfek.create_input_parameters()

    assert parameters.requested_task_definition is bsfek
    assert parameters.task_definition is bsfek

    parameters.ENERGY.NE = [1]

    assert parameters.requested_task_definition is bsfek
    assert parameters.task_definition is bsfkk


def test_task_can_construct_schema_lazily():
    created_for = []

    def create_schema(task):
        created_for.append(task)
        return InputParametersDefinition("lazy", [])

    task = TaskDefinition(
        name="LAZY",
        input_parameters_definition_factory=create_schema,
    )

    assert created_for == []
    schema = task.input_parameters_definition
    assert created_for == [task]
    assert schema.task_definition(SimpleNamespace()) is task
    assert task.input_parameters_definition is schema
    assert task._input_parameters_definition_factory is None
    assert created_for == [task]


def test_materialized_lazy_schema_with_local_factory_is_pickleable():
    def create_schema(_task):
        return InputParametersDefinition("lazy-pickle", [])

    task = TaskDefinition(
        name="LAZY-PICKLE",
        input_parameters_definition_factory=create_schema,
    )
    schema = task.input_parameters_definition

    restored = pickle.loads(pickle.dumps(schema))
    restored_task = restored.task_definition(SimpleNamespace())

    assert restored_task.name == "LAZY-PICKLE"
    assert restored_task.input_parameters_definition is restored
    assert restored_task._input_parameters_definition_factory is None


def test_task_definition_requires_exactly_one_schema_source():
    schema = InputParametersDefinition("test", [])

    with pytest.raises(ValueError, match="Exactly one"):
        TaskDefinition(name="NONE")
    with pytest.raises(ValueError, match="Exactly one"):
        TaskDefinition(
            schema,
            name="BOTH",
            input_parameters_definition_factory=lambda task: schema,
        )
    with pytest.raises(TypeError, match="must be callable"):
        TaskDefinition(
            name="NOT_CALLABLE",
            input_parameters_definition_factory=object(),
        )


def test_lazy_schema_factory_result_is_validated():
    task = TaskDefinition(
        name="INVALID",
        input_parameters_definition_factory=lambda task: object(),
    )

    with pytest.raises(TypeError, match="must return InputParametersDefinition"):
        task.input_parameters_definition


def test_task_policy_is_not_on_input_schema():
    task = TASKS.definition("PHAGEN")
    schema = task.input_parameters_definition

    assert task.executable == "kkrscf"
    assert task.mpi is False
    for attribute in (
        "executable",
        "mpi",
        "result_reader",
        "input_files",
        "output_files",
        "cache_files",
        "input_file_templates",
        "output_file_templates",
        "cache_file_templates",
    ):
        assert not hasattr(schema, attribute)


def test_dos_output_manifest_belongs_to_effective_task(tmp_path):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.output_files(str(tmp_path)) == (str(tmp_path / "Fe.dos"),)


def test_jxc_continuation_policy_belongs_to_effective_task():
    parameters = InputParameters.create("JXC")
    result = SimpleNamespace(
        last_iteration=SimpleNamespace(
            energy=SimpleNamespace(EMIN=lambda: -0.42),
        )
    )

    parameters.task_definition.prepare_parameters_from_result(parameters, result)

    assert parameters.ENERGY.EMIN() == -0.42
