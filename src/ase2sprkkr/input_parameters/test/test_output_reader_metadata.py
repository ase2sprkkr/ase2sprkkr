import pytest

from ..input_parameters import InputParameters
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition
from ...common.module_registry import RegistryEntryNotFound
from ...outputs.reader_registry import reader_registry
from ...outputs.readers.default import DefaultOutputReader
from ...outputs.task_result import KkrOutputReader
from ...sprkkr.executables import executables


def test_builtin_tasks_resolve_readers_and_cover_specialized_readers():
    claimed_readers = set()
    # Enumerating this lazy registry loads every task definition. Check both
    # sides of the reader mapping in the same traversal.
    for task in InputParameters.task_registry.values():
        assert task.executable in executables
        assert issubclass(reader_registry.class_for_task(task), KkrOutputReader)
        claimed_readers.add(
            task.output_reader_name
            if task.output_reader_name is not None
            else executables[task.executable].default_output_reader
        )

    specialized_readers = set(reader_registry.names(include_aliases=False)) - {
        "DEFAULT"
    }
    assert specialized_readers <= claimed_readers


def test_builtin_reader_mappings_are_declared_at_the_right_level():
    expected = {
        "SCF": "SCF",
        "DOS": "DOS",
        "JXC": "JXC",
        "BSFEK": "BSF",
        "BSFKK": "BSF",
        "ARPES": "ARPES",
        "AIPES": "SPEC",
        "BAND": "SPEC",
        "SPLEED": "SPEC",
        "PHAGEN": "DEFAULT",
        "TORQUE": "DEFAULT",
        "XAS": "DEFAULT",
        "GILBERT": "DEFAULT",
    }

    assert reader_registry["DEFAULT"] is DefaultOutputReader
    for task_name, reader_name in expected.items():
        assert reader_registry.class_for_task_name(task_name) is reader_registry[reader_name]


def test_explicit_output_reader_name_is_strict_even_for_custom_executable():
    task = TaskDefinition(
        "SCF",
        name="BROKEN-READER-MAPPING",
        executable="my-custom-sprkkr",
        output_reader="NO-SUCH-READER",
    )

    with pytest.raises(RegistryEntryNotFound):
        reader_registry.class_for_task(task)


def test_explicit_output_reader_does_not_require_registered_executable():
    task = TaskDefinition(
        "SCF",
        name="CUSTOM-SCF",
        executable="my-custom-sprkkr",
        output_reader="SCF",
    )

    assert reader_registry.class_for_task(task) is reader_registry["SCF"]


def test_unknown_executable_without_explicit_reader_is_not_guessed():
    task = TaskDefinition(
        "SCF",
        name="CUSTOM",
        executable="my-custom-sprkkr",
    )

    with pytest.raises(RegistryEntryNotFound, match="my-custom-sprkkr"):
        reader_registry.class_for_task(task)


@pytest.mark.parametrize("parameters_kind", ("matching", "detached", "missing"))
def test_task_output_reader_accepts_matching_detached_or_missing_parameters(
    tmp_path, parameters_kind,
):
    task = InputParameters.task_registry["SCF"]
    if parameters_kind == "matching":
        parameters = InputParameters.create("SCF")
        assert parameters.task_definition_required() is task
        reader = task.output_reader(parameters, directory=tmp_path)
        assert reader.directory == tmp_path
    elif parameters_kind == "detached":
        schema = InputParametersDefinition.definition_from_dict({"TEST": []})
        parameters = schema.create_object()
        assert parameters.task_definition is None
        reader = task.output_reader(parameters)
    else:
        parameters = None
        reader = task.output_reader()

    assert type(reader) is reader_registry["SCF"]
    assert reader.input_parameters is parameters
    assert reader.task_name == "SCF"


def test_task_output_reader_rejects_parameters_for_another_task():
    scf = InputParameters.task_registry["SCF"]
    parameters = InputParameters.create("DOS")

    with pytest.raises(ValueError, match="represent task 'DOS'.*'SCF'"):
        scf.output_reader(parameters)


def test_input_parameters_reader_helpers_remain_thin_task_facades(monkeypatch):
    parameters = InputParameters.create("SCF")
    task = parameters.task_definition_required()
    reader_marker = object()
    result_marker = object()

    def output_reader(
        self,
        input_parameters=None,
        calculator=None,
        print_info=False,
        directory=None,
    ):
        assert self is task
        assert input_parameters is parameters
        assert calculator == "calculator"
        assert print_info is True
        assert directory == "directory"
        return reader_marker

    def read_output_from_file(
        self,
        filename,
        *,
        input_parameters=None,
        calculator=None,
        directory=None,
    ):
        assert self is task
        assert filename == "result.out"
        assert input_parameters is parameters
        assert calculator is None
        assert directory == "directory"
        return result_marker

    monkeypatch.setattr(TaskDefinition, "output_reader", output_reader)
    assert parameters.output_reader(
        calculator="calculator",
        print_info=True,
        directory="directory",
    ) is reader_marker

    monkeypatch.setattr(TaskDefinition, "read_output_from_file", read_output_from_file)
    assert parameters.read_output_from_file(
        "result.out",
        directory="directory",
    ) is result_marker


def test_input_parameters_potential_convenience_is_preserved():
    assert isinstance(InputParameters.potential, property)
