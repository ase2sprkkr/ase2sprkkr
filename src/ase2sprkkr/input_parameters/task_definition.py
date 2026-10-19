"""Task-level metadata and registries for SPR-KKR calculations."""

from __future__ import annotations

from dataclasses import dataclass, replace
import os
import string
from typing import Callable, Iterable, Optional, Tuple, Union

from ..common.doc import process_input_parameters_definition
from ..common.module_registry import (
    ModuleRegistry,
    RegistryDefinitionError,
    RegistryEntryNotFound,
)


_FilePath = Union[str, os.PathLike[str]]
_FileTemplate = Union[_FilePath, Callable[[object], Optional[_FilePath]]]


def _normalize_identifier(value: str, field: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{field} must be a string")
    value = value.strip()
    if not value:
        raise ValueError(f"{field} must not be empty")
    return value.upper()


def _validate_executable(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("executable must be a string")
    if not value.strip():
        raise ValueError("executable must not be empty")
    return value


@dataclass(frozen=True, eq=False)
class TaskFile(os.PathLike):
    """One task file declaration, or its resolved path with retained metadata.

    A pathname may be a template or callback; a callback returning None means
    the file is not declared for that parameter set.  An optional declared
    file may be absent even after successful execution.
    """

    pathname: _FileTemplate
    type: Optional[str] = None
    info: Optional[str] = None
    optional: bool = False

    def __post_init__(self):
        path = self.pathname
        if not callable(path):
            path = os.fspath(path)
            if not isinstance(path, str):
                raise TypeError("TaskFile pathname must be a string")
            if not path:
                raise ValueError("TaskFile pathname must not be empty")
            object.__setattr__(self, "pathname", path)
        if self.type is not None:
            object.__setattr__(self, "type", _normalize_identifier(self.type, "type"))
        if self.info is not None and not isinstance(self.info, str):
            raise TypeError("TaskFile info must be a string or None")
        if type(self.optional) is not bool:
            raise TypeError("TaskFile optional must be a bool")

    def __fspath__(self):
        if callable(self.pathname):
            raise TypeError("Resolve the TaskFile callback before using its path")
        return self.pathname

    def __str__(self):
        return str(self.pathname)

    def __eq__(self, other):
        if isinstance(other, TaskFile):
            return (self.pathname, self.type, self.info, self.optional) == (
                other.pathname, other.type, other.info, other.optional
            )
        if isinstance(other, str):
            return (not self.optional and self.type is None and self.info is None
                    and self.pathname == other)
        return NotImplemented

    def __hash__(self):
        if not self.optional and self.type is None and self.info is None:
            return hash(self.pathname)
        return hash((self.pathname, self.type, self.info, self.optional))

    @property
    def effective_info(self):
        """Use the explicitly supplied info, or fall back to the result format."""
        from ..result_files.definition_registry import result_file_definition_registry
        return result_file_definition_registry.info_for_task_file(self)


class _MissingTaskFileValue(ValueError):
    """Optional filename template field is unset."""


class _TaskFileFormatter(string.Formatter):
    def __init__(self, parameters):
        self.parameters = parameters

    def get_field(self, field_name, args, kwargs):
        if field_name == "dataset":
            value = self.parameters.CONTROL.DATASET()
        else:
            value = self.parameters
            for key in field_name.split("."):
                value = value[key]
            if callable(value):
                value = value()
        if value is None:
            raise _MissingTaskFileValue(field_name)
        if value == "":
            raise ValueError(f"TaskFile template field {field_name!r} must not be empty")
        return value, field_name


def _normalize_task_files(files: Iterable[_FileTemplate], field: str) -> Tuple[TaskFile, ...]:
    if isinstance(files, (str, bytes, os.PathLike)) or callable(files):
        raise TypeError(f"{field} must be an iterable of path-like strings or TaskFile declarations")
    try:
        iterator = iter(files)
    except TypeError:
        raise TypeError(f"{field} must be an iterable of path-like strings or TaskFile declarations") from None

    normalized = []
    for index, item in enumerate(iterator):
        try:
            normalized.append(item if isinstance(item, TaskFile) else TaskFile(item))
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{field}[{index}]: {exc}") from None
    return tuple(normalized)


@dataclass(frozen=True, init=False)
class TaskDefinition:
    """Immutable declarative metadata identifying one SPR-KKR task.

    Convenience methods delegate task-bound operations to their owning service
    registries; resolution and parsing policy does not live in this value object.
    """

    input_parameters_definition_name: str
    name: str
    executable: str
    mpi: bool
    output_reader_name: Optional[str]
    input_files: Tuple[TaskFile, ...]
    output_files: Tuple[TaskFile, ...]
    cache_files: Tuple[TaskFile, ...]

    def __init__(
        self,
        input_parameters_definition_name: str,
        *,
        name: Optional[str] = None,
        executable: str = "kkrscf",
        mpi: bool = True,
        output_reader_name: Optional[str] = None,
        output_reader: Optional[str] = None,
        input_files: Iterable[_FileTemplate] = (),
        output_files: Iterable[_FileTemplate] = (),
        cache_files: Iterable[_FileTemplate] = (),
    ) -> None:
        schema_name = _normalize_identifier(
            input_parameters_definition_name, "input_parameters_definition_name"
        )
        task_name = (
            schema_name if name is None else _normalize_identifier(name, "name")
        )
        executable = _validate_executable(executable)
        if type(mpi) is not bool:
            raise TypeError("mpi must be a bool")
        if output_reader_name is not None and output_reader is not None:
            raise TypeError(
                "Specify only one of output_reader_name and output_reader"
            )
        if output_reader_name is None:
            output_reader_name = output_reader
        if output_reader_name is not None:
            if not isinstance(output_reader_name, str):
                raise TypeError("output_reader must be a string or None")
            output_reader_name = _normalize_identifier(
                output_reader_name, "output_reader"
            )

        object.__setattr__(self, "input_parameters_definition_name", schema_name)
        object.__setattr__(self, "name", task_name)
        object.__setattr__(self, "executable", executable)
        object.__setattr__(self, "mpi", mpi)
        object.__setattr__(self, "output_reader_name", output_reader_name)
        for kind, values in (
            ("input", input_files),
            ("output", output_files),
            ("cache", cache_files),
        ):
            object.__setattr__(
                self, f"{kind}_files",
                _normalize_task_files(values, f"{kind}_files"),
            )

    def output_reader(
        self,
        input_parameters=None,
        calculator=None,
        print_output=False,
        directory=None,
    ):
        """Create this task's output reader through the reader registry."""
        from ..outputs.reader_registry import reader_registry

        return reader_registry.reader_for_task(
            self,
            input_parameters=input_parameters,
            calculator=calculator,
            print_output=print_output,
            directory=directory,
        )

    def read_output_from_file(
        self,
        filename,
        *,
        input_parameters=None,
        calculator=None,
        directory=None,
    ):
        """Read a completed output file through the task's registered reader."""
        from ..outputs.reader_registry import reader_registry

        return reader_registry.read_for_task(
            self,
            filename,
            input_parameters=input_parameters,
            calculator=calculator,
            directory=directory,
        )


class TaskMismatchError(ValueError):
    """Parsed input values represent a different task than expected."""

    def __init__(self, expected_task_name, task_name):
        self.expected_task_name = expected_task_name
        self.task_name = task_name
        super().__init__(
            f"Input file read for task {expected_task_name!r} represents "
            f"task {task_name!r}"
        )


class InputParametersRegistry(ModuleRegistry):
    """Registry and factory service for input-parameter schemas and values."""

    def __init__(self, package):
        super().__init__(
            package,
            "input_parameters_definition",
            factory=True,
            include=lambda name: not name.startswith("_") and name != "sections",
        )

    def prepare_value(self, name, module, value):
        from .input_parameters_definitions import InputParametersDefinition

        if not isinstance(value, InputParametersDefinition):
            raise RegistryDefinitionError(
                f"{module.__name__}.input_parameters_definition must produce "
                f"InputParametersDefinition, got {type(value).__name__}"
            )
        if str(value.name).upper() != name:
            raise RegistryDefinitionError(
                f"Input-parameter definition in {module.__name__} is named "
                f"{value.name}, but the module defines {name}"
            )

        task = value._task_definition
        if task is not None:
            if not isinstance(task, TaskDefinition):
                raise RegistryDefinitionError(
                    f"{module.__name__}._task_definition must be TaskDefinition, "
                    f"got {type(task).__name__}"
                )
            if task.name != name:
                raise RegistryDefinitionError(
                    f"Canonical task {task.name} does not match input schema {name}"
                )
            if task.input_parameters_definition_name != name:
                raise RegistryDefinitionError(
                    f"Canonical task {task.name} references input schema "
                    f"{task.input_parameters_definition_name} instead of {name}"
                )

        process_input_parameters_definition(module, value)
        return value

    def definition_for_task(self, task):
        """Return the schema identified by declarative task metadata."""
        if not isinstance(task, TaskDefinition):
            raise TypeError("task must be a TaskDefinition")
        definition = self[task.input_parameters_definition_name]
        return definition

    def create_for_task(
        self,
        task,
        *,
        inputfile=None,
        outputfile=False,
    ):
        """Create InputParameters for a resolved TaskDefinition."""
        from .input_parameters import InputParameters

        definition = self.definition_for_task(task)
        return InputParameters(
            definition=definition,
            requested_task_name=task.name,
            inputfile=inputfile,
            outputfile=outputfile,
        )

    def read_for_task(
        self,
        task,
        file,
        *,
        input_definition=None,
        allow_dangerous=False,
    ):
        """Read input with an exact schema and require it to represent *task*.

        ``input_definition`` is primarily for historical reconstruction: when
        omitted, the task's current declared schema is used. Passing a schema
        name or definition selects that exact schema. An explicit definition
        object may intentionally be unregistered, which supports experimental
        or custom parsers while retaining the task-identity check.

        The expected task is deliberately *not* supplied as parsing context.
        Shared-schema defaults must therefore identify the file independently
        of what the caller expected. ``requested_task_name`` is provenance and
        is attached only after the parsed identity has been verified.
        """
        from .input_parameters_definitions import InputParametersDefinition

        if input_definition is None:
            definition = self.definition_for_task(task)
        elif isinstance(input_definition, str):
            definition = self[input_definition]
        elif isinstance(input_definition, InputParametersDefinition):
            definition = input_definition
        else:
            raise TypeError(
                "input_definition must be an InputParametersDefinition, string, or None"
            )

        out = definition.read_from_file(
            file,
            allow_dangerous=allow_dangerous,
            inputfile=file,
        )

        if task.executable == "kkrspec" and not out.TASK.TASK.is_set():
            raise TaskMismatchError(task.name, None)

        task_name = out.task_name
        if task_name != task.name:
            raise TaskMismatchError(task.name, task_name)

        out._requested_task_name = task.name
        return out

    def create_from_result(self, task, result):
        """Create InputParameters for *task* using snapshot defaults from *result*."""
        definition = self.definition_for_task(task)
        out = self.create_for_task(task)
        definition.initialize_from_result(out, result)
        return out

    @staticmethod
    def _resolve_file_templates(parameters, directory, templates):
        """Resolve templates to concrete TaskFile values without losing metadata."""
        formatter = _TaskFileFormatter(parameters)
        files = []
        for value in templates:
            declaration = value if isinstance(value, TaskFile) else TaskFile(value)
            pathname = declaration.pathname
            if callable(pathname):
                pathname = pathname(parameters)
            if pathname is None:
                continue
            if isinstance(pathname, TaskFile):
                declaration = pathname
                pathname = declaration.pathname
            path = os.fspath(pathname)
            if not isinstance(path, str):
                raise TypeError("File callback must return a string/path or None")
            if not path:
                raise ValueError("TaskFile pathname must not be empty")
            try:
                path = formatter.format(path)
            except _MissingTaskFileValue:
                continue
            resolved = path if os.path.isabs(path) else os.path.join(directory, path)
            files.append(replace(declaration, pathname=resolved))
        return tuple(files)

    def files_for_task(self, task, parameters, directory, kind):
        """Combine declarative task templates with schema-specific file rules."""
        definition = self.definition_for_task(task)
        declared = self._resolve_file_templates(
            parameters, directory, getattr(task, f"{kind}_files")
        )
        schema_files = getattr(definition, f"{kind}_files")(parameters, directory)
        # Keep schema extension hooks; shared resolution lives in this registry.
        # A TaskFile retains metadata, while plain paths remain valid shorthands.
        schema_files = tuple(
            value if isinstance(value, TaskFile) else TaskFile(value)
            for value in schema_files
        )
        return (*declared, *schema_files)


class TaskRegistry:
    """Authoritative namespace and resolver for SPR-KKR task identities.

    Task declarations intentionally have two sources. A schema with exactly
    one canonical task may co-locate its :class:`TaskDefinition` on the
    ``InputParametersDefinition`` as ``_task_definition``; this keeps the
    common one-schema/one-task case simple and declarative. Tasks sharing one
    input schema are declared separately in the extra task package and refer
    to that schema by name.

    The registry owns the task namespace and resolution rules, not necessarily
    the physical location of each declaration. Consumers should resolve task
    identities through this registry rather than reading ``_task_definition``
    directly.
    """

    normalize_name = staticmethod(ModuleRegistry.normalize_name)

    def __init__(self, input_parameters_registry, extra_tasks_package):
        self.input_parameters_registry = input_parameters_registry
        self.extra_tasks = ModuleRegistry(
            extra_tasks_package,
            "task_definition",
            factory=True,
        )

    def _input_task(self, name):
        try:
            definition = self.input_parameters_registry[name]
        except RegistryEntryNotFound:
            return None
        return definition._task_definition

    def _extra_task(self, name):
        try:
            task = self.extra_tasks[name]
        except RegistryEntryNotFound:
            return None
        if not isinstance(task, TaskDefinition):
            raise RegistryDefinitionError(
                f"Extra task {name} must be TaskDefinition, got {type(task).__name__}"
            )
        if task.name != name:
            raise RegistryDefinitionError(
                f"Extra task {task.name} is registered as {name}"
            )
        try:
            self.input_parameters_registry.definition_for_task(task)
        except RegistryEntryNotFound:
            raise RegistryDefinitionError(
                f"Task {task.name} references unknown input schema "
                f"{task.input_parameters_definition_name}"
            ) from None
        return task

    def __getitem__(self, name):
        name = self.normalize_name(name)
        input_task = self._input_task(name)
        extra_task = self._extra_task(name)
        if input_task is not None and extra_task is not None:
            raise RegistryDefinitionError(
                f"Task {name} is defined both by an input schema and explicitly"
            )
        if extra_task is not None:
            return extra_task
        if input_task is not None:
            return input_task
        raise RegistryEntryNotFound(name)

    def __contains__(self, name):
        try:
            self[name]
        except RegistryEntryNotFound:
            return False
        return True

    def names(self):
        names = set(self.extra_tasks.names())
        for schema_name in self.input_parameters_registry.names():
            task = self.input_parameters_registry[schema_name]._task_definition
            if task is None:
                continue
            name = task.name
            if name in names:
                raise RegistryDefinitionError(
                    f"Task {name} is defined both by an input schema and explicitly"
                )
            names.add(name)
        return tuple(sorted(names))

    def values(self):
        return tuple(self[name] for name in self.names())

    def task_for_parameters(self, parameters):
        """Resolve the task represented by concrete parameter values."""
        name = parameters.task_name
        return None if name is None else self[name]
