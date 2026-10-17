"""Task-level metadata and registries for SPR-KKR calculations."""

from __future__ import annotations

from dataclasses import dataclass
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
_FilePathResolver = Callable[[object], Optional[_FilePath]]


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


def _normalize_pathname(pathname: Union[_FilePath, _FilePathResolver]):
    if callable(pathname):
        try:
            hash(pathname)
        except TypeError:
            raise TypeError("pathname callback must be hashable") from None
        return pathname
    try:
        pathname = os.fspath(pathname)
    except TypeError:
        raise TypeError(
            "pathname must be a string, os.PathLike, or callable"
        ) from None
    if not isinstance(pathname, str):
        raise TypeError("pathname must resolve to a string")
    if not pathname:
        raise ValueError("pathname must not be empty")
    return pathname


@dataclass(frozen=True, init=False)
class TaskFile:
    """Declarative description of one file used or produced by a task.

    ``pathname`` may be a path template or a callback receiving concrete input
    parameters. Templates support ``{dataset}`` and dotted input-value paths,
    for example ``{CONTROL.POTFIL}``. A callback returns a path or ``None``;
    ``None`` means that the file is absent for those concrete parameters.
    Empty pathnames are invalid.

    ``type`` is an optional strict reference to a structured
    :class:`ResultFileDefinition`; parser/format metadata remains owned by the
    result-file registry. ``info`` is only a task-specific description override.
    """

    pathname: Union[str, _FilePathResolver]
    type: Optional[str]
    info: Optional[str]

    def __init__(
        self,
        pathname: Union[_FilePath, _FilePathResolver],
        *,
        type: Optional[str] = None,
        info: Optional[str] = None,
    ) -> None:
        pathname = _normalize_pathname(pathname)
        if type is not None:
            type = _normalize_identifier(type, "type")
        if info is not None and not isinstance(info, str):
            raise TypeError("info must be a string or None")
        object.__setattr__(self, "pathname", pathname)
        object.__setattr__(self, "type", type)
        object.__setattr__(self, "info", info)


def _normalize_files(
    files: Iterable[Union[TaskFile, _FilePath, _FilePathResolver]], field: str
) -> Tuple[TaskFile, ...]:
    if isinstance(files, (str, bytes, os.PathLike, TaskFile)) or callable(files):
        raise TypeError(
            f"{field} must be an iterable of TaskFile or path-like values, "
            "not a single file"
        )
    try:
        iterator = iter(files)
    except TypeError:
        raise TypeError(
            f"{field} must be an iterable of TaskFile or path-like values"
        ) from None

    normalized = []
    for index, item in enumerate(iterator):
        if isinstance(item, TaskFile):
            normalized.append(item)
            continue
        try:
            normalized.append(TaskFile(item))
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{field}[{index}]: {exc}") from None
    return tuple(normalized)


def _coalesce_file_declarations(files, templates, kind):
    if files is not None and templates:
        raise TypeError(
            f"Specify only one of {kind}_files and {kind}_file_templates"
        )
    if files is None:
        files = templates
    return _normalize_files(files, f"{kind}_files")


def _parameter_value(parameters, path):
    value = parameters
    for name in path.split("."):
        value = value[name]
    return value() if callable(value) else value


class _MissingTaskFileValue(ValueError):
    pass


class _TaskFileFormatter(string.Formatter):
    def __init__(self, parameters):
        self.parameters = parameters

    def get_field(self, field_name, args, kwargs):
        if field_name == "dataset":
            value = self.parameters.CONTROL.DATASET()
        else:
            value = _parameter_value(self.parameters, field_name)
        if value is None:
            raise _MissingTaskFileValue(field_name)
        if value == "":
            raise ValueError(
                f"TaskFile template field {field_name!r} must not be empty"
            )
        return value, field_name


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
        input_files: Optional[
            Iterable[Union[TaskFile, _FilePath, _FilePathResolver]]
        ] = None,
        output_files: Optional[
            Iterable[Union[TaskFile, _FilePath, _FilePathResolver]]
        ] = None,
        cache_files: Optional[
            Iterable[Union[TaskFile, _FilePath, _FilePathResolver]]
        ] = None,
        input_file_templates: Iterable[_FilePath] = (),
        output_file_templates: Iterable[_FilePath] = (),
        cache_file_templates: Iterable[_FilePath] = (),
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
        object.__setattr__(
            self,
            "input_files",
            _coalesce_file_declarations(input_files, input_file_templates, "input"),
        )
        object.__setattr__(
            self,
            "output_files",
            _coalesce_file_declarations(output_files, output_file_templates, "output"),
        )
        object.__setattr__(
            self,
            "cache_files",
            _coalesce_file_declarations(cache_files, cache_file_templates, "cache"),
        )

    @staticmethod
    def _template_view(files):
        return tuple(file.pathname for file in files)

    @property
    def input_file_templates(self):
        """Compatibility view of declared input pathnames."""
        return self._template_view(self.input_files)

    @property
    def output_file_templates(self):
        """Compatibility view of declared output pathnames."""
        return self._template_view(self.output_files)

    @property
    def cache_file_templates(self):
        """Compatibility view of declared cache pathnames."""
        return self._template_view(self.cache_files)

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


def _validate_task_file_types(task):
    from ..result_files.definition_registry import result_file_definition_registry

    for kind in ("input", "output", "cache"):
        for file in getattr(task, f"{kind}_files"):
            if (
                file.type is not None
                and file.type not in result_file_definition_registry
            ):
                raise RegistryDefinitionError(
                    f"Task {task.name} declares unknown result file type "
                    f"{file.type!r}"
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
            _validate_task_file_types(task)

        process_input_parameters_definition(module, value)
        return value

    def definition_for_task(self, task):
        """Return the registered schema referenced by task metadata.

        The task itself need not be registered or equal to the schema's canonical
        task. This intentionally allows callers to supply compatible ad-hoc task
        metadata, for example to work around an executable or file-manifest
        incompatibility. Built-in declarations are validated separately.
        """
        if not isinstance(task, TaskDefinition):
            raise TypeError("task must be a TaskDefinition")
        return self[task.input_parameters_definition_name]

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
                "input_definition must be an InputParametersDefinition, "
                "string, or None"
            )

        out = definition.read_from_file(
            file,
            allow_dangerous=allow_dangerous,
            inputfile=file,
        )

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
    def _resolve_declared_files(parameters, directory, files):
        formatter = _TaskFileFormatter(parameters)
        resolved = []
        for file in files:
            pathname = (
                file.pathname(parameters)
                if callable(file.pathname)
                else file.pathname
            )
            if pathname is None:
                continue
            try:
                pathname = os.fspath(pathname)
            except TypeError:
                raise TypeError(
                    "TaskFile pathname callback must return a string, "
                    "os.PathLike, or None"
                ) from None
            if not isinstance(pathname, str):
                raise TypeError("TaskFile pathname must resolve to a string")
            if not pathname:
                raise ValueError("TaskFile pathname must not be empty")
            try:
                pathname = formatter.format(pathname)
            except _MissingTaskFileValue:
                continue
            if not os.path.isabs(pathname):
                pathname = os.path.join(directory, pathname)
            resolved.append(pathname)
        return tuple(resolved)

    def files_for_task(self, task, parameters, directory, kind):
        """Resolve declarative file metadata for a concrete task invocation."""
        if kind not in ("input", "output", "cache"):
            raise ValueError(f"Unknown task file kind {kind!r}")
        self.definition_for_task(task)
        return self._resolve_declared_files(
            parameters,
            directory,
            getattr(task, f"{kind}_files"),
        )


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
        _validate_task_file_types(task)
        try:
            definition = self.input_parameters_registry.definition_for_task(task)
        except RegistryEntryNotFound:
            raise RegistryDefinitionError(
                f"Task {task.name} references unknown input schema "
                f"{task.input_parameters_definition_name}"
            ) from None
        canonical = definition._task_definition
        if canonical is not None and canonical != task:
            raise RegistryDefinitionError(
                f"Extra task {task.name} references input schema "
                f"{task.input_parameters_definition_name}, which is bound to "
                f"different canonical task metadata {canonical!r}"
            )
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
