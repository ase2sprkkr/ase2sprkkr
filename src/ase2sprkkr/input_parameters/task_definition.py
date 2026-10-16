"""Task-level metadata and registries for SPR-KKR calculations."""

from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Callable, Iterable, Optional, Tuple, Union

from ..common.doc import process_input_parameters_definition
from ..common.module_registry import (
    ModuleRegistry,
    RegistryDefinitionError,
    RegistryEntryNotFound,
)


_FileTemplate = Union[str, os.PathLike[str]]
_FileResolver = Callable[[object, str], Optional[_FileTemplate]]
_FileCondition = Callable[[object], bool]


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


@dataclass(frozen=True)
class FileSpec:
    """Declarative description of one task-owned file.

    ``path`` is either a static filename template or a resolver receiving the
    concrete input parameters and calculation directory.  ``condition`` can
    suppress the file for parameter combinations for which it is not used.
    Resolvers and conditions are declarations only; materialization is owned by
    :class:`InputParametersRegistry`.
    """

    path: Union[str, _FileResolver]
    condition: Optional[_FileCondition] = None

    def __post_init__(self):
        path = self.path
        if not callable(path):
            try:
                path = os.fspath(path)
            except TypeError:
                raise TypeError(
                    "FileSpec.path must be path-like or callable"
                ) from None
            if not isinstance(path, str):
                raise TypeError("FileSpec.path must resolve to a string")
            object.__setattr__(self, "path", path)
        if self.condition is not None and not callable(self.condition):
            raise TypeError("FileSpec.condition must be callable or None")


def _normalize_file_specs(values, field):
    if isinstance(values, (str, bytes, os.PathLike, FileSpec)):
        raise TypeError(
            f"{field} must be an iterable of FileSpec or path-like values, "
            "not a single file"
        )
    try:
        iterator = iter(values)
    except TypeError:
        raise TypeError(
            f"{field} must be an iterable of FileSpec or path-like values"
        ) from None

    out = []
    for index, value in enumerate(iterator):
        if isinstance(value, FileSpec):
            out.append(value)
            continue
        try:
            out.append(FileSpec(value))
        except TypeError as exc:
            raise TypeError(f"{field}[{index}]: {exc}") from None
    return tuple(out)


@dataclass(frozen=True, init=False)
class TaskFiles:
    """Input, output and reusable-cache declarations of one task."""

    input: Tuple[FileSpec, ...]
    output: Tuple[FileSpec, ...]
    cache: Tuple[FileSpec, ...]

    def __init__(self, *, input=(), output=(), cache=()):
        object.__setattr__(self, "input", _normalize_file_specs(input, "files.input"))
        object.__setattr__(self, "output", _normalize_file_specs(output, "files.output"))
        object.__setattr__(self, "cache", _normalize_file_specs(cache, "files.cache"))


def _legacy_templates(files, kind):
    """Return static unconditional declarations for compatibility/introspection."""
    return tuple(
        spec.path
        for spec in getattr(files, kind)
        if isinstance(spec.path, str) and spec.condition is None
    )


@dataclass(frozen=True, init=False)
class TaskDefinition:
    """Immutable declarative metadata identifying one SPR-KKR task."""

    input_parameters_definition_name: str
    name: str
    executable: str
    mpi: bool
    output_reader: Optional[str]
    files: TaskFiles

    def __init__(
        self,
        input_parameters_definition_name: str,
        *,
        name: Optional[str] = None,
        executable: str = "kkrscf",
        mpi: bool = True,
        output_reader: Optional[str] = None,
        files: Optional[TaskFiles] = None,
        input_file_templates: Iterable[_FileTemplate] = (),
        output_file_templates: Iterable[_FileTemplate] = (),
        cache_file_templates: Iterable[_FileTemplate] = (),
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
        if output_reader is not None:
            if not isinstance(output_reader, str):
                raise TypeError("output_reader must be a string or None")
            output_reader = _normalize_identifier(output_reader, "output_reader")

        legacy = (
            tuple(input_file_templates),
            tuple(output_file_templates),
            tuple(cache_file_templates),
        )
        if files is not None:
            if not isinstance(files, TaskFiles):
                raise TypeError("files must be a TaskFiles or None")
            if any(legacy):
                raise ValueError(
                    "files cannot be combined with legacy *_file_templates"
                )
        else:
            files = TaskFiles(
                input=legacy[0],
                output=legacy[1],
                cache=legacy[2],
            )

        object.__setattr__(self, "input_parameters_definition_name", schema_name)
        object.__setattr__(self, "name", task_name)
        object.__setattr__(self, "executable", executable)
        object.__setattr__(self, "mpi", mpi)
        object.__setattr__(self, "output_reader", output_reader)
        object.__setattr__(self, "files", files)

    @property
    def input_file_templates(self):
        return _legacy_templates(self.files, "input")

    @property
    def output_file_templates(self):
        return _legacy_templates(self.files, "output")

    @property
    def cache_file_templates(self):
        return _legacy_templates(self.files, "cache")


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
        canonical = definition._task_definition
        if canonical is not None and canonical != task:
            raise RegistryDefinitionError(
                f"Task {task.name} references input schema "
                f"{task.input_parameters_definition_name}, which is bound to "
                f"different canonical task metadata {canonical!r}"
            )
        return definition

    def create_for_task(self, task, *, inputfile=None, outputfile=False):
        """Create InputParameters for a resolved TaskDefinition."""
        from .input_parameters import InputParameters

        definition = self.definition_for_task(task)
        return InputParameters(
            definition=definition,
            requested_task_name=task.name,
            inputfile=inputfile,
            outputfile=outputfile,
        )

    def read_for_task(self, task, file, *, allow_dangerous=False):
        """Read parameters and require the parsed values to represent *task*."""
        out = self.create_for_task(task, inputfile=file)
        out.read_from_file(file, allow_dangerous=allow_dangerous)

        definition = out.definition
        if definition._task_definition is None:
            effective_name = definition.task_definition_name(out)
            if effective_name is not None:
                effective_name = self.normalize_name(effective_name)
            if effective_name != task.name:
                raise ValueError(
                    f"Input file read for task {task.name!r} represents "
                    f"task {effective_name!r}"
                )
        return out

    def create_from_result(self, task, result):
        """Create InputParameters for *task* initialized from *result*."""
        definition = self.definition_for_task(task)
        out = self.create_for_task(task)
        definition.update_from_result(out, result)
        return out

    @staticmethod
    def _materialize_file(spec, parameters, directory):
        if spec.condition is not None and not spec.condition(parameters):
            return None
        value = (
            spec.path(parameters, directory)
            if callable(spec.path)
            else spec.path.format(dataset=parameters.CONTROL.DATASET())
        )
        if value is None:
            return None
        try:
            path = os.fspath(value)
        except TypeError:
            raise TypeError(
                f"Task file resolver returned non-path value {value!r}"
            ) from None
        if not isinstance(path, str):
            raise TypeError(
                f"Task file resolver returned non-string path {type(path).__name__}"
            )
        return path if os.path.isabs(path) else os.path.join(directory, path)

    def files_for_task(self, task, parameters, directory, kind):
        """Materialize one class of files declared by the effective task."""
        self.definition_for_task(task)
        try:
            specs = getattr(task.files, kind)
        except AttributeError:
            raise ValueError(f"Unknown task file kind {kind!r}") from None
        return tuple(
            path
            for path in (
                self._materialize_file(spec, parameters, directory) for spec in specs
            )
            if path is not None
        )


class TaskRegistry:
    """Authoritative namespace and resolver for SPR-KKR task identities.

    Task declarations intentionally have two sources.  A schema with exactly
    one canonical task may co-locate its :class:`TaskDefinition` on the
    ``InputParametersDefinition`` as ``_task_definition``; this keeps the
    common one-schema/one-task case simple and declarative.  Tasks sharing one
    input schema are declared separately in the extra task package and refer
    to that schema by name.

    The registry owns the task namespace and resolution rules, not necessarily
    the physical location of each declaration.  Consumers should resolve task
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

    def effective_task(self, parameters):
        """Resolve the task identity represented by concrete parameter values."""
        definition = parameters._definition
        if definition._task_definition is not None:
            return definition._task_definition
        name = definition.task_definition_name(parameters)
        if name is None:
            return None
        return self[name]
