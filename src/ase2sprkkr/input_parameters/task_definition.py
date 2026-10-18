"""Task-level metadata and registries for SPR-KKR calculations."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
import os
from typing import Iterable, Optional, Tuple, Union

from ..common.doc import process_input_parameters_definition
from ..common.module_registry import (
    ModuleRegistry,
    RegistryDefinitionError,
    RegistryEntryNotFound,
)


_FileTemplate = Union[str, os.PathLike[str]]


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


def _normalize_file_templates(
    templates: Iterable[_FileTemplate], field: str
) -> Tuple[str, ...]:
    if isinstance(templates, (str, bytes, os.PathLike)):
        raise TypeError(
            f"{field} must be an iterable of path-like strings, not a single path"
        )
    try:
        iterator = iter(templates)
    except TypeError:
        raise TypeError(
            f"{field} must be an iterable of path-like strings"
        ) from None

    normalized = []
    for index, item in enumerate(iterator):
        try:
            path = os.fspath(item)
        except TypeError:
            raise TypeError(
                f"{field}[{index}] must be a string or os.PathLike, "
                f"got {type(item).__name__}"
            ) from None
        if not isinstance(path, str):
            raise TypeError(
                f"{field}[{index}] must resolve to a string, "
                f"got {type(path).__name__}"
            )
        normalized.append(path)
    return tuple(normalized)


@dataclass(frozen=True, init=False)
class TaskDefinition:
    """Immutable declarative metadata identifying one SPR-KKR task.

    Convenience methods delegate task-bound operations to their owning service
    registries; resolution and parsing policy does not live in this value object.
    """

    input_parameters_definition_name: str
    name: str
    aliases: Tuple[str, ...]
    executable: str
    mpi: bool
    output_reader_name: Optional[str]
    input_file_templates: Tuple[str, ...]
    output_file_templates: Tuple[str, ...]
    cache_file_templates: Tuple[str, ...]

    def __init__(
        self,
        input_parameters_definition_name: str,
        *,
        name: Optional[str] = None,
        aliases: Iterable[str] = (),
        executable: str = "kkrscf",
        mpi: bool = True,
        output_reader_name: Optional[str] = None,
        output_reader: Optional[str] = None,
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
        if isinstance(aliases, str):
            aliases = (aliases,)
        aliases = tuple(
            _normalize_identifier(alias, "aliases") for alias in aliases
        )
        if task_name in aliases or len(set(aliases)) != len(aliases):
            raise ValueError("Task aliases must be distinct from each other and the canonical name")
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
        object.__setattr__(self, "aliases", aliases)
        object.__setattr__(self, "executable", executable)
        object.__setattr__(self, "mpi", mpi)
        object.__setattr__(self, "output_reader_name", output_reader_name)
        object.__setattr__(
            self,
            "input_file_templates",
            _normalize_file_templates(input_file_templates, "input_file_templates"),
        )
        object.__setattr__(
            self,
            "output_file_templates",
            _normalize_file_templates(output_file_templates, "output_file_templates"),
        )
        object.__setattr__(
            self,
            "cache_file_templates",
            _normalize_file_templates(cache_file_templates, "cache_file_templates"),
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
        canonical = definition._task_definition
        if canonical is not None and canonical != task:
            raise RegistryDefinitionError(
                f"Task {task.name} references input schema "
                f"{task.input_parameters_definition_name}, which is bound to "
                f"different canonical task metadata {canonical!r}"
            )
        return definition

    def create_for_task(
        self,
        task,
        *,
        inputfile=None,
        outputfile=False,
        requested_task_name=None,
    ):
        """Create values for a task, preserving the invoked alias for defaults."""
        from .input_parameters import InputParameters

        definition = self.definition_for_task(task)
        requested_name = (
            task.name if requested_task_name is None
            else _normalize_identifier(requested_task_name, "requested_task_name")
        )
        if requested_name != task.name and requested_name not in task.aliases:
            raise ValueError(f"{requested_name} is not a name of task {task.name}")
        return InputParameters(
            definition=definition,
            requested_task_name=requested_name,
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

        task_name = out.task_name
        if out.task_definition != task:
            raise TaskMismatchError(task.name, task_name)

        out._requested_task_name = task_name
        return out

    def create_from_result(self, task, result):
        """Create InputParameters for *task* using snapshot defaults from *result*."""
        definition = self.definition_for_task(task)
        out = self.create_for_task(task)
        definition.initialize_from_result(out, result)
        return out

    @staticmethod
    def _resolve_file_templates(parameters, directory, templates):
        dataset = parameters.CONTROL.DATASET()
        return tuple(
            path if os.path.isabs(path) else os.path.join(directory, path)
            for path in (
                os.fspath(template).format(dataset=dataset)
                for template in templates
            )
        )

    def files_for_task(self, task, parameters, directory, kind):
        """Combine declarative task templates with schema-specific file rules."""
        definition = self.definition_for_task(task)
        templates = getattr(task, f"{kind}_file_templates")
        declared = self._resolve_file_templates(parameters, directory, templates)
        schema_files = getattr(definition, f"{kind}_files")(parameters, directory)
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

    @cached_property
    def _tasks_by_name(self):
        """All canonical task names and aliases; reject shadowing and collisions.

        The input and extra task registries are singleton-style declarations.
        Task identity and alias metadata must not change after first resolution.
        """
        canonical = {}
        for schema_name in self.input_parameters_registry.names():
            task = self.input_parameters_registry[schema_name]._task_definition
            if task is None:
                continue
            if task.name in canonical:
                raise RegistryDefinitionError(f"Duplicate task {task.name}")
            canonical[task.name] = task
        for name in self.extra_tasks.names(include_aliases=False):
            task = self._extra_task(name)
            if task.name in canonical:
                raise RegistryDefinitionError(
                    f"Task {task.name} is defined both by an input schema and explicitly"
                )
            canonical[task.name] = task

        by_name = dict(canonical)
        for task in canonical.values():
            for alias in task.aliases:
                if alias in by_name:
                    raise RegistryDefinitionError(
                        f"Task alias {alias!r} conflicts with another task name or alias"
                    )
                by_name[alias] = task
        return by_name

    def __getitem__(self, name):
        name = self.normalize_name(name)
        try:
            return self._tasks_by_name[name]
        except KeyError:
            raise RegistryEntryNotFound(name) from None

    def __contains__(self, name):
        return self.normalize_name(name) in self._tasks_by_name

    def names(self, *, include_aliases=True):
        names = self._tasks_by_name
        if not include_aliases:
            names = {task.name: task for task in names.values()}
        return tuple(sorted(names))

    def values(self):
        """Return each canonical task once, excluding duplicate alias entries."""
        return tuple(self[name] for name in self.names(include_aliases=False))

    def task_for_parameters(self, parameters):
        """Resolve the task represented by concrete parameter values."""
        name = parameters.task_name
        return None if name is None else self[name]
