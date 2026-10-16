"""Task-level metadata and registries for SPR-KKR calculations."""

from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Optional, Tuple

from ..common.doc import process_input_parameters_definition
from ..common.module_registry import (
    ModuleRegistry,
    RegistryDefinitionError,
    RegistryEntryNotFound,
)


@dataclass(frozen=True, init=False)
class TaskDefinition:
    """Immutable declarative metadata identifying one SPR-KKR task."""

    input_parameters_definition_name: str
    name: str
    executable: str
    mpi: bool
    output_reader: Optional[str]
    input_file_templates: Tuple[object, ...]
    output_file_templates: Tuple[object, ...]
    cache_file_templates: Tuple[object, ...]

    def __init__(
        self,
        input_parameters_definition_name,
        *,
        name=None,
        executable="kkrscf",
        mpi=True,
        output_reader=None,
        input_file_templates=(),
        output_file_templates=(),
        cache_file_templates=(),
    ):
        if not isinstance(input_parameters_definition_name, str):
            raise TypeError("input_parameters_definition_name must be a string")
        if output_reader is not None and not isinstance(output_reader, str):
            raise TypeError("output_reader must be a string or None")

        object.__setattr__(
            self,
            "input_parameters_definition_name",
            input_parameters_definition_name.upper(),
        )
        object.__setattr__(
            self,
            "name",
            str(name or input_parameters_definition_name).upper(),
        )
        object.__setattr__(self, "executable", executable)
        object.__setattr__(self, "mpi", mpi)
        object.__setattr__(
            self,
            "output_reader",
            output_reader.upper() if output_reader is not None else None,
        )
        object.__setattr__(self, "input_file_templates", tuple(input_file_templates))
        object.__setattr__(self, "output_file_templates", tuple(output_file_templates))
        object.__setattr__(self, "cache_file_templates", tuple(cache_file_templates))


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
        out = InputParameters(
            definition=definition,
            task_name=task.name,
            inputfile=inputfile,
            outputfile=outputfile,
        )
        definition.initialize_for_task(out, task.name)
        return out

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
