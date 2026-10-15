"""Task-level metadata and registries for SPR-KKR calculations."""

from __future__ import annotations

from functools import cached_property
import os

from ..common.doc import process_input_parameters_definition
from ..common.module_registry import (
    ModuleRegistry,
    RegistryDefinitionError,
    RegistryEntryNotFound,
)


class TaskDefinition:
    """Definition of one SPR-KKR task."""

    def __init__(
        self,
        input_parameters_definition_name,
        *,
        name=None,
        executable="kkrscf",
        mpi=True,
        result_reader=None,
        input_file_templates=(),
        output_file_templates=(),
        cache_file_templates=(),
        prepare_from_result=None,
    ):
        if not isinstance(input_parameters_definition_name, str):
            raise TypeError("input_parameters_definition_name must be a string")
        self.input_parameters_definition_name = input_parameters_definition_name.upper()
        self.name = str(name or self.input_parameters_definition_name).upper()
        self.executable = executable
        self.mpi = mpi
        self.result_reader = result_reader
        self.input_file_templates = tuple(input_file_templates)
        self.output_file_templates = tuple(output_file_templates)
        self.cache_file_templates = tuple(cache_file_templates)
        self._prepare_from_result = prepare_from_result
        self._input_parameters_registry = None

    def _bind_input_parameters_registry(self, registry):
        """Bind the registry used to resolve the task's input schema."""
        if (
            self._input_parameters_registry is not None
            and self._input_parameters_registry is not registry
        ):
            raise RegistryDefinitionError(
                f"Task {self.name} is already bound to another input-parameter registry"
            )
        self._input_parameters_registry = registry
        return self

    @cached_property
    def input_parameters_definition(self):
        """Return the input-file schema referenced by this task."""
        if self._input_parameters_registry is None:
            raise RuntimeError(
                f"Task {self.name} is not bound to an input-parameter registry"
            )
        return self._input_parameters_registry[self.input_parameters_definition_name]

    def create_input_parameters(self, *, inputfile=None, outputfile=False):
        from .input_parameters import InputParameters

        definition = self.input_parameters_definition
        out = InputParameters(
            definition=definition,
            task_name=self.name,
            inputfile=inputfile,
            outputfile=outputfile,
        )
        definition.initialize_for_task(out, self.name)
        return out

    def read_input_parameters(self, file, *, allow_dangerous=False):
        out = self.create_input_parameters(inputfile=file)
        out.read_from_file(file, allow_dangerous=allow_dangerous)
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

    def input_files(self, parameters, directory):
        return self._resolve_file_templates(
            parameters, directory, self.input_file_templates
        )

    def output_files(self, parameters, directory):
        return self._resolve_file_templates(
            parameters, directory, self.output_file_templates
        )

    def cache_files(self, parameters, directory):
        return self._resolve_file_templates(
            parameters, directory, self.cache_file_templates
        )

    def execution_task_name(self, parameters=None):
        if parameters is not None:
            try:
                value = parameters["TASK"]["TASK"]()
            except KeyError:
                value = None
            if value:
                return str(value).upper()
        return self.name

    def output_reader_class(self, parameters=None):
        if self.result_reader is not None:
            return self.result_reader
        from ..outputs.reader_registry import reader_registry

        return reader_registry.class_for_task(self.execution_task_name(parameters))

    def prepare_from_result(self, parameters, result):
        if self._prepare_from_result is not None:
            self._prepare_from_result(parameters, result)


class InputParametersRegistry(ModuleRegistry):
    """Registry of input-file schemas, independent of task identity."""

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
            task._bind_input_parameters_registry(self)

        process_input_parameters_definition(module, value)
        return value

    def task_definitions(self):
        """Return canonical tasks declared by standard input schemas."""
        return {
            name: definition._task_definition
            for name in self.names()
            if (definition := self[name])._task_definition is not None
        }


class TaskRegistry:
    """Registry combining colocated 1:1 tasks and explicit shared-schema tasks."""

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
        task._bind_input_parameters_registry(self.input_parameters_registry)
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
        for name in self.input_parameters_registry.task_definitions():
            if name in names:
                raise RegistryDefinitionError(
                    f"Task {name} is defined both by an input schema and explicitly"
                )
            names.add(name)
        return tuple(sorted(names))

    def values(self):
        return tuple(self[name] for name in self.names())
