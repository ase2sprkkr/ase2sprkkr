"""Task-level metadata and registries for SPR-KKR calculations."""

from __future__ import annotations

import os

from ..common.doc import process_input_parameters_definition
from ..common.module_registry import ModuleRegistry, RegistryEntryNotFound


class TaskDefinition:
    """Definition of one SPR-KKR task.

    The input-parameter schema can be supplied directly, lazily via a factory,
    or by registry name. Registry-name references are used by special tasks that
    share one input schema.
    """

    def __init__(
        self,
        input_parameters_definition,
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
        self._input_parameters_definition_source = input_parameters_definition
        self._input_parameters_definition = (
            None
            if callable(input_parameters_definition)
            or isinstance(input_parameters_definition, str)
            else input_parameters_definition
        )
        if name is None:
            if self._input_parameters_definition is None:
                raise TypeError("name is required for a lazy or registry input schema")
            name = self._input_parameters_definition.name
        self.name = str(name).upper()
        self.executable = executable
        self.mpi = mpi
        self.result_reader = result_reader
        self.input_file_templates = tuple(input_file_templates)
        self.output_file_templates = tuple(output_file_templates)
        self.cache_file_templates = tuple(cache_file_templates)
        self._prepare_from_result = prepare_from_result

    @property
    def input_parameters_definition(self):
        """Return the lazily resolved input-file schema used by this task."""
        if self._input_parameters_definition is None:
            source = self._input_parameters_definition_source
            if isinstance(source, str):
                from .input_parameters import InputParameters

                definition = InputParameters.input_parameters_registry[source]
            else:
                definition = source()
            self._input_parameters_definition = definition
        return self._input_parameters_definition

    def create_input_parameters(
        self, *, task_name=None, inputfile=None, outputfile=False
    ):
        from .input_parameters import InputParameters

        return InputParameters(
            definition=self.input_parameters_definition,
            task_name=task_name or self.name,
            inputfile=inputfile,
            outputfile=outputfile,
        )

    def read_input_parameters(
        self, file, *, task_name=None, allow_dangerous=False
    ):
        out = self.create_input_parameters(task_name=task_name, inputfile=file)
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
        return self._resolve_file_templates(parameters, directory, self.input_file_templates)

    def output_files(self, parameters, directory):
        return self._resolve_file_templates(parameters, directory, self.output_file_templates)

    def cache_files(self, parameters, directory):
        return self._resolve_file_templates(parameters, directory, self.cache_file_templates)

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
        from ..outputs.readers.default import DefaultOutputReader

        try:
            return reader_registry[self.execution_task_name(parameters)]
        except RegistryEntryNotFound:
            return DefaultOutputReader

    def prepare_from_result(self, parameters, result):
        if self._prepare_from_result is not None:
            self._prepare_from_result(parameters, result)


class InputParametersRegistry(ModuleRegistry):
    """Registry of input-file schemas, independent of task identity."""

    def __init__(self, package):
        super().__init__(package, "input_parameters_definition")

    def prepare_value(self, name, module, value):
        from .input_parameters_definitions import InputParametersDefinition

        if not isinstance(value, InputParametersDefinition):
            raise TypeError(
                f"{module.__name__}.input_parameters_definition must return "
                f"InputParametersDefinition, got {type(value).__name__}"
            )
        if str(value.name).upper() != name:
            raise ValueError(
                f"Input-parameter definition in {module.__name__} is named "
                f"{value.name}, but the module defines {name}"
            )
        task = value._task_definition
        if task is not None:
            if not isinstance(task, TaskDefinition):
                raise TypeError(
                    f"{module.__name__}._task_definition must be TaskDefinition, "
                    f"got {type(task).__name__}"
                )
            if task.name != name:
                raise ValueError(
                    f"Canonical task {task.name} does not match input schema {name}"
                )
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
    """Composite registry of canonical and extra task definitions.

    Standard tasks are exposed by their InputParametersDefinition through
    ``_task_definition``. Only tasks without a 1:1 schema relationship live in
    the extra module registry.
    """

    normalize_name = staticmethod(ModuleRegistry.normalize_name)

    def __init__(self, input_parameters_registry, extra_tasks_package):
        self.input_parameters_registry = input_parameters_registry
        self.extra_tasks = ModuleRegistry(extra_tasks_package, "task_definition")

    def _input_task(self, name):
        try:
            definition = self.input_parameters_registry[name]
        except RegistryEntryNotFound:
            return None
        task = definition._task_definition
        if task is None:
            return None
        if task.name != name:
            raise ValueError(
                f"Canonical task {task.name} does not match input schema {name}"
            )
        return task

    def _extra_task(self, name):
        try:
            task = self.extra_tasks[name]
        except RegistryEntryNotFound:
            return None
        if not isinstance(task, TaskDefinition):
            raise TypeError(
                f"Extra task {name} must be TaskDefinition, got {type(task).__name__}"
            )
        if task.name != name:
            raise ValueError(f"Extra task {task.name} is registered as {name}")
        return task

    def __getitem__(self, name):
        name = self.normalize_name(name)
        input_task = self._input_task(name)
        extra_task = self._extra_task(name)
        if input_task is not None and extra_task is not None:
            raise ValueError(f"Task {name} is defined both by an input schema and explicitly")
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
        for name, task in self.input_parameters_registry.task_definitions().items():
            if name in names:
                raise ValueError(
                    f"Task {name} is defined both by an input schema and explicitly"
                )
            if task is not None:
                names.add(name)
        return tuple(sorted(names))

    def values(self):
        return tuple(self[name] for name in self.names())
