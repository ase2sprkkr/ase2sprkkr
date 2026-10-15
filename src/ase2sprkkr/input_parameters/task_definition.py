"""Task-level metadata and discovery for SPR-KKR calculations."""

from __future__ import annotations

import os

from ..common.doc import process_input_parameters_definition
from ..common.module_registry import ModuleRegistry, RegistryEntryNotFound


class TaskDefinition:
    """Definition of one canonical SPR-KKR task.

    Task identity and execution/result semantics live here.  The input-file
    schema is a separate object and is materialized lazily; once materialized it
    points back to this task so runtime InputParameters do not need a second
    authoritative task reference.
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
            None if callable(input_parameters_definition) else input_parameters_definition
        )
        if name is None:
            if self._input_parameters_definition is None:
                raise TypeError("name is required for a lazy input-parameter schema")
            name = self._input_parameters_definition.name
        self.name = str(name).upper()
        self.executable = executable
        self.mpi = mpi
        self.result_reader = result_reader
        self.input_file_templates = tuple(input_file_templates)
        self.output_file_templates = tuple(output_file_templates)
        self.cache_file_templates = tuple(cache_file_templates)
        self._prepare_from_result = prepare_from_result
        self._definition_module = None
        self._definition_processed = False
        if self._input_parameters_definition is not None:
            self._bind_input_parameters_definition(self._input_parameters_definition)

    @property
    def input_parameters_definition(self):
        """Return this task's lazily materialized input-file schema."""
        if self._input_parameters_definition is None:
            definition = self._input_parameters_definition_source()
            self._bind_input_parameters_definition(definition)
            self._input_parameters_definition = definition
        self._process_input_parameters_definition()
        return self._input_parameters_definition

    def _bind_definition_module(self, module):
        self._definition_module = module
        self._process_input_parameters_definition()

    def _bind_input_parameters_definition(self, definition):
        bound = getattr(definition, "task_definition", None)
        if bound is not None and bound is not self:
            raise ValueError(
                f"Input-parameter definition {definition.name} is already bound "
                f"to task {bound.name}"
            )
        definition.task_definition = self

    def _process_input_parameters_definition(self):
        if (
            self._definition_processed
            or self._definition_module is None
            or self._input_parameters_definition is None
        ):
            return
        process_input_parameters_definition(
            self._definition_module, self._input_parameters_definition
        )
        self._definition_processed = True

    def create_input_parameters(
        self, *, task_name=None, inputfile=None, outputfile=False
    ):
        """Create parameters for this task.

        ``task_name`` is only the concrete lookup/view name (for example
        ``BSFEK``).  Canonical task identity remains this TaskDefinition.
        """
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
        """Return the low-level SPR-KKR TASK.TASK selected by parameters."""
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
        """Apply target-task-specific continuation setup."""
        if self._prepare_from_result is not None:
            self._prepare_from_result(parameters, result)


class TaskRegistry(ModuleRegistry):
    """Discover canonical TaskDefinition objects from task modules."""

    aliases = {"BSFEK": "BSF", "BSFKK": "BSF"}

    def __init__(self, package):
        super().__init__(package, "task_definition", aliases=self.aliases)

    def prepare_value(self, name, module, value):
        if not isinstance(value, TaskDefinition):
            raise TypeError(
                f"{module.__name__}.task_definition must return TaskDefinition, "
                f"got {type(value).__name__}"
            )
        if value.name != name:
            raise ValueError(
                f"Task definition in {module.__name__} is named {value.name}, "
                f"but the module defines task {name}"
            )
        value._bind_definition_module(module)
        return value

    def alias_value(self, alias, value):
        """Aliases select the same canonical task definition."""
        return value
