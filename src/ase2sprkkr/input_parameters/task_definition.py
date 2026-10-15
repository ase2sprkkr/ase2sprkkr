"""Task-level metadata and discovery for SPR-KKR calculations."""

from __future__ import annotations

import copy
import os

from ..common.doc import process_input_parameters_definition
from ..common.module_registry import ModuleRegistry, RegistryDefinitionError


class TaskDefinition:
    """Definition of one SPR-KKR task.

    The task owns the input-file schema together with execution/result metadata.
    The schema itself remains independent of the task that uses it.
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
        self.input_parameters_definition = input_parameters_definition
        self.name = str(name or input_parameters_definition.name).upper()
        self.executable = executable
        self.mpi = mpi
        self.result_reader = result_reader
        self.input_file_templates = tuple(input_file_templates)
        self.output_file_templates = tuple(output_file_templates)
        self.cache_file_templates = tuple(cache_file_templates)
        self._prepare_from_result = prepare_from_result

    def with_name(self, name):
        """Return the same task definition exposed under another task name."""
        out = copy.copy(self)
        out.name = str(name).upper()
        return out

    def create_input_parameters(self, *, inputfile=None, outputfile=False):
        """Create task-bound input parameters with default values."""
        from .input_parameters import InputParameters

        return InputParameters(
            task_definition=self,
            inputfile=inputfile,
            outputfile=outputfile,
        )

    def read_input_parameters(self, file, *, allow_dangerous=False):
        """Parse an input file into parameters bound to this task."""
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
        """Return task-specific required input files for one calculation."""
        return self._resolve_file_templates(
            parameters, directory, self.input_file_templates
        )

    def output_files(self, parameters, directory):
        """Return task-specific output files for one calculation."""
        return self._resolve_file_templates(
            parameters, directory, self.output_file_templates
        )

    def cache_files(self, parameters, directory):
        """Return optional reusable cache files for one calculation."""
        return self._resolve_file_templates(
            parameters, directory, self.cache_file_templates
        )

    def output_reader_class(self):
        """Return the reader class used for this task's textual output."""
        if self.result_reader is not None:
            return self.result_reader

        from ..outputs.task_result import KkrOutputReader

        return KkrOutputReader.class_for_task(
            self.input_parameters_definition.name
        )

    def prepare_from_result(self, parameters, result):
        """Prepare newly created task parameters from a previous result."""
        if self._prepare_from_result is not None:
            self._prepare_from_result(parameters, result)


class TaskRegistry(ModuleRegistry):
    """Discover :class:`TaskDefinition` objects from task definition modules."""

    aliases = {
        "BSFEK": "BSF",
        "BSFKK": "BSF",
    }

    def __init__(self, package):
        super().__init__(
            package,
            "task_definition",
            aliases=self.aliases,
            include=lambda name: not name.startswith("_") and name != "sections",
        )

    def value_from_module(self, name, module):
        value = super().value_from_module(name, module)
        return value() if callable(value) else value

    def prepare_value(self, name, module, value):
        if not isinstance(value, TaskDefinition):
            raise RegistryDefinitionError(
                f"{module.__name__}.task_definition must return TaskDefinition, "
                f"got {type(value).__name__}"
            )
        if value.name != name:
            raise RegistryDefinitionError(
                f"Task definition in {module.__name__} is named {value.name}, "
                f"but the module defines task {name}"
            )
        process_input_parameters_definition(
            module, value.input_parameters_definition
        )
        return value

    def alias_value(self, alias, value):
        return value.with_name(alias)
