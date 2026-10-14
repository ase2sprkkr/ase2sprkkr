"""SPR-KKR task definitions and discovery."""

from __future__ import annotations

import copy
import os
from typing import Callable, TYPE_CHECKING

from ..common.doc import process_input_parameters_definition
from ..common.module_registry import ModuleRegistry
from . import definitions

if TYPE_CHECKING:
    from .input_parameters import InputParameters
    from .input_parameters_definitions import InputParametersDefinition
    from ..outputs.task_result import TaskResult


PrepareFromResult = Callable[["InputParameters", "TaskResult"], None]


class TaskDefinition:
    """Definition of one SPR-KKR task.

    The task owns execution semantics, auxiliary-file policy and the input-file
    schema. The schema itself describes only the input format and intentionally
    has no back-reference to its task.
    """

    def __init__(
        self,
        input_parameters_definition: "InputParametersDefinition",
        *,
        name: str,
        executable: str = "kkrscf",
        mpi: bool = True,
        result_reader: type | None = None,
        input_file_templates=(),
        output_file_templates=(),
        cache_file_templates=(),
        prepare_from_result: PrepareFromResult | None = None,
    ):
        self.name = str(name).upper()
        self.input_parameters_definition = input_parameters_definition
        self.executable = executable
        self.mpi = mpi
        self.result_reader = result_reader
        self.input_file_templates = tuple(input_file_templates)
        self.output_file_templates = tuple(output_file_templates)
        self.cache_file_templates = tuple(cache_file_templates)
        self._prepare_from_result = prepare_from_result

    def with_name(self, name: str) -> "TaskDefinition":
        """Return this task definition exposed under an alias name."""
        out = copy.copy(self)
        out.name = str(name).upper()
        return out

    def create_input_parameters(self, *, inputfile=None, outputfile=False) -> "InputParameters":
        """Create default input parameters bound to this task."""
        from .input_parameters import InputParameters

        return InputParameters(
            task_definition=self,
            inputfile=inputfile,
            outputfile=outputfile,
        )

    def read_input_parameters(self, file, *, allow_dangerous: bool = False) -> "InputParameters":
        """Parse an input file and bind the result to this task."""
        out = self.create_input_parameters(inputfile=file)
        out.read_from_file(file, allow_dangerous=allow_dangerous)
        return out

    def read_input_parameters_from_string(
        self,
        string: str,
        *,
        allow_dangerous: bool = False,
    ) -> "InputParameters":
        """Parse input-parameter text and bind the result to this task."""
        from io import StringIO

        return self.read_input_parameters(
            StringIO(string),
            allow_dangerous=allow_dangerous,
        )

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
        """Return task-specific produced files for one calculation."""
        return self._resolve_file_templates(
            parameters, directory, self.output_file_templates
        )

    def cache_files(self, parameters, directory):
        """Return optional reusable files for one calculation."""
        return self._resolve_file_templates(
            parameters, directory, self.cache_file_templates
        )

    def prepare_parameters_from_result(
        self,
        parameters: "InputParameters",
        result: "TaskResult",
    ) -> None:
        """Apply target-task-specific defaults when continuing from a result."""
        if self._prepare_from_result is not None:
            self._prepare_from_result(parameters, result)


class TaskRegistry(ModuleRegistry):
    """Discover :class:`TaskDefinition` objects from task definition modules."""

    aliases = {
        "BSFEK": "BSF",
        "BSFKK": "BSF",
    }

    def __init__(self):
        super().__init__(definitions, "task_definition", aliases=self.aliases)

    def prepare_value(self, name, module, value):
        if not isinstance(value, TaskDefinition):
            raise TypeError(
                f"{module.__name__}.task_definition() must return TaskDefinition, "
                f"got {type(value).__name__}."
            )
        if value.name != name:
            raise ValueError(
                f"Task module {module.__name__} is registered as {name}, "
                f"but defines task {value.name}."
            )
        process_input_parameters_definition(
            module,
            value.input_parameters_definition,
        )
        return value

    def alias_value(self, alias, value):
        return value.with_name(alias)

    def definition(self, name: str) -> TaskDefinition:
        return self[name]

    def contains(self, name: str) -> bool:
        return name in self

    @property
    def definitions(self) -> dict[str, TaskDefinition]:
        return {
            name: self[name]
            for name in self.names(include_aliases=False)
        }


TASKS = TaskRegistry()
