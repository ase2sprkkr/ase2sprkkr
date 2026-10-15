"""SPR-KKR task definitions and current task-registry discovery."""

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
InputParametersDefinitionFactory = Callable[["TaskDefinition"], "InputParametersDefinition"]


class TaskDefinition:
    """Definition of one effective SPR-KKR task.

    A task owns execution/result semantics and references one input-file schema.
    Several distinct tasks may reference the same schema.  The schema resolves
    the effective task from concrete ``InputParameters`` values; registry lookup
    is supplied externally to those parameters rather than imported by the
    schema itself.
    """

    def __init__(
        self,
        input_parameters_definition: "InputParametersDefinition | None" = None,
        *,
        input_parameters_definition_factory: InputParametersDefinitionFactory | None = None,
        name: str,
        executable: str = "kkrscf",
        mpi: bool = True,
        result_reader: type | None = None,
        input_file_templates=(),
        output_file_templates=(),
        cache_file_templates=(),
        prepare_from_result: PrepareFromResult | None = None,
    ):
        if (input_parameters_definition is None) == (
            input_parameters_definition_factory is None
        ):
            raise ValueError(
                "Exactly one of input_parameters_definition and "
                "input_parameters_definition_factory must be supplied"
            )
        if (
            input_parameters_definition_factory is not None
            and not callable(input_parameters_definition_factory)
        ):
            raise TypeError("input_parameters_definition_factory must be callable")

        self.name = str(name).upper()
        self._input_parameters_definition = None
        self._input_parameters_definition_factory = input_parameters_definition_factory
        if input_parameters_definition is not None:
            self._set_input_parameters_definition(input_parameters_definition)
        self.executable = executable
        self.mpi = mpi
        self.result_reader = result_reader
        self.input_file_templates = tuple(input_file_templates)
        self.output_file_templates = tuple(output_file_templates)
        self.cache_file_templates = tuple(cache_file_templates)
        self._prepare_from_result = prepare_from_result

    def _set_input_parameters_definition(self, definition):
        from .input_parameters import InputParameters
        from .input_parameters_definitions import InputParametersDefinition

        if not isinstance(definition, InputParametersDefinition):
            raise TypeError(
                "Input-parameter definition factory must return "
                f"InputParametersDefinition, got {type(definition).__name__}"
            )
        if not issubclass(definition.result_class, InputParameters):
            raise TypeError(
                "InputParametersDefinition.result_class must be a subclass of "
                f"InputParameters, got {definition.result_class.__name__}"
            )

        definition._associate_task_definition(self)
        self._input_parameters_definition = definition
        return definition

    @property
    def input_parameters_definition(self) -> "InputParametersDefinition":
        """Return this task's input schema, constructing it on first use if needed."""
        definition = self._input_parameters_definition
        if definition is None:
            factory = self._input_parameters_definition_factory
            if factory is None:
                raise RuntimeError(f"Task {self.name} has no input-parameter definition")
            definition = self._set_input_parameters_definition(factory(self))
            self._input_parameters_definition_factory = None
        return definition

    def with_name(self, name: str) -> "TaskDefinition":
        """Return this task exposed under a registry alias.

        An alias is presentation/provenance, not another effective task of the
        schema.  Materialize the canonical schema first and share it without
        associating the alias back to the schema.
        """
        definition = self.input_parameters_definition
        out = copy.copy(self)
        out.name = str(name).upper()
        out._input_parameters_definition = definition
        out._input_parameters_definition_factory = None
        return out

    def create_input_parameters(
        self,
        *,
        inputfile=None,
        outputfile=False,
    ) -> "InputParameters":
        """Create parameters and remember that this task was explicitly requested."""
        return self.input_parameters_definition.create_object(
            inputfile=inputfile,
            outputfile=outputfile,
            requested_task_definition=self,
        )

    def read_input_parameters(
        self,
        file,
        *,
        allow_dangerous: bool = False,
    ) -> "InputParameters":
        """Parse an input file under an explicitly requested task context."""
        return self.input_parameters_definition.read_from_file(
            file,
            allow_dangerous=allow_dangerous,
            inputfile=file,
            requested_task_definition=self,
        )

    def read_input_parameters_from_string(
        self,
        string: str,
        *,
        allow_dangerous: bool = False,
    ) -> "InputParameters":
        """Parse input-parameter text under an explicitly requested task context."""
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
        return self._resolve_file_templates(parameters, directory, self.input_file_templates)

    def output_files(self, parameters, directory):
        return self._resolve_file_templates(parameters, directory, self.output_file_templates)

    def cache_files(self, parameters, directory):
        return self._resolve_file_templates(parameters, directory, self.cache_file_templates)

    def prepare_parameters_from_result(
        self,
        parameters: "InputParameters",
        result: "TaskResult",
    ) -> None:
        if self._prepare_from_result is not None:
            self._prepare_from_result(parameters, result)


class TaskRegistry(ModuleRegistry):
    """Current compatibility registry for task-definition modules.

    The independent registry refactor replaces this with the composite registry
    that combines canonical tasks exposed by the input-definition registry with
    additional shared-schema task modules.
    """

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
        process_input_parameters_definition(module, value.input_parameters_definition)
        return value

    def alias_value(self, alias, value):
        return value.with_name(alias)

    def definition(self, name: str) -> TaskDefinition:
        return self[name]

    def contains(self, name: str) -> bool:
        return name in self

    @property
    def definitions(self) -> dict[str, TaskDefinition]:
        return {name: self[name] for name in self.names(include_aliases=False)}


TASKS = TaskRegistry()

# Compatibility composition-root wiring.  InputParameters does not import this
# registry; the registry side supplies itself.  The independent registry work
# will replace this assignment with the composite TaskRegistry and will also
# provide InputParameters.definition_registry.
from .input_parameters import InputParameters  # NOQA: E402
InputParameters.task_registry = TASKS
