"""Discovery and selection of task-specific SPR-KKR output readers."""

import os
import re

from ..common.module_registry import (
    ModuleRegistry,
    RegistryDefinitionError,
    RegistryEntryNotFound,
)
from ..sprkkr.executables import executables as default_executables
from . import readers


_TASK_PATTERN = re.compile(rb"(?m)^\s*TASK\s*=\s*([A-Z]+)\s*$")


class OutputReaderRegistry(ModuleRegistry):
    """Load, select and detect ``<Task>OutputReader`` classes."""

    def __init__(self, package=readers, executables=default_executables):
        super().__init__(package)
        self._executables = executables

    def value_from_module(self, name, module):
        class_name = name.title() + "OutputReader"
        try:
            return vars(module)[class_name]
        except KeyError as exc:
            raise RegistryDefinitionError(
                f"Reader module {module.__name__!r} does not define expected "
                f"class {class_name!r}"
            ) from exc

    def prepare_value(self, name, module, value):
        from .task_result import KkrOutputReader

        if not isinstance(value, type) or not issubclass(value, KkrOutputReader):
            raise RegistryDefinitionError(
                f"{module.__name__}.{name.title()}OutputReader must be a "
                f"KkrOutputReader subclass, got {value!r}"
            )
        return value

    def _class_for_reader_name(self, name):
        """Resolve a reader identifier, with the generic reader as fallback."""
        from .readers.default import DefaultOutputReader

        if not name:
            return DefaultOutputReader
        try:
            return self[name]
        except RegistryEntryNotFound:
            return DefaultOutputReader

    def class_for_task_name(self, task_name):
        """Return the reader for a registered task name or guess by reader name."""
        if not task_name:
            return self._class_for_reader_name(task_name)

        # Registered tasks use declarative task/executable metadata. Unknown
        # task names retain the generic standalone/legacy reader fallback.
        from ..input_parameters.input_parameters import InputParameters

        try:
            task = InputParameters.task_registry[task_name]
        except RegistryEntryNotFound:
            return self._class_for_reader_name(task_name)
        return self.class_for_task(task)

    def class_for_task(self, task):
        """Return the reader explicitly declared by the task or its executable."""
        from ..input_parameters.task_definition import TaskDefinition

        if not isinstance(task, TaskDefinition):
            raise TypeError("task must be a TaskDefinition")
        if task.output_reader_name is not None:
            reader_name = task.output_reader_name
        else:
            executable = self._executables[task.executable]
            reader_name = executable.default_output_reader
        return self[reader_name]

    @staticmethod
    def _validate_input_parameters(task, input_parameters):
        if input_parameters is None:
            return
        parameter_task = input_parameters.task_definition
        if parameter_task is not None and parameter_task != task:
            raise ValueError(
                f"Input parameters represent task {parameter_task.name!r}, "
                f"not requested task {task.name!r}"
            )

    @staticmethod
    def _validate_directory(directory):
        if directory is None:
            return
        try:
            os.fspath(directory)
        except TypeError as exc:
            raise TypeError("directory must be a path-like object or None") from exc

    def reader_for_task(
        self,
        task,
        *,
        input_parameters=None,
        calculator=None,
        print_output=False,
        directory=None,
    ):
        """Construct the registered output reader for *task*."""
        self._validate_input_parameters(task, input_parameters)
        if directory is None and calculator is not None:
            directory = getattr(calculator, "directory", None)
        self._validate_directory(directory)
        reader = self.class_for_task(task)(
            input_parameters,
            calculator,
            directory,
            print_output,
        )
        if reader.task_name is None:
            reader.task_name = task.name
        return reader

    def read_for_task(
        self,
        task,
        filename,
        *,
        input_parameters=None,
        calculator=None,
        directory=None,
    ):
        """Read a completed output file using the registered reader for *task*."""
        if directory is None and isinstance(filename, (str, os.PathLike)):
            directory = os.path.dirname(os.fspath(filename))
        return self.reader_for_task(
            task,
            input_parameters=input_parameters,
            calculator=calculator,
            directory=directory,
        ).read_from_file(filename)

    def detect(self, file):
        """Detect the task represented by a standalone SPR-KKR output file.

        Modern outputs identify their task explicitly. Older outputs that do
        not do so are offered to task-specific reader fingerprints, keeping
        legacy format knowledge beside the reader that understands it.
        """
        with open(os.fspath(file), "rb") as output:
            raw_output = output.read()

        match = _TASK_PATTERN.search(raw_output)
        if match is not None:
            task = match.group(1).decode("ascii")
            if task != "NONE":
                return task

        detected = [
            task
            for task in self.names(include_aliases=False)
            if self[task].matches_output(raw_output)
        ]
        if len(detected) > 1:
            raise RegistryDefinitionError(
                "Output matches multiple task-specific readers: "
                + ", ".join(detected)
            )
        return detected[0] if detected else None


reader_registry = OutputReaderRegistry()
