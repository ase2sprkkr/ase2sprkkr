"""Discovery and selection of task-specific SPR-KKR output readers."""

import os
import re

from ..common.module_registry import (
    ModuleRegistry,
    RegistryDefinitionError,
    RegistryEntryNotFound,
)
from . import readers


_TASK_PATTERN = re.compile(rb"(?m)^\s*TASK\s*=\s*([A-Z]+)\s*$")


class OutputReaderRegistry(ModuleRegistry):
    """Load, select and detect ``<Task>OutputReader`` classes."""

    def __init__(self, package=readers):
        super().__init__(package)

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

    def class_for_task(self, task):
        """Return the specialized reader for *task*, or the default reader."""
        from .readers.default import DefaultOutputReader

        if not task:
            return DefaultOutputReader
        try:
            return self[task]
        except RegistryEntryNotFound:
            return DefaultOutputReader

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
