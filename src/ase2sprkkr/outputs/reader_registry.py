"""Discovery of task-specific SPR-KKR output readers."""

from ..common.module_registry import ModuleRegistry
from . import readers


class OutputReaderRegistry(ModuleRegistry):
    """Load ``<Task>OutputReader`` classes from ``outputs.readers.<task>``."""

    def __init__(self, package=readers):
        super().__init__(package)

    def value_from_module(self, name, module):
        class_name = name.title() + "OutputReader"
        try:
            return getattr(module, class_name)
        except AttributeError as exc:
            raise KeyError(name) from exc


reader_registry = OutputReaderRegistry()
