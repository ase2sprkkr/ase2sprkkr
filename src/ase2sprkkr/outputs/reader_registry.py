"""Discovery of task-specific SPR-KKR output readers."""

from ..common.module_registry import ModuleRegistry, RegistryDefinitionError
from . import readers


class OutputReaderRegistry(ModuleRegistry):
    """Load ``<Task>OutputReader`` classes from ``outputs.readers.<task>``."""

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


reader_registry = OutputReaderRegistry()
