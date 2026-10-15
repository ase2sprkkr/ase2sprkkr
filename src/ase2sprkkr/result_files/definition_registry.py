"""Discovery of structured SPR-KKR result-file definitions."""

from ..common.module_registry import ModuleRegistry
from . import definitions


class ResultFileDefinitionRegistry(ModuleRegistry):
    """Load result-file ``definition`` objects from definition modules."""

    def __init__(self, package=definitions):
        super().__init__(package, "definition")

    def by_extension(self):
        """Return loaded definitions keyed by their declared file extension."""
        return {
            definition.extension: definition
            for definition in self.values(include_aliases=False)
        }


result_file_definition_registry = ResultFileDefinitionRegistry()
