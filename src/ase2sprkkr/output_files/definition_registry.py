"""Discovery of structured SPR-KKR output-file definitions."""

from ..common.module_registry import ModuleRegistry
from . import definitions


class OutputFileDefinitionRegistry(ModuleRegistry):
    """Load output-file ``definition`` objects from definition modules."""

    def __init__(self, package=definitions):
        super().__init__(package, "definition")

    def by_extension(self):
        """Return loaded definitions keyed by their declared file extension."""
        return {
            definition.extension: definition
            for definition in self.values(include_aliases=False)
        }


output_file_definition_registry = OutputFileDefinitionRegistry()
