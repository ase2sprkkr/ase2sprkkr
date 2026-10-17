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

    def definition_for_task_file(self, task_file):
        """Return the structured result definition declared by *task_file*."""
        return None if task_file.type is None else self[task_file.type]

    def info_for_task_file(self, task_file):
        """Return task-specific info or the referenced result-format description."""
        if task_file.info is not None:
            return task_file.info
        definition = self.definition_for_task_file(task_file)
        return None if definition is None else definition.info


result_file_definition_registry = ResultFileDefinitionRegistry()
