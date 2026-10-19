"""Discovery of structured SPR-KKR result-file definitions."""

from ..common.module_registry import ModuleRegistry
from . import definitions


class ResultFileDefinitionRegistry(ModuleRegistry):
    """Load result-file ``definition`` objects from definition modules."""

    def __init__(self, package=definitions):
        super().__init__(package, "definition")

    def by_extension(self):
        """Map each extension to *all* definitions that can parse it."""
        grouped = {}
        for definition in self.values(include_aliases=False):
            extension = definition.extension.lower().lstrip(".")
            grouped.setdefault(extension, []).append(definition)
        return {extension: tuple(group) for extension, group in grouped.items()}

    def by_name(self):
        """Map unique format names to definitions, retaining unambiguous ext aliases.

        Unlike file extensions, module names identify a single parser. For
        compatibility, a unique extension such as `dat` remains a valid
        format name (the JXC parser), even though its module is named `jxc`.
        """
        definitions = {
            name.lower(): self[name]
            for name in self.names(include_aliases=False)
        }
        for extension, candidates in self.by_extension().items():
            if len(candidates) == 1:
                existing = definitions.get(extension)
                if existing is not None and existing is not candidates[0]:
                    from ..common.module_registry import RegistryDefinitionError

                    raise RegistryDefinitionError(
                        f"Format name {extension!r} conflicts with its file extension"
                    )
                definitions.setdefault(extension, candidates[0])
        return definitions

    def definition_for_task_file(self, task_file):
        """Return the known structured format associated with a task file."""
        return None if task_file.type is None else self[task_file.type]

    def info_for_task_file(self, task_file):
        """Task-specific info wins over its registered file-format description."""
        if task_file.info is not None:
            return task_file.info
        definition = self.definition_for_task_file(task_file)
        return None if definition is None else definition.info



result_file_definition_registry = ResultFileDefinitionRegistry()
