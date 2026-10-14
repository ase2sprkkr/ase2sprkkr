from ..output_files import OutputFile
from ..output_files_definitions import OutputFileDefinition


def test_output_file_registry_discovers_definitions():
    registry = OutputFile.output_file_registry

    assert "dos" in registry
    assert "DOS" in registry
    assert isinstance(registry["dos"], OutputFileDefinition)
    assert OutputFile.definition("DOS") is registry["dos"]
    assert OutputFile.definitions is registry


def test_output_file_registry_uses_module_names_as_keys():
    registry = OutputFile.output_file_registry

    assert {"bsf", "dos", "jxc", "rat", "sfn", "spc"}.issubset(registry.keys())
    assert all(definition.extension == name for name, definition in registry.items())
