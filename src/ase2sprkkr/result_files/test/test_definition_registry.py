"""Registry contracts for distinct file formats sharing one extension."""

from types import SimpleNamespace

from ase2sprkkr.result_files.definition_registry import ResultFileDefinitionRegistry
from ase2sprkkr.result_files import ResultFile


def test_multiple_definitions_share_extension_without_overwrite():
    photo = SimpleNamespace(extension="spc")
    band = SimpleNamespace(extension="spc")

    class TestRegistry(ResultFileDefinitionRegistry):
        def names(self, *, include_aliases=True):
            return ("SPC", "BAND_SPC")

        def __getitem__(self, name):
            return {"SPC": photo, "BAND_SPC": band}[name.upper()]

        def values(self, *, include_aliases=False):
            return (photo, band)

    registry = TestRegistry()
    assert registry.by_extension() == {"spc": (photo, band)}
    assert registry.by_name() == {"spc": photo, "band_spc": band}


def test_existing_unique_extension_alias_is_preserved():
    definitions = ResultFile.definitions
    assert definitions["dat"] is definitions["jxc"]
    assert definitions["spc"] is not None
