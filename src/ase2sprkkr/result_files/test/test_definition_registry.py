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

def test_format_autodetection_selects_a_unique_candidate(monkeypatch, tmp_path):
    from ase2sprkkr.result_files import definition_registry
    import pyparsing as pp

    filename = tmp_path / "example.spc"
    filename.write_text("BAND")

    class FakeDefinition:
        def __init__(self, name, valid):
            self.name = name
            self.extension = "spc"
            self.valid = valid

        def read_from_file(self, file, potential=None):
            if not self.valid:
                raise pp.ParseException("wrong file format")
            return self.name

    photo = FakeDefinition("SPC", False)
    band = FakeDefinition("BAND_SPC", True)
    monkeypatch.setattr(definition_registry.result_file_definition_registry, "by_extension", lambda: {"spc": (photo, band)})
    monkeypatch.setattr(definition_registry.result_file_definition_registry, "by_name", lambda: {"spc": photo, "band_spc": band})
    monkeypatch.setattr(definition_registry.result_file_definition_registry, "values", lambda include_aliases=False: (photo, band))

    # Cached definitions may already have been populated by other tests.
    monkeypatch.setattr(ResultFile, "definitions", {"spc": photo, "band_spc": band})
    assert ResultFile.from_file(filename, unknown=False) == "BAND_SPC"
    assert ResultFile.from_file(filename, try_only="band_spc") == "BAND_SPC"


def test_ambiguous_format_is_not_silently_selected(monkeypatch, tmp_path):
    import pytest
    from ase2sprkkr.result_files import definition_registry

    filename = tmp_path / "example.spc"
    filename.write_text("data")

    class FakeDefinition:
        extension = "spc"
        def __init__(self, name):
            self.name = name
        def read_from_file(self, file, potential=None):
            return self.name

    first, second = FakeDefinition("SPC"), FakeDefinition("BAND_SPC")
    monkeypatch.setattr(definition_registry.result_file_definition_registry, "by_extension", lambda: {"spc": (first, second)})
    monkeypatch.setattr(definition_registry.result_file_definition_registry, "by_name", lambda: {"spc": first, "band_spc": second})
    monkeypatch.setattr(ResultFile, "definitions", {"spc": first, "band_spc": second})
    with pytest.raises(ValueError, match="Ambiguous result file"):
        ResultFile.from_file(filename)

def test_first_try_extension_checks_all_shared_extension_formats(monkeypatch, tmp_path):
    from ase2sprkkr.result_files import definition_registry

    filename = tmp_path / "sample.spc"
    filename.write_text("data")

    class Candidate:
        extension = "spc"

        def __init__(self, name, success):
            self.name, self.success = name, success

        def read_from_file(self, file, potential=None):
            if not self.success:
                import pyparsing as pp
                raise pp.ParseException("not this format")
            return self.name

    photo, band = Candidate("SPC", False), Candidate("BAND_SPC", True)
    registry = definition_registry.result_file_definition_registry
    monkeypatch.setattr(registry, "by_extension", lambda: {"spc": (photo, band)})
    monkeypatch.setattr(registry, "by_name", lambda: {"spc": photo, "band_spc": band})
    monkeypatch.setattr(registry, "values", lambda include_aliases=False: (photo, band))
    monkeypatch.setattr(ResultFile, "definitions", {"spc": photo, "band_spc": band})

    assert ResultFile.from_file(filename, first_try="spc") == "BAND_SPC"
    assert ResultFile.from_file(filename, try_only="band_spc") == "BAND_SPC"


def test_programmer_errors_are_not_swallowed_in_format_detection(monkeypatch, tmp_path):
    import pytest
    from ase2sprkkr.result_files import definition_registry

    filename = tmp_path / "sample.spc"
    filename.write_text("data")

    class Broken:
        name, extension = "SPC", "spc"
        def read_from_file(self, file, potential=None):
            raise RuntimeError("broken parser implementation")

    broken = Broken()
    registry = definition_registry.result_file_definition_registry
    monkeypatch.setattr(registry, "by_extension", lambda: {"spc": (broken,)})
    monkeypatch.setattr(registry, "by_name", lambda: {"spc": broken})
    monkeypatch.setattr(ResultFile, "definitions", {"spc": broken})
    with pytest.raises(RuntimeError, match="broken parser implementation"):
        ResultFile.from_file(filename)
