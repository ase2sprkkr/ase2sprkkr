"""Several SPC parsers can coexist and results select one by task."""
from ase2sprkkr.result_files import ResultFile
from ase2sprkkr.result_files.definition_registry import result_file_definition_registry

def test_native_spc_definitions_are_distinct():
    names = ResultFile.definitions
    assert names["band_spc"].extension == "spc"
    assert names["spleed_spc"].extension == "spc"
    assert names["spc"].extension == "spc"
    assert len(result_file_definition_registry.by_extension()["spc"]) == 3

def test_band_and_spleed_results_register_typed_spectra(tmp_path):
    from ase2sprkkr.outputs.readers.spec import SpecResult
    (tmp_path / "Fe_data.spc").write_text("dummy")
    for task, kind in (("BAND", "band_spc"), ("SPLEED", "spleed_spc")):
        result = SpecResult(None, None, tmp_path, output_file="Fe_SPEC.out", task_name=task)
        assert result.output_values["SPC"].file_type == kind
