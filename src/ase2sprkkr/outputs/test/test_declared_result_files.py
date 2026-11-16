"""Reconcile task output declarations with parser-observed result files."""
import os
from pathlib import Path

import pytest

from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.outputs.readers.arpes import ArpesResult
from ase2sprkkr.outputs.readers.bsf import BsfResult
from ase2sprkkr.outputs.readers.dos import DosResult
from ase2sprkkr.outputs.readers.jxc import JxcResult
from ase2sprkkr.outputs.readers.spec import SpecResult
from ase2sprkkr.outputs.task_result import TaskResult


def _parameters(task, dataset="Fe"):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = dataset
    parameters.CONTROL.POTFIL = f"{dataset}.pot"
    return parameters


def _declare_file(parameters, directory, suffix):
    return next(
        Path(os.fspath(file))
        for file in parameters.output_files(str(directory))
        if os.fspath(file).endswith(suffix)
    )


def test_dos_is_registered_without_an_out_files_table_entry(tmp_path, monkeypatch):
    parameters = _parameters("DOS")
    output = tmp_path / "run.out"
    output.write_bytes(b"OUT")
    dos = _declare_file(parameters, tmp_path, ".dos")
    dos.write_bytes(b"DOS data")

    # Emulate a parsed OUT stream which does not list the declared DOS file.
    reader = parameters.output_reader(directory=str(tmp_path))
    monkeypatch.setattr(
        reader.parser,
        "read_from_file",
        lambda output_stream, error_stream, results, code: (None, None, code),
    )
    result = reader.read_from_file(str(output))

    assert isinstance(result, DosResult)
    assert result.files["DOS"].path() == str(dos)
    assert result.output_values["DOS"].file_type == "dos"
    assert result.dos_filename == str(dos)
    assert len([v for v in result.files.items().values() if v.path() == str(dos)]) == 1

    # Reconciliation is idempotent and must not duplicate the GUI entry.
    assert result.enrich_declared_output_files() == ()
    assert len([v for v in result.files.items().values() if v.path() == str(dos)]) == 1


def test_missing_mandatory_file_warns_but_is_not_a_result(tmp_path):
    parameters = _parameters("DOS")
    output = tmp_path / "run.out"
    output.write_bytes(b"OUT")
    result = DosResult(parameters, None, str(tmp_path), output_file="run.out", task_name="DOS")

    with pytest.warns(RuntimeWarning, match="missing required output"):
        result.complete(None, 0)

    assert "DOS" not in result.files
    assert "DOS" not in result.output_values
    assert "output" in result.output_values
    manifest = result.to_file_manifest()
    expected = _declare_file(parameters, tmp_path, ".dos")
    declared = next(x for x in manifest.outputs if x.local_path == str(expected))
    assert not declared.optional
    assert declared.availability_state == "failed"


def test_failed_run_does_not_warn_about_unproduced_outputs(tmp_path):
    parameters = _parameters("DOS")
    result = DosResult(parameters, None, str(tmp_path), task_name="DOS")
    import warnings

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result.complete("execution failed", 1)
    assert not caught
    assert "DOS" not in result.files


def test_bsf_main_and_optional_spin_output_keep_distinct_keys(tmp_path):
    parameters = _parameters("BSFEK")
    output_files = parameters.output_files(str(tmp_path))
    main = next(Path(os.fspath(x)) for x in output_files if not x.optional)
    spin = next(Path(os.fspath(x)) for x in output_files if x.optional)
    main.write_bytes(b"BSF data")

    result = BsfResult(parameters, None, str(tmp_path), task_name="BSFEK")
    result.complete(None, 0)
    assert result.files["Bloch-SF"].path() == str(main)
    assert spin.name not in result.files
    assert "Bloch-SF" in result.output_values

    spin.write_bytes(b"spin BSF")
    assert result.enrich_declared_output_files() == (spin.name,)
    assert result.files[spin.name].file_type == "bsf"
    assert list(result.files.items()).count("Bloch-SF") == 1


def test_jxc_registers_existing_dat_and_agr_without_overwriting_parser_files(tmp_path):
    parameters = _parameters("JXC")
    declared = parameters.output_files(str(tmp_path))
    jij = next(
        Path(os.fspath(x)) for x in declared
        if os.fspath(x).endswith("_XCPLTEN_Jij.dat")
    )
    plot_pattern = next(
        os.fspath(x) for x in declared
        if "_XCPLTEN_Jij_" in os.fspath(x) and ".agr" in os.fspath(x)
    )
    plot = Path(plot_pattern.replace("*", "Fe"))
    jij.write_bytes(b"couplings")
    plot.write_bytes(b"agr data")

    # A historical file manifest needs the primary OUT as its provenance root.
    (tmp_path / "run.out").write_bytes(b"OUT")
    result = JxcResult(
        parameters, None, str(tmp_path), output_file="run.out", task_name="JXC",
    )
    # A parser-provided entry wins, even if the task also declares its path.
    observed = result.files.add_file("jxc", str(jij), "jxc")
    with pytest.warns(RuntimeWarning, match="missing required output"):
        result.complete(None, 0)

    assert result.files["jxc"] is observed
    assert result.files[plot.name].path() == str(plot)
    assert plot.name in result.output_values
    assert "open" in result.output_values[plot.name].actions()
    assert len([v for v in result.files.items().values() if v.path() == str(jij)]) == 1

    manifest = result.to_file_manifest()
    assert str(plot) in {entry.local_path for entry in manifest.outputs}


def test_existing_unreported_jxc_data_keeps_public_accessor(tmp_path):
    parameters = _parameters("JXC")
    jij = next(
        Path(os.fspath(x)) for x in parameters.output_files(str(tmp_path))
        if os.fspath(x).endswith("_XCPLTEN_Jij.dat")
    )
    jij.write_bytes(b"couplings")

    result = JxcResult(parameters, None, str(tmp_path), task_name="JXC")
    result.enrich_declared_output_files(warn_missing=False)
    assert result.files["jxc"].path() == str(jij)
    assert result.output_values["jxc"].file_type == "jxc"


def test_existing_optional_scf_startup_file_becomes_result(tmp_path):
    parameters = _parameters("SCF")
    (tmp_path / "Fe.pot_new").write_bytes(b"potential")
    startup = tmp_path / "Fe_SCFSTART.dos"
    startup.write_bytes(b"DOS from SCF startup")

    result = TaskResult(parameters, None, str(tmp_path), task_name="SCF")
    result.complete(None, 0)
    assert startup.name in result.output_values
    assert result.files[startup.name].file_type == "dos"


def test_spectra_stay_owned_by_the_existing_spec_registration(tmp_path):
    parameters = _parameters("ARPES")
    spc = _declare_file(parameters, tmp_path, ".spc")
    spc.write_bytes(b"SPC result")

    result = ArpesResult(parameters, None, str(tmp_path), task_name="ARPES")
    with pytest.warns(RuntimeWarning, match="missing required output"):
        result.complete(None, 0)
    # SPC is registered by ArpesResult, not added under an unrelated filename.
    assert spc.name not in result.files
    assert "SPC" in result.output_values
    assert result.files["SPC"].path() == str(spc)
    assert sum(x.path() == str(spc) for x in result.files.items().values()) == 1


def test_spec_header_only_spc_is_not_promoted_to_spectrum(tmp_path):
    parameters = _parameters("SPLEED")
    parameters.SPEC_EL.TYP = 0  # The i(E) diagram is embedded in SPEC.out.
    spc = _declare_file(parameters, tmp_path, ".spc")
    spc.write_bytes(b"header only")

    result = SpecResult(parameters, None, str(tmp_path), task_name="SPLEED")
    result._file_content = ""
    with pytest.warns(RuntimeWarning, match="missing required output"):
        result.complete(None, 0)
    assert spc.name not in result.files
    assert "SPC" not in result.output_values


def test_only_existing_files_are_exposed_even_with_wildcard_declarations(tmp_path):
    parameters = _parameters("GILBERT")
    declared = parameters.output_files(str(tmp_path))
    plot = next(Path(os.fspath(x)) for x in declared if os.fspath(x).endswith(".agr"))
    plot.write_bytes(b"agr test")

    result = TaskResult(parameters, None, str(tmp_path), task_name="GILBERT")
    with pytest.warns(RuntimeWarning, match="missing required output"):
        result.complete(None, 0)

    assert plot.name in result.output_values
    assert "open" in result.output_values[plot.name].actions()
    assert all(
        os.path.isfile(option.path()) for key, option in result.output_values.items()
        if key != "output"
    )


def test_historical_mismatched_schema_skips_declarations_but_keeps_result(tmp_path):
    parameters = _parameters("ARPES")
    parameters.TASK.TASK = "SPLEED"
    output = tmp_path / "run.out"
    output.write_text("historical output")

    # Execution validation stays strict: this is a historical, not runnable,
    # combination of effective task and input schema.
    with pytest.raises(RuntimeError, match="different input schema"):
        parameters.output_files(str(tmp_path))

    result = TaskResult(
        parameters, None, str(tmp_path), output_file="run.out", task_name="SPLEED",
    )
    result.complete(None, 0)
    assert result.task_name == "spleed"
    assert result.input_parameters.definition_name == "ARPES"
    assert tuple(result.output_values) == ("output",)

    # Do not infer SPLEED outputs using the incompatible ARPES schema.
    manifest = result.to_file_manifest()
    assert manifest.stdout.local_path == str(output)
    assert not manifest.outputs
    assert any(item.local_path.endswith("Fe.pot") for item in manifest.inputs)

