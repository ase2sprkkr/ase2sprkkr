"""Filename contracts for native DATSET0 versus DATSET in SPR-KKR 9.7."""

import io
import pickle

import pytest

from ..input_parameters import InputParameters
from ..task_definition import InputParametersRegistry, TaskFile


@pytest.mark.parametrize("adsi, effective", [("", "Fe"), ("CUSTOM", "Fe_CUSTOM")])
def test_dataset_template_fields_and_compatibility_alias(tmp_path, adsi, effective):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "Fe"
    parameters.CONTROL.ADSI = adsi
    declared = InputParametersRegistry._resolve_file_templates(
        parameters, tmp_path,
        ("{DATASET}.raw", "{dataset}.legacy", TaskFile("{DATASET_ADSI}.dos", type="DOS")),
    )
    assert tuple(map(str, declared)) == (
        str(tmp_path / "Fe.raw"), str(tmp_path / "Fe.legacy"),
        str(tmp_path / f"{effective}.dos"),
    )
    assert declared[-1].type == "DOS"


# These names come from the routines' DATSET/DOSFIL/TAUDIR expressions,
# not from a task-name convention (BSF and DMI have distinct ADSI defaults).
@pytest.mark.parametrize("task, adsi_default, kind, suffix", [
    ("DOS", "DOS", "output", ".dos"),
    ("TDOBS", "TDOBS", "output", ".dos"),
    ("BSFEK", "BLOCHSF", "output", ".bsf"),
    ("BSFKK", "BLOCHSF", "output", "_spol.bsf"),
    ("ARPES", "ARPES", "output", "_data.spc"),
    ("AIPES", "AIPES", "output", "_SPEC.out"),
    ("BAND", "BAND", "output", "_data.spc"),
    ("SPLEED", "SPLEED", "output", "_SPEC.out"),
    ("CLPES", "CLPES", "output", "_corewf"),
    ("XAS", "XAS", "output", "_*_L23.rat"),
    ("XES", "XES", "cache", ".tau.dir"),
    ("XRS", "XRS", "cache", ".tau.dir"),
    ("AES", "AES", "input", ".tau.dir"),
    ("CLXPS", "CLXPS", "output", "*.L23.xps"),
    ("APS", "APS", "input", ".dos"),
    ("NRAES", "NRAES", "output", ".*.L23.aes"),
    ("COMPTON", "COMPTON", "output", "_TOT.mcp"),
    ("POSANI", "POSANI", "input", "_pos_dat.pan"),
    ("POSANI", "POSANI", "output", "_TOT.mpa"),
    ("POSANIPREP", "POSANIPREP", "output", "_pos_dat.pan"),
    ("POSANITAU", "POSANITAU", "cache", "_positron.rho"),
    ("BSDOSIJ", "BSDOSIJ", "output", "_XCPLTEN_Jij.dat"),
    ("JXC", "JXC", "output", "_XCPLTEN_Jij.dat"),
    ("DMI", "JXC", "output", "_DMIVEC_Dij.dat"),
    ("CHI", "CHI", "cache", ".tau"),
    ("CHIXAS", "CHIXAS", "output", "_*_L23.rat"),
    ("EKREL", "EKREL", "output", ".bnd"),
    ("FMAG", "FMAG", "output", "_fmag_*.agr"),
    ("FSCAT", "FSCAT", "output", "_GEO?_*.agr"),
    ("GILBERT", "GILBERT", "output", "alfa_T.dat"),
    ("MEPLOT", "MEPLOT", "output", "_VB_*.agr"),
    ("NUCSSCPL", "NUCSSCPL", "output", "_A_ij.dat"),
    ("PSHIFT", "PSHIFT", "output", "_pshift_*.agr"),
    ("SOCPAR", "SOCPAR", "output", "_soc_*.agr"),
    ("STONER", "STONER", "output", "_*_g-spin.dat"),
    ("VBPES", "VBPES", "input", ".tau.dir"),
    ("WFPLOT", "WFPLOT", "output", "_WF_match_l?_*.agr"),
])
@pytest.mark.parametrize("adsi", [None, "", "CUSTOM"])
def test_native_task_prefixes(tmp_path, task, adsi_default, kind, suffix, adsi):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = "Fe"
    assert parameters.CONTROL.ADSI() == adsi_default
    if adsi is not None:
        parameters.CONTROL.ADSI = adsi
    effective_adsi = adsi_default if adsi is None else adsi
    prefix = "Fe" + ("_" + effective_adsi if effective_adsi else "")
    assert str(tmp_path / (prefix + suffix)) in tuple(
        map(str, getattr(parameters, f"{kind}_files")(tmp_path))
    )


@pytest.mark.parametrize("task, suffix", [
    ("JALF", "_JALF.dat"), ("JORB", "_JORB.dat"),
    ("PLOTPOT", "_V_*.agr"), ("PLOTRHO", "_RHO_*.agr"),
    ("PLOTSFN", "_SFN_M*.agr"), ("PLOTPRS", "_SFN_M*.agr"),
])
def test_native_datset0_outputs_ignore_adsi(tmp_path, task, suffix):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = "Fe"
    parameters.CONTROL.ADSI = "CUSTOM"
    assert str(tmp_path / ("Fe" + suffix)) in parameters.output_files(tmp_path)
    assert parameters.CONTROL.SFNFIL() == "Fe.sfn"


@pytest.mark.parametrize("task", ["SCF", "DOS", "BSFEK", "TDOBS", "ARPES"])
def test_empty_adsi_roundtrip_and_pickle(task):
    parameters = InputParameters.create(task)
    parameters.CONTROL.POTFIL = "Fe.pot"
    parameters.CONTROL.ADSI = ""
    if task == "BSFEK":
        parameters.TASK.KPATH = 1
    text = parameters.to_string()
    assert "ADSI" not in text
    parsed = parameters._definition.read_from_string(text)
    assert parsed.CONTROL.ADSI() == ""
    assert "ADSI" not in parsed.to_string()
    assert pickle.loads(pickle.dumps(parsed)).CONTROL.ADSI() == ""
    # Reading into an existing object must replace the old suffix atomically.
    parameters.CONTROL.ADSI = "OLD"
    parameters.read_from_file(io.StringIO(text))
    assert parameters.CONTROL.ADSI() == ""


@pytest.mark.parametrize("npol", [0, 3])
@pytest.mark.parametrize("adsi, prefix", [("", "Fe"), ("CUSTOM", "Fe_CUSTOM")])
def test_spec_polarized_output_names(tmp_path, npol, adsi, prefix):
    parameters = InputParameters.create("ARPES")
    parameters.CONTROL.DATASET = "Fe"
    parameters.CONTROL.ADSI = adsi
    parameters.SPEC_PH.NPOL = npol
    assert tuple(map(str, parameters.output_files(tmp_path))) == (
        str(tmp_path / f"{prefix}_SPEC.out"),
        str(tmp_path / f"{prefix}_data_pol1.spc"),
        str(tmp_path / f"{prefix}_data_pol2.spc"),
    )


@pytest.mark.parametrize("adsi, filename", [(None, "Fe_DOS.dos"), ("", "Fe.dos"), ("CUSTOM", "Fe_CUSTOM.dos")])
def test_dos_result_filename_matches_manifest(tmp_path, adsi, filename):
    from ...outputs.readers.dos import DosResult

    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "Fe"
    if adsi is not None:
        parameters.CONTROL.ADSI = adsi
    result = DosResult(parameters, None, tmp_path)
    assert result.dos_filename == str(tmp_path / filename)
    assert result.dos_filename == str(parameters.output_files(tmp_path)[0])
