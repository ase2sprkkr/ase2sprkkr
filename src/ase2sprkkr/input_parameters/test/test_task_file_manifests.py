import os

import pytest

from ase2sprkkr.input_parameters.input_parameters import InputParameters


DATASET = "manifest_test"


def parameters(task):
    out = InputParameters.create(task)
    out.CONTROL.DATASET = DATASET
    return out


def basenames(paths):
    return {os.path.basename(os.fspath(path)) for path in paths}


def manifest(task, kind, tmp_path):
    p = parameters(task)
    return basenames(getattr(p._definition, f"{kind}_files")(p, tmp_path))


@pytest.mark.parametrize("task", ["AES", "COMPTON", "POSANI", "VBPES"])
def test_required_tau_is_input(task, tmp_path):
    p = parameters(task)
    tau = f"{DATASET}.tau.dir"
    assert tau in manifest(task, "input", tmp_path)
    assert tau not in basenames(p._definition.cache_files(p, tmp_path))


@pytest.mark.parametrize("task", ["BSF", "STONER", "TDOBS", "XAS", "XES", "XRS"])
def test_reusable_tau_is_cache(task, tmp_path):
    p = parameters(task)
    tau = f"{DATASET}.tau.dir"
    assert tau in manifest(task, "cache", tmp_path)
    assert tau not in basenames(p._definition.input_files(p, tmp_path))


def test_non_tau_cache_roles(tmp_path):
    aes = manifest("AES", "cache", tmp_path)
    assert {"ME", *(f"ME{i}" for i in range(1, 13))} <= aes
    assert f"{DATASET}_positron.rho" in manifest("POSANITAU", "cache", tmp_path)


@pytest.mark.parametrize(
    "task, expected",
    [
        ("APS", {f"{DATASET}*.aps"}),
        ("NRAES", {f"{DATASET}*.aes"}),
        ("CLXPS", {f"{DATASET}.pes"}),
        ("FMAG", {f"{DATASET}_fmag_*.agr", f"{DATASET}_dens_*.agr"}),
        ("SCF", {"*_SCFSTART.dos"}),
        ("SIGMA", {"sigma?????.dat"}),
        ("XAS", {f"{DATASET}*.rat"}),
        ("XES", {f"{DATASET}*.rat"}),
        ("XRS", {f"{DATASET}*.rat"}),
    ],
)
def test_audited_output_manifests(task, expected, tmp_path):
    assert expected <= manifest(task, "output", tmp_path)


def test_nr_spectroscopy_conditional_outputs(tmp_path):
    aps = parameters("APS")
    names = basenames(aps._definition.output_files(aps, tmp_path))
    assert f"{DATASET}.dos.brd" not in names
    assert not any(name.endswith(".aps.??") for name in names)
    aps.TASK.GAMMA = 0.2
    aps.TASK.CALCLS = True
    aps.TASK.OUTGNU = True
    names = basenames(aps._definition.output_files(aps, tmp_path))
    assert f"{DATASET}.dos.brd" in names
    assert f"{DATASET}*.aps.??" in names
    assert {"D??[UD][UD]", "E??[UD][UD]", "DE??[UD][UD]", "W??[UD][UD]"} <= names

    nraes = parameters("NRAES")
    nraes.TASK.CALCLS = True
    nraes.TASK.OUTGNU = True
    names = basenames(nraes._definition.output_files(nraes, tmp_path))
    assert f"{DATASET}*.aes.??" in names
    assert {"D??[UD][UD]", "E??[UD][UD]", "DE??[UD][UD]", "W??[UD][UD]"} <= names


def test_posani_previous_result_is_required_input(tmp_path):
    assert f"{DATASET}_pos_dat.pan" in manifest("POSANI", "input", tmp_path)


def test_posani_2d_output_is_conditional(tmp_path):
    p = parameters("POSANI")
    assert not any("pan2D" in name for name in basenames(p._definition.output_files(p, tmp_path)))
    p.TASK.MD2D = True
    assert f"{DATASET}*_pan2D*" in basenames(p._definition.output_files(p, tmp_path))
