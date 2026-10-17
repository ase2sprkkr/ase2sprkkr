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
    assert {f"{DATASET}.tau", "TKTK_tables"} <= manifest("CHI", "cache", tmp_path)


@pytest.mark.parametrize(
    "task, expected",
    [
        ("AES", {f"{DATASET}*.L23.aes"}),
        ("APS", {f"{DATASET}.*.L23.aps"}),
        ("NRAES", {f"{DATASET}.*.L23.aes"}),
        ("CLXPS", {f"{DATASET}*.L23.xps"}),
        ("COMPTON", {f"{DATASET}_BND.mcp", f"{DATASET}_COR.mcp", f"{DATASET}_TOT.mcp"}),
        ("EKREL", {f"{DATASET}.bnd"}),
        ("FMAG", {f"{DATASET}_fmag_*.agr", f"{DATASET}_dens_*.agr"}),
        ("MEPLOT", {f"{DATASET}_VB_*.agr"}),
        ("POSANI", {f"{DATASET}_BND.mpa", f"{DATASET}_COR.mpa", f"{DATASET}_TOT.mpa"}),
        ("SCF", {"*_SCFSTART.dos"}),
        ("SIGMA", {"sigma?????.dat"}),
        ("SOCPAR", {f"{DATASET}_soc_*.agr", f"{DATASET}_exc_*.agr"}),
        ("WFPLOT", {f"{DATASET}_?_*.agr", f"{DATASET}_WF_match_l?_*.agr"}),
        ("XAS", {f"{DATASET}_*_L23.rat"}),
        ("XES", {f"{DATASET}_*_L23.rat"}),
        ("XRS", {f"{DATASET}_*_L23.rat"}),
    ],
)
def test_audited_output_manifests(task, expected, tmp_path):
    assert expected <= manifest(task, "output", tmp_path)


@pytest.mark.parametrize(
    "task, core_level, expected",
    [
        ("AES", "1S", f"{DATASET}*.K.aes"),
        ("CLXPS", "3D", f"{DATASET}*.M45.xps"),
        ("APS", "2S", f"{DATASET}.*.L1.aps"),
        ("NRAES", "5F", f"{DATASET}.*.O67.aes"),
        ("XAS", "1S", f"{DATASET}_*_K.rat"),
        ("XES", "3D", f"{DATASET}_*_M45.rat"),
        ("XRS", "5F", f"{DATASET}_*_O67.rat"),
    ],
)
def test_core_level_output_manifest_tracks_selected_shell(task, core_level, expected, tmp_path):
    p = parameters(task)
    p.TASK.CL = core_level
    assert expected in basenames(p._definition.output_files(p, tmp_path))


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
    assert f"{DATASET}.*.L23.aps.??" in names
    assert {"D??[UD][UD]", "E??[UD][UD]", "DE??", "W??[UD][UD]"} <= names

    nraes = parameters("NRAES")
    nraes.TASK.CALCLS = True
    nraes.TASK.OUTGNU = True
    names = basenames(nraes._definition.output_files(nraes, tmp_path))
    assert f"{DATASET}.*.L23.aes.??" in names
    assert {"D??[UD][UD]", "E??[UD][UD]", "DE??", "W??[UD][UD]"} <= names


def test_momentum_density_output_variants(tmp_path):
    compton = parameters("COMPTON")
    compton.TASK.MD2D = True
    names = basenames(compton._definition.output_files(compton, tmp_path))
    assert {f"{DATASET}_{part}_mcp2d.dat" for part in ("BND", "COR", "TOT")} <= names
    assert f"{DATASET}_BND.mcp" not in names

    compton.TASK.MD2D = False
    compton.TASK.EMD2D = True
    names = basenames(compton._definition.output_files(compton, tmp_path))
    assert {f"{DATASET}_{part}_emd2d.dat" for part in ("BND", "COR", "TOT")} <= names

    posani = parameters("POSANI")
    posani.TASK.MD2D = True
    names = basenames(posani._definition.output_files(posani, tmp_path))
    assert {f"{DATASET}_{part}_pan2d.dat" for part in ("BND", "COR", "TOT")} <= names
    assert f"{DATASET}_BND.mpa" not in names


def test_momentum_density_debug_output(tmp_path):
    for task, extension in (("COMPTON", "mcp"), ("POSANI", "mpa")):
        p = parameters(task)
        p.CONTROL.PRINT = 1
        assert f"{DATASET}_data.{extension}" in basenames(p._definition.output_files(p, tmp_path))


@pytest.mark.parametrize("task", ["ARPES", "AIPES", "BAND", "CLPES", "SPLEED"])
def test_kkrspec_outputs_include_task_log(task, tmp_path):
    assert f"{DATASET}_SPEC.out" in manifest(task, "output", tmp_path)


@pytest.mark.parametrize("task", ["ARPES", "AIPES", "BAND", "CLPES", "SPLEED"])
def test_kkrspec_default_data_output(task, tmp_path):
    names = manifest(task, "output", tmp_path)
    assert f"{DATASET}_data.spc" in names
    assert f"{DATASET}_ARPES_data.spc" not in names


@pytest.mark.parametrize("task", ["ARPES", "AIPES", "CLPES"])
def test_kkrspec_polarization_outputs_follow_pol_p(task, tmp_path):
    p = parameters(task)
    p.SPEC_PH.POL_P = "S"
    names = basenames(p._definition.output_files(p, tmp_path))
    assert f"{DATASET}_data.spc" not in names
    assert {f"{DATASET}_data_pol1.spc", f"{DATASET}_data_pol2.spc"} <= names


def test_mecheck_enabled_outputs(tmp_path):
    p = parameters("MECHECK")
    assert not basenames(p._definition.output_files(p, tmp_path))

    p.TASK.INITIALSTATE = "CORE"
    p.ENERGY.NE = 3
    names = basenames(p._definition.output_files(p, tmp_path))
    assert {f"MPZAZB_{polarization}_*-*-*.dat" for polarization in ("-", "z", "+")} <= names
    assert not any("_check_" in name for name in names)


def test_pshift_branch_outputs(tmp_path):
    p = parameters("PSHIFT")
    names = basenames(p._definition.output_files(p, tmp_path))
    assert f"{DATASET}_pshift_*.agr" in names
    assert f"{DATASET}_logdrv_*.agr" in names
    assert f"{DATASET}_t-matrix_*.agr" not in names

    p.MODE.MODE = "SREL"
    names = basenames(p._definition.output_files(p, tmp_path))
    assert f"{DATASET}_pshift_*.agr" in names
    assert f"{DATASET}_t-matrix_*.agr" in names
    assert f"{DATASET}_logdrv_*.agr" not in names


def test_meplot_output_variants(tmp_path):
    p = parameters("MEPLOT")
    p.TASK.L = 2
    assert f"{DATASET}_d_*.agr" in basenames(p._definition.output_files(p, tmp_path))

    p.TASK.SEFIN = 2
    names = basenames(p._definition.output_files(p, tmp_path))
    assert f"{DATASET}_d_SE2_*.agr" in names
    assert f"{DATASET}_SEfin_E_SE2_*.agr" in names
    assert f"{DATASET}_SEfin_r_SE2_*.agr" in names

    p.TASK.INITIALSTATE = "CORE"
    p.TASK.CL = "3D"
    p.TASK.SEFIN = 0
    assert f"{DATASET}_3D_*.agr" in basenames(p._definition.output_files(p, tmp_path))


def test_wfplot_output_variants(tmp_path):
    p = parameters("WFPLOT")
    p.TASK.L = 3
    names = basenames(p._definition.output_files(p, tmp_path))
    assert f"{DATASET}_f_*.agr" in names
    assert f"{DATASET}_WF_match_lf_*.agr" in names

    p.TASK.STATE = "CORE"
    p.TASK.CL = "3D"
    names = basenames(p._definition.output_files(p, tmp_path))
    assert f"{DATASET}_3D_*.agr" in names
    assert not any("WF_match" in name for name in names)


def test_vbpes_final_state_tau_cache(tmp_path):
    p = parameters("VBPES")
    assert not any(".tau.dir.F_" in name for name in basenames(p._definition.cache_files(p, tmp_path)))
    assert not any(".dos.F_" in name for name in basenames(p._definition.output_files(p, tmp_path)))

    p.TASK.USETAUNN = True
    assert f"{DATASET}.tau.dir.F_1253" in basenames(p._definition.cache_files(p, tmp_path))
    assert f"{DATASET}.dos.F_1253" in basenames(p._definition.output_files(p, tmp_path))

    p.TASK.EPHOT = 750.9
    assert f"{DATASET}.tau.dir.F_750" in basenames(p._definition.cache_files(p, tmp_path))
    assert f"{DATASET}.dos.F_750" in basenames(p._definition.output_files(p, tmp_path))

    p.TASK.EPHOT = 1.0
    assert f"{DATASET}.tau.dir.F_1253" in basenames(p._definition.cache_files(p, tmp_path))

    p.TASK.USETFIN = True
    assert not any(".tau.dir.F_" in name for name in basenames(p._definition.cache_files(p, tmp_path)))
    assert not any(".dos.F_" in name for name in basenames(p._definition.output_files(p, tmp_path)))


def test_posani_previous_result_is_required_input(tmp_path):
    assert f"{DATASET}_pos_dat.pan" in manifest("POSANI", "input", tmp_path)
