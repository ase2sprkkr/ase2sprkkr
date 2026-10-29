"""Preparation paths must describe the same files in Python and in SPR-KKR."""

import os
from pathlib import Path

import pytest
from ase.build import bulk

from ...input_parameters.input_parameters import InputParameters
from ..calculator import SPRKKR


@pytest.fixture(scope="module")
def li_bulk():
    return bulk("Li")


def _calculator(directory, atoms=None):
    return SPRKKR(
        atoms=atoms.copy() if atoms is not None else None,
        directory=directory,
        empty_spheres=False,
        mpi=False,
    )


def test_nested_files_use_calculation_directory_for_serialized_references(
    tmp_path, monkeypatch, li_bulk
):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / "runs" / "case"
    work.mkdir(parents=True)
    shape = tmp_path / "shared" / "shape.sfn"
    calculator = _calculator("runs/case", li_bulk)

    prepared = calculator.prepare(
        input_file="inputs/main.inp",
        potential_file="potentials/main.pot",
        output_file="logs/main.out",
        options={"SFNFIL": str(shape)},
        create_subdirs=True,
    )
    try:
        assert Path(prepared.input_file.name) == work / "inputs/main.inp"
        assert Path(prepared.potential_file) == work / "potentials/main.pot"
        assert Path(prepared.output_file.name) == work / "logs/main.out"
        assert prepared.input_parameters.directory == str(work)
        assert calculator.input_parameters.directory is None

        saved = InputParameters.from_file(prepared.input_file.name)
        assert saved.CONTROL.POTFIL() == "potentials/main.pot"
        assert saved.CONTROL.SFNFIL() == os.path.relpath(shape, work)
        assert (work / saved.CONTROL.POTFIL()).is_file()
        assert (work / saved.CONTROL.SFNFIL()).resolve() == shape
    finally:
        prepared.close()


def test_reused_nested_input_resolves_potfil_against_calculation_directory(
    tmp_path, li_bulk
):
    source = _calculator(tmp_path, li_bulk)
    _, input_name, potential_name, _ = source.save_input(
        input_file="inputs/source.inp",
        potential_file="potentials/source.pot",
        output_file=False,
        create_subdirs=True,
    )
    original = Path(input_name).read_bytes()
    prepared = _calculator(tmp_path).prepare(
        input_parameters="inputs/source.inp", potential=False
    )
    try:
        assert prepared.input_file.name == input_name
        assert prepared.potential_file == potential_name
        assert Path(prepared.potential_file).is_file()
        assert Path(input_name).read_bytes() == original
    finally:
        prepared.close()


@pytest.mark.parametrize("absolute", [False, True])
def test_equivalent_potential_paths_do_not_rewrite_existing_input(
    tmp_path, li_bulk, absolute
):
    source = _calculator(tmp_path, li_bulk)
    _, input_name, potential_name, _ = source.save_input(
        input_file="source.inp", potential_file="source.pot", output_file=False
    )
    original = Path(input_name).read_bytes()
    original_mtime = Path(input_name).stat().st_mtime_ns
    selected = Path(potential_name) if absolute else Path("subdir/../source.pot")

    prepared = _calculator(tmp_path).prepare(
        input_parameters="source.inp", potential=selected
    )
    try:
        assert prepared.input_file.name == input_name
        assert prepared.potential_file == potential_name
        assert Path(input_name).read_bytes() == original
        assert Path(input_name).stat().st_mtime_ns == original_mtime
    finally:
        prepared.close()


def test_potfil_override_is_resolved_after_options_are_applied(tmp_path, li_bulk):
    source = _calculator(tmp_path, li_bulk)
    _, input_name, _, _ = source.save_input(
        input_file="source.inp", potential_file="source.pot", output_file=False
    )
    original = Path(input_name).read_bytes()
    replacement = tmp_path / "replacement.pot"
    replacement.write_text("existing potential, not to be rewritten\n")

    prepared = _calculator(tmp_path).prepare(
        input_parameters="source.inp",
        potential=False,
        input_file="modified.inp",
        options={"POTFIL": "replacement.pot"},
    )
    try:
        assert prepared.potential_file == str(replacement)
        saved = InputParameters.from_file(prepared.input_file.name)
        assert saved.CONTROL.POTFIL() == "replacement.pot"
        assert Path(input_name).read_bytes() == original
        assert replacement.read_text() == "existing potential, not to be rewritten\n"
    finally:
        prepared.close()


def test_caller_owned_file_names_keep_their_process_cwd_base(
    tmp_path, monkeypatch, li_bulk
):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / "work"
    work.mkdir()

    with open("external.inp", "w+") as input_file, open("external.pot", "w+") as potential_file:
        prepared = _calculator(work, li_bulk).prepare(
            input_file=input_file, potential_file=potential_file
        )
        try:
            assert prepared.input_file is input_file
            assert prepared.potential_file == str(tmp_path / "external.pot")
            saved = InputParameters.from_file("external.inp")
            assert saved.CONTROL.POTFIL() == "../external.pot"
        finally:
            prepared.close()


def test_filename_templates_share_one_counter_per_preparation(tmp_path, li_bulk):
    calculator = _calculator(tmp_path, li_bulk)
    for number in (1, 2):
        prepared = calculator.prepare(
            input_file="input_%c.inp",
            potential_file="potential_%c.pot",
            output_file="output_%c.out",
        )
        try:
            assert Path(prepared.input_file.name).name == f"input_{number}.inp"
            assert Path(prepared.potential_file).name == f"potential_{number}.pot"
            assert Path(prepared.output_file.name).name == f"output_{number}.out"
            assert prepared.counter == number
            saved = InputParameters.from_file(prepared.input_file.name)
            assert saved.CONTROL.POTFIL() == f"potential_{number}.pot"
            assert saved.CONTROL.DATASET() == f"input_{number}"
        finally:
            prepared.close()


@pytest.mark.parametrize(
    ("input_name", "output_name"),
    [("case.inp", "case.out"), ("case.in", "case.out"), ("case.data", "case.data.out")],
)
def test_automatic_output_name_follows_reused_input_basename(
    tmp_path, li_bulk, input_name, output_name
):
    source = _calculator(tmp_path, li_bulk)
    source.save_input(
        input_file="inputs/" + input_name,
        potential_file="case.pot",
        output_file=False,
        create_subdirs=True,
    )
    prepared = _calculator(tmp_path).prepare(
        input_parameters="inputs/" + input_name,
        potential=False,
        output_file=True,
    )
    try:
        assert Path(prepared.output_file.name) == tmp_path / output_name
        assert Path(prepared.potential_file) == tmp_path / "case.pot"
    finally:
        prepared.close()


def test_fluctuation_file_is_serialized_relative_to_run_directory(tmp_path):
    parameters = InputParameters.create("GILBERT")
    parameters.CONTROL.POTFIL = "existing.pot"
    parameters.TASK.SETFLUCT = "M_T"
    fluctuation_file = tmp_path / "shared" / "magnetization.dat"
    parameters.TASK.FLUCTFIL = str(fluctuation_file)
    work = tmp_path / "work"
    work.mkdir()

    prepared = _calculator(work).prepare(input_parameters=parameters, potential=False)
    try:
        saved = InputParameters.from_file(prepared.input_file.name)
        assert saved.TASK.FLUCTFIL() == "../shared/magnetization.dat"
        assert tuple(map(str, prepared.input_parameters.input_files(work))) == (
            str(fluctuation_file),
        )
        assert parameters.TASK.FLUCTFIL() == str(fluctuation_file)
    finally:
        prepared.close()


@pytest.mark.parametrize("task", ["TDOBS", "POSANI", "COMPTON"])
@pytest.mark.parametrize("parent", [None, ".", "relative", "absolute"])
@pytest.mark.parametrize("from_dataset", [False, True])
@pytest.mark.parametrize("adsi", [None, "", "CUSTOM"])
def test_tau_path_pair_keeps_native_and_declared_paths_equal(
    tmp_path, task, parent, from_dataset, adsi
):
    work = tmp_path / "work"
    work.mkdir()
    parameters = InputParameters.create(task)
    parameters.CONTROL.POTFIL = "existing.pot"
    if adsi is not None:
        parameters.CONTROL.ADSI = adsi
    dataset = tmp_path / "shared" / "existing"
    suffix = parameters.CONTROL.ADSI() if from_dataset else ""
    tau_directory = Path(str(dataset) + (f"_{suffix}" if suffix else "") + ".tau.dir")
    tau_directory.mkdir(parents=True)
    if from_dataset:
        parameters.CONTROL.DATASET = str(dataset)
    else:
        parameters.TAU.TAUDIR = str(tau_directory)
    if parent == "relative":
        (work / "cache").mkdir()
        parameters.TAU.TAUPATH = "cache"
    elif parent == "absolute":
        parameters.TAU.TAUPATH = str(tmp_path / "shared")
    elif parent is not None:
        parameters.TAU.TAUPATH = parent

    prepared = _calculator(work).prepare(input_parameters=parameters, potential=False)
    try:
        saved = InputParameters.from_file(prepared.input_file.name)
        taupath = saved.TAU.TAUPATH() or "."
        taudir = saved.TAU.TAUDIR()
        assert not os.path.isabs(taupath)
        assert not os.path.isabs(taudir)
        # Native kkrmain.f concatenates TAUPATH, '/', and TAUDIR.
        assert (work / taupath / taudir).resolve() == tau_directory
        assert (work / taupath / taudir).is_dir()
        declared = (
            prepared.input_parameters.cache_files(work)
            if task == "TDOBS"
            else prepared.input_parameters.input_files(work)
        )
        assert str(tau_directory) in tuple(map(str, declared))
        assert parameters.TAU.TAUDIR.is_set() is not from_dataset
    finally:
        prepared.close()


def test_absolute_dataset_prefix_is_serialized_relative_to_run_directory(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    dataset = tmp_path / "results" / "Fe"
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.POTFIL = "existing.pot"
    parameters.CONTROL.DATASET = str(dataset)

    prepared = _calculator(work).prepare(input_parameters=parameters, potential=False)
    try:
        saved = InputParameters.from_file(prepared.input_file.name)
        assert saved.CONTROL.DATASET() == "../results/Fe"
        assert tuple(map(str, prepared.input_parameters.output_files(work))) == (
            str(dataset) + "_DOS.dos",
        )
        assert parameters.CONTROL.DATASET() == str(dataset)
    finally:
        prepared.close()


@pytest.mark.parametrize("task", ["TDOBS", "POSANI", "COMPTON"])
def test_relative_taudir_is_preserved_relative_to_taupath(tmp_path, task):
    parameters = InputParameters.create(task)
    parameters.CONTROL.POTFIL = "existing.pot"
    parameters.CONTROL.DATASET = str(tmp_path / "results" / "Fe")
    parameters.CONTROL.ADSI = ""
    parameters.TAU.TAUDIR = "explicit.tau.dir"
    parameters.TAU.TAUPATH = str(tmp_path / "shared")
    work = tmp_path / "work"
    work.mkdir()
    prepared = _calculator(work).prepare(input_parameters=parameters, potential=False)
    try:
        saved = InputParameters.from_file(prepared.input_file.name)
        assert saved.TAU.TAUDIR() == "explicit.tau.dir"
        assert saved.TAU.TAUPATH() == "../shared"
        assert saved.CONTROL.ADSI() == ""
        declared = (
            prepared.input_parameters.cache_files(work)
            if task == "TDOBS" else prepared.input_parameters.input_files(work)
        )
        assert str(tmp_path / "shared" / "explicit.tau.dir") in tuple(map(str, declared))
    finally:
        prepared.close()


def test_serialization_without_run_directory_keeps_parent_paths(tmp_path, monkeypatch):
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    parameters = InputParameters.create("SCF")
    potential_file = tmp_path / "shared" / "existing.pot"
    shape_file = tmp_path / "shared" / "existing.sfn"
    parameters.CONTROL.POTFIL = str(potential_file)
    parameters.CONTROL.SFNFIL = str(shape_file)

    assert parameters.directory is None
    saved = parameters._definition.read_from_string(parameters.to_string())
    assert saved.CONTROL.POTFIL() == "../shared/existing.pot"
    assert saved.CONTROL.SFNFIL() == "../shared/existing.sfn"
    assert parameters.CONTROL.POTFIL() == str(potential_file)
