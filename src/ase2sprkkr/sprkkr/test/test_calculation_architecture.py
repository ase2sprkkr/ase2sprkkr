import inspect
import os
from pathlib import Path

from ase.build import bulk
import pytest

from ...execution.base import RunState
from ...input_parameters.input_parameters import InputParameters
from ...potentials.potentials import Potential
from .. import calculation
from ..calculator import SPRKKR
from ..profile import SPRKKRProfile


class _UnstartedRunner:
    def __init__(self, spec):
        self.spec = spec
        self.state = RunState.NEW
        self.result = None
        self.reference = None
        self.exception = None

    def add_output_listener(self, callback):
        pass


class _Executor:
    def resolve_mpi(self, mpi):
        return mpi

    def create_run(self, reader, spec, callback=None):
        return _UnstartedRunner(spec)


def test_input_parameters_do_not_own_execution():
    assert not hasattr(InputParameters, "run_process")
    assert not hasattr(InputParameters, "calculate")
    assert InputParameters.create("SCF").task_definition.executable


def test_preparation_is_independent_of_execution_profile():
    parameters = inspect.signature(calculation.PreparedCalculation.prepare).parameters
    assert "profile" not in parameters
    assert "executor" not in parameters
    assert "persistent" not in parameters
    assert callable(calculation.PreparedCalculation.create_run)


def test_calculator_exposes_preparation_as_public_boundary():
    assert callable(SPRKKR.prepare)
    assert not hasattr(SPRKKR, "_prepare_calculation")
    assert "return_run" not in inspect.signature(SPRKKR.calculate).parameters


def test_unstarted_run_owns_temporary_calculation_directory():
    calculator = SPRKKR(
        atoms=bulk("Fe"),
        profile=SPRKKRProfile(mpi=False),
        empty_spheres=False,
    )

    runner = calculator.create_run(
        directory=False,
        executor=_Executor(),
        empty_spheres=False,
    )
    directory = runner.spec.directory

    assert os.path.isdir(directory)
    assert os.path.exists(runner.spec.input_file.name)

    runner.close()

    assert not os.path.exists(directory)


def test_temporary_input_and_potential_are_path_backed_until_runner_close(tmp_path):
    calculator = SPRKKR(
        atoms=bulk("Fe"),
        directory=tmp_path,
        profile=SPRKKRProfile(mpi=False),
        empty_spheres=False,
    )

    runner = calculator.create_run(
        input_file=False,
        potential_file=False,
        output_file=False,
        executor=_Executor(),
    )

    input_name = runner.spec.input_file.name
    potential_names = [
        value for value in runner.spec.input_files if isinstance(value, str)
    ]

    assert isinstance(input_name, str)
    assert Path(input_name).is_file()
    assert len(potential_names) == 1
    assert Path(potential_names[0]).is_file()

    runner.close()

    assert not Path(input_name).exists()
    assert not Path(potential_names[0]).exists()


def test_preparation_snapshots_input_parameters(tmp_path):
    calculator = SPRKKR(
        atoms=bulk("Fe"),
        directory=tmp_path,
        profile=SPRKKRProfile(mpi=False),
        empty_spheres=False,
    )
    defaults = calculator.input_parameters
    original_potfil = defaults.CONTROL.POTFIL.result

    prepared = calculator.prepare(
        output_file=False,
        allow_temporary_files=True,
    )

    try:
        assert prepared.input_parameters is not defaults
        assert prepared.input_parameters.CONTROL.POTFIL.result
        assert defaults.CONTROL.POTFIL.result == original_potfil

        prepared_ne = prepared.input_parameters.ENERGY.NE()
        defaults.ENERGY.NE = prepared_ne + 1
        assert prepared.input_parameters.ENERGY.NE() == prepared_ne
    finally:
        prepared.close()


def test_effective_task_survives_preparation_snapshot(tmp_path):
    parameters = InputParameters.create("ARPES")
    parameters.TASK.TASK = "SPLEED"
    calculator = SPRKKR(
        atoms=bulk("Fe"),
        directory=tmp_path,
        profile=SPRKKRProfile(mpi=False),
        empty_spheres=False,
    )

    prepared = calculator.prepare(
        input_parameters=parameters,
        output_file=False,
        allow_temporary_files=True,
    )

    try:
        assert prepared.input_parameters.task_name == "ARPES"
        assert prepared.input_parameters.task_definition.name == "SPLEED"
        assert prepared.input_parameters.execution_task_name == "SPLEED"
    finally:
        prepared.close()


def test_prepared_calculation_constructs_its_runner(tmp_path):
    calculator = SPRKKR(
        atoms=bulk("Fe"),
        directory=tmp_path,
        profile=SPRKKRProfile(mpi=False),
        empty_spheres=False,
    )
    prepared = calculator.prepare(output_file=False)

    runner = prepared.create_run(
        result_context=calculator,
        profile=calculator.profile,
        executor=_Executor(),
    )

    assert runner.spec.directory == str(tmp_path)
    assert runner.spec.input_file is prepared.input_file

    runner.close()
    assert prepared.closed


def test_failed_potential_serialization_preserves_primary_error_and_cleans_files(
    monkeypatch, tmp_path
):
    calculator = SPRKKR(
        atoms=bulk("Fe"),
        directory=tmp_path,
        profile=SPRKKRProfile(mpi=False),
        empty_spheres=False,
    )
    opened = []
    real_open_file = calculation._open_file

    def tracked_open_file(*args, **kwargs):
        file = real_open_file(*args, **kwargs)
        opened.append(file)
        return file

    def fail_save(self, file, atoms=None, **kwargs):
        raise RuntimeError("potential serialization failed")

    monkeypatch.setattr(calculation, "_open_file", tracked_open_file)
    monkeypatch.setattr(Potential, "save_to_file", fail_save)

    with pytest.raises(RuntimeError, match="potential serialization failed"):
        calculator.save_input()

    assert len(opened) == 2
    assert all(file.closed for file in opened)


def test_reader_file_not_found_is_not_relabelled_as_missing_executable(
    monkeypatch, tmp_path
):
    calculator = SPRKKR(
        atoms=bulk("Fe"),
        directory=tmp_path,
        profile=SPRKKRProfile(mpi=False),
        empty_spheres=False,
    )
    prepared = calculator.prepare(
        allow_temporary_files=True,
        output_file=False,
        empty_spheres=False,
    )

    def missing_reader(self, *args, **kwargs):
        raise FileNotFoundError("reader dependency missing")

    monkeypatch.setattr(InputParameters, "output_reader", missing_reader)

    with pytest.raises(FileNotFoundError, match="reader dependency missing"):
        prepared.create_run(
            result_context=calculator,
            profile=calculator.profile,
            executor=_Executor(),
        )

    assert prepared.closed
    assert prepared.input_file.closed
