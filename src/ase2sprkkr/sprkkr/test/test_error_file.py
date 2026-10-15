import os
from pathlib import Path

from ase.build import bulk

from ...execution.base import RunState
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


def _calculator(tmp_path):
    return SPRKKR(
        atoms=bulk("Fe"),
        directory=tmp_path,
        profile=SPRKKRProfile(mpi=False),
        empty_spheres=False,
    )


def test_default_error_file_is_output_file_plus_err(tmp_path):
    runner = _calculator(tmp_path).create_run(executor=_Executor())
    try:
        output_name = runner.spec.output_file.name
        assert runner.spec.error_file.name == output_name + ".err"
        assert runner.spec.error_file is not runner.spec.output_file
        assert runner.spec.error_file in runner.spec.output_files
    finally:
        runner.close()


def test_error_file_true_explicitly_merges_with_output(tmp_path):
    runner = _calculator(tmp_path).create_run(
        executor=_Executor(),
        error_file=True,
    )
    try:
        assert runner.spec.error_file is runner.spec.output_file
        assert sum(
            item is runner.spec.output_file for item in runner.spec.output_files
        ) == 1
    finally:
        runner.close()


def test_explicit_error_file_path_is_resolved_in_calculation_directory(tmp_path):
    runner = _calculator(tmp_path).create_run(
        executor=_Executor(),
        error_file="logs/stderr.txt",
        create_subdirs=True,
    )
    try:
        assert runner.spec.error_file.name == os.fspath(
            tmp_path / "logs" / "stderr.txt"
        )
        assert Path(runner.spec.error_file.name).is_file()
    finally:
        runner.close()


def test_error_file_false_requests_managed_temporary_file(tmp_path):
    runner = _calculator(tmp_path).create_run(
        executor=_Executor(),
        error_file=False,
    )
    error_name = runner.spec.error_file.name
    try:
        assert runner.spec.error_file is not runner.spec.output_file
        assert Path(error_name).is_file()
    finally:
        runner.close()
    assert not Path(error_name).exists()
