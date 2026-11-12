"""Prepared SPR-KKR resources belong to the generic Run lifecycle."""

import sys
from threading import Event
from types import SimpleNamespace

import pytest

from kkrunner import RunState
from kkrunner.backends.contracts import SubmissionNotStarted
from kkrunner.client import EndpointClient
from kkrunner.runs import Run

from ase2sprkkr.sprkkr.calculation import PreparedCalculation
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


class _Task:
    name = "SCF"
    # Python reads the prepared input file from stdin, so this exercises the
    # real direct runner without requiring an installed SPR-KKR executable.
    executable = sys.executable
    mpi = False

    def output_reader(self, **kwargs):
        return SimpleNamespace(
            read_from_file=lambda *args, **options: SimpleNamespace(
                return_code=options.get("return_code"),
            ),
        )


class _Parameters:
    CONTROL = SimpleNamespace(
        POTFIL=SimpleNamespace(result=None),
        DATASET=SimpleNamespace(result=None),
    )

    def task_definition_required(self):
        return _Task()

    def input_files(self, directory):
        return ()

    def output_files(self, directory):
        return ()

    def cache_files(self, directory):
        return ()


class _Prepared:
    create_run = PreparedCalculation.create_run

    def __init__(self, input_file, output_file, directory):
        self.input_parameters = _Parameters()
        self.input_file = input_file
        self.potential_file = None
        self.output_file = output_file
        self.directory = str(directory)
        self.closed = False

    def close(self):
        if not self.closed:
            self.closed = True
            self.input_file.close()
            self.output_file.close()


def _prepared(tmp_path):
    input_path = tmp_path / "run.inp"
    input_path.write_text("print('ok')\n")
    input_file = input_path.open("rb")
    output_file = (tmp_path / "run.out").open("w+b")
    return _Prepared(input_file, output_file, tmp_path)


def _create_run(prepared, *, executor=None, callback=None):
    return prepared.create_run(
        result_context=None,
        profile=SPRKKRProfile(mpi=False),
        executor=executor or EndpointClient(),
        callback=callback,
    )


def test_files_are_closed_when_runner_creation_fails(tmp_path, monkeypatch):
    prepared = _prepared(tmp_path)
    executor = EndpointClient()

    def fail_factory(*args, **kwargs):
        raise RuntimeError("factory failed")

    monkeypatch.setattr(executor, "create_run_parts", fail_factory)
    with pytest.raises(RuntimeError, match="factory failed"):
        _create_run(prepared, executor=executor)
    assert prepared.closed
    assert prepared.input_file.closed
    assert prepared.output_file.closed


def test_files_are_closed_when_submission_fails_before_start(tmp_path, monkeypatch):
    prepared = _prepared(tmp_path)
    runner = _create_run(prepared)
    assert isinstance(runner, Run)
    assert not prepared.closed

    def fail_submission():
        raise SubmissionNotStarted("staging failed")

    monkeypatch.setattr(runner._job, "submit", fail_submission)
    with pytest.raises(SubmissionNotStarted, match="staging failed"):
        runner.submit()
    assert runner.state is RunState.FAILED
    assert prepared.closed
    assert prepared.input_file.closed
    assert prepared.output_file.closed


def test_completion_releases_files_without_wait(tmp_path):
    prepared = _prepared(tmp_path)
    completed = Event()
    runner = _create_run(prepared, callback=lambda result: completed.set())
    runner.submit()
    assert completed.wait(10), "completion callback was not called"
    assert runner.state is RunState.FINISHED
    assert prepared.closed
    assert prepared.input_file.closed
    assert prepared.output_file.closed


def test_unstarted_runner_can_release_owned_files_explicitly(tmp_path):
    prepared = _prepared(tmp_path)
    runner = _create_run(prepared)
    runner.close()
    runner.close()
    assert prepared.closed
    assert prepared.input_file.closed
    assert prepared.output_file.closed
