import asyncio
import threading
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution.base import ExecutionSpec, LaunchSpec, RunState
from ase2sprkkr.execution import (
    RunRecord, Execution, WorkingDirectory, ExecutorEndpoint, FileManifest, FileEntry,
)
from ase2sprkkr.execution import psij as psij_backend


class _Result:
    def complete(self, error, return_code):
        self.error = error
        self.return_code = return_code


class _Parser:
    print_output = False

    async def read_output(self, stdout, result):
        result.output = await stdout.read()
        return result.output

    async def read_error(self, stderr, result):
        result.error_output = await stderr.read()
        return result.error_output

    def result(self, output, error, return_code):
        if return_code:
            raise ValueError(f"Process returned {return_code}")
        return output, error, return_code


class _Reader:
    print_output = False

    def __init__(self):
        self.parser = _Parser()

    def create_result(self, output_file, input_file):
        return _Result()

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(
            self.parser.read_output(stdout, result),
            self.parser.read_error(stderr, result),
        )

    def complete_result(self, result, output, error, return_code):
        _, error, return_code = self.parser.result(output, error, return_code)
        result.complete(error, return_code)
        return result


class _Status:
    def __init__(self, name="NEW", exit_code=None):
        self.state = SimpleNamespace(name=name)
        self.exit_code = exit_code


class _JobSpec:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class _JobAttributes:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class _ResourceSpecV1:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class _Job:
    def __init__(self, spec=None):
        self.spec = spec
        self.native_id = None
        self.status = _Status()
        self._done = threading.Event()

    def wait(self):
        self._done.wait(2)
        return self.status

    def cancel(self):
        self.status = _Status("CANCELED", 1)
        self._done.set()


class _Backend:
    def __init__(self):
        self.submitted = None
        self.attached = None

    def submit(self, job):
        self.submitted = job
        job.native_id = "123"
        job.status = _Status("QUEUED")

        def run():
            job.status = _Status("ACTIVE")
            with open(job.spec.stdout_path, "wb") as out:
                out.write(b"hello\n")
            with open(job.spec.stderr_path, "wb") as err:
                err.write(b"warning\n")
            job.status = _Status("COMPLETED", 0)
            job._done.set()

        threading.Thread(target=run, daemon=True).start()

    def attach(self, job, native_id):
        self.attached = native_id
        job.native_id = native_id
        job.status = _Status("COMPLETED", 0)
        job._done.set()


class _JobExecutor:
    @staticmethod
    def get_instance(*args, **kwargs):
        raise AssertionError("A supplied backend should be used")


_FAKE_PSIJ = SimpleNamespace(
    Job=_Job,
    JobSpec=_JobSpec,
    JobAttributes=_JobAttributes,
    JobExecutor=_JobExecutor,
    ResourceSpecV1=_ResourceSpecV1,
)


def _executor(monkeypatch, backend):
    monkeypatch.setattr(psij_backend, "_load_psij", lambda: _FAKE_PSIJ)
    return psij_backend.PsijExecutor(
        "slurm", job_executor=backend, poll_interval=0.001
    )


def _psij_run(directory, backend_id=None):
    return RunRecord(
        id="psij-run", task="DOS", state=RunState.RUNNING,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            ExecutorEndpoint(executor="psij", psij_backend="slurm"),
            WorkingDirectory(str(directory)), backend_id,
        ),
        files=FileManifest(
            inputs=(FileEntry(local_path=str(directory / "run.inp")),),
            stdout=FileEntry(local_path=str(directory / "run.out")),
            stderr=FileEntry(local_path=str(directory / "run.err")),
        ),
    )


def test_psij_runner_submits_and_streams_output(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    input_file = open(tmp_path / "run.inp", "wb+")
    input_file.write(b"input")
    output_file = open(tmp_path / "run.out", "wb")
    spec = ExecutionSpec(
        argv=("program", "arg"),
        input_file=input_file,
        output_file=output_file,
        directory=str(tmp_path),
        stdin=input_file,
    )
    runner = executor.create_run(spec)
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))

    runner.submit()
    backend_id = runner.backend_id
    result = runner.wait()

    assert backend.submitted.spec.executable == "program"
    assert backend.submitted.spec.arguments == ["arg"]
    assert backend.submitted.spec.stdin_path == str(tmp_path / "run.inp")
    assert backend.submitted.spec.stdout_path == str(tmp_path / "run.out")
    assert backend.submitted.spec.resources is None
    assert backend.submitted.spec.launcher is None
    assert backend_id == "123"
    assert runner.state == RunState.FINISHED
    assert (tmp_path / "run.out").read_bytes() == b"hello\n"
    assert (tmp_path / "run.out.err").read_bytes() == b"warning\n"
    assert result.return_code == 0
    assert ("out", "hello\n") in seen
    assert ("err", "warning\n") in seen
    input_file.close()


def test_psij_mpi_process_count_becomes_resource_spec(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    input_file = open(tmp_path / "run.inp", "wb+")
    output_file = open(tmp_path / "run.out", "wb")
    spec = ExecutionSpec(
        argv=("program",),
        mpi_argv=("programMPI", str(tmp_path / "run.inp")),
        launch=LaunchSpec.from_mpi(16),
        input_file=input_file,
        output_file=output_file,
        directory=str(tmp_path),
        stdin=input_file,
    )

    runner = executor.create_run(spec)
    job_spec = runner._job.spec

    assert job_spec.executable == "programMPI"
    assert job_spec.arguments == [str(tmp_path / "run.inp")]
    assert job_spec.resources.process_count == 16
    assert job_spec.launcher == "srun"
    assert job_spec.stdin_path is None
    input_file.close()
    output_file.close()


def test_psij_auto_uses_mpi_without_fixed_process_count(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    input_file = open(tmp_path / "run.inp", "wb+")
    output_file = open(tmp_path / "run.out", "wb")
    spec = ExecutionSpec(
        argv=("program",),
        mpi_argv=("programMPI", str(tmp_path / "run.inp")),
        launch=LaunchSpec.from_mpi("auto"),
        input_file=input_file,
        output_file=output_file,
        directory=str(tmp_path),
        stdin=input_file,
    )

    runner = executor.create_run(spec)
    job_spec = runner._job.spec

    assert job_spec.executable == "programMPI"
    assert job_spec.resources is None
    assert job_spec.launcher == "srun"
    input_file.close()
    output_file.close()


def test_psij_project_run_exposes_recovery_metadata(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    input_file = open(tmp_path / "run.inp", "wb+")
    output_file = open(tmp_path / "run.out", "wb")
    spec = ExecutionSpec(
        argv=("program",),
        input_file=input_file,
        output_file=output_file,
        directory=str(tmp_path),
        stdin=input_file,
    )

    runner = executor.create_project_run(
        spec, run_id="persistent-run"
    )

    assert runner._job.spec.environment == {
        "ASE2SPRKKR_RUN_ID": "persistent-run"
    }
    assert runner._job.spec.attributes.custom_attributes == {
        "slurm.comment": "ase2sprkkr:persistent-run"
    }
    input_file.close()
    output_file.close()



def test_psij_record_status_uses_native_status(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    expected = psij_backend.RunStatus(
        RunState.FINISHED, backend_id="456", return_code=0,
    )
    seen = []
    def status_native(native_id):
        seen.append(native_id)
        return expected

    monkeypatch.setattr(executor, "status_native", status_native)
    assert executor.status(_psij_run(tmp_path, "456")) == (expected,)
    assert seen == ["456"]


def test_psij_submission_without_native_identity_is_not_started(monkeypatch, tmp_path):
    from ase2sprkkr.execution.base import FailurePhase, SubmissionNotStarted

    class FailingBackend(_Backend):
        def submit(self, job):
            raise OSError("scheduler unavailable")

    executor = _executor(monkeypatch, FailingBackend())
    runner = executor.create_project_run(
        ExecutionSpec(
            argv=("program",),
            directory=str(tmp_path),
            output_file=str(tmp_path / "run.out"),
        ),
        run_id="run-123",
    )
    with pytest.raises(SubmissionNotStarted, match="scheduler unavailable"):
        runner.submit()
    assert runner.failure_phase == FailurePhase.SUBMISSION

def test_psij_executor_attaches_to_native_job(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    output = tmp_path / "run.out"
    error = tmp_path / "run.err"
    output.write_bytes(b"old output\n")
    error.write_bytes(b"")
    run = _psij_run(tmp_path, "456")
    runner = executor.attach(run)
    result = runner.wait()

    assert backend.attached == "456"
    assert runner.backend_id == "456"
    assert runner.state == RunState.FINISHED
    assert (tmp_path / "run.out").read_bytes() == b"old output\n"


def test_psij_recovers_live_slurm_submission_by_run_id(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    calls = []

    def scheduler_output(argv):
        calls.append(tuple(argv))
        if argv[0] == "squeue":
            return f"123|ase2sprkkr:run-1|RUNNING|{tmp_path}\n"
        if argv[:4] == ("scontrol", "--oneliner", "show", "job"):
            return (
                f"JobId=123 WorkDir={tmp_path} "
                f"StdOut={tmp_path / 'run.out'} StdErr={tmp_path / 'run.err'}\n"
            )
        raise AssertionError(argv)

    monkeypatch.setattr(psij_backend, "_scheduler_output", scheduler_output)

    status = executor.recover_submission(
        "run-1", _psij_run(tmp_path)
    )

    assert status.state == RunState.RUNNING
    assert status.backend_id == "123"
    assert status.execution.directory.local == str(tmp_path)
    assert status.files.stdout.remote_path == str(tmp_path / "run.out")
    assert status.files.stderr.remote_path == str(tmp_path / "run.err")
    assert calls[0][0] == "squeue"


def test_psij_recovers_completed_slurm_submission_from_accounting(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())

    def scheduler_output(argv):
        if argv[0] == "squeue":
            return ""
        if argv[:3] == ("scontrol", "show", "config"):
            return "AccountingStoreFlags = job_comment,job_env\n"
        if argv[0] == "sacct":
            return (
                f"456|ase2sprkkr:run-2|COMPLETED|{tmp_path}|"
                f"{tmp_path / 'run.out'}|{tmp_path / 'run.err'}|0:0|\n"
            )
        raise AssertionError(argv)

    monkeypatch.setattr(psij_backend, "_scheduler_output", scheduler_output)

    status = executor.recover_submission(
        "run-2",
        _psij_run(tmp_path),
        submitted_at="2026-10-13T08:00:00+02:00",
    )

    assert status.state == RunState.FINISHED
    assert status.return_code == 0
    assert status.backend_id == "456"
    assert status.files.stdout.remote_path == str(tmp_path / "run.out")


def test_psij_accounting_confirms_missing_submission(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())

    def scheduler_output(argv):
        if argv[0] == "squeue":
            return ""
        if argv[:3] == ("scontrol", "show", "config"):
            return "AccountingStoreFlags = job_comment\n"
        if argv[0] == "sacct":
            return ""
        raise AssertionError(argv)

    monkeypatch.setattr(psij_backend, "_scheduler_output", scheduler_output)

    assert executor.recover_submission(
        "missing",
        _psij_run(tmp_path),
        submitted_at="2026-10-13T08:00:00+02:00",
    ) is None


def test_psij_missing_accounting_history_keeps_submission_uncertain(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())

    def scheduler_output(argv):
        if argv[0] == "squeue":
            return ""
        if argv[:3] == ("scontrol", "show", "config"):
            return "AccountingStoreFlags = job_env\n"
        raise AssertionError(argv)

    monkeypatch.setattr(psij_backend, "_scheduler_output", scheduler_output)

    with pytest.raises(NotImplementedError, match="comments are unavailable"):
        executor.recover_submission(
            "uncertain",
            _psij_run(tmp_path),
            submitted_at="2026-10-13T08:00:00+02:00",
        )


def test_psij_recovery_rejects_ambiguous_live_jobs(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    monkeypatch.setattr(
        psij_backend,
        "_scheduler_output",
        lambda argv: (
            "123|ase2sprkkr:duplicate|RUNNING|/one\n"
            "124|ase2sprkkr:duplicate|PENDING|/two\n"
        ),
    )

    with pytest.raises(RuntimeError, match="Multiple live Slurm jobs"):
        executor.recover_submission(
            "duplicate", _psij_run(tmp_path)
        )


def test_psij_submit_is_serialized_across_threads(monkeypatch, tmp_path):
    class BlockingBackend(_Backend):
        def __init__(self):
            super().__init__()
            self.calls = 0
            self.first_entered = threading.Event()
            self.second_entered = threading.Event()
            self.release = threading.Event()

        def submit(self, job):
            self.calls += 1
            (self.first_entered if self.calls == 1 else self.second_entered).set()
            self.release.wait(1)
            return super().submit(job)

    backend = BlockingBackend()
    executor = _executor(monkeypatch, backend)
    input_file = open(tmp_path / "run.inp", "wb+")
    input_file.write(b"input")
    output_file = open(tmp_path / "run.out", "wb")
    spec = ExecutionSpec(
        argv=("program",),
        input_file=input_file,
        output_file=output_file,
        directory=str(tmp_path),
        stdin=input_file,
    )
    runner = executor.create_run(spec)
    errors = []

    def submit():
        try:
            runner.submit()
        except RuntimeError as exc:
            errors.append(exc)

    first = threading.Thread(target=submit)
    second = threading.Thread(target=submit)
    first.start()
    assert backend.first_entered.wait(1)
    second.start()
    assert not backend.second_entered.wait(0.2)
    backend.release.set()
    first.join(2)
    second.join(2)

    assert not first.is_alive()
    assert not second.is_alive()
    assert backend.calls == 1
    assert len(errors) == 1
    runner.wait()
    input_file.close()
