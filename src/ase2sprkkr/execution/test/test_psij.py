import asyncio
import subprocess
import threading
from dataclasses import replace
from datetime import timedelta
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
        self._status = _Status()
        self._status_callback = None
        self._done = threading.Event()

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, status):
        self._status = status
        self.emit(status)

    def emit(self, status):
        if self._status_callback is not None:
            self._status_callback(self, status)

    def set_job_status_callback(self, callback):
        self._status_callback = callback

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
            if job.spec.stdout_path is not None:
                with open(job.spec.stdout_path, "wb") as out:
                    out.write(b"hello\n")
            if job.spec.stderr_path is not None:
                mode = "ab" if job.spec.stderr_path == job.spec.stdout_path else "wb"
                with open(job.spec.stderr_path, mode) as err:
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
        error_file=str(tmp_path / "run.err"),
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
    assert backend.submitted.spec.stderr_path == str(tmp_path / "run.err")
    assert (tmp_path / "run.err").read_bytes() == b"warning\n"
    assert result.return_code == 0
    assert ("out", "hello\n") in seen
    assert ("err", "warning\n") in seen
    input_file.close()



def test_psij_without_output_files_does_not_create_any(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    spec = ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=None, output_file=None, error_file=None,
    )
    seen = []
    runner = executor.create_run(spec)
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))
    runner.submit()
    result = runner.wait()

    assert backend.submitted.spec.stdout_path is None
    assert backend.submitted.spec.stderr_path is None
    assert result.output_file is None
    assert result.error_file is None
    assert result.return_code == 0
    assert not seen
    assert list(tmp_path.iterdir()) == []


def test_psij_stderr_only_does_not_create_stdout(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=None, output_file=None,
        error_file=str(tmp_path / "run.err"),
    ))
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))
    runner.submit()
    result = runner.wait()

    assert backend.submitted.spec.stdout_path is None
    assert backend.submitted.spec.stderr_path == str(tmp_path / "run.err")
    assert result.output_file is None
    assert (tmp_path / "run.err").read_bytes() == b"warning\n"
    assert ("err", "warning\n") in seen
    assert not (tmp_path / "run.out").exists()


def test_psij_merged_output_is_followed_once(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    output = str(tmp_path / "merged.out")
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=None, output_file=output, error_file=output,
    ))
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))
    runner.submit()
    runner.wait()

    assert (tmp_path / "merged.out").read_bytes() == b"hello\nwarning\n"
    assert "".join(value for kind, value in seen if kind == "out") == "hello\nwarning\n"
    assert not any(kind == "err" for kind, value in seen)



@pytest.mark.parametrize(
    "native_name, return_code, attached, expected",
    [
        ("COMPLETED", 7, False, RunState.FAILED),
        ("COMPLETED", 0, False, RunState.FINISHED),
        ("COMPLETED", None, False, RunState.FINISHED),
        ("FAILED", 0, False, RunState.FAILED),
        ("CANCELED", 1, False, RunState.CANCELLED),
        ("NEW", None, False, RunState.NEW),
        ("NEW", None, True, RunState.UNKNOWN),
        ("ACTIVE", None, True, RunState.RUNNING),
        (None, None, True, RunState.UNKNOWN),
    ],
)
def test_native_psij_state_normalization(
    native_name, return_code, attached, expected,
):
    native = None if native_name is None else _Status(native_name, return_code)
    assert psij_backend._psij_state(native, attached=attached) == expected


def test_psij_runner_prefers_observed_return_code_over_cached_job_status(
    monkeypatch, tmp_path,
):
    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=None, output_file=None,
    ))
    runner._job.status = _Status("COMPLETED", 0)
    runner._return_code = 7
    assert runner.state == RunState.FAILED


def test_psij_non_slurm_status_uses_shared_normalization(monkeypatch):
    executor = _executor(monkeypatch, _Backend())
    executor.backend = "pbs"
    executor.url = "remote-pbs"

    class Job:
        def __init__(self):
            self.status = _Status("COMPLETED", 9)

    class Backend:
        def attach(self, job, native_id):
            assert native_id == "pbs-123"

    executor._psij = SimpleNamespace(Job=Job)
    executor._job_executor = Backend()
    result = executor.status_native("pbs-123")
    assert result.state == RunState.FAILED
    assert result.return_code == 9
    assert result.failure_phase == psij_backend.FailurePhase.EXECUTION


def test_psij_nonzero_completed_job_is_failed(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path), input_file=None,
        output_file=str(tmp_path / "run.out"),
    ))
    runner._job.status = _Status("COMPLETED", 7)
    assert runner.state == RunState.FAILED


def test_psij_slurm_nonzero_accounted_job_is_failed(monkeypatch):
    executor = _executor(monkeypatch, _Backend())
    monkeypatch.setattr(psij_backend, "_slurm_details", lambda native_id: {})
    monkeypatch.setattr(
        psij_backend, "_scheduler_output",
        lambda argv: "123|COMPLETED|7:0|\\n".replace("\\n", "\n"),
    )
    status = executor.status_native("123")
    assert status.state == RunState.FAILED
    assert status.return_code == 7
    assert status.failure_phase == psij_backend.FailurePhase.EXECUTION


def test_psij_observation_accepts_non_utf8_output(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), input_file=None, directory=str(tmp_path),
        output_file=str(tmp_path / "run.out"),
    ))
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))
    runner._output_observers.emit_bytes(b"raw:\xff", "out")
    assert seen == [("out", "raw:\ufffd")]
    assert runner.exception is None


def test_psij_project_manifest_marks_concrete_and_wildcard_files(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    (tmp_path / "data-1.bin").write_bytes(b"cache")
    files = FileManifest(
        stdout=FileEntry(local_path=str(tmp_path / "run.out"), availability_state="pending"),
        outputs=(FileEntry(local_path=str(tmp_path / "data-*.bin"), availability_state="pending"),),
    )
    runner = executor.create_project_run(ExecutionSpec(
        argv=("program",), input_file=None, directory=str(tmp_path),
        output_file=str(tmp_path / "run.out"),
    ), run_id="manifest")
    runner.set_manifest(files)
    runner.submit()
    runner.wait()
    assert runner.files.stdout.availability_state == "available"
    assert runner.files.outputs[0].availability_state == "available"
    matches = [entry for entry in runner.files.outputs if entry.source_pattern is not None]
    assert len(matches) == 1
    assert matches[0].availability_state == "available"


def test_psij_missing_optional_output_is_skipped(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    files = FileManifest(outputs=(
        FileEntry(local_path=str(tmp_path / "optional.dat"), availability_state="pending"),
    ))
    runner = executor.create_project_run(ExecutionSpec(
        argv=("program",), input_file=None, directory=str(tmp_path),
        output_file=None,
    ), run_id="optional")
    runner.set_manifest(files)
    runner.submit()
    runner.wait()
    assert runner.files.outputs[0].availability_state == "skipped"


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
    assert job_spec.stdin_path == str(tmp_path / "run.inp")
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

def test_psij_status_is_nonblocking_and_explicit_wait_synchronizes(monkeypatch):
    executor = _executor(monkeypatch, _Backend())
    executor.backend = "pbs"
    executor.url = "remote-pbs"
    observed = []

    class Job:
        def __init__(self):
            self.status = _Status("NEW")
        def wait(self, *, timeout, target_states):
            observed.append((timeout, target_states))
            self.status = _Status("ACTIVE")
            return self.status

    class Backend:
        def attach(self, job, native_id):
            observed.append(native_id)

    executor._job_executor = Backend()
    executor._psij = SimpleNamespace(
        Job=Job, JobState=SimpleNamespace(QUEUED="QUEUED"),
    )
    status = executor.status_native("native-42")
    assert status.state == RunState.UNKNOWN
    assert observed == ["native-42"]

    status = executor.wait_status_native("native-42", timeout=2.5)
    assert status.state == RunState.RUNNING
    assert observed == ["native-42", (timedelta(seconds=2.5), ["QUEUED"])]
    assert executor.status_native("native-42").state == RunState.RUNNING
    assert len(observed) == 2


def test_psij_explicit_wait_timeout_keeps_cached_new_unknown(monkeypatch):
    executor = _executor(monkeypatch, _Backend())
    executor.backend = "pbs"
    executor.url = "remote-pbs"
    observed = []

    class Job:
        def __init__(self):
            self.status = _Status("NEW")

        def wait(self, *, timeout, target_states):
            observed.append((timeout, target_states))
            return None

    class Backend:
        def attach(self, job, native_id):
            observed.append(native_id)

    executor._job_executor = Backend()
    executor._psij = SimpleNamespace(
        Job=Job, JobState=SimpleNamespace(QUEUED="QUEUED"),
    )
    first = executor.status_native("native-42")
    assert first.state == RunState.UNKNOWN
    assert observed == ["native-42"]
    status = executor.wait_status_native("native-42", timeout=0.05)
    assert status.state == RunState.UNKNOWN
    assert status.backend_id == "native-42"
    assert executor.status_native("native-42").state == RunState.UNKNOWN
    assert observed == ["native-42", (timedelta(seconds=0.05), ["QUEUED"])]
    assert executor.wait_status_native("native-42", timeout=0).state == RunState.UNKNOWN
    assert len(observed) == 2
    with pytest.raises(ValueError, match="negative"):
        executor.wait_status_native("native-42", timeout=-1)


def test_psij_local_slurm_explicit_wait_uses_native_job_wait(monkeypatch):
    executor = _executor(monkeypatch, _Backend())
    observed = []

    class Job:
        def wait(self, **kwargs):
            observed.append(("wait", kwargs))

    class Backend:
        def attach(self, job, native_id):
            observed.append(("attach", native_id))

    executor._psij = SimpleNamespace(
        Job=Job, JobState=SimpleNamespace(QUEUED="QUEUED")
    )
    executor._job_executor = Backend()
    # First immediate probe is UNKNOWN; a second UNKNOWN triggers a native
    # wait and the final probe observes the job's newly reported state.
    states = [RunState.UNKNOWN, RunState.UNKNOWN, RunState.RUNNING]
    monkeypatch.setattr(executor, "status_native",
                        lambda native_id: psij_backend.RunStatus(
                            states.pop(0), backend_id=native_id
                        ))
    assert executor.wait_status_native("slurm-123", timeout=0).state == RunState.UNKNOWN
    assert executor.wait_status_native("slurm-123", timeout=0.1).state == RunState.RUNNING
    assert observed[0] == ("attach", "slurm-123")
    assert observed[1][0] == "wait"



def test_psij_status_and_attach_check_complete_endpoint(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    executor.url = "ssh://cluster-a"
    run = _psij_run(tmp_path, "456")
    foreign = replace(
        run,
        execution=replace(
            run.execution,
            endpoint=ExecutorEndpoint(
                executor="psij", psij_backend="slurm", psij_url="ssh://cluster-b"
            ),
        ),
    )
    with pytest.raises(ValueError, match="PSI/J endpoint"):
        executor.status(foreign)
    with pytest.raises(ValueError, match="PSI/J endpoint"):
        executor.attach(foreign, start=False)

    matching = replace(
        run,
        execution=replace(
            run.execution,
            endpoint=ExecutorEndpoint(
                executor="psij", psij_backend="slurm", psij_url=executor.url
            ),
        ),
    )
    seen = []
    state = psij_backend.RunStatus(RunState.RUNNING, backend_id="456")

    def status_native(native_id):
        seen.append(native_id)
        return state

    monkeypatch.setattr(executor, "status_native", status_native)
    assert executor.status(matching) == (state,)
    attached = executor.attach(matching, start=False)
    assert attached.backend_id == "456"
    assert seen == ["456"]

    # A completed native result must also propagate its return code.
    state = psij_backend.RunStatus(
        RunState.FINISHED, backend_id="456", return_code=0,
    )
    assert executor.status(matching) == (state,)
    assert seen == ["456", "456"]


def test_psij_uncertain_submit_without_native_id_stays_recoverable(monkeypatch, tmp_path):
    from ase2sprkkr.execution.base import SubmissionNotStarted

    class FailingBackend(_Backend):
        def submit(self, job):
            raise OSError("scheduler reply lost")

    executor = _executor(monkeypatch, FailingBackend())
    runner = executor.create_project_run(
        ExecutionSpec(
            argv=("program",), directory=str(tmp_path), input_file=None,
            output_file=str(tmp_path / "run.out"),
        ), run_id="run-123",
    )
    with pytest.raises(OSError, match="scheduler reply lost") as caught:
        runner.submit()
    assert not isinstance(caught.value, SubmissionNotStarted)
    assert runner.backend_id is None
    assert runner.failure_phase is None


def test_psij_preparation_failure_is_definitively_not_submitted(monkeypatch, tmp_path):
    from ase2sprkkr.execution.base import FailurePhase, SubmissionNotStarted

    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_project_run(
        ExecutionSpec(
            argv=("program",), directory=str(tmp_path), input_file=None,
            output_file=str(tmp_path / "run.out"),
        ), run_id="run-123",
    )
    def fail_preparation():
        raise OSError("cannot prepare output")

    monkeypatch.setattr(runner, "_prepare_io", fail_preparation)
    with pytest.raises(SubmissionNotStarted, match="cannot prepare output"):
        runner.submit()
    assert executor._job_executor.submitted is None
    assert runner.failure_phase == FailurePhase.SUBMISSION

def test_psij_native_attachment_and_status_share_one_job(monkeypatch):
    class Backend:
        def __init__(self):
            self.attaches = []

        def attach(self, job, native_id):
            self.attaches.append((job, native_id))
            job.native_id = native_id
            job.status = _Status("ACTIVE")

    backend = Backend()
    executor = _executor(monkeypatch, backend)
    executor.backend = "pbs"
    executor.url = "remote-pbs"

    first = executor.attach_native_job("native-42")
    second = executor.attach_native_job("native-42")
    assert first is second
    assert executor.status_native("native-42").state == RunState.RUNNING
    assert executor._status_job("native-42") is first
    assert backend.attaches == [(first, "native-42")]


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



def test_psij_attach_to_stderr_only(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    from dataclasses import replace
    run = _psij_run(tmp_path, "456")
    run = replace(run, files=replace(run.files, stdout=None))
    (tmp_path / "run.err").write_bytes(b"existing error\n")

    runner = executor.attach(run)
    assert runner.wait().output_file is None
    assert runner.result.error_file == str(tmp_path / "run.err")


def test_psij_attach_without_recorded_streams(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    from dataclasses import replace
    run = _psij_run(tmp_path, "456")
    run = replace(run, files=replace(run.files, stdout=None, stderr=None))
    attached = executor.attach(run, start=False)
    assert attached.spec.output_file is None
    assert attached.spec.error_file is None
    assert attached.files == run.files


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


@pytest.mark.parametrize(
    "exit_code, accounting_flags, use_workdir_paths, expected_state",
    [
        (0, "job_comment,job_env", True, RunState.FINISHED),
        (7, "job_comment", False, RunState.FAILED),
    ],
)
def test_psij_recovers_completed_slurm_submission_from_accounting(
    monkeypatch, tmp_path, exit_code, accounting_flags, use_workdir_paths,
    expected_state,
):
    executor = _executor(monkeypatch, _Backend())

    def scheduler_output(argv):
        if argv[0] == "squeue":
            return ""
        if argv[:3] == ("scontrol", "show", "config"):
            return f"AccountingStoreFlags = {accounting_flags}\n"
        if argv[0] == "sacct":
            directory, stdout, stderr = (
                (str(tmp_path), str(tmp_path / "run.out"), str(tmp_path / "run.err"))
                if use_workdir_paths else ("/tmp", "out", "err")
            )
            return (
                f"456|ase2sprkkr:run-2|COMPLETED|{directory}|"
                f"{stdout}|{stderr}|{exit_code}:0|\n"
            )
        raise AssertionError(argv)

    monkeypatch.setattr(psij_backend, "_scheduler_output", scheduler_output)
    status = executor.recover_submission(
        "run-2", _psij_run(tmp_path),
        submitted_at="2026-10-13T08:00:00+02:00",
    )
    assert status.state == expected_state
    assert status.return_code == exit_code
    assert status.backend_id == "456"
    if use_workdir_paths:
        assert status.files.stdout.remote_path == str(tmp_path / "run.out")


@pytest.mark.parametrize(
    "accounting_flags, accounting_available, run_id, expected_error",
    [
        ("job_comment", True, "missing", "still uncertain"),
        ("job_env", False, "uncertain", "comments are unavailable"),
    ],
)
def test_psij_unavailable_accounting_keeps_submission_uncertain(
    monkeypatch, tmp_path, accounting_flags, accounting_available,
    run_id, expected_error,
):
    executor = _executor(monkeypatch, _Backend())

    def scheduler_output(argv):
        if argv[0] == "squeue":
            return ""
        if argv[:3] == ("scontrol", "show", "config"):
            return f"AccountingStoreFlags = {accounting_flags}\n"
        if argv[0] == "sacct" and accounting_available:
            return ""
        raise AssertionError(argv)

    monkeypatch.setattr(psij_backend, "_scheduler_output", scheduler_output)
    with pytest.raises(NotImplementedError, match=expected_error):
        executor.recover_submission(
            run_id, _psij_run(tmp_path),
            submitted_at="2026-10-13T08:00:00+02:00",
        )


def test_psij_remote_url_does_not_search_local_scheduler(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    executor.url = "ssh://other-cluster"
    monkeypatch.setattr(
        psij_backend, "_scheduler_output",
        lambda argv: pytest.fail("must not search local Slurm"),
    )
    with pytest.raises(NotImplementedError, match="local Slurm"):
        executor.recover_native("run-remote", str(tmp_path))


def test_psij_recovery_handles_job_disappearing_between_queries(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())

    def scheduler_output(argv):
        if argv[0] == "squeue":
            return f"123|ase2sprkkr:run-race|RUNNING|{tmp_path}\n"
        if argv[:4] == ("scontrol", "--oneliner", "show", "job"):
            raise subprocess.CalledProcessError(1, argv)
        if argv[:3] == ("scontrol", "show", "config"):
            return "AccountingStoreFlags = job_comment\n"
        if argv[0] == "sacct":
            return (
                f"123|ase2sprkkr:run-race|COMPLETED|{tmp_path}|"
                f"{tmp_path / 'run.out'}|{tmp_path / 'run.err'}|0:0|\n"
            )
        raise AssertionError(argv)

    monkeypatch.setattr(psij_backend, "_scheduler_output", scheduler_output)
    facts = executor.recover_native(
        "run-race", str(tmp_path),
        submitted_at="2026-10-13T08:00:00+02:00",
    )
    assert facts["native_id"] == "123"
    assert facts["state"] == RunState.FINISHED
    assert facts["return_code"] == 0


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

def test_psij_mpi_respects_explicit_stdin_policy(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    spec = ExecutionSpec(
        argv=("program",), mpi_argv=("parallel",),
        launch=LaunchSpec.from_mpi(2),
        input_file=None, output_file=None, directory=str(tmp_path),
        stdin=str(tmp_path / "stdin.dat"),
        mpi_uses_stdin=False,
    )
    runner = executor.create_run(spec)
    assert runner._job.spec.stdin_path is None
    assert replace(spec, mpi_uses_stdin=True).mpi_uses_stdin


@pytest.mark.parametrize("merged", (False, True))
def test_psij_preparation_closes_each_output_handle_once(monkeypatch, tmp_path, merged):
    class NamedOutput:
        def __init__(self, name):
            self.name = str(name)
            self.close_calls = 0
            self.flush_calls = 0

        def flush(self):
            self.flush_calls += 1

        def close(self):
            self.close_calls += 1
            if self.close_calls != 1:
                raise AssertionError("Output handle was closed more than once")

    executor = _executor(monkeypatch, _Backend())
    output = NamedOutput(tmp_path / "run.out")
    error = output if merged else NamedOutput(tmp_path / "run.err")
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=None, output_file=output, error_file=error,
    ))
    runner._prepare_io()
    assert output.close_calls == 1
    assert error.close_calls == 1
    assert output.flush_calls >= 1
    assert error.flush_calls >= 1
    assert (tmp_path / "run.out").exists()
    if not merged:
        assert (tmp_path / "run.err").exists()

def test_psij_active_job_stays_active_without_blocking_wait(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path), input_file=None,
        output_file=None,
    ))
    runner._job.status = _Status("ACTIVE")
    assert runner.state is RunState.RUNNING

@pytest.mark.parametrize(
    "native_state, exit_code, expected",
    (("COMPLETED", 0, RunState.FINISHED),
     ("COMPLETED", 7, RunState.FAILED),
     ("FAILED", 9, RunState.FAILED)),
)
def test_psij_wait_terminal_status_overrides_stale_job_status(
    monkeypatch, tmp_path, native_state, exit_code, expected,
):
    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=None, output_file=None,
    ))
    runner._job.status = _Status("ACTIVE")
    runner._job.status = _Status(native_state, exit_code)
    asyncio.run(runner._monitor())
    assert runner.state is expected
    assert runner.return_code == exit_code
    assert runner.result.return_code == exit_code


def test_psij_attached_new_status_does_not_resubmit(monkeypatch, tmp_path):
    from ase2sprkkr.execution.manager import ExecutionManager, RunStore

    class PendingAttachBackend(_Backend):
        def submit(self, job):
            raise AssertionError("Attached job must not be resubmitted")

        def attach(self, job, native_id):
            self.attached = native_id
            job.native_id = native_id
            job.status = _Status("NEW")

    backend = PendingAttachBackend()
    executor = _executor(monkeypatch, backend)
    output = tmp_path / "run.out"
    output.write_bytes(b"existing output")
    record = _psij_run(tmp_path, "456")
    manager = ExecutionManager(RunStore(tmp_path / "history"))
    manager.store.save(record)
    run = manager.attach(
        record, executor=executor, process_result=lambda outcome: outcome,
    )
    assert run._runner._runner.state is RunState.UNKNOWN
    assert output.read_bytes() == b"existing output"
    job = run._runner._runner._job
    job.status = _Status("COMPLETED", 0)
    run.wait()
    assert output.read_bytes() == b"existing output"


def test_psij_output_permission_error_is_unknown_not_execution_failure(
    monkeypatch, tmp_path,
):
    from ase2sprkkr.execution.growing_file_tail import GrowingFileTail

    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path), input_file=None,
        output_file=str(tmp_path / "run.out"),
    ))
    def denied(self, *args, **kwargs):
        raise PermissionError("cannot read scheduler output")

    monkeypatch.setattr(GrowingFileTail, "read_available", denied)
    asyncio.run(runner._monitor())
    assert runner.state is RunState.UNKNOWN
    assert runner.exception is None
    assert runner.failure_phase is None
    assert "cannot read" in runner.failure_message


@pytest.mark.parametrize("native_state", ["COMPLETED", "FAILED", "CANCELED"])
def test_psij_does_not_invent_missing_exit_code(monkeypatch, tmp_path, native_state):
    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=None, output_file=None,
    ))
    runner._job.status = _Status(native_state)
    asyncio.run(runner._monitor())
    assert runner.return_code is None
    assert runner.result.return_code is None
    assert runner.state is {
        "COMPLETED": RunState.FINISHED,
        "FAILED": RunState.FAILED,
        "CANCELED": RunState.CANCELLED,
    }[native_state]



def test_psij_terminal_wait_without_exit_does_not_use_stale_status(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path), input_file=None, output_file=None,
    ))
    runner._job.status = _Status("FAILED", 7)
    runner._job.status = _Status("COMPLETED")
    asyncio.run(runner._monitor())
    assert runner.state is RunState.FINISHED
    assert runner.return_code is None
    assert runner.failure_phase is None


def test_psij_confirmed_completion_is_not_undone_by_unknown_refresh(
    monkeypatch, tmp_path,
):
    from ase2sprkkr.execution.manager import ExecutionManager, RunStore
    from ase2sprkkr.execution.base import RunStatus
    executor = _executor(monkeypatch, _Backend())
    manager = ExecutionManager(RunStore(tmp_path / "history"))
    record = manager.new_record(task="DOS", directory=tmp_path, endpoint=None)
    manager.store.save(record)
    record = manager.transition(record, state=RunState.RUNNING)
    record = manager.update(record, RunStatus(RunState.FINISHED))
    assert record.state is RunState.POSTPROCESSING
    older = manager.update(record, RunStatus(RunState.UNKNOWN, message="disconnected"))
    assert older.state is RunState.POSTPROCESSING
    assert older.message is None


def test_psij_callback_finishes_monitor_and_unsubscribes(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=None, output_file=None,
    ))
    runner._job.status = _Status("ACTIVE")
    def forbidden_wait():
        raise AssertionError("job.wait must not be used")
    runner._job.wait = forbidden_wait

    async def check():
        monitor = asyncio.create_task(runner._monitor())
        await asyncio.sleep(0)
        assert runner._job._status_callback is not None
        runner._job.status = _Status("COMPLETED", 0)
        await asyncio.wait_for(monitor, 2)
    asyncio.run(check())
    assert runner.state == RunState.FINISHED
    assert runner._job._status_callback is psij_backend._ignore_job_status


def test_psij_subscriptions_share_job_and_detach_independently(monkeypatch):
    executor = _executor(monkeypatch, _Backend())
    job = _Job()
    job.native_id = "shared"
    executor._cache_submitted_job(job)
    calls = []
    a = executor._observe_job(job, lambda status: calls.append(("a", status)))
    b = executor._observe_job(job, lambda status: calls.append(("b", status)))
    assert executor.attach_native_job("shared") is job
    assert job._status_callback is not None
    a()
    job.status = _Status("ACTIVE")
    assert [tag for tag, _ in calls] == ["b"]
    b()
    assert job._status_callback is psij_backend._ignore_job_status
    assert executor.attach_native_job("shared") is job



def test_psij_simultaneous_attaches_share_callback_and_finish(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    job = _Job()
    job.native_id = "shared"
    job.status = _Status("ACTIVE")
    executor._cache_submitted_job(job)

    run = _psij_run(tmp_path, "shared")
    first = executor.attach(run, start=False)
    second = executor.attach(run, start=False)
    assert first._job is second._job is job
    assert len(executor._jobs) == 1

    async def check():
        monitors = [asyncio.create_task(runner._monitor())
                    for runner in (first, second)]
        await asyncio.sleep(0)
        assert len(executor._observers[job]) == 2
        job.status = _Status("COMPLETED", 0)
        await asyncio.wait_for(asyncio.gather(*monitors), 3)

    asyncio.run(check())
    assert first.state is second.state is RunState.FINISHED
    assert job._status_callback is psij_backend._ignore_job_status
    assert job not in executor._observers
    assert "shared" not in executor._jobs


def test_psij_monitor_cancel_detaches_without_cancelling_job(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=None, output_file=None,
    ))
    runner._job.status = _Status("ACTIVE")

    async def check():
        monitor = asyncio.create_task(runner._monitor())
        await asyncio.sleep(0)
        assert runner._job._status_callback is not None
        monitor.cancel()
        with pytest.raises(asyncio.CancelledError):
            await monitor

    asyncio.run(check())
    assert runner._job._status_callback is psij_backend._ignore_job_status
    assert runner._job.status.state.name == "ACTIVE"



def test_psij_stop_observation_ends_monitor_without_cancel(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path), input_file=None,
        output_file=None,
    ))
    runner._job.status = _Status("ACTIVE")
    async def exercise():
        task = asyncio.create_task(runner._monitor())
        await asyncio.sleep(0)
        runner.stop_observation()
        await asyncio.wait_for(task, 2)
    asyncio.run(exercise())
    assert runner.state is RunState.UNKNOWN
    assert runner._job.status.state.name == "ACTIVE"
    assert runner._job._status_callback is psij_backend._ignore_job_status


def test_psij_status_only_cache_is_bounded(monkeypatch):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    executor.backend = "pbs"
    executor.url = "remote-pbs"
    for number in range(90):
        executor.attach_native_job(str(number))
    assert len(executor._jobs) <= executor._max_idle_jobs


def test_project_observer_detaches_without_stopping_other_attached_client(
    monkeypatch, tmp_path,
):
    from ase2sprkkr.execution.manager import ExecutionManager, RunStore

    class ActiveBackend(_Backend):
        def attach(self, job, native_id):
            self.attached = native_id
            job.native_id = native_id
            job.status = _Status("ACTIVE")

    executor = _executor(monkeypatch, ActiveBackend())
    record = _psij_run(tmp_path, "running-job")
    store = RunStore(tmp_path / "history")
    manager = ExecutionManager(store)
    store.save(record)

    subscribed = threading.Event()
    num_subscribers = 0
    subscriptions_lock = threading.Lock()
    register = executor._observe_job

    def observe(job, listener):
        nonlocal num_subscribers
        detach = register(job, listener)
        with subscriptions_lock:
            num_subscribers += 1
            if num_subscribers == 2:
                subscribed.set()
        return detach

    monkeypatch.setattr(executor, "_observe_job", observe)
    first = manager.attach(
        record, executor=executor, process_result=lambda outcome: outcome,
    )
    second = manager.attach(
        record, executor=executor, process_result=lambda outcome: outcome,
    )
    try:
        assert subscribed.wait(3)
        job = first._runner._runner._job
        assert second._runner._runner._job is job
        first.stop_observation()
        first.wait()
        assert first._runner.state is RunState.UNKNOWN
        assert job.status.state.name == "ACTIVE"
        assert not second._runner._done.is_set()
        assert len(executor._observers[job]) == 1
        # A scheduler status callback from another thread finishes the
        # remaining client without reattaching or resubmitting the job.
        change = threading.Thread(
            target=lambda: setattr(job, "status", _Status("COMPLETED", 0)),
        )
        change.start()
        change.join(3)
        assert not change.is_alive()
        second.wait()
        assert second._runner.state is RunState.FINISHED
        assert job not in executor._observers
    finally:
        first.stop_observation()
        second.stop_observation()


def test_psij_failed_callback_registration_allows_later_attach(monkeypatch):
    executor = _executor(monkeypatch, _Backend())
    job = _Job()
    attempts = []
    original = job.set_job_status_callback

    def flaky_set_callback(callback):
        attempts.append(callback)
        if len(attempts) == 1:
            raise RuntimeError("native callback registration failed")
        original(callback)

    job.set_job_status_callback = flaky_set_callback
    got = []
    with pytest.raises(RuntimeError, match="registration failed"):
        executor._observe_job(job, lambda status: got.append(("failed", status)))
    assert job not in executor._observers

    detach = executor._observe_job(job, lambda status: got.append(("retry", status)))
    job.status = _Status("COMPLETED", 0)
    assert [kind for kind, _ in got] == ["retry"]
    detach()
    assert job not in executor._observers


def test_psij_detach_callback_failure_preserves_original_error(monkeypatch, caplog):
    executor = _executor(monkeypatch, _Backend())
    job = _Job()
    detach = executor._observe_job(job, lambda status: None)
    real_set = job.set_job_status_callback

    def cannot_disable(callback):
        if callback is psij_backend._ignore_job_status:
            raise RuntimeError("native callback cannot be removed")
        return real_set(callback)

    monkeypatch.setattr(job, "set_job_status_callback", cannot_disable)

    with pytest.raises(ValueError, match="primary monitor error"):
        try:
            raise ValueError("primary monitor error")
        finally:
            detach()
    assert job not in executor._observers
    assert "native callback cannot be removed" in caplog.text


def test_psij_detach_callback_failure_reports_problem_and_does_not_leak_listener(
    monkeypatch, caplog,
):
    executor = _executor(monkeypatch, _Backend())
    job = _Job()
    received = []
    detach = executor._observe_job(job, lambda status: received.append(status))
    original_set = job.set_job_status_callback

    def cannot_disable(callback):
        if callback is psij_backend._ignore_job_status:
            raise RuntimeError("native unsubscribe unavailable")
        return original_set(callback)

    monkeypatch.setattr(job, "set_job_status_callback", cannot_disable)
    detach()
    assert job not in executor._observers
    assert "native unsubscribe unavailable" in caplog.text

    # The native callback may remain installed but must not notify observers
    # that are no longer attached.
    job.status = _Status("COMPLETED", 0)
    assert received == []


def test_psij_preparation_reuses_declared_outputs_for_check_and_cleanup(
    monkeypatch, tmp_path,
):
    """Streams and declared wildcard outputs are cleaned, inputs/caches kept."""
    input_file = tmp_path / "run.inp"
    cache_file = tmp_path / "restart.pot"
    output_file = tmp_path / "run.out"
    error_file = tmp_path / "run.err"
    extra_output = tmp_path / "result-1.dat"
    for path in (input_file, cache_file, output_file, error_file, extra_output):
        path.write_bytes(b"previous contents")

    executor = _executor(monkeypatch, _Backend())
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(tmp_path),
        input_file=str(input_file), output_file=str(output_file),
        error_file=str(error_file),
        output_files=("result-*.dat", "run.out"),
        cache_files=("restart.pot",),
    ))
    runner._prepare_io()

    assert input_file.read_bytes() == b"previous contents"
    assert cache_file.read_bytes() == b"previous contents"
    assert output_file.read_bytes() == b""
    assert error_file.read_bytes() == b""
    assert not extra_output.exists()
