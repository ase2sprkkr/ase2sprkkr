import asyncio
import threading
from types import SimpleNamespace

from ase2sprkkr.execution.base import ExecutionSpec, RunReference, RunState
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


_FAKE_PSIJ = SimpleNamespace(Job=_Job, JobSpec=_JobSpec, JobExecutor=_JobExecutor)


def _executor(monkeypatch, backend):
    monkeypatch.setattr(psij_backend, "_load_psij", lambda: _FAKE_PSIJ)
    return psij_backend.PsijExecutor(
        "slurm", job_executor=backend, poll_interval=0.001
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
    runner = executor.create_run(_Reader(), spec)
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))

    runner.submit()
    reference = runner.reference
    result = runner.wait()

    assert backend.submitted.spec.executable == "program"
    assert backend.submitted.spec.arguments == ["arg"]
    assert backend.submitted.spec.stdin_path == str(tmp_path / "run.inp")
    assert backend.submitted.spec.stdout_path == str(tmp_path / "run.out")
    assert reference.executor == "psij:slurm"
    assert reference.native_id == "123"
    assert runner.state == RunState.FINISHED
    assert result.output == b"hello\n"
    assert result.error_output == b"warning\n"
    assert result.return_code == 0
    assert ("out", "hello\n") in seen
    assert ("err", "warning\n") in seen
    input_file.close()


def test_psij_executor_attaches_to_native_job(monkeypatch, tmp_path):
    backend = _Backend()
    executor = _executor(monkeypatch, backend)
    output = tmp_path / "run.out"
    error = tmp_path / "run.err"
    output.write_bytes(b"old output\n")
    error.write_bytes(b"")
    reference = RunReference(
        executor="psij:slurm",
        native_id="456",
        directory=str(tmp_path),
        input_file=str(tmp_path / "run.inp"),
        output_file=str(output),
        error_file=str(error),
    )

    runner = executor.attach(reference, _Reader())
    result = runner.wait()

    assert backend.attached == "456"
    assert runner.reference.native_id == "456"
    assert runner.state == RunState.FINISHED
    assert result.output == b"old output\n"
