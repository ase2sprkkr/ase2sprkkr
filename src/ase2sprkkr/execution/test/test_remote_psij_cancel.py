from types import SimpleNamespace
from threading import Event, Thread

from ..base import RunReference, RunState
from ..remote import _run_directory, _write_state, status_remote
from ..remote_psij import cancel_remote_psij, run_remote_psij_worker


def test_remote_psij_cancel_waits_for_submit_handoff(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    run_id = "psij-cancel-race"
    run_dir = _run_directory(workspace, run_id)
    run_dir.mkdir(parents=True)
    _write_state(
        workspace,
        run_id,
        id=run_id,
        state="queued",
        directory=str(workspace),
        pid=None,
        backend_reference=None,
        finished_at=None,
    )

    request = {
        "id": run_id,
        "executor": "psij",
        "executor_options": {"backend": "slurm", "url": None},
        "directory": str(workspace),
        "argv": ["true"],
        "input_file": None,
        "stdin": None,
        "output_file": str(workspace / "out"),
        "error_file": str(workspace / "err"),
    }

    submit_entered = Event()
    allow_submit = Event()
    cancelled = Event()
    cancel_done = Event()
    worker_errors = []
    killed = []

    class FakeRunner:
        result = SimpleNamespace(return_code=1)
        reference = RunReference(
            executor="psij:slurm",
            native_id="job-1",
            directory=str(workspace),
            output_file=str(workspace / "out"),
            error_file=str(workspace / "err"),
        )

        @property
        def state(self):
            return RunState.CANCELLED if cancelled.is_set() else RunState.QUEUED

        def submit(self):
            submit_entered.set()
            assert allow_submit.wait(2)

        def cancel(self):
            cancelled.set()

        def wait(self):
            assert cancelled.wait(2)
            return self.result

    runner = FakeRunner()

    class FakeExecutor:
        def create_run(self, reader, spec):
            return runner

        def attach(self, reference, reader):
            assert reference.native_id == "job-1"
            return runner

    executor = FakeExecutor()
    monkeypatch.setattr("ase2sprkkr.execution.remote_psij._executor", lambda request: executor)
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij.os.killpg",
        lambda pid, sig: killed.append((pid, sig)),
    )
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij.os.kill",
        lambda pid, sig: (_ for _ in ()).throw(ProcessLookupError()),
    )

    def worker():
        try:
            run_remote_psij_worker(workspace, run_id, request)
        except BaseException as exc:
            worker_errors.append(exc)

    worker_thread = Thread(target=worker)
    worker_thread.start()
    assert submit_entered.wait(2)

    cancel_result = []

    def cancel():
        cancel_result.append(cancel_remote_psij(workspace, run_id, request))
        cancel_done.set()

    cancel_thread = Thread(target=cancel)
    cancel_thread.start()

    assert not cancel_done.wait(0.1)
    allow_submit.set()

    cancel_thread.join(2)
    worker_thread.join(2)

    assert not cancel_thread.is_alive()
    assert not worker_thread.is_alive()
    assert worker_errors == []
    assert cancel_result[0]["state"] == "cancelled"
    assert status_remote(workspace, run_id)["state"] == "cancelled"
    assert killed
