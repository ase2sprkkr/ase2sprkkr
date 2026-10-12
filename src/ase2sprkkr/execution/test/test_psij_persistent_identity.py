from types import SimpleNamespace

from ase2sprkkr.execution import ExecutionSpec
from ase2sprkkr.execution import psij as psij_backend


class _JobSpec:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class _JobAttributes:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class _Job:
    def __init__(self, spec=None):
        self.spec = spec
        self.native_id = None
        self.status = SimpleNamespace(state=SimpleNamespace(name="NEW"))


class _JobExecutor:
    @staticmethod
    def get_instance(*args, **kwargs):
        raise AssertionError("test supplies the job executor")


_FAKE_PSIJ = SimpleNamespace(
    Job=_Job,
    JobSpec=_JobSpec,
    JobAttributes=_JobAttributes,
    JobExecutor=_JobExecutor,
)


class _Reader:
    print_output = False

    @staticmethod
    def create_result(output_file, input_file):
        return object()


def _executor(monkeypatch, backend="slurm"):
    monkeypatch.setattr(psij_backend, "_load_psij", lambda: _FAKE_PSIJ)
    return psij_backend.PsijExecutor(backend, job_executor=object())


def _spec(tmp_path):
    output = tmp_path / "run.out"
    return ExecutionSpec(
        argv=("program",),
        input_file=None,
        output_file=str(output),
        directory=str(tmp_path),
    )


def test_persistent_slurm_job_carries_ase2sprkkr_identity(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, "slurm")

    runner = executor.create_persistent_run(
        _Reader(), _spec(tmp_path), run_id="run-123"
    )
    job_spec = runner._job.spec

    assert job_spec.environment == {"ASE2SPRKKR_RUN_ID": "run-123"}
    assert job_spec.attributes.custom_attributes == {
        "slurm.comment": "ase2sprkkr:run-123"
    }


def test_persistent_non_slurm_job_uses_environment_only(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, "pbs")

    runner = executor.create_persistent_run(
        _Reader(), _spec(tmp_path), run_id="run-123"
    )
    job_spec = runner._job.spec

    assert job_spec.environment == {"ASE2SPRKKR_RUN_ID": "run-123"}
    assert job_spec.attributes is None


def test_nonpersistent_psij_job_has_no_recovery_metadata(monkeypatch, tmp_path):
    executor = _executor(monkeypatch, "slurm")

    runner = executor.create_run(_Reader(), _spec(tmp_path))
    job_spec = runner._job.spec

    assert job_spec.environment is None
    assert job_spec.attributes is None
