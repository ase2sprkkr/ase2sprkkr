from ase2sprkkr.execution import (
    ExecutionManager,
    ExecutorEndpoint,
    ExecutorPool,
    LocalExecutor,
    RemoteEndpointExecutor,
    RunReference,
    RunState,
    RunStore,
    SshTransport,
)
from ase2sprkkr.execution.endpoint import endpoint_for_executor
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


def test_profile_uses_execution_endpoint_value():
    profile = SPRKKRProfile(
        executor="ssh",
        host="cluster",
        ssh_options=("-p", "2222"),
        remote_command=("ase2sprkkr", "remote"),
        psij_backend="slurm",
        psij_url="local://slurm",
    )

    assert profile.endpoint == ExecutorEndpoint(
        executor="ssh",
        host="cluster",
        ssh_options=("-p", "2222"),
        remote_command=("ase2sprkkr", "remote"),
        psij_backend="slurm",
        psij_url="local://slurm",
    )


def test_run_store_roundtrips_endpoint(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store, executor_pool=ExecutorPool())
    endpoint = ExecutorEndpoint(executor="ssh", host="cluster", ssh_options=("-v",))
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        dataset="case",
        endpoint=endpoint,
    )
    stored = store.save(record)

    loaded = store.get(stored.id)

    assert loaded.endpoint == endpoint
    assert isinstance(loaded.endpoint, ExecutorEndpoint)


def test_manager_reconstructs_local_executor_after_restart(tmp_path):
    store = RunStore(tmp_path / "runs")
    endpoint = ExecutorEndpoint()
    first = ExecutionManager(store, executor_pool=ExecutorPool())
    record = first.new_record(
        task="DOS",
        directory=tmp_path,
        endpoint=endpoint,
    )
    store.save(record)

    restarted = ExecutionManager(store, executor_pool=ExecutorPool())
    executor = restarted.executor_for(store.get(record.id))

    assert isinstance(executor, LocalExecutor)


def test_attach_reconstructs_persisted_endpoint_and_reader_after_restart(tmp_path):
    class AttachedRunner:
        def __init__(self, reference):
            self.reference = reference
            self.state = RunState.RUNNING
            self.result = None

    class ReconstructedExecutor:
        def __init__(self):
            self.calls = []

        def attach(self, reference, reader, callback=None):
            self.calls.append((reference, reader, callback))
            return AttachedRunner(reference)

    class RecordingPool:
        def __init__(self, executor):
            self.executor = executor
            self.endpoints = []

        def for_endpoint(self, endpoint):
            self.endpoints.append(endpoint)
            return self.executor

    store = RunStore(tmp_path / "runs")
    endpoint = ExecutorEndpoint(executor="ssh", host="cluster")
    reference = RunReference(executor="ssh", native_id="42", computer="cluster")
    first = ExecutionManager(store, executor_pool=ExecutorPool())
    record = first.new_record(
        task="DOS",
        directory=tmp_path,
        endpoint=endpoint,
    )
    store.save(record)
    record = first.transition(record, state=RunState.RUNNING, reference=reference)

    executor = ReconstructedExecutor()
    pool = RecordingPool(executor)
    restarted = ExecutionManager(store, executor_pool=pool)

    attached = restarted.attach(record.id)

    assert pool.endpoints == [endpoint]
    assert attached.executor is executor
    assert executor.calls[0][0] == reference
    assert executor.calls[0][1] is not None
    assert executor.calls[0][1].directory == str(tmp_path)


def test_pool_reconstructs_and_reuses_ssh_endpoint():
    endpoint = ExecutorEndpoint(
        executor="ssh",
        host="cluster",
        ssh_options=("-p", "2222"),
        remote_command=("ase2sprkkr", "remote"),
        psij_backend="slurm",
        psij_url="local://slurm",
    )
    pool = ExecutorPool()

    first = pool.for_endpoint(endpoint)
    second = pool.for_endpoint(endpoint)

    assert first is second
    assert isinstance(first, RemoteEndpointExecutor)
    assert isinstance(first.transport, SshTransport)
    assert first.transport.host == "cluster"
    assert first.transport.options == ("-p", "2222")
    assert first.backend == "psij"
    assert first.backend_options == {"backend": "slurm", "url": "local://slurm"}


def test_default_remote_psij_endpoint_survives_restart():
    live = RemoteEndpointExecutor(
        SshTransport("cluster"),
        backend="psij",
    )

    endpoint = endpoint_for_executor(live)
    restored = ExecutorPool().for_endpoint(endpoint)

    assert endpoint.psij_backend == "slurm"
    assert isinstance(restored, RemoteEndpointExecutor)
    assert restored.backend == "psij"
    assert restored.backend_options == {"backend": "slurm", "url": None}
