from threading import Event
from kkrunner.backends import ExecutorEndpoint, ExecutorPool, LocalExecutor
from kkrunner.remote import RemoteEndpointExecutor, SshTransport
from kkrunner import RunState, RunStore
from kkrunner.runs.manager import RunManager
from kkrunner.backends.endpoint import endpoint_for_executor
from ase2sprkkr.sprkkr.profile import SPRKKRProfile

def test_profile_uses_execution_endpoint_value():
    profile = SPRKKRProfile(executor='ssh', host='cluster', ssh_options=('-p', '2222'), ssh_zstd=True, remote_command=('kkrunner', 'remote'), psij_backend='slurm', psij_url='local://slurm')
    assert profile.endpoint == ExecutorEndpoint(executor='ssh', host='cluster', ssh_options=('-p', '2222'), ssh_zstd=True, remote_command=('kkrunner', 'remote'), psij_backend='slurm', psij_url='local://slurm')

def test_run_store_roundtrips_endpoint(tmp_path):
    store = RunStore(tmp_path / 'runs')
    manager = RunManager(store, executor_pool=ExecutorPool())
    endpoint = ExecutorEndpoint(executor='ssh', host='cluster', ssh_options=('-v',), ssh_zstd=True)
    record = manager.new_record(directory=tmp_path, endpoint=endpoint, metadata={'ase2sprkkr': {'task': str('DOS').upper(), 'dataset': 'case'}}, processing_required=True)
    stored = store.save(record)
    loaded = store.get(stored.id)
    assert loaded.execution.endpoint == endpoint
    assert isinstance(loaded.execution.endpoint, ExecutorEndpoint)

def test_manager_reconstructs_local_executor_after_restart(tmp_path):
    store = RunStore(tmp_path / 'runs')
    endpoint = ExecutorEndpoint()
    first = RunManager(store, executor_pool=ExecutorPool())
    record = first.new_record(directory=tmp_path, endpoint=endpoint, metadata={'ase2sprkkr': {'task': str('DOS').upper()}}, processing_required=True)
    store.save(record)
    restarted = RunManager(store, executor_pool=ExecutorPool())
    executor = restarted.executor_for(store.get(record.id))
    assert isinstance(executor, LocalExecutor)

def test_attach_reconstructs_persisted_endpoint_without_task_reader(tmp_path):

    class AttachedRunner:

        def __init__(self, run):
            self.backend_id = run.execution.backend_id
            self.state = RunState.RUNNING
            self.result = None
            self.exception = None
            self.failure_phase = None
            self.failure_message = None
            self._done = Event()

        def start(self):
            pass

        def wait(self):
            self._done.wait()

        def cancel(self):
            self.state = RunState.CANCELLED
            self._done.set()

        def add_output_listener(self, callback):
            pass

        def add_output_error_listener(self, callback):
            pass

    class ReconstructedExecutor:

        def __init__(self):
            self.calls = []

        def attach(self, run, callback=None, output_listener=None, start=True):
            self.calls.append((run, callback, output_listener, start))
            runner = AttachedRunner(run)
            if start:
                runner.start()
            return runner

    class RecordingPool:

        def __init__(self, executor):
            self.executor = executor
            self.endpoints = []

        def for_endpoint(self, endpoint):
            self.endpoints.append(endpoint)
            return self.executor
    store = RunStore(tmp_path / 'runs')
    endpoint = ExecutorEndpoint(executor='ssh', host='cluster')
    first = RunManager(store, executor_pool=ExecutorPool())
    record = first.new_record(directory=tmp_path, endpoint=endpoint, metadata={'ase2sprkkr': {'task': str('DOS').upper()}}, processing_required=True)
    store.save(record)
    record = first.transition(record, state=RunState.RUNNING, execution=record.execution.submitted('42'))
    executor = ReconstructedExecutor()
    pool = RecordingPool(executor)
    restarted = RunManager(store, executor_pool=pool)
    attached = restarted.attach(record.id, process_result=lambda outcome: outcome)
    assert pool.endpoints == [endpoint]
    assert attached.executor is executor
    assert executor.calls == [(record, None, None, False)]
    attached.cancel()

def test_pool_reconstructs_and_reuses_ssh_endpoint():
    endpoint = ExecutorEndpoint(executor='ssh', host='cluster', ssh_options=('-p', '2222'), ssh_zstd=True, remote_command=('kkrunner', 'remote'), psij_backend='slurm', psij_url='local://slurm')
    pool = ExecutorPool()
    first = pool.for_endpoint(endpoint)
    second = pool.for_endpoint(endpoint)
    assert first is second
    assert isinstance(first, RemoteEndpointExecutor)
    assert isinstance(first.transport, SshTransport)
    assert first.transport.host == 'cluster'
    assert first.transport.options == ('-p', '2222')
    assert first.transport.zstd is True
    assert first.backend == 'psij'
    assert first.backend_options == {'backend': 'slurm', 'url': 'local://slurm'}

def test_default_remote_psij_endpoint_survives_restart():
    live = RemoteEndpointExecutor(SshTransport('cluster'), backend='psij')
    endpoint = endpoint_for_executor(live)
    restored = ExecutorPool().for_endpoint(endpoint)
    assert endpoint.psij_backend == 'slurm'
    assert isinstance(restored, RemoteEndpointExecutor)
    assert restored.backend == 'psij'
    assert restored.backend_options == {'backend': 'slurm', 'url': None}
