from types import SimpleNamespace

from ase2sprkkr.execution import ExecutionManager, Executor, RunReference, RunState, RunStore
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


class _Runner:
    def __init__(self, callback=None):
        self._callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=0)

    def submit(self):
        self.reference = RunReference(executor="fake", native_id="job-1")
        self.state = RunState.QUEUED

    def wait(self):
        self.state = RunState.FINISHED
        if self._callback:
            self._callback(self.result)
        return self.result

    def run(self):
        self.submit()
        return self.wait()

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass


class _Executor(Executor):
    name = "fake"

    def __init__(self):
        self.runners = []

    def resolve_mpi(self, mpi=None):
        return False

    def create_run(self, reader, spec, callback=None):
        runner = _Runner(callback)
        self.runners.append(runner)
        return runner


def _parameters(monkeypatch):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "sample"
    parameters.CONTROL.POTFIL = "sample.pot"
    monkeypatch.setattr(
        InputParameters,
        "output_reader",
        lambda self, *args, **kwargs: SimpleNamespace(print_output=False),
    )
    monkeypatch.setattr(parameters, "input_files", lambda directory: ())
    monkeypatch.setattr(parameters, "output_files", lambda directory: ())
    monkeypatch.setattr(parameters, "cache_files", lambda directory: ())
    return parameters


def _calculator(profile, *, executor=None, persistent=False):
    return SimpleNamespace(
        profile=profile,
        executor=executor,
        persistent=persistent,
        value_or_default=lambda name, value: value,
    )


def test_sprkkr_constructor_stores_persistent_policy(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))

    calculator = SPRKKR(persistent=manager)

    assert calculator.persistent is manager


def test_explicit_manager_persists_canonical_task_and_dataset(tmp_path, monkeypatch):
    parameters = _parameters(monkeypatch)
    executor = _Executor()
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    input_file = (tmp_path / "test.inp").open("w+")
    output_file = (tmp_path / "test.out").open("wb")

    runner = parameters.run_process(
        _calculator(SPRKKRProfile(), executor=executor),
        input_file,
        output_file,
        directory=str(tmp_path),
        executor=executor,
        persistent=manager,
        run_async=True,
    )

    assert runner.id
    runner.run()
    record = manager.store.get(runner.id)
    assert record.task == parameters.task_name.upper()
    assert record.dataset == "sample"
    assert record.state == RunState.FINISHED


def test_persistence_defaults_to_calculator_setting(tmp_path, monkeypatch):
    parameters = _parameters(monkeypatch)
    executor = _Executor()
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    input_file = (tmp_path / "test.inp").open("w+")
    output_file = (tmp_path / "test.out").open("wb")

    runner = parameters.run_process(
        _calculator(SPRKKRProfile(), executor=executor, persistent=manager),
        input_file,
        output_file,
        directory=str(tmp_path),
        executor=executor,
        run_async=True,
    )

    runner.submit()
    assert manager.store.get(runner.id).state == RunState.QUEUED


def test_run_level_profile_selects_executor_instead_of_calculator_executor(
    tmp_path, monkeypatch
):
    parameters = _parameters(monkeypatch)
    calculator_executor = _Executor()
    profile_executor = _Executor()
    profile = SPRKKRProfile(executor="ssh", host="cluster")
    monkeypatch.setattr(
        "ase2sprkkr.input_parameters.input_parameters.executor_from_profile",
        lambda selected: profile_executor,
    )
    input_file = (tmp_path / "test.inp").open("w+")
    output_file = (tmp_path / "test.out").open("wb")

    parameters.run_process(
        _calculator(SPRKKRProfile(), executor=calculator_executor),
        input_file,
        output_file,
        directory=str(tmp_path),
        profile=profile,
        run_async=True,
    )

    assert len(profile_executor.runners) == 1
    assert calculator_executor.runners == []
