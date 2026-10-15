from types import SimpleNamespace

from ase2sprkkr.execution import ExecutionManager, Executor, RunReference, RunState, RunStore
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.sprkkr import calculation as calculation_module
from ase2sprkkr.sprkkr import calculator as calculator_module
from ase2sprkkr.sprkkr.calculation import create_sprkkr_run
from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


class _Runner:
    def __init__(self, callback=None):
        self._callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=0)
        self.exception = None

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


class _Prepared:
    def __init__(self, parameters, input_file, output_file, directory):
        self.input_parameters = parameters
        self.input_file = input_file
        self.potential_file = None
        self.output_file = output_file
        self.directory = str(directory)
        self.closed = False

    def close(self):
        if self.closed:
            return
        self.closed = True
        self.input_file.close()
        self.output_file.close()


def _parameters(monkeypatch):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "sample"
    parameters.CONTROL.POTFIL = "sample.pot"
    monkeypatch.setattr(
        InputParameters,
        "output_reader",
        lambda self, *args, **kwargs: SimpleNamespace(print_output=False),
    )
    monkeypatch.setattr(InputParameters, "input_files", lambda self, directory: ())
    monkeypatch.setattr(InputParameters, "output_files", lambda self, directory: ())
    monkeypatch.setattr(InputParameters, "cache_files", lambda self, directory: ())
    return parameters


def test_explicit_manager_persists_canonical_task_and_dataset(tmp_path, monkeypatch):
    parameters = _parameters(monkeypatch)
    executor = _Executor()
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    input_file = (tmp_path / "test.inp").open("w+")
    output_file = (tmp_path / "test.out").open("wb")

    runner = create_sprkkr_run(
        _Prepared(parameters, input_file, output_file, tmp_path),
        result_context=SimpleNamespace(directory=str(tmp_path)),
        profile=SPRKKRProfile(),
        executor=executor,
        persistent=manager,
    )

    assert runner.id
    runner.run()
    record = manager.store.get(runner.id)
    assert record.task == parameters.execution_task_name.upper()
    assert record.dataset == "sample"
    assert record.state == RunState.FINISHED


def test_persistence_defaults_to_calculator_setting(tmp_path, monkeypatch):
    executor = _Executor()
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    calculator = SPRKKR(directory=tmp_path, executor=executor, persistent=manager)
    prepared = object()
    captured = {}

    monkeypatch.setattr(calculator, "_prepare_calculation", lambda **kwargs: prepared)

    def create_run(value, **kwargs):
        assert value is prepared
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(calculator_module, "_create_sprkkr_run", create_run)

    calculator.create_run()

    assert captured["persistent"] is manager
    assert captured["executor"] is executor


def test_run_level_profile_selects_executor_instead_of_calculator_executor(
    tmp_path, monkeypatch
):
    parameters = _parameters(monkeypatch)
    calculator_executor = _Executor()
    profile_executor = _Executor()
    profile = SPRKKRProfile(executor="ssh", host="cluster")
    calculator = SPRKKR(directory=tmp_path, executor=calculator_executor)
    input_file = (tmp_path / "test.inp").open("w+")
    output_file = (tmp_path / "test.out").open("wb")
    prepared = _Prepared(parameters, input_file, output_file, tmp_path)

    monkeypatch.setattr(calculator, "_prepare_calculation", lambda **kwargs: prepared)
    monkeypatch.setattr(
        calculation_module,
        "executor_from_profile",
        lambda selected: profile_executor,
    )

    runner = calculator.create_run(profile=profile)
    try:
        assert len(profile_executor.runners) == 1
        assert calculator_executor.runners == []
    finally:
        runner.close()
