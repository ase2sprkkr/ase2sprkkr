from types import SimpleNamespace

from ase2sprkkr.execution import (
    ExecutionManager,
    Executor,
    RunReference,
    RunState,
    RunStore,
)


class _Runner:
    def __init__(self, callback=None):
        self.callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=None)

    def submit(self):
        self.reference = RunReference(executor="history", native_id="job-1")
        self.state = RunState.QUEUED

    def wait(self):
        return self.result

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass

    def finish(self):
        self.result.return_code = 0
        self.state = RunState.FINISHED
        self.callback(self.result)


class _HistoryExecutor(Executor):
    name = "history"

    def __init__(self):
        self.runner = None
        self.fetch_calls = []

    def create_run(self, reader, spec, callback=None):
        self.runner = _Runner(callback)
        return self.runner

    def fetch_result_files(self, reference, result_files):
        self.fetch_calls.append((reference, result_files))
        with open(result_files.output_file, "w") as output:
            output.write("restored output")


def _spec(tmp_path):
    input_file = tmp_path / "run.inp"
    output_file = tmp_path / "run.out"
    input_file.write_text("input")
    return SimpleNamespace(
        directory=str(tmp_path),
        input_file=str(input_file),
        output_file=str(output_file),
        input_files=(str(input_file),),
        output_files=(str(output_file),),
        cache_files=(),
    )


def test_history_listing_does_not_fetch_missing_result_files(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _HistoryExecutor()
    spec = _spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("original output")
    executor.runner.finish()
    (tmp_path / "run.out").unlink()

    records = store.runs()

    assert records[0].id == run.id
    assert executor.fetch_calls == []


def test_load_result_fetches_missing_result_files_on_demand(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _HistoryExecutor()
    spec = _spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("original output")
    executor.runner.finish()
    (tmp_path / "run.out").unlink()

    parsed = []

    class Reader:
        def read_from_file(self, output, **kwargs):
            parsed.append((output, kwargs))
            return SimpleNamespace(return_code=kwargs["return_code"])

    result = manager.load_result(run.id, reader=Reader())

    assert len(executor.fetch_calls) == 1
    assert (tmp_path / "run.out").read_text() == "restored output"
    assert parsed[0][0] == str(tmp_path / "run.out")
    assert result.return_code == 0
    assert result.run_record.state == RunState.FINISHED


def test_failed_history_without_return_code_is_not_parsed_as_success(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    spec = _spec(tmp_path)
    (tmp_path / "run.out").write_text("failed output")
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    store.save(record)
    record = manager.finalize(record, state=RunState.FAILED)
    parsed = []

    class Reader:
        def read_from_file(self, output, **kwargs):
            parsed.append(kwargs)
            return SimpleNamespace(return_code=kwargs["return_code"])

    result = manager.load_result(record.id, reader=Reader())

    assert parsed[0]["return_code"] != 0
    assert result.return_code != 0
    assert result.run_record.state == RunState.FAILED


def test_missing_failed_history_uses_inferred_return_code(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    spec = _spec(tmp_path)
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    store.save(record)
    record = manager.finalize(record, state=RunState.FAILED)

    class Result:
        return_code = None

        def complete(self, error, return_code):
            self.return_code = return_code

    class Reader:
        def create_result(self, output, input_file):
            return Result()

    result = manager.load_result(record.id, reader=Reader())

    assert result.return_code == 1
    assert result.run_record.state == RunState.FAILED
