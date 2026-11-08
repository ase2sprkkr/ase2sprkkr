from types import SimpleNamespace

from kkrunner import RunState
from kkrunner.backends import LocalExecutor
from kkrunner.runs import LaunchSpec

from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.input_parameters.task_definition import TaskDefinition
from ase2sprkkr.sprkkr import calculation as calculation_module
from ase2sprkkr.sprkkr.calculation import PreparedCalculation
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


class _Prepared:
    create_run = PreparedCalculation.create_run

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


def test_input_parameters_preserve_mpi_launch_intent(tmp_path, monkeypatch):
    reader = object()
    parameters = InputParameters.create('SCF')
    monkeypatch.setattr(TaskDefinition, 'output_reader', lambda self, **kwargs: reader)
    seen = {}
    runner = SimpleNamespace(state=RunState.NEW, result=object(), reference=None, exception=None)

    class Executor:

        def resolve_mpi(self, mpi):
            seen['executor_mpi'] = mpi
            return LaunchSpec.from_mpi(mpi)

        def create_run(self, spec, callback=None):
            seen['spec'] = spec
            return runner
    profile = SPRKKRProfile(mpi='auto')
    calculator = SimpleNamespace(directory=str(tmp_path))
    with (tmp_path / 'run.inp').open('w+') as input_file, (tmp_path / 'run.out').open('wb') as output_file:
        created = _Prepared(parameters, input_file, output_file, tmp_path).create_run(result_context=calculator, profile=profile, executor=Executor())
        assert created.result is None
        assert seen['executor_mpi'] == 'auto'
        assert seen['spec'].argv == ('kkrscf',)
        assert seen['spec'].mpi_argv == ('kkrscfMPI', input_file.name)
        assert seen['spec'].launch == LaunchSpec.from_mpi('auto')
        assert seen['spec'].stdin is input_file
        assert seen['spec'].mpi_uses_stdin is False
        created.close()


def test_explicit_local_executor_overrides_remote_profile(tmp_path, monkeypatch):
    reader = object()
    parameters = InputParameters.create('SCF')
    monkeypatch.setattr(TaskDefinition, 'output_reader', lambda self, **kwargs: reader)
    runner = SimpleNamespace(state=RunState.NEW, result=object(), reference=None, exception=None)
    seen = {}

    class ExplicitLocalExecutor(LocalExecutor):

        def resolve_mpi(self, mpi):
            seen['mpi'] = mpi
            return LaunchSpec()

        def create_run(self, spec, callback=None):
            seen['spec'] = spec
            return runner

    def must_not_resolve_profile(profile):
        raise AssertionError('explicit executor was replaced by the profile backend')
    monkeypatch.setattr(calculation_module, 'executor_from_profile', must_not_resolve_profile)
    profile = SPRKKRProfile(executor='ssh', host='cluster', mpi='auto')
    executor = ExplicitLocalExecutor()
    calculator = SimpleNamespace(directory=str(tmp_path))
    with (tmp_path / 'run.inp').open('w+') as input_file, (tmp_path / 'run.out').open('wb') as output_file:
        created = _Prepared(parameters, input_file, output_file, tmp_path).create_run(result_context=calculator, profile=profile, executor=executor)
        assert created.result is None
        assert seen['mpi'] == 'auto'
        assert seen['spec'].argv == ('kkrscf',)
        assert seen['spec'].mpi_argv is None
        assert seen['spec'].launch == LaunchSpec()
        created.close()
