"""Inner workload for the PSI/J Slurm Docker integration test."""

from tempfile import TemporaryDirectory

from ase.build import bulk
from ase2sprkkr.execution import PsijExecutor, RunState
from ase2sprkkr.sprkkr.calculator import SPRKKR


def main():
    executor = PsijExecutor("slurm")

    with TemporaryDirectory(prefix="ase2sprkkr-psij-") as directory:
        calculator = SPRKKR(
            atoms=bulk("Al"),
            directory=directory,
            executor=executor,
            executable_dir="/opt/sprkkr",
            executable_suffix="",
            mpi=False,
            empty_spheres=False,
            print_output=False,
        )
        runner = calculator.create_run(
            options={"NITER": 2, "NKTAB": 5, "NE": 5},
        )
        assert runner.executor is executor
        assert runner.state == RunState.NEW

        runner.submit()
        assert runner.backend_id is not None

        result = runner.wait()
        assert runner.state == RunState.FINISHED
        assert result.return_code == 0
        assert len(result.iterations) == 2

    print("PSI/J -> Slurm -> SPR-KKR Al SCF integration test passed")


if __name__ == "__main__":
    main()
