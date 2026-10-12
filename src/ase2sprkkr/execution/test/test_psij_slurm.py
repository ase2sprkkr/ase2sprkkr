"""Explicit Docker/Slurm integration test for the PSI/J executor."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("ASE2SPRKKR_RUN_SLURM_TEST") != "1",
    reason="set ASE2SPRKKR_RUN_SLURM_TEST=1 to run the Docker/Slurm integration test",
)


def test_psij_slurm_docker():
    from ase2sprkkr.configuration import config

    docker = shutil.which("docker")
    if docker is None:
        pytest.fail("Docker is required for the explicit PSI/J Slurm integration test")

    executable_dir = config.executables.dir()
    if not executable_dir:
        pytest.fail("config.executables.dir must point to the SPR-KKR executable directory")
    executable_dir = Path(executable_dir).expanduser().resolve()
    suffix = config.executables.suffix() or ""
    executable = executable_dir / f"kkrscf{suffix}"
    if not executable.is_file():
        pytest.fail(f"SPR-KKR executable not found: {executable}")

    repository = Path(__file__).resolve().parents[4]
    docker_context = Path(__file__).with_name("slurm")
    image = "ase2sprkkr-psij-slurm-test"

    subprocess.run(
        [docker, "build", "-t", image, str(docker_context)],
        check=True,
    )
    subprocess.run(
        [
            docker,
            "run",
            "--rm",
            "--privileged",
            "--cgroupns=private",
            "--hostname",
            "node1",
            "--mount",
            f"type=bind,src={repository},dst=/workspace/ase2sprkkr,readonly",
            "--mount",
            f"type=bind,src={executable_dir},dst=/opt/sprkkr,readonly",
            "-e",
            "PYTHONPATH=/workspace/ase2sprkkr/src",
            "-e",
            f"SPRKKR_EXECUTABLE_SUFFIX={suffix}",
            image,
            "python",
            "/workspace/ase2sprkkr/src/ase2sprkkr/execution/test/slurm/run_test.py",
        ],
        check=True,
    )
