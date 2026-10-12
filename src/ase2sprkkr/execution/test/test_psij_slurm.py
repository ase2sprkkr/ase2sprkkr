"""Explicit Docker/Slurm integration test for the PSI/J executor."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("ASE2SPRKKR_RUN_DOCKER_TESTS") != "1",
    reason="set ASE2SPRKKR_RUN_DOCKER_TESTS=1 to run Docker integration tests",
)


def test_psij_slurm_docker():
    docker = shutil.which("docker")
    if docker is None:
        pytest.fail("Docker is required for the explicit PSI/J Slurm integration test")

    repository = Path(__file__).resolve().parents[4]
    docker_context = Path(__file__).with_name("docker")
    image = "ase2sprkkr-execution-test"

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
            "-e",
            "ASE2SPRKKR_TEST_SERVICES=slurm",
            "-e",
            "PYTHONPATH=/workspace/ase2sprkkr/src",
            image,
            "python",
            "/workspace/ase2sprkkr/src/ase2sprkkr/execution/test/docker/"
            "psij_slurm_test_in_docker.py",
        ],
        check=True,
    )
