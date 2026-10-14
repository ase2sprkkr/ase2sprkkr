"""Shared fixtures for execution integration tests."""

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess

import pytest

from ase2sprkkr.execution import pool as execution_pool


@dataclass(frozen=True)
class DockerExecutionImage:
    """Built Docker image and paths shared by execution integration tests."""

    docker: str
    repository: Path
    image: str


def _package_mount_root() -> Path:
    """Return a host path importable inside both source and wheel test layouts."""
    package_dir = Path(__file__).resolve().parents[2]
    package_parent = package_dir.parent
    if package_parent.name == "src":
        return package_parent.parent
    return package_parent


@pytest.fixture(autouse=True)
def isolated_default_executor_pool():
    """Keep the process-wide executor pool from leaking between test cases."""
    previous = execution_pool._DEFAULT_EXECUTOR_POOL
    execution_pool._DEFAULT_EXECUTOR_POOL = execution_pool.ExecutorPool()
    try:
        yield
    finally:
        execution_pool._DEFAULT_EXECUTOR_POOL.close()
        execution_pool._DEFAULT_EXECUTOR_POOL = previous


@pytest.fixture(scope="session")
def docker_execution_image():
    """Build the common execution-test image once per pytest session."""
    docker = shutil.which("docker")
    if docker is None:
        pytest.skip("Docker is required for the execution integration tests")

    repository = _package_mount_root()
    context = Path(__file__).with_name("docker")
    image = "ase2sprkkr-execution-test"
    subprocess.run([docker, "build", "-t", image, str(context)], check=True)
    return DockerExecutionImage(docker=docker, repository=repository, image=image)
