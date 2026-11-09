"""Shared fixtures for execution integration tests."""

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import pytest
from kkrunner.backends import pool as execution_pool


@dataclass(frozen=True)
class DockerExecutionImage:
    """Built Docker image and paths shared by execution integration tests."""

    docker: str
    repository: Path
    image: str


def _package_mount_root() -> Path:
    """Return a host path importable inside both source and wheel test layouts."""
    import ase2sprkkr

    package_dir = Path(ase2sprkkr.__file__).resolve().parent
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
    source_context = Path(__file__).with_name("docker")
    image = "ase2sprkkr-execution-test"
    kkrunner_version = version("kkrunner")
    with tempfile.TemporaryDirectory(prefix="ase2sprkkr-docker-") as temporary:
        context = Path(temporary)
        shutil.copytree(source_context, context, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns("__pycache__"))
        subprocess.run(
            [docker, "build", "--build-arg", f"KKRUNNER_VERSION={kkrunner_version}", "-t", image, str(context)],
            check=True,
        )
    return DockerExecutionImage(docker=docker, repository=repository, image=image)
