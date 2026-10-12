"""Shared fixtures for execution integration tests."""

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess

import pytest


@dataclass(frozen=True)
class DockerExecutionImage:
    """Built Docker image and paths shared by execution integration tests."""

    docker: str
    repository: Path
    image: str


@pytest.fixture(scope="session")
def docker_execution_image():
    """Build the common execution-test image once per pytest session."""
    docker = shutil.which("docker")
    if docker is None:
        pytest.skip("Docker is required for the execution integration tests")

    repository = Path(__file__).resolve().parents[4]
    context = Path(__file__).with_name("docker")
    image = "ase2sprkkr-execution-test"
    subprocess.run([docker, "build", "-t", image, str(context)], check=True)
    return DockerExecutionImage(docker=docker, repository=repository, image=image)
