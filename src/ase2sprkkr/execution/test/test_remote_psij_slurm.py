"""Explicit SSH -> remote PSI/J -> Slurm -> SPR-KKR integration test."""

import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from uuid import uuid4

import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("ASE2SPRKKR_RUN_SSH_TEST") != "1"
    or os.environ.get("ASE2SPRKKR_RUN_SLURM_TEST") != "1",
    reason=(
        "set ASE2SPRKKR_RUN_SSH_TEST=1 and ASE2SPRKKR_RUN_SLURM_TEST=1 "
        "to run the SSH/PSI-J/Slurm integration test"
    ),
)


def _wait_for_ssh(ssh, attempts=80):
    for _ in range(attempts):
        if subprocess.run(
            ssh + ["true"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode == 0:
            return
        time.sleep(0.1)
    pytest.fail("SSH server did not become ready")


def test_remote_executor_runs_psij_slurm_over_ssh():
    from ase.build import bulk

    from ase2sprkkr.configuration import config
    from ase2sprkkr.execution import RemoteRunner
    from ase2sprkkr.sprkkr.calculator import SPRKKR
    from ase2sprkkr.sprkkr.profile import SPRKKRProfile

    docker = shutil.which("docker")
    ssh_keygen = shutil.which("ssh-keygen")
    if docker is None:
        pytest.fail("Docker is required for the explicit SSH/PSI-J/Slurm test")
    if ssh_keygen is None:
        pytest.fail("ssh-keygen is required for the explicit SSH/PSI-J/Slurm test")

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
    container = f"ase2sprkkr-psij-ssh-test-{uuid4().hex[:12]}"

    subprocess.run(
        [docker, "build", "-t", image, str(docker_context)],
        check=True,
    )

    with tempfile.TemporaryDirectory(prefix="ase2sprkkr-psij-ssh-") as tmp:
        tmp = Path(tmp)
        key = tmp / "id_ed25519"
        subprocess.run(
            [
                ssh_keygen,
                "-q",
                "-t",
                "ed25519",
                "-N",
                "",
                "-f",
                str(key),
            ],
            check=True,
        )

        try:
            subprocess.run(
                [
                    docker,
                    "run",
                    "--rm",
                    "-d",
                    "--name",
                    container,
                    "--privileged",
                    "--cgroupns=private",
                    "--hostname",
                    "node1",
                    "-p",
                    "127.0.0.1::22",
                    "--mount",
                    f"type=bind,src={repository},dst=/workspace/ase2sprkkr,readonly",
                    "--mount",
                    f"type=bind,src={executable_dir},dst=/opt/sprkkr,readonly",
                    "--mount",
                    f"type=bind,src={key}.pub,dst=/authorized_keys,readonly",
                    "-e",
                    f"SPRKKR_EXECUTABLE_SUFFIX={suffix}",
                    image,
                    "sleep",
                    "infinity",
                ],
                check=True,
                stdout=subprocess.DEVNULL,
            )

            port = subprocess.check_output(
                [docker, "port", container, "22/tcp"], text=True
            ).strip().splitlines()[0].rsplit(":", 1)[1]
            ssh_options = (
                "-i",
                str(key),
                "-p",
                port,
                "-o",
                "StrictHostKeyChecking=no",
                "-o",
                "UserKnownHostsFile=/dev/null",
            )
            ssh = ["ssh", *ssh_options, "tester@127.0.0.1"]
            _wait_for_ssh(ssh)

            workdir = tmp / "calculation"
            workdir.mkdir()
            profile = SPRKKRProfile(
                executable_dir="/opt/sprkkr",
                executable_suffix=suffix,
                mpi=False,
                executor="ssh",
                host="tester@127.0.0.1",
                local_root=str(tmp),
                remote_root="/home/tester/ase2sprkkr-runs",
                ssh_options=ssh_options,
                psij_backend="slurm",
            )
            calculator = SPRKKR(
                atoms=bulk("Al"),
                directory=workdir,
                profile=profile,
                empty_spheres=False,
                print_output=False,
            )
            runner = calculator.create_run(
                options={"NITER": 2, "NKTAB": 5, "NE": 5},
            )

            assert isinstance(runner, RemoteRunner)
            assert runner.executor.backend == "psij"
            assert runner.executor.backend_options == {
                "backend": "slurm",
                "url": None,
            }

            runner.submit()
            reference = runner.reference
            assert reference is not None
            assert reference.executor == "remote:psij"
            assert reference.computer == "tester@127.0.0.1"
            assert reference.directory == "/home/tester/ase2sprkkr-runs/calculation"

            result = runner.wait()
            assert runner.state.value == "finished"
            assert result.return_code == 0
            assert len(result.iterations) == 2

            remote_state = json.loads(
                subprocess.check_output(
                    ssh
                    + [
                        "ase2sprkkr",
                        "remote",
                        "status",
                        reference.directory,
                        reference.native_id,
                    ],
                    text=True,
                )
            )
            backend_reference = remote_state["backend_reference"]
            assert backend_reference["executor"] == "psij:slurm"
            assert backend_reference["native_id"]
            assert remote_state["state"] == "finished"
            assert remote_state["return_code"] == 0
        finally:
            subprocess.run(
                [docker, "rm", "-f", container],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
