"""Explicit Docker/SSH integration test for RemoteExecutor."""

import os
import posixpath
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from threading import Event
from uuid import uuid4

import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("ASE2SPRKKR_RUN_DOCKER_TESTS") != "1",
    reason="set ASE2SPRKKR_RUN_DOCKER_TESTS=1 to run Docker integration tests",
)


def _wait_for_ssh(ssh, attempts=50):
    for _ in range(attempts):
        if subprocess.run(
            ssh + ["true"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode == 0:
            return
        time.sleep(0.1)
    pytest.fail("SSH server did not become ready")


def _remote_error(exc):
    stderr = exc.stderr
    if isinstance(stderr, bytes):
        stderr = stderr.decode("utf8", errors="replace")
    pytest.fail(f"{exc}\n\nRemote stderr:\n{stderr or '<empty>'}")


def _wait_for_live_output(runner, output_seen, output, timeout=30):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if output_seen.wait(0.1):
            return
        state = runner.state
        if state.value in {"finished", "failed", "cancelled"}:
            try:
                result = runner.wait()
            except Exception as exc:
                pytest.fail(
                    f"remote run ended before live output: state={state.value}, "
                    f"exception={exc!r}, output={output!r}"
                )
            pytest.fail(
                f"remote run ended before live output: state={state.value}, "
                f"return_code={result.return_code}, output={output!r}"
            )
    pytest.fail(
        f"no live SPR-KKR output arrived over SSH; "
        f"state={runner.state.value}, output={output!r}"
    )


def test_remote_executor_over_ssh(docker_execution_image):
    from ase.build import bulk
    from ase2sprkkr.execution import (
        ExecutionManager,
        ExecutorPool,
        RemoteExecutor,
        ProjectRepository,
        SshTransport,
    )
    from ase2sprkkr.outputs.project_results import load_project_result
    from ase2sprkkr.sprkkr.calculator import SPRKKR
    from ase2sprkkr.sprkkr.profile import SPRKKRProfile

    docker = docker_execution_image.docker
    repository = docker_execution_image.repository
    image = docker_execution_image.image
    ssh_keygen = shutil.which("ssh-keygen")
    if not ssh_keygen:
        pytest.fail("ssh-keygen is required")

    container = f"ase2sprkkr-ssh-test-{uuid4().hex[:12]}"

    with tempfile.TemporaryDirectory(prefix="ase2sprkkr-ssh-") as tmp:
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
                    "-p",
                    "127.0.0.1::22",
                    "--mount",
                    f"type=bind,src={repository},dst=/workspace/ase2sprkkr,readonly",
                    "--mount",
                    f"type=bind,src={key}.pub,dst=/authorized_keys,readonly",
                    "-e",
                    "ASE2SPRKKR_TEST_SERVICES=ssh",
                    image,
                ],
                check=True,
                stdout=subprocess.DEVNULL,
            )

            port = subprocess.check_output(
                [docker, "port", container, "22/tcp"], text=True
            ).strip().splitlines()[0].rsplit(":", 1)[1]
            options = (
                "-i",
                str(key),
                "-p",
                port,
                "-o",
                "StrictHostKeyChecking=no",
                "-o",
                "UserKnownHostsFile=/dev/null",
            )
            ssh = ["ssh", *options, "tester@127.0.0.1"]
            _wait_for_ssh(ssh)

            executor = RemoteExecutor(
                SshTransport("tester@127.0.0.1", options=options),
            )
            profile = SPRKKRProfile(
                executable_dir="/opt/sprkkr",
                executable_suffix="",
                mpi=False,
                local_root=tmp,
                remote_root="/home/tester/ase2sprkkr-runs",
            )
            projects = ProjectRepository(tmp / "projects")
            project = projects.create("SSH integration")
            manager = ExecutionManager(
                executor_pool=ExecutorPool(), project_repository=projects
            ).for_project(project)
            output_seen = Event()
            output = []
            first_output_at = []

            def read_output(text, kind):
                output.append((kind, text))
                if kind == "out" and text:
                    if not first_output_at:
                        first_output_at.append(time.monotonic())
                    output_seen.set()

            workdir = tmp / "calculation"
            workdir.mkdir()
            calculator = SPRKKR(
                atoms=bulk("Al"),
                directory=workdir,
                executor=executor,
                profile=profile,
                project=manager,
                empty_spheres=False,
                print_output=False,
            )
            runner = calculator.create_run(
                options={"NITER": 2, "NKTAB": 5, "NE": 5},
                read_callback=read_output,
            )
            try:
                runner.submit()
            except subprocess.CalledProcessError as exc:
                _remote_error(exc)

            assert runner.backend_id is not None
            assert runner.execution.endpoint.host == "tester@127.0.0.1"
            assert (
                runner.execution.directory.remote
                == "/home/tester/ase2sprkkr-runs/calculation"
            )
            _wait_for_live_output(runner, output_seen, output)

            result = runner.wait()
            finished_at = time.monotonic()
            assert result.return_code == 0
            assert len(result.iterations) > 0
            assert any(kind == "out" and text for kind, text in output)
            assert first_output_at
            assert finished_at - first_output_at[0] >= 0.5

            potfil = result.input_parameters.CONTROL.POTFIL()
            remote_potential = (
                potfil
                if posixpath.isabs(potfil)
                else posixpath.join(runner.execution.directory.remote, potfil)
            )
            remote_pot_new = remote_potential + "_new"
            subprocess.run(ssh + ["test", "-s", remote_pot_new], check=True)
            subprocess.run(
                ssh + ["cmp", "--", remote_potential, remote_pot_new],
                check=True,
            )

            # Simulate a client restart. A fresh manager/pool knows only the
            # persisted run id. It must reconstruct endpoint identity from the
            # client-owned record, while terminal result replay is handled by
            # the outputs layer rather than by attaching a live runner.
            run_id = runner.id
            restarted = ExecutionManager(
                executor_pool=ExecutorPool(), project_repository=projects
            ).for_project(project)
            record = restarted.store.get(run_id)
            restarted_executor = restarted.executor_for(record)
            replayed = load_project_result(restarted, run_id)

            assert restarted_executor is not executor
            assert restarted_executor.transport.host == "tester@127.0.0.1"
            assert restarted_executor.transport.options == options
            assert replayed.return_code == 0
            assert len(replayed.iterations) > 0
        finally:
            subprocess.run(
                [docker, "rm", "-f", container],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
