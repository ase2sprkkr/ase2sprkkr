"""Explicit Docker/SSH integration test for RemoteExecutor."""

import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from threading import Event
from uuid import uuid4

import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("ASE2SPRKKR_RUN_SSH_TEST") != "1",
    reason="set ASE2SPRKKR_RUN_SSH_TEST=1 to run this integration test",
)


def _wait_for_ssh(ssh, attempts=50):
    for _ in range(attempts):
        if subprocess.run(ssh + ["true"], stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL).returncode == 0:
            return
        time.sleep(0.1)
    pytest.fail("SSH server did not become ready")


def _remote_error(exc):
    stderr = exc.stderr
    if isinstance(stderr, bytes):
        stderr = stderr.decode("utf8", errors="replace")
    pytest.fail(f"{exc}\n\nRemote stderr:\n{stderr or '<empty>'}")


def test_remote_executor_over_ssh():
    from ase.build import bulk
    from ase2sprkkr.configuration import config
    from ase2sprkkr.execution import RemoteExecutor, SshTransport
    from ase2sprkkr.sprkkr.calculator import SPRKKR

    docker = shutil.which("docker")
    ssh_keygen = shutil.which("ssh-keygen")
    if not docker:
        pytest.fail("Docker is required")
    if not ssh_keygen:
        pytest.fail("ssh-keygen is required")

    executable_dir = config.executables.dir()
    if not executable_dir:
        pytest.fail("config.executables.dir must point to SPR-KKR executables")
    executable_dir = Path(executable_dir).expanduser().resolve()
    suffix = config.executables.suffix() or ""
    executable = executable_dir / f"kkrscf{suffix}"
    if not executable.is_file():
        pytest.fail(f"SPR-KKR executable not found: {executable}")

    repository = Path(__file__).resolve().parents[4]
    context = Path(__file__).with_name("ssh")
    image = "ase2sprkkr-ssh-test"
    container = f"ase2sprkkr-ssh-test-{uuid4().hex[:12]}"
    subprocess.run([docker, "build", "-t", image, str(context)], check=True)

    with tempfile.TemporaryDirectory(prefix="ase2sprkkr-ssh-") as tmp:
        tmp = Path(tmp)
        key = tmp / "id_ed25519"
        subprocess.run([ssh_keygen, "-q", "-t", "ed25519", "-N", "",
                        "-f", str(key)], check=True)
        try:
            subprocess.run([
                docker, "run", "--rm", "-d", "--name", container,
                "-p", "127.0.0.1::22",
                "--mount", f"type=bind,src={repository},dst=/workspace/ase2sprkkr,readonly",
                "--mount", f"type=bind,src={executable_dir},dst=/opt/sprkkr,readonly",
                "--mount", f"type=bind,src={key}.pub,dst=/authorized_keys,readonly",
                image,
            ], check=True, stdout=subprocess.DEVNULL)

            port = subprocess.check_output(
                [docker, "port", container, "22/tcp"], text=True
            ).strip().splitlines()[0].rsplit(":", 1)[1]
            options = (
                "-i", str(key), "-p", port,
                "-o", "StrictHostKeyChecking=no",
                "-o", "UserKnownHostsFile=/dev/null",
            )
            ssh = ["ssh", *options, "tester@127.0.0.1"]
            _wait_for_ssh(ssh)

            executor = RemoteExecutor(
                SshTransport("tester@127.0.0.1", options=options),
                local_root=tmp,
                remote_root="/home/tester/ase2sprkkr-runs",
            )
            output_seen = Event()
            output = []

            def read_output(text, kind):
                output.append((kind, text))
                if kind == "out" and text:
                    output_seen.set()

            workdir = tmp / "calculation"
            workdir.mkdir()
            calculator = SPRKKR(
                atoms=bulk("Al"), directory=workdir, executor=executor,
                executable_dir="/opt/sprkkr", executable_suffix=suffix,
                mpi=False, empty_spheres=False, print_output=False,
            )
            runner = calculator.create_run(
                options={"NITER": 2, "NKTAB": 5, "NE": 5},
                read_callback=read_output,
            )
            try:
                runner.submit()
            except subprocess.CalledProcessError as exc:
                _remote_error(exc)

            assert runner.reference is not None
            assert runner.reference.computer == "tester@127.0.0.1"
            assert runner.reference.directory == "/home/tester/ase2sprkkr-runs/calculation"
            assert output_seen.wait(30), "no live SPR-KKR output arrived over SSH"
            result = runner.wait()
            assert result.return_code == 0
            assert len(result.iterations) == 2
            assert any(kind == "out" and text for kind, text in output)
        finally:
            subprocess.run([docker, "rm", "-f", container],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
