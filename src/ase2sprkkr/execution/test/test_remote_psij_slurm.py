"""Explicit SSH -> remote PSI/J -> Slurm -> SPR-KKR integration test."""

import base64
import json
import os
import struct
import subprocess
import tempfile
import time
from pathlib import Path
from uuid import uuid4

import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("ASE2SPRKKR_RUN_DOCKER_TESTS") != "1",
    reason="set ASE2SPRKKR_RUN_DOCKER_TESTS=1 to run Docker integration tests",
)

# Static, public test-only Ed25519 key material. Keeping it in this module makes
# it part of the installed test package while avoiding key generation on every
# integration-test run. It must never be used outside this isolated test.
_TEST_SSH_SEED = bytes(range(1, 33))
_TEST_SSH_PUBLIC = bytes.fromhex(
    "79b5562e8fe654f94078b112e8a98ba7901f853ae695bed7e0e3910bad049664"
)
_TEST_SSH_PUBLIC_LINE = (
    "ssh-ed25519 "
    "AAAAC3NzaC1lZDI1NTE5AAAAIHm1Vi6P5lT5QHixEuipi6eQH4U65pW+1+DjkQutBJZk "
    "ase2sprkkr test-only key\n"
)


def _ssh_string(value):
    return struct.pack(">I", len(value)) + value


def _write_test_ssh_key(directory):
    """Serialize the packaged static test key into OpenSSH files."""
    public_blob = _ssh_string(b"ssh-ed25519") + _ssh_string(_TEST_SSH_PUBLIC)
    comment = b"ase2sprkkr test-only key"
    private_blob = (
        struct.pack(">II", 0x12345678, 0x12345678)
        + _ssh_string(b"ssh-ed25519")
        + _ssh_string(_TEST_SSH_PUBLIC)
        + _ssh_string(_TEST_SSH_SEED + _TEST_SSH_PUBLIC)
        + _ssh_string(comment)
    )
    padding = 1
    while len(private_blob) % 8:
        private_blob += bytes((padding,))
        padding += 1

    key_blob = (
        b"openssh-key-v1\0"
        + _ssh_string(b"none")
        + _ssh_string(b"none")
        + _ssh_string(b"")
        + struct.pack(">I", 1)
        + _ssh_string(public_blob)
        + _ssh_string(private_blob)
    )
    encoded = base64.b64encode(key_blob).decode("ascii")
    pem = (
        "-----BEGIN OPENSSH PRIVATE KEY-----\n"
        + "\n".join(encoded[i : i + 70] for i in range(0, len(encoded), 70))
        + "\n-----END OPENSSH PRIVATE KEY-----\n"
    )

    key = Path(directory) / "id_ed25519"
    key.write_text(pem, encoding="ascii")
    key.chmod(0o600)
    Path(str(key) + ".pub").write_text(_TEST_SSH_PUBLIC_LINE, encoding="ascii")
    return key


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


def _remote_error(exc):
    stderr = exc.stderr
    if isinstance(stderr, bytes):
        stderr = stderr.decode("utf8", errors="replace")
    pytest.fail(f"{exc}\n\nRemote stderr:\n{stderr or '<empty>'}")


def _ssh_capture(ssh, *command):
    return subprocess.run(
        ssh + list(command),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _diagnostic_command(ssh, title, *command):
    result = _ssh_capture(ssh, *command)
    parts = [f"{title} (exit {result.returncode}):"]
    if result.stdout:
        parts.append(result.stdout.rstrip())
    if result.stderr:
        parts.append("stderr:\n" + result.stderr.rstrip())
    if len(parts) == 1:
        parts.append("<empty>")
    return "\n".join(parts)


def _remote_run_diagnostics(ssh, reference):
    """Return stateless endpoint and scheduler diagnostics."""
    status = _ssh_capture(
        ssh,
        "ase2sprkkr",
        "remote",
        "status",
        reference.native_id,
    )
    sections = []
    if status.returncode:
        sections.append(
            f"Remote status command failed with exit {status.returncode}.\n"
            f"stdout:\n{status.stdout or '<empty>'}\n"
            f"stderr:\n{status.stderr or '<empty>'}"
        )
    else:
        try:
            state = json.loads(status.stdout)
            sections.append("Remote status:\n" + json.dumps(state, indent=2, sort_keys=True))
            error_file = state.get("error_file")
            if error_file:
                sections.append(
                    _diagnostic_command(
                        ssh, "Scheduler stderr", "cat", "--", error_file
                    )
                )
        except json.JSONDecodeError as exc:
            sections.append(
                f"Could not parse remote status: {exc}\n"
                f"stdout:\n{status.stdout or '<empty>'}\n"
                f"stderr:\n{status.stderr or '<empty>'}"
            )

    sections.append(_diagnostic_command(ssh, "Slurm queue", "squeue", "-l"))
    sections.append(
        _diagnostic_command(
            ssh,
            "Slurm accounting",
            "sacct",
            "-X",
            "-n",
            "-P",
            "--format=JobIDRaw,Comment,State,ExitCode",
        )
    )
    return "\n\n".join(sections)


def test_remote_executor_runs_psij_slurm_over_ssh(docker_execution_image):
    from ase.build import bulk

    from ase2sprkkr.sprkkr.calculator import SPRKKR
    from ase2sprkkr.sprkkr.profile import SPRKKRProfile

    docker = docker_execution_image.docker
    repository = docker_execution_image.repository
    image = docker_execution_image.image
    container = f"ase2sprkkr-psij-ssh-test-{uuid4().hex[:12]}"

    with tempfile.TemporaryDirectory(prefix="ase2sprkkr-psij-ssh-") as tmp:
        tmp = Path(tmp)
        key = _write_test_ssh_key(tmp)

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
                    f"type=bind,src={key}.pub,dst=/authorized_keys,readonly",
                    "-e",
                    "ASE2SPRKKR_TEST_SERVICES=ssh,slurm",
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
                executable_suffix="",
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

            assert runner.executor.backend == "psij"
            assert runner.executor.backend_options == {
                "backend": "slurm",
                "url": None,
            }

            try:
                runner.submit()
            except subprocess.CalledProcessError as exc:
                _remote_error(exc)
            reference = runner.reference
            assert reference is not None
            assert reference.executor == "remote:psij"
            assert reference.computer == "tester@127.0.0.1"
            assert reference.directory == "/home/tester/ase2sprkkr-runs/calculation"
            assert reference.native_id

            try:
                result = runner.wait()
            except Exception as exc:
                diagnostics = _remote_run_diagnostics(ssh, reference)
                raise RuntimeError(
                    f"Remote PSI/J/Slurm run failed: {exc}\n\n{diagnostics}"
                ) from exc

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
                        reference.native_id,
                    ],
                    text=True,
                )
            )
            assert remote_state["state"] == "finished"
            assert remote_state["return_code"] == 0

            lifecycle_dir = f"{reference.directory}/.ase2sprkkr/runs"
            assert subprocess.run(
                ssh + ["test", "!", "-e", lifecycle_dir]
            ).returncode == 0
        finally:
            subprocess.run(
                [docker, "rm", "-f", container],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
