"""Low-level command and file transports used by remote execution."""

from __future__ import annotations

import os
import posixpath
import shlex
import shutil
import subprocess
from pathlib import Path


class Transport:
    """Run commands and move files on one target machine."""

    def exec(self, argv, *, cwd=None, input=None, check=True):
        raise NotImplementedError

    def exec_stream(self, argv, *, cwd=None, stdin=None):
        raise NotImplementedError

    def put_file(self, source, destination):
        raise NotImplementedError

    def get_file(self, source, destination):
        raise NotImplementedError

    def read_file(self, path):
        raise NotImplementedError


class LocalTransport(Transport):
    """Transport backed by the local subprocess and filesystem APIs."""

    def exec(self, argv, *, cwd=None, input=None, check=True):
        return subprocess.run(
            tuple(os.fspath(arg) for arg in argv),
            cwd=cwd,
            input=input,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=check,
        )

    def exec_stream(self, argv, *, cwd=None, stdin=None):
        return subprocess.Popen(
            tuple(os.fspath(arg) for arg in argv),
            cwd=cwd,
            stdin=stdin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def put_file(self, source, destination):
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)

    def get_file(self, source, destination):
        self.put_file(source, destination)

    def read_file(self, path):
        return Path(path).read_bytes()


class SshTransport(Transport):
    """Transport using the system ``ssh`` client.

    The remote side only needs a POSIX shell and common file utilities.  No
    Python SSH library is required on the controlling machine.
    """

    def __init__(self, host, *, ssh="ssh", options=()):
        self.host = host
        self.ssh = ssh
        self.options = tuple(options)

    def _command(self, argv, cwd=None):
        command = shlex.join([os.fspath(arg) for arg in argv])
        if cwd is not None:
            command = "cd {} && exec {}".format(
                shlex.quote(os.fspath(cwd)), command
            )
        return [self.ssh, *self.options, self.host, command]

    def exec(self, argv, *, cwd=None, input=None, check=True):
        return subprocess.run(
            self._command(argv, cwd),
            input=input,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=check,
        )

    def exec_stream(self, argv, *, cwd=None, stdin=None):
        return subprocess.Popen(
            self._command(argv, cwd),
            stdin=stdin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def put_file(self, source, destination):
        destination = os.fspath(destination)
        parent = posixpath.dirname(destination)
        if parent:
            self.exec(("mkdir", "-p", "--", parent))

        # Stream through ssh instead of relying on scp path parsing/quoting.
        command = self._command(("sh", "-c", 'cat > "$1"', "sh", destination))
        with open(source, "rb") as source_file:
            subprocess.run(
                command,
                stdin=source_file,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )

    def get_file(self, source, destination):
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        command = self._command(("cat", "--", os.fspath(source)))
        with destination.open("wb") as destination_file:
            subprocess.run(
                command,
                stdout=destination_file,
                stderr=subprocess.PIPE,
                check=True,
            )

    def read_file(self, path):
        return self.exec(("cat", "--", os.fspath(path))).stdout
