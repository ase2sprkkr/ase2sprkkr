"""Minimal direct-process supervisor used by the stateless remote service.

It owns no ASE2SPRKKR lifecycle state.  Its only durable output beyond the
calculation streams is the immutable process return code, which a bare POSIX
PID cannot provide after the process has been reaped.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from uuid import uuid4


_SPEC_ENV = "ASE2SPRKKR_REMOTE_PROCESS_SPEC"


def _write_return_code(path, return_code):
    path = Path(path)
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        temporary.write_text(str(int(return_code)))
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def main():
    raw = os.environ.pop(_SPEC_ENV, None)
    if raw is None:
        raise RuntimeError(f"Missing {_SPEC_ENV}")
    spec = json.loads(raw)

    opened = []
    stdin = subprocess.DEVNULL
    try:
        if spec.get("stdin"):
            stdin = open(spec["stdin"], "rb")
            opened.append(stdin)

        output = open(spec["output_file"], "wb")
        opened.append(output)
        if spec["error_file"] == spec["output_file"]:
            error = subprocess.STDOUT
        else:
            error = open(spec["error_file"], "wb")
            opened.append(error)

        process = subprocess.Popen(
            spec["argv"],
            cwd=spec["directory"],
            stdin=stdin,
            stdout=output,
            stderr=error,
            close_fds=True,
        )
        return_code = process.wait()
    finally:
        for file in opened:
            file.close()

    _write_return_code(spec["return_code_file"], return_code)
    return return_code


if __name__ == "__main__":
    sys.exit(main())
