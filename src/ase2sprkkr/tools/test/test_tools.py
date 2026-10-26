import shlex
import subprocess
from pathlib import Path

import pytest


TOOLS_DIR = Path(__file__).parent.parent
RUN_SCRIPTS = sorted((TOOLS_DIR / "examples").glob("*/run.sh"))


@pytest.mark.parametrize(
    "script",
    RUN_SCRIPTS,
    ids=[path.parent.name for path in RUN_SCRIPTS],
)
def test_tool_example(script):
    lines = script.read_text().splitlines()
    assert len(lines) >= 2, f"{script} does not contain a command"

    # Keep the example command unchanged, but suppress its GUI in automated
    # runs; examples should still open the viewer when run normally.
    command = [*shlex.split(lines[1]), "--no-view"]
    result = subprocess.run(
        command,
        cwd=script.parent,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"Command {command!r} failed in {script.parent}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
