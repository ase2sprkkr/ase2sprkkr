import os
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

    command = shlex.split(lines[1])
    if command[1:2] in (["show-structure"], ["view"]) and os.environ.get("DISPLAY"):
        command.extend(["--close_after", "0.1"])
    else:
        # Headless test workers cannot create Tk windows.
        command.append("--no-view")
    result = subprocess.run(
        command,
        cwd=script.parent,
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert result.returncode == 0, (
        f"Command {command!r} failed in {script.parent}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
