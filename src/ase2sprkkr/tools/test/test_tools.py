import shlex
import subprocess
from pathlib import Path

import pytest


TOOLS_DIR = Path(__file__).parent.parent
RUN_SCRIPTS = sorted((TOOLS_DIR / "examples").glob("*/run.sh"))


def test_tool_examples_exist():
    assert RUN_SCRIPTS, "No CLI example run.sh files were found"


@pytest.mark.parametrize("script", RUN_SCRIPTS, ids=lambda path: path.parent.name)
def test_tool_example(script):
    lines = script.read_text().splitlines()
    assert len(lines) >= 2, f"{script} does not contain a command"

    command = [argument for argument in shlex.split(lines[1]) if argument != "-a"]
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
