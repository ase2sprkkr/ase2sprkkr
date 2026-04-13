""" Classes, for opening ASE2SPRKKR in shell or Jupyter notebook """
import code
import codeop
import contextlib
import os
import sys
import ase2sprkkr
import re
from ase2sprkkr.gui.examples import Example
import subprocess
from pathlib import Path

try:
  import nbformat
  from nbformat.v4 import new_notebook, new_code_cell
  from nbconvert.preprocessors import ExecutePreprocessor
except ImportError:
  nbformat = None



def split_top_level_chunks(code: str):
    """
    Split Python code into top-level chunks.
    Only splits at blank lines that are preceded by a line with no indentation.
    """
    chunks = []
    current_chunk = []
    empty = False

    for line in code.splitlines():
        # Check if the current line is blank
        if line.strip() == "":
            empty=True
        elif empty:
            if not line.startswith((" ", "\t")) and current_chunk:
                chunks.append("\n".join(current_chunk))
                current_chunk = []
            empty=False
        current_chunk.append(line)

    if current_chunk:
        chunks.append("\n".join(current_chunk))

    return chunks


@contextlib.contextmanager
def chdir(path):
    if path is None:
        yield
        return
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def format_as_python_shell(code_str: str) -> str:
    """
    Format Python code as it would appear in a Python shell.
    Prepend >>> to new statements and ... to continuation lines.

    Handles multi-line blocks like classes, functions, loops, and conditionals.
    """
    lines = code_str.splitlines()
    formatted_lines = []
    buffer = []

    for line in lines:
        buffer.append(line)
        # Join buffer and check if it's a complete statement
        joined = "\n".join(buffer)
        try:
            compiled = codeop.compile_command(joined, symbol="exec")
        except (SyntaxError, ValueError):
            # Invalid code; treat as complete to avoid infinite loop
            compiled = True

        if compiled is None:
            # Incomplete -> next line is continuation
            continue
        else:
            # Complete statement; format lines
            for i, l in enumerate(buffer):
                prefix = ">>> " if i == 0 else "... "
                formatted_lines.append(f"{prefix}{l}")
            buffer = []  # reset buffer for next statement

    # If any remaining lines (rare), treat them as continuation
    for i, l in enumerate(buffer):
        prefix = ">>> " if i == 0 else "... "
        formatted_lines.append(f"{prefix}{l}")

    return "\n".join(formatted_lines)


class Shell:
    """ A base class for ASE2SPRKKR interactive computation using ase2sprkkr tool"""

    """ Prefix code, do it always """
    initialization = [
    """from ase2sprkkr import (
        SPRKKRAtoms,
        InputParameters,
        OutputFile,
        Potential
    )"""
    ]

    def __init__(self):
        self.code = Shell.initialization.copy()
        self.output_file_counter=""
        self.dir = None

    def change_working_dir(self, dir):
        self.dir = dir

    def add_output_file(self, file):
        of = 'out' + str(self.output_file_counter)
        self.output_file_counter = \
            self.output_file_counter + 1 if \
            self.output_file_counter else 1
        self.code.append(
            f"{of} = OutputFile.from_file({repr(str(file))})\n"
            f"{of}.plot()"
        )

    def run_example(self, example, dir):
        ex=Example.by_number(example)
        ex.copy(dir)

        content = ex.body_of_main()
        content = split_top_level_chunks(content)
        self.code += content

class CommandlineShell(Shell):
    """ Run the given code in a commandline shell """
    def save(self, filename):
        with open(filename, "w") as f:
            f.write("\n\n".join(self.code))


class Python(CommandlineShell):
    """ Run the given code in a Python shell """

    def open(self, temp=False):

        print("\n Opening python shell...\n")

        with chdir(self.dir):
            local_ns = {}
            for cmd in self.code:
                print(format_as_python_shell(cmd))
                exec(cmd, globals(), local_ns)
                print("\n")
            code.interact(local=local_ns)


class Pdb(CommandlineShell):
    """ Run the given code in a Python shell """
    def open(self, temp=False):

        import pdb
        import linecache

        print("\n Opening pdb session...\n")

        source = "\n".join(self.code)
        filename = "<user_code>"

        lines = source.splitlines(keepends=True)
        linecache.cache[filename] = (
            len(source),   # size
            None,          # mtime (None is fine)
            lines,         # list of lines WITH newline chars
            filename,
        )
        code = compile(source, filename, "exec")

        pdb.runctx(code, {}, {})


if nbformat is None:
  JupyterLab = None
else:
  class JupyterLab(Shell):
    """Run the given code in a JupyterLab notebook."""

    def __init__(self, run=False):
        self.saved = False
        self.run = run
        super().__init__()

    def save(self, filename):
        """Save the notebook to disk."""
        filename = Path(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)

        nb = new_notebook()
        for cmd in self.code:
            nb.cells.append(new_code_cell(cmd))

        if self.run:
            ep = ExecutePreprocessor(timeout=600)
            ep.preprocess(nb, {'metadata': {'path': self.dir}})

        with open(filename, "w", encoding="utf-8") as f:
            nbformat.write(nb, f)
        self.saved = filename

    def open(self):
        print("\n Opening JupyterLab... Press CTRL+C to stop the Jupyter computational kernel.\n")
        dr = Path(self.dir or '.')

        if not self.saved:
            self.save(dr / "notebook.ipynb")


        try:
            import jupyterlab
        except ImportError:
            raise Exception("Jupyter lab not installed, please install (e.g. using pip install jupyter-lab")

        with chdir(self.dir):
            cmd = ["jupyter-lab", "--NotebookApp.shutdown_no_activity_timeout=60", self.saved]
            # Launch JupyterLab
            try:
                proc=subprocess.Popen(
                    cmd,
                    cwd=self.dir,
                       stdin=None,  # inherit parent's stdin
                       #stdout=None, # inherit parent's stdout
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL,
                       #stderr=None, # inherit parent's stderr
                        start_new_session=False  # important: do not detach
                )
                proc.wait()
            except KeyboardInterrupt:
                print("Terminating the kernel")
                proc.terminate()
                try:
                  proc.wait(timeout=5)  # wait gracefully
                except subprocess.TimeoutExpired:
                  proc.kill()       # force kill
                  proc.wait()  # wait gracefully
"""
                start_new_session=True,
                stdout=subprocess.DEVNULL,
                stdin=subprocess.DEVNULL
            )"""
