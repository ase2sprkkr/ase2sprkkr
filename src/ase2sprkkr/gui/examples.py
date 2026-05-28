"""Function related to ASE2SPRKKR examples"""

from pathlib import Path
import ast
import textwrap
import re
import shutil

from ase2sprkkr.common.decorators import cached_property
import ase2sprkkr


def examples_dir():
    """Base dir for all examples"""
    return Path(ase2sprkkr.__file__).parent / "examples"


class Example:
    @staticmethod
    def by_number(example: int):
        """Directory for the given example"""
        # Find the matching example directory
        prefix = f"A{example:02d}_"
        for d in examples_dir().iterdir():
            if d.is_dir() and d.name.startswith(prefix):
                return Example(d)
        raise FileNotFoundError(f"No example directory matching {prefix}* found.")

    def __init__(self, path: Path):
        self.dir = path

    def copy(self, dest_dir: str):
        """
        Copy contents of ase2sprkkr.examples.A{n:02d}_* directory into dest_dir
        and return the path to the .py script inside that example directory
        (excluding __init__.py).

        :param n: Example number (integer)
        :param dest_dir: Destination directory where files will be copied
        :return: Path to the .py example script
        """
        # Create destination directory
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)

        # Copy contents of example dir → destination
        for item in self.dir.iterdir():
            dst = dest_dir / item.name
            if item.name.startswith("__"):
                continue
            if item.is_dir():
                shutil.copytree(item, dst, dirs_exist_ok=True)
            else:
                shutil.copy2(item, dst)

    @property
    def name(self):
        return self.dir.name

    @cached_property
    def main_script(self):
        """Find Python script inside example directory (excluding __init__.py).
        Now we suspose, that the first such one is the main.
        """
        for p in self.dir.glob("*.py"):
            if p.name != "__init__.py":
                return p

        raise FileNotFoundError("No .py script found inside example directory.")

    @cached_property
    def docstring(self):
        # Read the file content
        source = self.source()
        # Parse the AST
        try:
            tree = ast.parse(source, filename=str(self.main_script))
            # Return the module-level docstring
            return ast.get_docstring(tree)
        except SyntaxError:
            return None  # in case the file is invalid Python

    @property
    def short_docstring(self):

        def first_sentence(text: str) -> str:
            match = re.match(r"^\s*([^.]*\.)", text)
            if match:
                return match.group(1)
            return text.strip()

        if self.docstring is None:
            return None
        out = first_sentence(self.docstring)
        out = re.sub(r"\s*\n\s*", " ", out)
        return out

    def source(self):
        return self.main_script.read_text(encoding="utf-8")

    def body_of_main(self):
        pattern = re.compile(
            r"^def\s+main\s*\([^)]*\)\s*:\s*\n"  # match 'def main(...):'
            r"((?:[ \t]+.*\n|\n)+)",  # capture indented block
            re.MULTILINE,
        )
        code = self.source()

        match = pattern.search(code)
        if not match:
            return None

        body = match.group(1)
        if not body:
            return code

        # Dedent the body
        return textwrap.dedent(body)

    def satisfy_regex(self, regex):
        regex = re.compile(regex)
        return regex.search(self.name) or regex.search(self.source())


def list_of_examples(regex=None):
    """
    Return list of tuples: (subdir_name, docstring of example_main_script)
    """
    base_dir = examples_dir()
    pattern = re.compile(r"^A\d{2}_")  # matches A{number*2}_
    out = [
        Example(d) for d in base_dir.iterdir() if d.is_dir() and pattern.match(d.name)
    ]
    if regex:
        regex = re.compile(regex)
        out = [i for i in out if i.satisfy_regex(regex)]

    return out
