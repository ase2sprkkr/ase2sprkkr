import contextlib
import tempfile
import os
from io import BufferedIOBase

_RAISE = object()

def filename_from_file(file, default=_RAISE, realpath=False):
    """Return a filesystem name for a path or a file backed by one.

    If the object has no filename, return ``default`` or raise ``TypeError``
    when no default was supplied.
    """
    if isinstance(file, BufferedIOBase):
        file = getattr(file, "name", file)
    try:
        out = os.fspath(file)
    except TypeError:
        if default is _RAISE:
            raise
        return default

    if realpath:
        out = os.path.realpath(out)

    return out

class Directory:
    @staticmethod
    def new(dir, default="."):
        if isinstance(dir, Directory):
            return dir
        return Directory(dir, default)

    def __init__(self, dir, default="."):
        if dir is None:
            dir = default
        if dir is None:
            raise ValueError("No directory has been specified")

        self.enters = 0
        self.dir = dir
        self.path = dir
        if self.dir == ".":
            self.chdir = contextlib.suppress

    def __fspath__(self):
        return self.path

    def __str__(self):
        return str(self.path) if self.path is not False else "<tempdir>"

    def __repr__(self):
        return f"<Directory {str(self)}>"

    def __enter__(self):
        self.enters += 1
        if self.enters == 1:
            if self.dir is False:
                self.handler = tempfile.TemporaryDirectory()
                self.path = self.handler.name
            else:
                self.handler = None

            if self.handler:
                self.handler.__enter__()
        return self

    def __exit__(self, type, value, traceback):
        self.enters -= 1
        if self.enters == 0 and self.handler:
            self.handler.__exit__(type, value, traceback)
            self.path = self.dir

    @contextlib.contextmanager
    def chdir(self):
        cwd = os.getcwd()
        if self.path:
            os.chdir(self.path)
            try:
                yield self.path
            finally:
                os.chdir(cwd)
        else:
            yield cwd
