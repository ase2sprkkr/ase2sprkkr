import contextlib
import tempfile
import os


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

    def __str__(self):
        return str(self.path) if self.path is not False else '<tempdir>'

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
