import contextlib
import tempfile
import os

class Directory:

    def new(dir, default):
        if isinstance(dir, Directory):
            return dir
        return Directory(dir, default)

    def __init__(self, dir, default):
        if dir is None:
           dir = default
        if dir is None:
           raise ValueError('No directory has been specified')

        self.enters = 0
        self.dir = dir
        if self.dir == '.':
            self.chdir=contextlib.suppress

    def __str__(self):
        return str(self.dir) if self.dir is not False else self.path

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
                self.path = self.dir

            if self.handler:
               self.handler.__enter__()
        return self

    def __exit__(self, type, value, traceback):
         self.enters -= 1
         if self.enters == 0 and self.handler:
           self.handler.__exit__(type, value, traceback)

    @contextlib.contextmanager
    def chdir(self):
        try:
          cwd=os.getcwd()
          os.chdir(self.path)
          yield self.path
        finally:
          os.chdir(cwd)
