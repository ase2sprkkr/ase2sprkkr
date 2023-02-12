""" Context manager, that suppress output to the given descriptor.
The common usage is to suppress output of C/Fortran functions """

import os, sys
from typing import Union
from collections.abc import Iterable


class NoOutput:
  """ Context manager, that suppress output to the given descriptor.
  The common usage is to suppress output of C/Fortran functions """

  def __init__(self, fd:Union[int, Iterable]=1, suppress=True):
      """
      Parameters
      ----------
      fd
        File descriptor to suppress, the default is 1 (stdout).
      """
      self.fd=fd if isinstance(fd, Iterable) else [ fd ]
      self.suppress = suppress
      self.special = {
          1: sys.stdout,
          2: sys.stderr
      }

  def __enter__(self):
      if not self.suppress: return
      self.null = open(os.devnull, "wb")
      null_fd = self.null.fileno()
      def silence(fd):
          try:
            os.fsync(fd)
          except OSError:
            pass
          dup = os.dup(fd)
          s_dup = s_fd = None
          if fd in self.special:
             stream = self.special[fd]
             stream.flush()
             #pytest replaces the standard output
             if stream.fileno() != fd:
                s_fd = stream.fileno()
                try:
                  os.fsync(s_fd)
                except OSError:
                  pass
                s_dup = os.dup(s_fd)
                os.dup2(null_fd, s_fd)
          os.dup2(null_fd, fd)
          return dup, s_fd, s_dup

      self.dups = [ silence(fd) for fd in self.fd ]

  def __exit__(self, type, value, traceback):
      if not self.suppress: return

      itr=list(zip(self.fd, self.dups))
      itr.reverse()

      for fd, i in itr:
          dup, s_fd, s_dup = i
          if fd in self.special:
             stream = self.special[fd]
             stream.flush()
             if s_fd:
                try:
                    os.fsync(s_fd)
                except OSError:
                    pass
                os.dup2(s_dup, s_fd)
                os.close(s_dup)

          try:
              os.fsync(fd)
          except OSError:
              pass
          os.dup2(dup, fd)
          os.close(dup)
      self.null.close()
