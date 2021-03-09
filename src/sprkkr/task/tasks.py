from collections import OrderedDict
import os
import tempfile
from ..common.conf_containers import RootConfContainer

class Task(RootConfContainer):
  def __init__(self, definition, inputfile=None, outputfile=False):
      super().__init__(definition)
      self._inputfile = inputfile
      self._outputfile = inputfile if outputfile is False else outputfile

  @property
  def task_name(self):
      return self._definition.name

  def executable_params(self, directory=None, ):
      """
      Return
      ------
      ([command], stdin, stdout) params for process.Popen
      """

      def filename(self, filename, dir=None):
         if filename is None:
            return tempfile.TemporaryFile()
         if dir:
            filename = os.path.join(filename, dir)
         return open(filename, "w+b")

      return self._definition.command, filename(self._inputfile), filename(self._outputfile)


  def __str__(self):
      from io import StringIO
      s = StringIO()
      self.save_to_file(s)
      return s.getvalue()

  def set_option(name, value):
      for i in self._members:
          if name in i:
             self._members[i][name] = value
      else:
          raise KeyError("No option with name {} in any of the members".format(name))
