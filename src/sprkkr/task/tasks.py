from collections import OrderedDict
import os
import tempfile

class Task:
  def __init__(self, definition, inputfile=None, outputfile=False):
      self._sections = OrderedDict()
      self._definition = definition
      self._inputfile = inputfile
      self._outputfile = inputfile if outputfile is False else outputfile
      for d in definition.sections:
          self._sections[d] = Section(definition.sections[d])

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

  def __getattr__(self, name):
      return self._sections[name]

  def __getitem__(self, name):
      return self._sections[name]

  def __iter__(self):
      yield from self._sections.values()

  def clear(self):
      for i in self.sections:
          i.clear()

  def save_to_file(self, file):
      if not hasattr(file, 'write'):
         file = open(file, "w")
      for i in self:
         i.save_to_file(file)

  def read_from_file(self, file, clear_first=True):
      values = self._definitions.grammar().parseFile(file, True)
      if clear_first:
         self.clear()
      self.set(values)

  def set(self, values, allow_add=True):
      for i in values:
          if not i in self._sections:
             if not allow_add:
                raise KeyError("No section with name {} in task {}".format(i, self._defintion.task_name))
          self[i].set(values[i], allow_add)

  def __str__(self):
      from io import StringIO
      s = StringIO()
      self.save_to_file(s)
      return s.getvalue()

  def add_section(self, name):
      self._sections = CustomSection(self, name)

  def remove_section(self, name):
      if not isinstance(self._sections[name], CustomSection):
         raise KeyError("No custom section with name {} to remove".format(name))
      del self._sections[name]

  def set_option(name, value):
      for i in self._sections:
          if name in i:
             self._sections[i][name] = value
      else:
          raise KeyError("No option with name {} in any of the sections".format(name))
