""" This module describe common header, that appears in the output
files of SPRKKR """

from ...sprkkr.configuration import ConfigurationFile, ConfigurationValue
import pyparsing as pp
import numpy as np
import pkgutil
import sys
from ...common.decorators import cached_class_property
from ...common.grammar_types.data import RestOfTheFile
import io

class UnknownDataValue(ConfigurationValue):

  def as_array(self):
      out = self()
      if not isinstance(out, np.ndarray):
          out = np.genfromtxt( io.StringIO(val) )
      return out


class OutputFile(ConfigurationFile):
  """ Objects of this class holds datas of an output file """

  @cached_class_property
  def unknown_output_file_definition(cls):
      """ a definition of unwnown (not yet known) output file that can hold any data in the rest
      of the file """
      V = output_files_definitions.OutputFileValueDefinition
      return output_files_definitions.create_output_file_definition(
          V('KEYWORD', str),
          [ V('DATA', RestOfTheFile(), name_in_grammar=False, result_class=UnknownDataValue) ]
      )

  @cached_class_property
  def definitions(cls):
      """ Return all known definitions of the SPR-KKR output files """
      out = {}
      for imp, module, ispackage in pkgutil.iter_modules(path=__path__, prefix=__name__+'.'):
           __import__(module)
           mod = sys.modules[module]
           out[mod.__name__.rsplit('.',1)[1]] = mod
      return out

  @classmethod
  def from_file(cls, filename):
      for i in cls.definitions.values():
          #try:
             out = i.definition.read_from_file(filename)
             return out
          #except Exception as e:
             print(e)
             breakpoint()
             pass
      try:
          return cls.unknown_output_file_definition.read_from_file(filename)
      except pp.ParseBaseException as e:
          raise Exception(f'Can not parse file: {filename}') from e

#at last, import this file that need this module
from .. import output_files_definitions
