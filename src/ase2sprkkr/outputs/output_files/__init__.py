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
  def from_file(cls, filename, first_try=None, try_only=None):
      if first_try is None and not try_only:
         fname = filename
         if hasattr(filename, 'name'):
             fname = filename.name
         if isinstance(filename, str):
             first_try = filename.rsplit('.',1)[1].lower()
         else:
             first_try = ''
      if isinstance(first_try,str):
          first_try=[ first_try ]

      if first_try:
         out = cls.from_file(filename, first_try=False, try_only=first_try)
         if out: return out

      for ext, i in cls.definitions.items():
          if try_only and ext not in try_only:
             continue
          if first_try and ext in first_try:
             continue
          try:
             out = i.definition.read_from_file(filename)
             return out
          except Exception as e:
             pass
      try:
          return cls.unknown_output_file_definition.read_from_file(filename)
      except pp.ParseBaseException as e:
          raise Exception(f'Can not parse file: {filename}') from e

class CommonOutputFile(ConfigurationFile):

    def n_atoms(self):
        return len(self.ORBITALS)

    def n_types(self):
        return len(self.TYPES)

    def site_type_index(self, type):
        if isinstance(type, int):
            return type
        for i,t in enumerate(self.TYPES):
            if t[0] == type:
               return i
        raise ValueError(f'There is no {type} atom in the output file')

    def n_orbitals(self, type):
        type = self.site_type_index(type)
        return self.ORBITALS[self.TYPES[type]['IQAT'][0]]

class Arithmetic(ConfigurationFile):

    def _check_arithmetic(self, other):
        pass

    def __add__(self, other):
        out = self.copy(copy_values=True)
        out+=other
        return out

    def __sub__(self, other):
        out = self.copy(copy_values=True)
        out+=other
        return out

    def __mul__(self, other):
        out = self.copy(copy_values=True)
        out*=other
        return out

    def __div__(self, other):
        out = self.copy(copy_values=True)
        out/=other
        return out

    def __rmul__(self, other):
        self._check_arithmetic(other)
        out = self.copy()
        out*=other
        return out

    def _do_arithmetic(self, func, other):
        """ Run given function for all "summable/subtractable/etc... data"""
        for val, selector in self._arithmetic_values:
            getattr(self[val]()[selector],func)(other[val]()[selector])

    def __iadd__(self, other):
        self._check_arithmetic(other)
        self._do_arithmetic('__iadd__', other)
        return self

    def __isub__(self, other):
        self._check_arithmetic(other)
        self._do_arithmetic('__isub__', other)
        return self

    def __imul__(self, other):
        self._do_arithmetic('__imul__', other)
        return self

    def __idiv__(self, other):
        self._do_arithmetic('__idiv__', other)
        return self


#at last, import this file that need this module
from .. import output_files_definitions
