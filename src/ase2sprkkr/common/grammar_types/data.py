""" This module contains special GrammarTypes used for large data in output files """

from .grammar_type import GrammarType, compare_numpy_values
from ..decorators import add_to_signature
import pyparsing as pp
import re
import io
import numpy as np
import copy

class RestOfTheFile(GrammarType):
    """ Match anything up to the end of the file """

    datatype = str
    datatype_name = 'string'

    _grammar = pp.Regex('.*$', re.M | re.S).setParseAction(lambda x:x[0])

    def grammar_name(self):
      return '<the rest of the file>'

class NumpyArray(GrammarType):
    """ Match anything up to the end of the file, as numpy array """

    array_access = True

    @add_to_signature(GrammarType.__init__)
    def __init__(self, *args, delimiter=None, written_shape=None, **kwargs):
        self.delimiter=delimiter
        self.written_shape=written_shape
        super().__init__(*args, **kwargs)

    def _validate(self, value, why='set'):
       return isinstance(value, np.ndarray)

    def convert(self, value):
       return np.asarray(value)

    def _string(self, value):
       out = io.StringIO()
       delimiter = self.delimiter
       if isinstance(delimiter, int):
           delimiter = ''
       if self.written_shape:
          out=out.reshape(self.written_shape)
       np.savetxt(out, value, delimiter=delimiter)
       return np.savetxt(io)

    is_the_same_value = staticmethod(compare_numpy_values)

    def _grammar(self, param_name=False):
         out = RestOfTheFile._grammar.copy()
         out.setParseAction(
                lambda v: np.genfromtxt( io.StringIO(v[0]), delimiter=self.delimiter )
         )
         return out

    def copy_value(self, value):
        return copy.deepcopy(value)
