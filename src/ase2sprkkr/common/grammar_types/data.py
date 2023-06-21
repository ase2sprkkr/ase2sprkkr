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
    def __init__(self, *args, delimiter=None, written_shape=None, item_format=None, indented=False, **kwargs):
        """
        Parameters
        ----------

        delimiter
          None - default behavior.
          int  - the number will take given fixed number of chars

        written_shape
          resize to given shape before writing

        indented
          If the file has the following structure::

             .......................................
                  ....rest of the splitted line.....
                  ....rest of the splitted line.
             ....... The second line................
                  ..... the rest of the ............
                  ............ second line ...
             ............

          Pass a tuple with two integers into this argument.
          The first number of tuple is the number of characters on on line max,
          longer lines will be splitted.
          The second number is the number of spaces placed on the begining of the
          new lines created by splitting the old.

        format
          Format string used to write, e.g.`%.4e`

        """
        self.delimiter=delimiter
        self.written_shape=written_shape
        self.item_format=item_format
        self.indented=indented
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
       np.savetxt(out, value, delimiter=delimiter, format=self.item_format)
       np.savetxt(io)
       out=io.getvalue()
       if indented:
          first = indented[0]
          nexts = first - indented[1]
          prefix = ' '*indented[1]
          def g():
              for i in '\n'.split(''):
                  yield i[:first]
                  s=firsts
                  l=len(i)
                  while s<l:
                    e=firsts+nexts
                    yield prefix+i[s:e]
                    s=e
          out = '\n'.join(g)
       return out


    is_the_same_value = staticmethod(compare_numpy_values)

    def _grammar(self, param_name=False):
         out = RestOfTheFile._grammar.copy()
         def parse(v):
             if self.indented:
                v=v.replace('\n'+' '*self.indented[1], '')
             v=np.genfromtxt( io.StringIO(v), delimiter=self.delimiter )
             return v

         out.setParseAction(
                lambda v: parse(v[0])
         )
         return out

    def copy_value(self, value):
        return copy.deepcopy(value)
