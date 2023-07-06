""" This module contains special GrammarTypes used for large data in output files """

from .grammar_type import GrammarType, compare_numpy_values
from ..decorators import add_to_signature, cached_property
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

class Prefixed(GrammarType):
    """ This value consists from a few lines, each prefixed with a given prefix """

    @add_to_signature(GrammarType.__init__, prepend=True)
    def __init__(self, data_prefix, allow_empty=True, *args, **kwargs):
        self.data_prefix = data_prefix
        self.allow_empty=allow_empty
        super().__init__(*args, **kwargs)

    @cached_property
    def _grammar(self):
        pref = re.escape(self.data_prefix)
        out = f'({pref}[^\n]*)(\n{pref}[^\n]*)*'
        if self.allow_empty:
            out=f'({out})?'
        return pp.Regex(out)

    def _string(self, value):
        return re.replace('^|\n',f'{self.data_prefix}\\1', value)


class NumpyArray(GrammarType):
    """ Match anything up to the end of the file, as numpy array """

    array_access = True

    @add_to_signature(GrammarType.__init__)
    def __init__(self, *args, delimiter=None, shape=None, written_shape=None,
                              lines=None, item_format=None, indented=False,
                              dtype=None,
                              **kwargs):
        """
        Parameters
        ----------

        delimiter
          None - default behavior.
          int  - the number will take given fixed number of chars

        shape
          Resize to given shape after read

        written_shape
          Resize to given shape before writing

        lines
          Number of lines to read. Can be given as string - then
          the value of given option determine the number of lines.

        item_format
          Output format of the array (just for writing).

        indented
          If there are <n> spaces before data, pass n to this arg.

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

        dtype
          Type of the resulting data. Pass ``'line'`` to get array of whole lines

        **kwargs
          Any other arguments are passed to the :meth:`GrammarType constructor<GrammarType.__init__>`
        """
        self.delimiter=delimiter
        self.written_shape=written_shape
        self.item_format=item_format
        self.indented=' '*indented if isinstance(indented, int) else indented
        self.lines=lines
        self.shape=shape
        self.dtype=dtype
        self.remove_forward=None
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
       indented = self.indented
       if indented:
          if isinstance(indented, tuple):
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
          else:
            out=re.sub('(^|\n)',r'\1' + indented, out)
       return out


    is_the_same_value = staticmethod(compare_numpy_values)

    def _n_lines_grammar(self, lines):
         """ return a grammar for n lines of text """
         out=pp.Regex(f"([^\n]*\n){{{lines-1}}}[^\n]*(?=\n|$)", re.S)
         out.leaveWhitespace()
         #out.addParseAction(lambda x: breakpoint() or x)
         out=self._parse_numpy_array_grammar(out)
         return out

    def _parse_numpy_array_grammar(self, grammar):
         """ Change a parse action of given grammar such that it returns
         numpy array """
         def parse(v):
             if self.indented:
                if isinstance(self.indented, tuple):
                  v=v.replace('\n'+' '*self.indented[1], '')
                else:
                  v=re.sub(f'(^|\n){self.indented}',r'\1', v)
             if self.dtype=='line':
                v=np.array([ i.rstrip() for i in v.split('\n')], dtype=object)
             else:
                v=np.genfromtxt( io.StringIO(v), delimiter=self.delimiter, dtype=self.dtype )
             if self.shape:
                v.shape=self.shape
             return v

         grammar.setParseAction(
                lambda v: parse(v[0])
         )
         return grammar

    def _grammar(self, param_name=False):
         if self.lines:
             if isinstance(self.lines, int):
                return self._n_lines_grammar(self.lines)
             else:
                return self.forward
         else:
             out = RestOfTheFile._grammar.copy()
             return self._parse_numpy_array_grammar(out)

    def copy_value(self, value):
         return copy.deepcopy(value)

    def added_to_container(self, container):
        if not self.lines or isinstance(self.lines, int):
            return
        if self.remove_forward:
           self.remove_forward
        if container:
           self.forward=pp.Forward()
           obj = container[self.lines]
           def paction(parsed):
               self.forward << self._n_lines_grammar(parsed[0][1])
               return parsed
           hook = lambda grammar: grammar.addParseAction(paction)
           obj.add_grammar_hook(hook)
           self.remove_forward = lambda: obj.remove_grammar_hook(hook)
        else:
           self.remove_forward=None
        super().added_to_container(container)

    def __del__(self):
        if self.remove_forward:
           self.remove_forward
        self.remove_forward=None
