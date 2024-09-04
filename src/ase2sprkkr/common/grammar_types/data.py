""" This module contains special GrammarTypes used for large data in output files """

from .grammar_type import GrammarType, compare_numpy_values
from ..decorators import add_to_signature, cached_property
import pyparsing as pp
import re
from ..grammar import SkipToRegex
import io
import numpy as np
import copy
import os
from typing import Union


class RestOfTheFile(GrammarType):
    """ Match anything up to the end of the file """

    datatype = str
    datatype_name = 'string'

    _grammar = pp.Regex('.*$', re.M | re.S).setParseAction(lambda x:x[0])
    _grammar.skipWhitespace=False

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


class RawData(GrammarType):
    """ Match anything up to the end of the file or to the given delimiter """

    @add_to_signature(GrammarType.__init__)
    def __init__(self, *args, lines=None, indented=False,
                              line_length=None,
                              ends_with:Union[str,re.Pattern]=None,
                              ends_with_str=None,
                              include_ends_with=False,
                              **kwargs):
        """
        Parameters
        ----------

        lines
          Number of lines to read. Can be given as string - then
          the value of the given option determines the number of lines.

        indented
If there are <n> spaces before data, pass n to this arg.

          If the file has the following structure:
          .. code-block:: text

             .......................................
                  ....rest of the splitted line.....
                  ....rest of the splitted line.
             ....... The second line................
                  ..... the rest of the ............
                  ............ second line ...
             ............

          Pass a tuple with two integers into this argument.
          The first number of tuple is the max. number of characters on a line,
          longer lines will be splitted.
          The second number is the number of spaces placed on the begining of the
          new lines created by splitting the old.

        line_length
          Wrap the lines longer than a given number

        ends_with
          The data ends with a given string.

        ends_with_str
          If ends_with is regex, print this on end of the data

        include_ends_with
          Include the ending delimiter to the data.

        **kwargs
          Any other arguments are passed to the :meth:`GrammarType constructor<GrammarType.__init__>`
        """
        self.ends_with = ends_with
        self.ends_with_str = self.ends_with if ends_with_str is None else ends_with_str
        self.indented=' ' * indented if isinstance(indented, int) else indented
        self.lines=lines
        self.line_length = line_length
        self.include_ends_with = include_ends_with
        self.remove_forward=None
        super().__init__(*args, **kwargs)

    def _n_lines_grammar(self, lines):
         """ return a grammar for n lines of text """
         out=pp.Regex(f"([^\n]*\n){{{lines-1}}}[^\n]*(?=\n|$)", re.S)
         out.leaveWhitespace()
         return out

    def _grammar(self, param_name=False):
        if self.lines:
            if isinstance(self.lines, int):
               out = self._n_lines_grammar(self.lines)
            else:
               out = self.forward
        elif self.ends_with:
            if isinstance(self.ends_with, re.Pattern):
                out=SkipToRegex(self.ends_with, include_pattern=self.include_ends_with)
            else:
                out=pp.SkipTo(pp.Suppress(self.ends_with), include=self.include_ends_with)
                out.setParseAction(lambda x: x[0])
        else:
            out = RestOfTheFile._grammar.copy()

        def parse(v):
            v=v[0]
            if self.indented:
                if isinstance(self.indented, tuple):
                  v=v.replace('\n' + ' ' * self.indented[1], '')
                else:
                  v=re.sub(f'(^|\n){self.indented}',r'\1', v)
            if self.line_length:
                v=re.sub(f'([^\n]{{{self.line_length}}}[^{self.written_delimiter}\n]*)\n',f'\\1{self.written_delimiter}', v)
            return v

        if self.indented or self.line_length:
            out.addParseAction(parse)
        return out

    def _string(self, val):
        out = str(val)
        if self.line_length:
            out=re.sub(f'([^\n]{{{self.line_length}}}[^{self.written_delimiter}\n]*){self.written_delimiter}','\\1\n', out)

        indented = self.indented
        if indented:
          if isinstance(indented, tuple):
            first = indented[0]
            nexts = first - indented[1]
            prefix = ' ' * indented[1]

            def g():
                for i in out.split('\n'):
                    yield i[:first]
                    s=first
                    ln=len(i)
                    while s<ln:
                      e=s + nexts
                      yield prefix + i[s:e]
                      s=e
            out = '\n'.join(g())
          else:
            out=re.sub('(^|\n(?!$))',r'\1' + indented, out)
        if self.ends_with_str and self.include_ends_with:
            out+=self.ends_with_str
        return out

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

    def convert(self, val):
        return str(val)


class NumpyArray(RawData):
    """ Match anything up to the end of the file or to the given delimiter, as numpy array """

    array_access = True

    @add_to_signature(GrammarType.__init__)
    def __init__(self, *args, delimiter=None, shape=None, written_shape=None,
                              item_format='% .18e', dtype=None, no_newline_at_end=True,
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

        item_format
          Output format of the array (just for writing).

        dtype
          Type of the resulting data. Pass ``'line'`` to get array of whole lines

        **kwargs
          Any other arguments are passed to the :meth:`GrammarType constructor<GrammarType.__init__>`
        """
        self.delimiter=delimiter
        self.written_delimiter = delimiter or ' '
        self.written_shape=written_shape
        self.item_format=item_format
        self.shape=shape
        self.no_newline_at_end=no_newline_at_end
        self.dtype=dtype
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
          value=value.reshape(self.written_shape)
       np.savetxt(out, value, delimiter=self.written_delimiter, fmt=self.item_format)
       if self.no_newline_at_end:
           out.seek(out.tell() - 1 , os.SEEK_SET)
           out.truncate()
       out=out.getvalue()
       return super()._string(out)

    is_the_same_value = staticmethod(compare_numpy_values)

    def _grammar(self, param_name=False):
         grammar = super()._grammar(param_name)

         def parse(v):
             v=v[0]
             if self.dtype=='line':
                v=np.array([ i.rstrip() for i in v.split('\n')], dtype=object)
             else:
                v=np.genfromtxt( io.StringIO(v), delimiter=self.delimiter, dtype=self.dtype )
             if self.shape:
                v.shape=self.shape
             return v

         grammar.addParseAction(
             parse
         )
         return grammar

    def copy_value(self, value):
         return copy.deepcopy(value)
