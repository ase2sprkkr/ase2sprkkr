"""This module contains special GrammarTypes used for large data in output files"""

import pyparsing as pp
import re
import io
import math
from numbers import Integral
import numpy as np
import copy
import os
from typing import Union
from functools import partial

from .grammar_type import GrammarType, compare_numpy_values
from ..dependencies import DependentValue
from ..decorators import add_to_signature, cached_property
from ..grammar import SkipToRegex, Forward

def _shape_lines(shape, *, items_per_line):
    size = math.prod(shape)
    if any(x < 0 for x in shape):
        raise ValueError("Invalid shape")
    return (size + items_per_line - 1) // items_per_line


class RestOfTheFile(GrammarType):
    """Match anything up to the end of the file"""

    datatype = str
    datatype_name = "string"

    _grammar = pp.Regex(".*$", re.M | re.S).set_parse_action(lambda x: x[0])
    _grammar.skipWhitespace = False

    def grammar_name(self):
        return "<the rest of the file>"


class Prefixed(GrammarType):
    """This value consists from a few lines, each prefixed with a given prefix"""

    @add_to_signature(GrammarType.__init__, prepend=True)
    def __init__(self, data_prefix, allow_empty=True, *args, **kwargs):
        self.data_prefix = data_prefix
        self.allow_empty = allow_empty
        super().__init__(*args, **kwargs)

    @cached_property
    def _grammar(self):
        pref = re.escape(self.data_prefix)
        out = f"({pref}[^\n]*)(\n{pref}[^\n]*)*"
        if self.allow_empty:
            out = f"({out})?"
        return pp.Regex(out)

    def _string(self, value):
        return re.replace("^|\n", f"{self.data_prefix}\\1", value)


class RawData(GrammarType):
    """Match anything up to the end of the file or to the given delimiter"""

    @add_to_signature(GrammarType.__init__)
    def __init__(
        self,
        *args,
        lines=None,
        indented=False,
        line_length=None,
        ends_with: Union[str, re.Pattern] = None,
        ends_with_str=None,
        include_ends_with=False,
        **kwargs,
    ):
        """
        Parameters
        ----------

        lines
          Number of lines to read. A string denotes the exact path of an
          option containing the count. A callable derives the count from the
          options named by its parameters; :class:`DependentValue` supports
          explicit paths such as ``..N``.

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
        self.indented = " " * indented if isinstance(indented, int) else indented
        self.lines = lines
        self.line_length = line_length
        self.include_ends_with = include_ends_with
        self._line_dependency = (
            DependentValue.create(lines)
            if lines is not None and not isinstance(lines, Integral)
            else None
        )
        self.forward = None
        super().__init__(*args, **kwargs)

    def _n_lines_grammar(self, lines):
        """return a grammar for n lines of text"""
        lines = int(lines)
        if lines < 0:
            raise ValueError("The number of lines cannot be negative")
        if lines == 0:
            return pp.Empty().set_parse_action(lambda: [""])
        out = pp.Regex(f"([^\n]*\n){{{lines - 1}}}[^\n]*(?=\n|$)", re.S)
        out.leave_whitespace()
        return out

    def _grammar(self, param_name=False):
        if self.lines is not None:
            if isinstance(self.lines, Integral):
                out = self._n_lines_grammar(self.lines)
            else:
                if self.forward is None or not self._line_dependency.hooks:
                    missing = ", ".join(self._line_dependency.paths)
                    raise KeyError(
                        f"No line-count dependencies found: {missing}"
                    )
                out = self.forward
        elif self.ends_with:
            if isinstance(self.ends_with, re.Pattern):
                out = SkipToRegex(self.ends_with, include_pattern=self.include_ends_with)
            else:
                out = pp.SkipTo(pp.Suppress(self.ends_with), include=self.include_ends_with)
                out.set_parse_action(lambda x: x[0])
        else:
            out = RestOfTheFile._grammar.copy()

        def parse(v):
            v = v[0]
            if self.indented:
                if isinstance(self.indented, tuple):
                    v = v.replace("\n" + " " * self.indented[1], "")
                else:
                    v = re.sub(f"(^|\n){self.indented}", r"\1", v)
            if self.line_length:
                v = re.sub(
                    f"([^\n]{{{self.line_length}}}[^{self.written_delimiter}\n]*)\n", f"\\1{self.written_delimiter}", v
                )
            return v

        if self.indented or self.line_length:
            out.add_parse_action(parse)
        return out

    def _string(self, val):
        out = str(val)
        if self.line_length:
            out = re.sub(
                f"([^\n]{{{self.line_length}}}[^{self.written_delimiter}\n]*){self.written_delimiter}", "\\1\n", out
            )

        indented = self.indented
        if indented:
            if isinstance(indented, tuple):
                first = indented[0]
                nexts = first - indented[1]
                prefix = " " * indented[1]

                def g():
                    for i in out.split("\n"):
                        yield i[:first]
                        s = first
                        ln = len(i)
                        while s < ln:
                            e = s + nexts
                            yield prefix + i[s:e]
                            s = e

                out = "\n".join(g())
            else:
                out = re.sub("(^|\n(?!$))", r"\1" + indented, out)
        if self.ends_with_str and self.include_ends_with:
            out += self.ends_with_str
        return out

    def added_to_container(self, container):
        if self._line_dependency is not None:
            self.forward = Forward()
            self._line_dependency.bind(container, self._set_number_of_lines)
        super().added_to_container(container)

    def _set_number_of_lines(self, lines):
        self.forward << self._n_lines_grammar(lines)

    def __del__(self):
        if self._line_dependency is not None:
            sef._line_dependency.unbind()

    def copy(self):
        out = super().copy()
        if self._line_dependency is not None:
            out._line_dependency = self._line_dependency.copy()
            out.forward = None
        return out

    def convert(self, val):
        return str(val)


class NumpyArray(RawData):
    """Match anything up to the end of the file or to the given delimiter, as numpy array"""

    array_access = True

    @add_to_signature(GrammarType.__init__)
    def __init__(
        self,
        *args,
        delimiter=None,
        written_delimiter=None,
        shape=None,
        written_shape=None,
        items_per_line=None,
        item_format="% .18e",
        dtype=None,
        dtypes=None,
        no_newline_at_end=True,
        **kwargs,
    ):
        """
        Parameters
        ----------

        delimiter
          None - default behavior.
          int  - the number will take given fixed number of chars

        written_delimiter
          Delimiter used only for writing. By default, use ``delimiter`` or
          one space when the parsing delimiter is not specified.

        shape
          Resize to given shape after read. A string, callable or
          :class:`DependentValue` derives the shape from other options.

        written_shape
          Resize to given shape before writing

        items_per_line
          Flatten the array and wrap its output after this many scalar items.
          If ``shape`` is dynamic, its value also determines the number of
          input lines.

        item_format
          Output format of the array (just for writing).

        dtype
          Type of the resulting data. Pass ``'line'`` to get array of whole lines

        dtypes
          More dtypes can be given. Then, the first, that match the data, is used

        **kwargs
          Any other arguments are passed to the :meth:`GrammarType constructor<GrammarType.__init__>`
        """
        self.delimiter = delimiter
        self.written_delimiter = (
            ("" if isinstance(delimiter, int) else delimiter or " ")
            if written_delimiter is None
            else written_delimiter
        )
        self.written_shape = written_shape
        self.items_per_line = items_per_line
        self.item_format = item_format
        self.shape = shape
        self._shape_dependency = (
            DependentValue.create(shape)
            if isinstance(shape, (str, DependentValue)) or callable(shape)
            else None
        )
        self._parsed_shape = None
        self.no_newline_at_end = no_newline_at_end
        if dtypes is None:
            if dtype == "line":
                dtypes = dtype
            else:
                dtypes = [dtype]
        else:
            if dtype is not None:
                raise ValueError("Use either dtype or dtypes, but not both")
        self.dtypes = dtypes
        lines = kwargs.pop("lines", None)
        if items_per_line is not None:
            if not isinstance(items_per_line, int) or items_per_line <= 0:
                raise ValueError("items_per_line has to be a positive integer")
            if lines is None:
                if self._shape_dependency is not None:
                    lines = self._shape_dependency.mapped(
                        partial(_shape_lines, items_per_line=items_per_line)
                    )
                elif shape is not None and -1 not in shape:
                    lines = _shape_lines(shape, items_per_line=items_per_line)
        super().__init__(*args, lines=lines, **kwargs)

    def _validate(self, value, why="set"):
        return isinstance(value, np.ndarray)

    def validate(
        self,
        value,
        param_name="<Unknown>",
        why="set",
        option=None,
    ):
        """Validate the array and, when possible, its configured shape."""
        super().validate(
            value,
            param_name=param_name,
            why=why,
            option=option,
        )
        if option is None:
            return True

        expected = self.expected_shape(option)
        if expected is None:
            return True

        actual = value.shape
        matches = len(actual) == len(expected) and all(
            wanted == -1 or current == wanted
            for current, wanted in zip(actual, expected)
        )
        if not matches:
            self._valueError(
                value,
                f"array has shape {actual}, expected {expected}",
                option._get_path,
            )
        return True

    def convert(self, value):
        return np.asarray(value)

    def _string(self, value):
        value = np.asarray(value)
        if self.written_shape:
            value = value.reshape(self.written_shape)
        if self.items_per_line is not None:
            values = value.reshape(-1)
            lines = []
            for start in range(0, len(values), self.items_per_line):
                chunk = values[start : start + self.items_per_line]
                lines.append(
                    self.written_delimiter.join(
                        self.item_format % item for item in chunk
                    )
                )
            output = "\n".join(lines)
            if lines and not self.no_newline_at_end:
                output += "\n"
            return super()._string(output)

        out = io.StringIO()
        np.savetxt(out, value, delimiter=self.written_delimiter, fmt=self.item_format)
        if self.no_newline_at_end and out.tell():
            out.seek(out.tell() - 1, os.SEEK_SET)
            out.truncate()
        out = out.getvalue()
        return super()._string(out)

    is_the_same_value = staticmethod(compare_numpy_values)

    def _grammar(self, param_name=False):
        grammar = super()._grammar(param_name)

        def parse(v):
            v = v[0]
            if self.dtypes == "line":
                v = np.array([i.rstrip() for i in v.split("\n")], dtype=object)
            else:
                if self.items_per_line is not None:
                    if isinstance(self.delimiter, int):
                        shape = (
                            self._parsed_shape
                            if self._shape_dependency is not None
                            else self.shape
                        )
                        if shape and all(i >= 0 for i in shape):
                            remaining = math.prod(shape)
                            lines = []
                            for line in v.splitlines():
                                count = min(self.items_per_line, remaining)
                                # Pyparsing may skip padding before the first
                                # fixed-width field of a value. Restore it
                                # from the known record width before parsing.
                                lines.append(line.rjust(count * self.delimiter))
                                remaining -= count
                            v = "".join(lines)
                        else:
                            v = v.replace("\n", "")
                    else:
                        v = v.replace("\n", " ")
                last_error = None
                for dt in self.dtypes:
                    try:
                        v = np.genfromtxt(
                              io.StringIO(v),
                              delimiter=self.delimiter,
                              dtype=dt,
                        )
                        break
                    except Exception as exc:
                        last_error = exc
                else:
                    if not last_error:
                        raise ValueError("No dtype specified")
                    raise last_error
            shape = (
                self._parsed_shape
                if self._shape_dependency is not None
                else self.shape
            )
            if shape:
                v.shape = tuple(shape)
            return v

        grammar.add_parse_action(parse)
        return grammar

    def copy_value(self, value):
        return copy.deepcopy(value)

    def added_to_container(self, container):
        if self._shape_dependency is not None:
            self._shape_dependency.bind(container, self._set_parsed_shape)
        super().added_to_container(container)

    def __del__(self):
        dependency = getattr(self, "_shape_dependency", None)
        if dependency is not None:
            dependency.unbind()
        super().__del__()

    def _set_parsed_shape(self, shape):
        self._parsed_shape = tuple(int(item) for item in shape)

    def expected_shape(self, item):
        """Return the shape expected for a runtime option."""

        if self._shape_dependency is not None:
            shape = self._shape_dependency.runtime_value(item)
        else:
            shape = self.shape
        return None if shape is None else tuple(int(value) for value in shape)

    def copy(self):
        out = super().copy()
        if self._shape_dependency is not None:
            out._shape_dependency = self._shape_dependency.copy()
            out._parsed_shape = None
        return out
