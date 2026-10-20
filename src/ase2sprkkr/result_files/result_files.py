"""Common parsing support for structured SPR-KKR result files."""

from ..sprkkr.configuration import ConfigurationFile, ConfigurationValue
from ..potentials.potentials import Potential
import pyparsing as pp
import numpy as np
from pathlib import Path
from ..common.decorators import cached_class_property, cached_property
from ..common.grammar_types.data import RestOfTheFile
import io
import os
import re
from numbers import Number


class UnknownDataValue(ConfigurationValue):
    def as_array(self):
        out = self()
        if not isinstance(out, np.ndarray):
            out = np.genfromtxt(io.StringIO(out))
        return out


class ResultFile(ConfigurationFile):
    """Representation of a structured SPR-KKR result file."""

    plot_parameters = {}
    additional_actions = ()

    def __init__(self, definition=None, container=None, potential=None):
        super().__init__(definition, container)
        self._potential = None
        self._potential_filename = None
        if potential is not None:
            self.potential = potential

    def read_from_file(
        self,
        file,
        clear_first=True,
        allow_dangerous=False,
        potential=None,
    ):
        """Read result data and optionally associate its potential."""

        if potential is not None:
            self.potential = potential
        return super().read_from_file(
            file,
            clear_first=clear_first,
            allow_dangerous=allow_dangerous,
        )

    @classmethod
    def can_be_plotted(cls):
        return callable(getattr(cls, "plot", None))

    def set_potential_filename(self, potential_filename):
        self.__dict__["potential_filename"] = potential_filename
        self.__dict__.pop("potential", None)

    @cached_property
    def potential_filename(self):
        if not getattr(self, "_filename", None):
            raise ValueError("The result file has no filename. Set the potential filename manually.")

        path = Path(self._filename)
        directory = path.parent
        stem = path.stem

        while True:
            for suffix in (".pot_new", ".pot"):
                candidate = Path(directory) / (stem + suffix)
                if candidate.is_file():
                    return str(candidate)
            stem = stem.rsplit("_", 1)
            if len(stem) == 1:
                break
            stem = stem[0]

        raise ValueError(
            f"No potential file found for {self._filename}. Set it manually using set_potential_filename()."
        )

    @potential_filename.setter
    def potential_filename(self, value):
        self.set_potential_filename(value)

    @cached_property
    def potential(self):
        return Potential.from_file(self.potential_filename)

    @potential.setter
    def potential(self, value):
        if isinstance(value, Potential):
            self.__dict__["potential"] = value
            if hasattr(value, "_filename"):
                self.__dict__["potential_filename"] = value._filename
            else:
                self.__dict__.pop("potential_filename", None)
        else:
            self.__dict__["potential_filename"] = value
            self.__dict__.pop("potential", None)

    @property
    def atoms(self):
        return self.potential.atoms

    @cached_class_property
    def unknown_result_file_definition(cls):
        """Definition used for an unknown result file that stores the remaining data verbatim."""
        V = result_files_definitions.ResultFileValueDefinition
        return result_files_definitions.create_result_file_definition(
            V("KEYWORD", str), [V("DATA", RestOfTheFile(), name_in_grammar=False, result_class=UnknownDataValue)]
        )

    @cached_class_property
    def definitions(cls):
        """Return all known definitions of the SPR-KKR result files."""
        # Discovery imports individual definition modules. Defer construction
        # until ResultFile and result_files_definitions are fully initialized so
        # definition modules can safely import this class.
        from .definition_registry import result_file_definition_registry

        return result_file_definition_registry.by_extension()

    @classmethod
    def definition(cls, self):
        return cls.definitions[self]

    @classmethod
    def from_file(
        cls,
        filename,
        first_try=None,
        try_only=None,
        unknown=None,
        potential=None,
    ):
        """
        Read SPRKKR result file (DOS, BSF....). The type of content of the
        result file is guessed from the content, however you can get hint what
        to try first or force to read only certain file type(s).

        Parameters
        ----------
        filename
          The file to read

        first_try
          Result file type(s) to try first. If omitted, the type is guessed
          from the file extension.

        try_only
          Result file type(s) allowed to be read.
          None means any known result file type.

        unknown
          If True, fall back to an unknown result file: its header is parsed and
          the remaining content is stored as text.
          If False, raise an exception if no known (and allowed) file type is
          recognized.
          None means True if try_only is None, False otherwise.

        potential
          A :class:`Potential` or its filename. If omitted, it is looked up
          lazily from the result filename when it is needed.
        """
        fname = filename
        if hasattr(filename, "name"):
            fname = filename.name
        elif isinstance(filename, str):
            if not os.path.exists(filename):
                raise ValueError(f"File {filename} to be read does not exist")

        if first_try is None and not try_only:
            if isinstance(fname, str):
                try:
                    first_try = fname.rsplit(".", 1)[1].lower()
                except IndexError:
                    first_try = None
                special = re.match(r"^.+_([^/\\]+\.[^/\\]+)$", fname)
                if special:
                    first_try = [special.groups(1)[0], first_try]
            else:
                first_try = ""

        first = None
        if first_try:
            out = None
            try:
                out = cls.from_file(
                    filename,
                    first_try=False,
                    try_only=first_try,
                    potential=potential,
                )
            except Exception as e:
                first = e
            if out:
                return out
        else:
            last = None

        if try_only:
            if isinstance(try_only, str):
                try_only = [try_only]
            for i in try_only:
                if i in cls.definitions:
                    try:
                        return cls.definitions[i].read_from_file(
                            filename, potential=potential
                        )
                    except Exception as e:
                        last = e
        else:
            for ext, i in cls.definitions.items():
                if first_try and ext in first_try:
                    continue
                try:
                    return i.read_from_file(filename, potential=potential)
                except Exception as e:
                    last = e
        if unknown is None:
            unknown = try_only is None
        if unknown:
            try:
                return cls.unknown_result_file_definition.read_from_file(
                    filename, potential=potential
                )
            except pp.ParseBaseException as e:
                raise Exception(f"Can not parse file: {filename}") from e
        raise (first or last or ValueError(f"File {filename} is not recognized as any known file type"))


class CommonResultFile(ResultFile):
    def n_atoms(self):
        return len(self.ORBITALS)

    def n_types(self):
        return len(self.TYPES)

    def site_type_index(self, type):
        if isinstance(type, int):
            return type
        for i, t in enumerate(self.TYPES):
            if t["TXT_T"] == type:
                return i
        raise ValueError(f"There is no {type} atom in the result file")

    def n_orbitals(self, type):
        type = self.site_type_index(type)
        return self.ORBITALS[self.TYPES[type]["IQAT"][0]]["NLQ"]


class Arithmetic:
    def _check_arithmetic(self, other):
        if isinstance(other, Number):
            return
        if not isinstance(other, type(self)):
            raise TypeError(
                f"Cannot combine {type(self).__name__} with {type(other).__name__}"
            )
        check = getattr(self, "_assert_arithmetic", None)
        if check is not None:
            check(other)

    def __add__(self, other):
        out = self.copy(copy_values=True)
        out += other
        return out

    def __sub__(self, other):
        out = self.copy(copy_values=True)
        out -= other
        return out

    def __mul__(self, other):
        out = self.copy(copy_values=True)
        out *= other
        return out

    def __truediv__(self, other):
        out = self.copy(copy_values=True)
        out /= other
        return out

    def __rmul__(self, other):
        out = self.copy(copy_values=True)
        out *= other
        return out

    def _do_arithmetic(self, func, other):
        """Run given function for all "summable/subtractable/etc... data"""
        vals = self._arithmetic_values
        if callable(vals):
            vals = vals()
        for val, selector in vals:
            operand = other if isinstance(other, Number) else other[val]()[selector]
            getattr(self[val]()[selector], func)(operand)

    def __iadd__(self, other):
        self._check_arithmetic(other)
        self._do_arithmetic("__iadd__", other)
        return self

    def __isub__(self, other):
        self._check_arithmetic(other)
        self._do_arithmetic("__isub__", other)
        return self

    def __imul__(self, other):
        self._check_arithmetic(other)
        self._do_arithmetic("__imul__", other)
        return self

    def __itruediv__(self, other):
        self._check_arithmetic(other)
        self._do_arithmetic("__itruediv__", other)
        return self


from . import result_files_definitions  # NOQA: E402
