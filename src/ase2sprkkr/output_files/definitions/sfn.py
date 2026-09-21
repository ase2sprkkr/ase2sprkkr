"""Definition of the SPR-KKR shape-function (``.sfn``) file."""

from collections.abc import Sequence
import math

import numpy as np
import pyparsing as pp

from ..output_files import OutputFile
from ..output_files_definitions import (
    OutputFileDefinition,
    OutputFileSectionDefinition,
    OutputFileValueDefinition,
)
from ...common.configuration_definitions import gather
from ...common.generated_configuration_definitions import GeneratedValueDefinition
from ...common.grammar import TokenConverter, White, line_end
from ...common.grammar_types import GrammarType, Integer, Real, compare_numpy_values
from ...common.parsing_results import dict_from_parsed
from ...common.repetition import (
    Repeated,
    RepeatedItemGrammar,
    VariableRepeatedItemGrammar,
)


class _SFNValueDefinition(OutputFileValueDefinition):
    """An output value without the leading blank used by standard outputs."""

    prefix = ""


class _SFNSectionDefinition(OutputFileSectionDefinition):
    """A section whose surrounding grammar owns its trailing newline."""

    write_last_delimiter = False
    def _create_grammar(self, allow_dangerous=False):
        values = self._grammar_of_values(allow_dangerous, self.delimiter)
        grammar = TokenConverter(values)

        def tupled(tokens):
            parsed = tokens.as_list()
            while len(parsed) == 1 and isinstance(parsed[0], list):
                parsed = parsed[0]
            return self.is_repeated.key_type(self.name), parsed

        grammar.add_parse_action(tupled)
        grammar.set_name(self.name)
        return grammar

    def _grammar_of_values(self, allow_dangerous=False, delimiter=None):
        """Build the strictly positional SFN record without boundary delimiters."""

        delimiter = delimiter or self.delimiter
        grammars = []
        for item in self._members.values():
            grammar = item._grammar and item._grammar(allow_dangerous)
            if not grammar:
                continue
            if item.is_repeated and item.repeated_with_name:
                repeated_delimiter = (
                    item.repeated_delimiter
                    if item.repeated_delimiter is not None
                    else delimiter
                )
                grammar = item.repeated_count.grammar(grammar, repeated_delimiter)
            if item.is_optional:
                grammar = pp.Optional(grammar)
            grammars.append(grammar)

        values = pp.And(
            [grammars[0]]
            + [delimiter + grammar for grammar in grammars[1:]]
        )
        values.set_parse_action(lambda tokens: dict_from_parsed(tokens.as_list()))
        if self.is_repeated and not self.repeated_with_name:
            values = self.repeated_count.grammar(values, self.repeated_delimiter)
            values.set_name(f"<{self.name}[]>")
        return values


class _FortranArray(GrammarType):
    """A counted NumPy array written as fixed-width Fortran records.

    ``count`` has the same meaning and path syntax as ``repeated_count``.  A
    record can contain more than one value (the radial mesh stores ``r, dr``
    pairs, for example), while ``records_per_line`` only controls writing.
    """

    array_access = True

    def __init__(
        self,
        count,
        value_type,
        item_format,
        *,
        record_size=1,
        records_per_line=1,
    ):
        super().__init__()
        self.count = count
        self.value_type = value_type
        self.item_format = item_format
        self.record_size = record_size
        self.records_per_line = records_per_line
        self._repetition = RepeatedItemGrammar.create(Repeated.REPEATED, count)

    @property
    def _dtype(self):
        return self.value_type.numpy_dtype()[0]

    def added_to_container(self, container):
        self._repetition.apply_hooks(container)

    def _grammar(self, param_name=False):
        value = self.value_type.grammar(param_name)
        whitespace = White(" \t\r\n").suppress()
        if self.record_size > 1:
            value = pp.Group(
                value + (whitespace + value) * (self.record_size - 1)
            )
        grammar = self._repetition.grammar(value, whitespace)
        grammar.add_parse_action(lambda tokens: self.convert(tokens.as_list()))
        return grammar

    def convert(self, value):
        value = np.asarray(value, dtype=self._dtype)
        if self.record_size == 1:
            return value.reshape(-1)
        return value.reshape(-1, self.record_size)

    def _validate(self, value, why="set"):
        if not isinstance(value, np.ndarray):
            return "A NumPy array is required"
        if self.record_size == 1:
            if value.ndim != 1:
                return "A one-dimensional array is required"
        elif value.ndim != 2 or value.shape[1] != self.record_size:
            return f"An array with shape (n, {self.record_size}) is required"
        if isinstance(self.count, int) and len(value) != self.count:
            return f"The array has to contain exactly {self.count} records"
        return True

    def _string(self, value):
        value = self.convert(value)
        if self.record_size == 1:
            value = value.reshape(-1, 1)
        lines = []
        for start in range(0, len(value), self.records_per_line):
            records = value[start : start + self.records_per_line]
            lines.append(
                "".join(self.item_format % item for item in records.flat)
            )
        return "\n".join(lines)

    def copy(self):
        return type(self)(
            self._repetition.copy(),
            self.value_type.copy(),
            self.item_format,
            record_size=self.record_size,
            records_per_line=self.records_per_line,
        )

    def copy_value(self, value):
        return value.copy()

    def numpy_dtype(self):
        shape = () if self.record_size == 1 else (self.record_size,)
        return self._dtype, shape

    def grammar_name(self):
        return f"{self.count} Fortran array records"

    is_the_same_value = staticmethod(compare_numpy_values)


class _RecordItemCount(VariableRepeatedItemGrammar):
    """Use one item of an array-valued option as a repetition count."""

    def __init__(self, source, index):
        super().__init__(source)
        self.index = index

    def call(self, value):
        value = int(value[self.index])
        return value, value

    def copy(self):
        return type(self)(self.count, self.index)


def _fortran_record(*items):
    """Gather scalar options into one fixed-width record without separators."""

    items = gather(*items)
    items[0].output_definition.value_delimiter = ""
    return items


class ShapeFunctionMesh:
    """Generated combined view of one radial and one Voronoi mesh block."""

    def __init__(self, output, index):
        self._output = output
        self._index = index

    @property
    def radial(self):
        return self._output.RADIAL_MESHES[self._index]

    @property
    def voronoi(self):
        return self._output.VORONOI_MESHES[self._index]

    @property
    def idx(self):
        return self.voronoi.IM()

    @property
    def npan(self):
        return self.radial.NPAN()

    @property
    def nr(self):
        return self.radial.NR()

    @property
    def jrcut(self):
        return self.radial.JRCUT()

    @property
    def radial_mesh(self):
        """The native, unscaled ``(r, dr)`` array stored in the file."""

        return self.radial.RMESH()

    @property
    def rmesh(self):
        return self.radial_mesh[:, 0] * self._output.alat

    @property
    def drmesh(self):
        return self.radial_mesh[:, 1] * self._output.alat

    @property
    def nsfn(self):
        return self.radial.NSFN()

    @property
    def sfn_lm(self):
        return np.fromiter(
            (item.LM() for item in self.radial.SHAPE_FUNCTIONS.values()),
            dtype=int,
            count=self.nsfn,
        )

    @property
    def sfn(self):
        if not self.nsfn:
            return np.empty((0, self.nr))
        return np.stack(
            [item.VALUES() for item in self.radial.SHAPE_FUNCTIONS.values()]
        )

    @property
    def rmt(self):
        return self.voronoi.RMTRED0() * self._output.alat

    @property
    def rmtfill(self):
        return self.voronoi.RMTFILL()

    @property
    def vol(self):
        return self.voronoi.VOL()

    @property
    def nface(self):
        return self.voronoi.NFACE()

    @property
    def faces(self):
        return tuple(self.voronoi.FACES.values())

    @staticmethod
    def _index2lm(index):
        l = math.isqrt(index - 1)
        return l, index - l * l - l - 1

    @staticmethod
    def _real_sph_harm(m, l, phi, theta):
        from scipy.special import sph_harm

        if m == 0:
            return np.real(sph_harm(0, l, phi, theta))
        if m > 0:
            return np.real(
                (sph_harm(-m, l, phi, theta) + (-1) ** m * sph_harm(m, l, phi, theta))
                / np.sqrt(2)
            )
        return np.real(
            1j
            * (sph_harm(-m, l, phi, theta) - (-1) ** m * sph_harm(m, l, phi, theta))
            / np.sqrt(2)
        )

    def to_3d_grid(self, n=80):
        """Rebuild the shape function on a cubic grid."""

        rmax = self.rmesh[-1] * 1.02
        linear = np.linspace(-rmax, rmax, n)
        x, y, z = np.meshgrid(linear, linear, linear, indexing="ij")
        radius = np.sqrt(x**2 + y**2 + z**2)
        with np.errstate(invalid="ignore", divide="ignore"):
            theta = np.where(
                radius > 0,
                np.arccos(np.clip(z / radius, -1, 1)),
                0.0,
            )
            phi = np.where(
                radius > 0,
                np.sign(y)
                * np.arccos(
                    np.clip(x / np.sqrt(np.maximum(x**2 + y**2, 1e-30)), -1, 1)
                ),
                0.0,
            )

        grid = np.zeros((n, n, n), dtype=float)
        for values, lm_index in zip(self.sfn, self.sfn_lm):
            l, m = self._index2lm(lm_index)
            radial = np.interp(radius, self.rmesh, values, left=0.0, right=0.0)
            grid += radial * self._real_sph_harm(m, l, phi, theta)
        return linear, grid

    def __repr__(self):
        return (
            f"ShapeFunctionMesh(idx={self.idx}, npan={self.npan}, "
            f"nr={self.nr}, nsfn={self.nsfn}, nface={self.nface})"
        )


class _MeshesView(Sequence):
    def __init__(self, output):
        self.output = output

    def __len__(self):
        return self.output.NM()

    def __getitem__(self, index):
        if isinstance(index, slice):
            return tuple(ShapeFunctionMesh(self.output, i) for i in range(*index.indices(len(self))))
        if index < 0:
            index += len(self)
        if index < 0 or index >= len(self):
            raise IndexError(index)
        return ShapeFunctionMesh(self.output, index)


class SFNOutputFile(OutputFile):
    """Grammar-backed representation of an SPR-KKR shape-function file."""

    def __init__(self, definition=None, container=None, alat=1.0):
        super().__init__(definition, container)
        object.__setattr__(self, "alat", alat)

    @property
    def nm(self):
        return self.NM()

    @property
    def meshes(self):
        return _MeshesView(self)

    def mesh_for_idx(self, index):
        return next((mesh for mesh in self.meshes if mesh.idx == index), None)

    def __repr__(self):
        return f"SFNOutputFile(nm={self.nm}, alat={self.alat:.6g})"


class SFNDefinition(OutputFileDefinition):
    result_class = SFNOutputFile


def create_definition():
    V = _SFNValueDefinition
    S = _SFNSectionDefinition
    GV = GeneratedValueDefinition

    header = V(
        "HEADER",
        _FortranArray(3, Integer(), "%5d", records_per_line=3),
        name_in_grammar=False,
    )

    radial_mesh = S(
        "RADIAL_MESHES",
        [
            V(
                "MESH_COUNTS",
                _FortranArray(2, Integer(), "%5d", records_per_line=2),
                name_in_grammar=False,
            ),
            GV("NPAN", lambda section: int(section.MESH_COUNTS()[0])),
            GV("NR", lambda section: int(section.MESH_COUNTS()[1])),
            V(
                "JRCUT",
                _FortranArray(
                    _RecordItemCount("MESH_COUNTS", 0),
                    Integer(),
                    "%5d",
                    records_per_line=16,
                ),
                name_in_grammar=False,
            ),
            V(
                "RMESH",
                _FortranArray(
                    _RecordItemCount("MESH_COUNTS", 1),
                    Real(),
                    "%22.14E",
                    record_size=2,
                    records_per_line=2,
                ),
                name_in_grammar=False,
            ),
            V("NSFN", Integer(format="5d"), name_in_grammar=False),
            S(
                "SHAPE_FUNCTIONS",
                [
                    V("LM", Integer(format="5d"), name_in_grammar=False),
                    V(
                        "VALUES",
                        _FortranArray(
                            _RecordItemCount("..MESH_COUNTS", 1),
                            Real(),
                            "%22.14E",
                            records_per_line=4,
                        ),
                        name_in_grammar=False,
                    ),
                ],
                name_in_grammar=False,
                is_repeated=True,
                repeated_count="NSFN",
                repeated_with_name=False,
                repeated_delimiter=line_end,
            ),
        ],
        name_in_grammar=False,
        is_repeated=True,
        repeated_count=_RecordItemCount("HEADER", 0),
        repeated_with_name=False,
        repeated_delimiter=line_end,
    )

    face = S(
        "FACES",
        [
            V(
                "FACE_COUNTS",
                _FortranArray(2, Integer(), "%10d", records_per_line=2),
                written_name="IFC NVERT_FCM(IFC,IM)",
                delimiter=line_end,
            ),
            GV("IFC", lambda section: int(section.FACE_COUNTS()[0])),
            GV("NVERT", lambda section: int(section.FACE_COUNTS()[1])),
            V(
                "VERTICES",
                _FortranArray(
                    _RecordItemCount("FACE_COUNTS", 1),
                    Real(),
                    "%25.17E",
                    record_size=3,
                    records_per_line=1,
                ),
                name_in_grammar=False,
            ),
            V(
                "ACOEF",
                _FortranArray(4, Real(), "%25.17E", records_per_line=4),
                written_name="ACOEF_FM",
                delimiter=line_end,
            ),
        ],
        name_in_grammar=False,
        is_repeated=True,
        repeated_count="NFACE",
        repeated_with_name=False,
        repeated_delimiter=line_end,
    )

    voronoi_mesh = S(
        "VORONOI_MESHES",
        [
            *_fortran_record(
                V(
                    "IM",
                    Integer(format="10d"),
                    written_name="IM RMTRED0_M(IM) RMTFILL_M(IM) VOL_M(IM)",
                    delimiter=line_end,
                ),
                V("RMTRED0", Real(format="25.17E"), name_in_grammar=False),
                V("RMTFILL", Real(format="25.17E"), name_in_grammar=False),
                V("VOL", Real(format="25.17E"), name_in_grammar=False),
            ),
            V(
                "NFACE",
                Integer(format="10d"),
                written_name="NFACE_M(IM)",
                delimiter=line_end,
            ),
            face,
        ],
        name_in_grammar=False,
        is_repeated=True,
        repeated_count=_RecordItemCount("HEADER", 0),
        repeated_with_name=False,
        repeated_delimiter=line_end,
    )

    return SFNDefinition(
        "SFN",
        [
            header,
            GV("NM", lambda section: int(section.HEADER()[0])),
            GV("IFMTSFN", lambda section: int(section.HEADER()[1])),
            GV("NL", lambda section: int(section.HEADER()[2])),
            radial_mesh,
            voronoi_mesh,
        ],
    )


definition = create_definition()


__all__ = ["SFNOutputFile", "ShapeFunctionMesh", "SFNDefinition", "definition"]
