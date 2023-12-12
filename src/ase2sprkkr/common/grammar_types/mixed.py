""" Grammar types for variant types (that can accepts values of more types) and the types derived from it (e.g. the Range type is in fact
Variant type that can accepts either number or two-item array) """

import pyparsing as pp
import copy

from .grammar_type import GrammarType, compare_numpy_values, recognized_set_types, type_from_type, type_from_value
from ..decorators import cached_property, add_to_signature

from .basic import Energy, Real, Integer, QString, LineString, Flag, Bool
from .arrays import SetOf, set_of_integers, set_of_reals


class BaseMixed(GrammarType):
  """
  A variant type - it can hold "anything".
  """

  type = None
  """ The types, that the value can hold. To be redefined in the descendants. """

  string_type = None
  """ Type of string grammar_type to be used.  To be redefined in the descendants. """

  def _grammar(self, param_name=False):
      return pp.MatchFirst((
        i.grammar(param_name) for i in self.types
      ))

  def get_type(self, value):
      """ Return the type of the value.
      Actualy, this implementation is a simple implementation that suits for the common
      Mixed types, so if you make a custom Mixed type, redefine it.
      """
      return self.string_type if isinstance(value, str) else type_from_value(value)

  def _validate(self, value, why='set'):
      if value is None:
          return True
      type = self.get_type(value)
      if type is value:
          return 'Can not determine the type of value {}'.format(value)
      return type.validate(value, why)

  def write(self, f, val):
      """ Output the value to the stream (in the propper format). """
      type = self.get_type(val)
      if type is val:
          return super().write(f, val)
      type.write(f, val)

  def grammar_name(self):
      return '<mixed>'

  def convert(self, value):
      if value is None:
          return None
      return self.get_type(value).convert(value)

  def copy_value(self, value):
      return copy.deepcopy(value)


class Range(BaseMixed):
  """ A range type - it accepts either one value or range of two values of a given type."""

  @add_to_signature(BaseMixed.__init__, prepend=True)
  def __init__(self, type, *args, **kwargs):
      self._type = type_from_type(type)
      super().__init__(*args, **kwargs)

  @cached_property
  def types(self):
      return [
          self._type,
          SetOf(self._type, min_length=2, max_length=2)
      ]

  def get_type(self, value):
      return self.types[1 if isinstance(value, recognized_set_types) else 0]


class Mixed(BaseMixed):
  """ A variant value to be used in input files (in unknown - custom - options) """

  types = [
        Energy.I,
        Real.I,
        Integer.I,
        set_of_integers,
        set_of_reals,
        QString.I,
        Flag.I,
  ]
  """ Possible types of the value """

  string_type = QString.I
  """ Input files use quoted strings. """

  def missing_value(self):
    return True, True, False

  is_the_same_value = staticmethod(compare_numpy_values)


class PotMixed(BaseMixed):
  """ A variant value to be used in potential files (in unknown - custom - options) """

  types = [
        Energy.I,
        Real.I,
        Integer.I,
        Bool.I,
        set_of_integers,
        set_of_reals,
        LineString.I,
  ]
  """ Possible types of the value """

  string_type = LineString.I
  """ Potential files use line strings. """

  def _string(self, val):
    if isinstance(val, bool):
       return Bool._string(self, val)
    else:
       return super()._string(val)

  is_the_same_value = staticmethod(compare_numpy_values)


mixed = Mixed.I = Mixed()   # NOQA: E741
""" A standard grammar type instance for variant (mixed) in input files """
pot_mixed = PotMixed.I = PotMixed()    # NOQA: E741
""" A standard grammar type instance for variant (mixed) values in potential files """
