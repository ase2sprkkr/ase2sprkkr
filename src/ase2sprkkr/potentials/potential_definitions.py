"""
This file contains the classes for definitions of Potentials:
the list of sections and their allowed (or standard) options
and their value formats. InputParameterss and potentials have a simlilar
structure, so they share common functionalities from
sprkkr.common.configuration_definitions
"""

import pyparsing as pp
from ..common.grammar import line_end
from ..common.grammar_types import separator, pot_mixed
from ..common.configuration_definitions import SeparatorDefinition
from ..sprkkr.configuration import \
         ConfigurationValueDefinition, ConfigurationSectionDefinition, \
         ConfigurationFileDefinition, CustomConfigurationValue
from .custom_potential_section import CustomPotentialSection, SectionString
from .potentials import Potential
from .potential_sections import PotentialSection, ASEArraySection
from ..common.decorators import cached_class_property, cache
from ..common.decorators import add_to_signature


class PotValueDefinition(ConfigurationValueDefinition):
  """
  Definition of a configuration option in a potential
  """
  @cached_class_property
  def grammar_of_delimiter():
    return pp.Empty().setName(' ')

  prefix = ''
  name_value_delimiter = '\t'

  def __init__(self, *args, required=None, **kwargs):
      super().__init__(*args, required=required, **kwargs)

  type_of_dangerous = pot_mixed


class Separator(SeparatorDefinition):
  """
  A special (hidden) value, that appears in a potential header section.

  The separator is a line of doublecross
  """
  separator_type = separator


class PotSectionDefinition(ConfigurationSectionDefinition):
  """ This class describes the format of one
  value of a standard potential section """

  force_order = True
  """ The order of items in potential file is fixed """

  value_name_format = '<12'

  child_class = PotValueDefinition
  """ standard child class """

  custom_class = staticmethod(CustomConfigurationValue.factory(PotValueDefinition, pot_mixed))
  """ Adding a custom values is allowed """

  delimiter = '\n'
  """ options are delimited by newline in ouptut. """

  grammar_of_delimiter = line_end
  """ items are separated by newline """

  def depends_on(self):
      """ The order of processing of sections during reading can be different than the order during a write. So, if the function should not be processed before given named sections, name then.

      Return
      ------
      prerequisites: [ str, str, ... ]
      """
      return []

  result_class = PotentialSection

  def __init__(self, *args, mandatory:bool=True, **kwargs):
      """
      For the documentation of the other parameters, see
      :meth:`ase2sprkkr.common.SectionDefinition`

      Parameters
      ----------
      mandatory
        Is the section mandatory? I.e. the potential file is required to
        contain this sections.
      """
      self.mandatory = mandatory
      super().__init__(*args, **kwargs)


class ASEArraySectionDefinition(PotSectionDefinition):
  """
  A definition of a section, that contains an ASE datas (Atoms.setArray)
  """

  @add_to_signature(PotSectionDefinition.__init__)
  def __init__(self, *args, array_name, **kwargs):
      """
      For the documentation of the other parameters, see
      :meth:`ase2sprkkr.potential_definitions.PotSectionDefinition`
      and its predecessor
      :meth:`ase2sprkkr.common.SectionDefinition`

      Parameters
      ----------
      array_name: str
        The name of the ASE array that contains the section's data
      """
      super().__init__(*args, **kwargs)
      self.array_name = array_name

  def depends_on(self):
      """ Array size is required """
      return [ 'SITES' ]

  result_class = ASEArraySection


class PotentialDefinition(ConfigurationFileDefinition):
  """ This class describes the format of a potential file """

  child_class = PotSectionDefinition
  """ Definition of the standard child class: """

  result_class = Potential
  """ The parsing of a potential file results in an instance of Potential. """

  configuration_type_name = 'POTENTIAL FILE'
  """ Name of the container type in the runtime documentation """

  force_order = True
  """ The order of items in potential file is fixed """

  delimiter="*" * 79 + "\n"
  """ Sections delimiter """

  @cached_class_property
  def grammar_of_delimiter():
      """ Grammar of the sections delimiter """
      return SectionString.grammar_of_delimiter()

  custom_class = CustomPotentialSection
  """ Unknown sections will be of this type """

  def _generic_info(self):
      return "Definition of the format of SPRKKR potential file"

  @classmethod
  @cache
  def custom_value_grammar(cls):
      """ Unknown sections are parsed by this grammar """
      return SectionString._grammar

  custom_name_characters = ConfigurationFileDefinition.custom_name_characters + ' '
  """ There can be space in a potential-section name """
