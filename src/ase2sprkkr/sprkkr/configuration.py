""" In this module the base classes for SPR-KKR configuration sections and files are present.
"""


from ..common.configuration_containers import RootConfigurationContainer, Section, CustomSection
from ..common.configuration_definitions import ConfigurationRootDefinition, SectionDefinition, ValueDefinition
from ..common.options import Option, CustomOption

from .io_data import WriteIoData
from ase import Atoms
from typing import Optional, Union

from .sprkkr_atoms import SPRKKRAtoms

class ConfigurationSectionTrait:

  def set_from_atoms(self, atoms:Atoms, io_data:Optional[WriteIoData]=None):
      """ Set the sections' values of the potential according to the given ASE atoms object.

      Parameters
      ----------
      atoms
        The atoms object, from which the data will be set.

      io_data
        The additional (in the time of the creation frozen) state of the atoms object,
        that contains e.g. numbering of the sites, atomic types etc.
        If it is not set, it is created from the atoms.

      validate: bool or str
        If False, do not validate the values
        String value: accepts any of values accepted by ``why`` argument of
        :func:`GrammarType.validate<ase2sprkkr.common.grammar_types.GrammarType.validate>`
        True (default) means full validation: i.e. the same as ``save``.
      """
      atoms = SPRKKRAtoms.promote_ase_atoms(atoms)
      io_data = io_data or WriteIoData(atoms)
      self._set_from_atoms(atoms, io_data)

  def _set_from_atoms(self, atoms:SPRKKRAtoms, io_data:WriteIoData):
      """ Set the sections' values of the potential according to the given ASE atoms object.
          Unlike the non_underscored routine, this one requires the io_data to be set.
      """
      for i in self:
          i._set_from_atoms(atoms, io_data)

# Containers and values

class ConfigurationFile(RootConfigurationContainer, ConfigurationSectionTrait):
  """ 'Root' configuration container for SPRKKR configuration file """

  def save_to_file(self, file, atoms=None, *, validate='save'):
      if atoms is not False:
         self.set_from_atoms(atoms)
      if validate:
         self.validate(validate)
      super().save_to_file(file)

class ConfigurationSection(Section, ConfigurationSectionTrait):
  """ Configuration section to be used in SPRKKR configuration files """


class CustomConfigurationSection(CustomSection, ConfigurationSectionTrait):
  """ Custom configuration section to be used in SPRKKR configuration files """


class ConfigurationValue(Option):
  """ Value (option) in a SPRKKR configuration file. """

  def _set_from_atoms(self, atoms:Atoms, io_data:WriteIoData):
      """ Some types can/should be updated by the data in :class:`atoms<Atoms>` object
      (and/or :class:`io_data<WriteIoData>`) during save.
      """
      if hasattr(self._definition.type, 'set_from_atoms'):
          self._definition.type.set_from_atoms(self, atoms, io_data)


class CustomConfigurationValue(CustomOption):
  """ Custom value (option) in a SPRKKR configuration file. """

  def _set_from_atoms(self, atoms:Atoms, io_data:WriteIoData):
      pass


# Definitions


class ConfigurationValueDefinition(ValueDefinition):
  """ Definition o a configuration value, used in SPRKKR configuration file."""

  result_class = ConfigurationValue
  """ The values (or more precisely the objects holding the values) of SPR-KKR configuration should be of this class. """


class ConfigurationSectionDefinition(SectionDefinition):
  """ Default class for definitions of SPRKKR configuration files' sections """

  result_class = ConfigurationSection
  """ The sections should have this class (if not stated otherwise) """


class ConfigurationFileDefinition(ConfigurationRootDefinition):
  """ The class for the configuration files (the root nodes of the configuration)."""