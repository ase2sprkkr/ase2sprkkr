""" In this module, the :class:`GrammarTypes<ase2sprkkr.common.grammar_types.GrammarType>`,
that are specific for SPR-KKR, are present. """

from ..common.grammar_types import ObjectNumber
from ..common.alternative_types import allowed_types
import numpy as np
from .sites import Site as _Site
from .atomic_types import AtomicType as _AtomicType

class Site(ObjectNumber):

  type = _Site

  def set_from_atoms(self, option, atoms, io_data):
      value = option()
      if isinstance(value, Site):
          option.result = io_data.sites[value]
Site.I = Site()

class AtomicType(ObjectNumber):

  type = _AtomicType

  def set_from_atoms(self, option, atoms, io_data):
      value = option()
      if isinstance(value, Site):
          option.result = io_data.types[value]
AtomicType.I = AtomicType()
