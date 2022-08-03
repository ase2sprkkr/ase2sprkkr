""" More atomic types can be on one site. """


import numpy as np
from .atomic_types import AtomicType
from ..common.misc import OrderedDict


class Occupation:
  """ Occupation of the atomic site, given by AtomicType : value

   The value determine the probability, that a given atomic type will
   be found on a given place (which can be used e.g. for computing
   alloys).
  """

  def __init__(self, dct, site=None):
      self._site = site
      self.set(dct)

  def copy(self):
      return Occupation({ a.copy(): v for a,v in self._occupation.items() })

  def set(self, dct):
      if isinstance(dct, (int, str, AtomicType)):
         dct = {dct : 1.0}
      if hasattr(dct, 'items'):
         iterator = dct.items() #dict
      else:
         iterator = dct
      self._occupation = OrderedDict( (AtomicType.to_atomic_type(i), j) for i,j in iterator )
      self._normalize()
      self._update_atoms()

  def items(self):
      """ dict.items() like enumeration """
      return self._occupation.items()

  def __repr__(self):
      return f"Occupation {self._occupation}"

  def __str__(self):
      return f"Occupation {self._occupation}"

  def __iter__(self):
      return iter(self._occupation)

  def _update_atoms(self):
      if self._site:
         self._site.update_atoms()

  def __iter__(self):
      return iter(self._occupation)

  def _find_key(self, name):
      if isinstance(name, AtomicType):
          return name
      for i in self:
          if i.symbol == name: return i
      raise KeyError(f"No {name} in the occupation")

  def __getitem__(self, name):
      name = self._find_key(name)
      return self._occupation[name]

  def __setitem__(self, name, value):
      try:
          name = self._find_key(name)
          self._normalize(1. - value, name)
          self._occupation[name] = value
      except KeyError:
          self.add(name, value)
          return
      self._update_atoms()

  @property
  def primary_atomic_type(self):
      """ Return the atomic type.
          If there are more atoms on the site, return the one
          with the largest occupation """
      m = 0.
      prim = None
      for at,occ in self._occupation.items():
          if occ > m and at.atomic_number > 0:
             m = occ
             prim = at
      return prim

  @property
  def primary_atomic_number(self):
      """ Return the atomic number of the atom at the site.
          If there are more atoms on the site, return the "main one". """
      primary = self.primary_atomic_type
      return primary.atomic_number if primary else 0

  @property
  def primary_symbol(self):
      """ Return the chemical symbol of the atom at the site.
          If there are more atoms on the site, return the "main one". """
      primary = self.primary_atomic_type
      return primary.symbol if primary else 'X'

  def __len__(self):
      return len(self._occupation)

  def add(self, name, value=None):
      """ Add atom to the site """
      if value is None:
         value = 1. / (len(self) + 1.)
      self._normalize(1. - value)
      self._occupation[AtomicType.to_atomic_type(name)] = value
      self._update_atoms()

  def __delitem__(self, name):
      name = self._find_key(name)
      del self._occupation[name]
      self._normalize()

  def _normalize(self, to=1., except_from=None):
      """
      Normalizes occupation so the sum will be equal to value

      Parameters
      ----------
      value : float
        Desired value

      except_from: AtomicType
        Skip given atom during normalizing
      """
      suma = self.total_occupation
      if except_from:
         suma -= self._occupation[except_from]
      if suma == to:
         return
      ratio = to / suma
      for i in self._occupation:
         if i != except_from:
            self._occupation[i] *= ratio

  @property
  def as_dict(self):
      occ = {}
      for at in self:
         if at.atomic_number == 0:
            continue
         occ[at.symbol] = occ.get(at.symbol, 0) + self[at]
      return occ

  @as_dict.setter
  def as_dict(self, x):
      self.set(x)

  @property
  def total_occupation(self):
      return sum(self._occupation.values())

  @staticmethod
  def to_occupation(occupation):
      if not isinstance(occupation, Occupation):
         occupation = Occupation(occupation)
      return occupation

  def check(self):
      if not np.isclose(sum(self),1.0):
         raise ValueError("Total occupation of the site should be equal to one")

  def to_tuple(self):
      return zip(self.keys(), self.values())

  def atomic_types(self):
      return self._occupation.keys()
