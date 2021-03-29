import numpy as np
from .atomic_types import AtomicType

class Occupation:
  """ Occupation of the atomic site, given by 
      AtomicSite : value
  """

  def __init__(self, dct, site=None):
      self._site = site
      self.set(dct)

  def set(self, dct):
      if isinstance(dct, (int, str, AtomicType)):
         dct = {dct : 1.0}
      if hasattr(dct, 'items'):
         iterator = dct.items() #dict 
      else:
         iterator = dct
      self._occupation = { AtomicType.to_atomic_type(i): j for i,j in iterator }
      self._normalize()

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
      raise KeyError("No {name} in occupation")

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

  def just_one_symbol(self):
      """ Return the chemical symbol of the atom.
          If there are more atoms on the site, then return X. """
      if len(self._occupation) != 1:
         return 'X'
      return next(iter(self._occupation)).symbol

  def __len__(self):
      return len(self._occupation)

  def add(self, name, value=None):
      """ Add atom to the site """
      if value is None:
         value = 1. / (len(self) + 1.)
      self._normalize(1. - value)
      self._occupation[AtomicType.to_atomic_type(i)] = value
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
  def total_occupation(self):
      return sum(self._occupation.values())

  @staticmethod
  def to_occupation(self, occupation):
      if not isinstance(occupation, Occupation):
         occupation = Occupation(occupation)
      return occupation

  def check(self):
      if not np.isclose(sum(self),1.0):
         raise ValueError("Total occupation of the site should be equal to one")

  def to_tuple(self):
      return zip(self.keys(), self.values())

  @classmethod
  def to_occupation(cls, occupation):
      if isinstance(occupation, (int, str)):
        occupation = AtomicType(occupation)
      if isinstance(occupation, AtomicType):
        occupation = cls({ occupation : 1. })
      return occupation

  def atomic_types(self):
      return self._occupation.keys()
