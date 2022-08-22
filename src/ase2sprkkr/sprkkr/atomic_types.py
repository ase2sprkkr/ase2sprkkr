""" Atomic could be placed at the atomic site. """

from ..common.misc import cached_property
import copy

class AtomicType:
    """ Atomic type represent a type of atom, that could be placed at the atomic site.

    It can be either a real chemical element, or vacuum pseudoelement.
    It also determine the number of electrons and valence electrons.
    """

    _mendeleev_module = None

    @cached_property
    def mendeleev(self):
        if self.atomic_number == 0:
           raise Exception("Vaccuum (pseudo)atom is not in Mendeleev package database")
        if not AtomicType._mendeleev_module:
          try:
            import mendeleev
          except ImportError as e:
            raise ImportError("Cannot import Mendeleev package to guess the atomic type properties") from e
          AtomicType._mendeleev_module = mendeleev
        return AtomicType._mendeleev_module.element(self.atomic_number or self.symbol)

    def __init__(self, symbol, atomic_number=None, n_core=None, n_valence=None, n_semicore=None, n_electrons=None):
        """
        Parameters
        ----------
        symbol: str or int
          Chemical symbol, e.g. Fe, N, ....
          If int is given, attempt to guess the symbol using third party library
          is done

        atomic_number: int
          Atomic number
          If it's zero, an attempt to guess it from the chemical symbol is done.

        n_core: int
          Number of core electrons
          If it's zero, an attempt to guess it from the chemical symbol/atomic number is done.

        n_valence: int
          Number of valence electrons
          If it's zero, an attempt to guess it from the chemical symbol/atomic_number is done.

        n_semicore: int
          Number of semicore electrons
          If it's zero, an attempt to guess it from the chemical symbol/atomic_number is done.

        n_electrons: int
        """
        if isinstance(symbol, int):
           if atomic_number is not None and atomic_number != symbol:
              raise ValueError(f'Number of electrons in symbol ({symbol}) and atomic_number ({atomic_number})differs')
           self.atomic_number = symbol
           self.symbol = None
        else:
           self.symbol = symbol
           self.atomic_number = atomic_number

        if self.symbol == 'Vc' or self.atomic_number == 0 or \
           (self.symbol == 'X' and self.atomic_number is None):
           self.symbol = 'Vc'
           self.atomic_number = self.n_electrons = self.n_core = self.n_valence = self.n_semicore = 0
        else:
          if not self.symbol: self.symbol = self.mendeleev.symbol
          if self.atomic_number is None: self.atomic_number = self.mendeleev.atomic_number

          if n_electrons is not None:
             self.n_electrons = n_electrons
          elif n_valence is not None and n_core is not None and n_semicore is not None:
             self.n_electrons = n_valence+n_core+n_semicore
          else:
             self.n_electrons = self.atomic_number
          self.n_valence = self.mendeleev.nvalence() if n_valence is None else n_valence
          self.n_core = self.n_electrons - self.n_valence if n_core is None else n_core;
          self.n_semicore = self.n_electrons - self.n_core - self.n_valence  if n_semicore is None else n_semicore

    def copy(self):
        return copy.copy(self)

    def __repr__(self):
        return f"({self.atomic_number})" if self.symbol == 'X' else self.symbol

    def __str__(self):
        return self.__repr__()

    def to_tuple(self):
        return (self.symbol, self.n_electrons, self.n_core, self.n_valence, self.n_semicore)

    @classmethod
    def to_atomic_type(cls, value):
        if isinstance(value, cls):
           return value
        return cls(value)

    def is_vacuum(self):
        return self.atomic_number == 0
