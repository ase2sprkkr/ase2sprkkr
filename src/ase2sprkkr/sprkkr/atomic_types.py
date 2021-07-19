from ..common.misc import cached_property
import copy

class AtomicType:

    _mendeleev_module = None

    @cached_property
    def mendeleev(self):
        if self.n_electrons == 0:
           raise Exception("Vaccuum (pseudo)atom is not in Mendeleev package database")
        if not AtomicType._mendeleev_module:
          try:
            import mendeleev
          except ImportError as e:
            raise ImportError("Cannot import Mendeleev package to guess the atomic type properties") from e
          _mendeleev_module = mendeleev
        return _mendeleev_module.element(self.n_electrons or self.symbol)

    def __init__(self, symbol, n_electrons=None, n_core=None, n_valence=None, n_semicore=None):
        """
        Parameters
        ----------
        symbol: str or int
          Chemical symbol, e.g. Fe, N, ....
          If int is given, attempt to guess the symbol using third party library
          is done

        n_electrons: int
          Number of core electrons
          If it's zero, an attempt to guess it from the chemical symbol is done.

        n_valence: int
          Number of core electrons
          If it's zero, an attempt to guess it from the chemical symbol is done.
        """
        self._mendeleev = None
        if isinstance(symbol, int):
           if n_electrons is not None and n_electrons != symbol:
              raise ValueError(f'Number of electrons in symbol ({symbol}) and n_electrons ({n_electrons})differs')
           self.n_electrons = symbol
           self.symbol = None
        else:
           self.symbol = symbol
           self.n_electrons = n_electrons

        if self.symbol == 'Vc' or self.n_electrons == 0:
           self.symbol = 'Vc'
           self.n_electrons = self.n_core = self.n_valence = self.n_semicore = 0
        else:
          if not self.symbol: self.symbol = self.mendeleev.symbol
          if not self.n_electrons: self.n_electrons = self.mendeleev.atomic_number
          self.n_valence = self.mendeleev.nvalence() if n_valence is None else n_valence
          self.n_core = self.n_electrons - self.n_valence if n_core is None else n_core;
          self.n_semicore = self.n_electrons - self.n_core - self.n_valence  if n_semicore is None else n_semicore

    def copy(self):
        return copy.copy(self)

    def __repr__(self):
        return f"({self.n_electrons})" if self.symbol == 'X' else self.symbol

    def __str__(self):
        return self.__repr__()

    def to_tuple(self):
        return (self.symbol, self.n_electrons, self.n_core, self.n_valence, self.n_semicore)

    @property
    def atomic_number(self):
        return self.n_electrons

    @classmethod
    def to_atomic_type(cls, value):
        if isinstance(value, cls):
           return value
        return cls(value)
