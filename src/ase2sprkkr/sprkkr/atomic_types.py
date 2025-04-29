""" Atomic could be placed at the atomic site. """

from ..common.decorators import cached_property
import copy
from .radial import RadialPotential, RadialCharge
from .moments import Moments
from typing import Dict


class AtomicType:
    """ Atomic type represent a type of atom, that could be placed at the atomic site.

    It can be either a real chemical element, or vacuum pseudoelement.
    It also determine the number of electrons and valence electrons.
    """

    _mendeleev_module = None

    @classmethod
    def to_atomic_type(cls, value, for_mesh=None):
        if isinstance(value, cls):
            return value.for_mesh(for_mesh)
        return cls(value, mesh=for_mesh)

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

    def __init__(self, symbol, atomic_number=None, n_core=None, n_valence=None, n_semicore=None, n_electrons=None, mesh=None):
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
              raise ValueError(f'Number of electrons in symbol ({symbol}) and atomic_number ({atomic_number}) differ')
           atomic_number = symbol
           symbol = None
        else:
           symbol = symbol
           atomic_number = atomic_number

        if atomic_number is None and symbol is None:
           raise ValueError("Unknown atomic type")

        if symbol == 'Vc' or atomic_number == 0 or \
           (symbol == 'X' and atomic_number is None):
           self._symbol = 'Vc'
           self._atomic_number = 0
        else:
           if symbol is None:
              self._atomic_number = atomic_number
              self._symbol = symbol if symbol is not None else self.mendeleev.symbol
           else:
              self._symbol = symbol
              self._atomic_number = None
              self._atomic_number = atomic_number if atomic_number is not None else self.mendeleev.atomic_number

        self._n_electrons = n_electrons
        self._n_valence = n_valence
        self._n_core = n_core
        self._n_semicore = n_semicore

        self._potential = None
        self._charge = None
        self._moments = None
        self._mesh = mesh

        self._check_n_electrons()

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        if self._mesh is mesh:
            return
        self._mesh = mesh
        if self.potential:
            self.potential.mesh = mesh
        if self.charge:
            self.charge.mesh = mesh

    def _check_n_electrons(self):
        if (
            self._n_core is not None and
            self._n_valence is not None and
            self._n_semicore is not None and
            self._n_electrons is not None and
            self._n_core + self._n_valence + self._n_semicore != self._n_electrons
           ):
              raise ValueError(f"""The following atom setup is inconsistent:
n_electrons: {self._n_electrons},
n_core: {self._n_core},
n_valence: {self._n_valence},
n_semicore: {self._n_semicore}""")

    @property
    def atomic_number(self):
        return self._atomic_number

    def _clear_symbol_cache(self):
        try:
          del self.mendeleev
        except (AttributeError, KeyError):
          pass

        try:
          del self.atomic_number
        except (AttributeError, KeyError):
          pass

        try:
          del self.symbol
        except (AttributeError, KeyError):
          pass

    @atomic_number.setter
    def atomic_number(self, v):
        self._atomic_number = v
        self._clear_symbol_cache()
        self._symbol = self.mendeleev.symbol if v else 'Vc'

    @property
    def symbol(self):
        return self._symbol

    @symbol.setter
    def symbol(self, v):
        self._symbol = v
        self._atomic_number = None
        self._clear_symbol_cache()
        self._atomic_number = self.mendeleev.atomic_number if v not in ['X', 'Vc'] else 0

    @property
    def n_electrons(self):
        if self._n_electrons is not None:
           return self._n_electrons
        if self._n_valence is not None and self._n_core is not None and self._n_semicore is not None:
           return self._n_valence + self._n_core + self._n_semicore
        return self.atomic_number

    @n_electrons.setter
    def n_electrons(self, v):
        self._n_electrons = v
        if self._potential:
            self._potential.z = v
        self.__check_n_electrons()

    @property
    def n_valence(self):
        if self._n_valence is not None:
           return self._n_valence
        if self._n_core is not None:
           return self.n_electrons - self._n_core
        if self.n_electrons == 0:
           return 0
        return self.mendeleev.nvalence()

    @n_valence.setter
    def n_valence(self, v):
        self._n_valence = v
        self.__check_n_electrons()

    @property
    def n_core(self):
        if self._n_core is not None:
           self._n_core
        return self.n_electrons - self.n_valence - self.n_semicore

    @n_core.setter
    def n_core(self, v):
        self._n_core = v
        self.__check_n_electrons()

    @property
    def n_semicore(self):
        if self._n_semicore is not None:
           return self._n_semicore
        if self._n_core is not None and self._n_valence is not None:
           return self.n_electrons - self._n_core - self._n_valence
        return 0

    @n_semicore.setter
    def n_semicore(self, v):
        self._n_semicore = v
        self.__check_n_electrons()

    def __repr__(self):
        return f"({self.atomic_number})" if self.symbol == 'X' else self.symbol

    def __str__(self):
        return self.__repr__()

    def to_tuple(self):
        return (self.symbol, self.n_electrons, self.n_core, self.n_valence, self.n_semicore)

    def is_vacuum(self):
        return self.atomic_number == 0

    @property
    def potential(self):
        """ The radial potential data of the site """
        return self._potential

    @potential.setter
    def potential(self, value):
        if value is None:
            self._potential = None
            return
        if not isinstance(value, RadialPotential):
            self._potential = RadialPotential(value, self._mesh, self.n_electrons)
        elif value.z is not None:
            self._potential = value.copy(self.mesh)
            self._potential.z = self._n_electrons
        else:
            self._potential = value
            self._potential.z = self._n_electrons

    @property
    def charge(self):
        """ The radial charge data of the site """
        return self._charge

    @charge.setter
    def charge(self, value):
        if value is None:
            self._charge = None
            return
        if not isinstance(value, RadialCharge):
            self._charge = RadialCharge(value, self.mesh)
        else:
            self._charge = value.for_mesh(self.mesh)

    @property
    def moments(self):
        """ The moments data of the site """
        return self._moments

    @moments.setter
    def moments(self, value):
        if not isinstance(value, Moments):
            value = Moments(*value)
        self._moments = value

    def for_mesh(self, mesh):
        """ Return copy for a given mesh """
        if mesh is None:
            return self
        smesh = self.mesh
        if smesh is mesh:
            return self
        if smesh is None:
            self.mesh = mesh
            return self
        return self.copy(mesh)

    def copy(self, for_mesh=None):
        out = copy.copy(self)
        if out.potential is not None:
            out.potential = out.potential.copy(for_mesh)
        if out.charge is not None:
            out.charge = out.charge.copy(for_mesh)
        if out.moments is not None:
            out.moments = out.moments.copy()
        if for_mesh is not None:
            out._mesh = for_mesh
        return out

    def remesh(self, mesh, map:'Dict[AtomicType,AtomicType]'):
        """ Return the copy for a given mesh, updating the informations
        in the map of the old: new objects, so the remeshed items will
        share the meshes. """
        if self in map:
            return map[self]
        out = map[self] = self.for_mesh(mesh)
        return out

    def has_converged_data(self):
        return self.potential and self.charge and self.moments
