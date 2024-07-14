""" IoData classes serves are intermediate object used as storage during
reading or writing of potential/input_parameters files."""

from ..common.decorators import cached_property
from ..common.unique_values import UniqueValuesMapping
from functools import wraps


def unique_mapping(fce):
    """
    Create function (with a cached result), that returns UniqueValuesMapping from a given iterator
    """
    @wraps(fce)
    def get_unique_mapping(self):
        lst = [ i for i in fce(self) ]
        return UniqueValuesMapping.from_values(lst)
    return cached_property(get_unique_mapping)


class BaseIoData(dict):
    """ Base class for object used during reading/writing.
        For a potential future use. """
    pass


class WriteIoData(BaseIoData):
    """ During writing of potential file, some lists are needed in more sections,
        typically one section contains the list itself, while another(s) reference(s)
        the items of the list.
        So the lists are created on demand and cached for its further use.
    """

    def __init__(self, atoms):
        self.atoms = atoms
        self._has_converged_data = None

    @unique_mapping
    def sites(self):
        return self.atoms.sites

    @unique_mapping
    def types(self):
        for s in self.sites.iter_unique():
            yield from s.occupation.atomic_types()

    @unique_mapping
    def reference_systems(self):
        return (s.reference_system for s in self.sites.iter_unique())

    @unique_mapping
    def meshes(self):
        return (s.mesh for s in self.sites.iter_unique())

    def has_converged_data(self, potential):

        def compute():
            if potential.SCF_INFO.FULLPOT():
                return False
            for site, i in self.sites.unique_items():
              if not site.charge or not site.potential or not site.moments:
                  return False

            return True

        if self._has_converged_data is None:
            self._has_converged_data = compute()
        return self._has_converged_data


class ReadIoData(BaseIoData):
    """ Object to store data during reading of a potential file. """

    def __init__(self):
        self._postponed = []

    def apply_on_atoms(self, handler, atoms):
        """ Apply the given function on the atoms object, if not-None is given.
        However, the store the handler, to be executed (or repeated) on a newly
        given atoms - for the case, that a new atoms object will be created later
        """
        if atoms:
           handler(atoms)
        self._postponed.append(handler)

    def update_atoms(self, atoms):
        """ Replay are the stored handlers on the given atoms object """
        if atoms.are_sites_inited():
            for i in atoms.sites:
                i._clear_data()
        for handler in self._postponed:
            handler(atoms)
