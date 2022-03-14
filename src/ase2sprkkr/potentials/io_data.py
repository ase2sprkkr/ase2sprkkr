""" IoData classes serves are intermediate object used as storage during
reading or writing a potential file."""

from ..common.misc import cached_property
from ..common.unique_values import UniqueValuesMapping


def unique_mapping(fce):
    """
    Create function (with a cached result), that returns UniqueValuesMapping from a given iterator
    """
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

    @unique_mapping
    def sites(self):
        return self.atoms.sites

    @unique_mapping
    def types(self):
        for s in self.sites.iter_unique():
            yield from s.occupation.atomic_types()

    @unique_mapping
    def reference_systems(self):
        return (s.reference_system for s in  self.sites.iter_unique())

    @unique_mapping
    def meshes(self):
        return (s.mesh for s in  self.sites.iter_unique())


class ReadIoData(BaseIoData):
    """ Object to store data during reading of a potential file.
        Nothing special is here yet, it is just a dict """
