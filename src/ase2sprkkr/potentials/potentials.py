"""Potential is object that holds data form SPR-KKR potential file"""

from os import PathLike
from pathlib import Path
from io import StringIO
from typing import Union

from ..sprkkr.configuration import ConfigurationFile
from ..sprkkr.io_data import ReadIoData, WriteIoData
from ..common.decorators import class_property, cache
from ..common.file_fingerprints import file_fingerprint
from ..sprkkr.sprkkr_atoms import SPRKKRAtoms


class PotentialFileChanged(RuntimeError):
    """A potential file changed while it was being parsed."""


class Potential(ConfigurationFile):
    """It holds data form SPR-KKR potential file

    It, in addition to being a containers for their sections, can read/write
    its properties from/to an ASE atoms object."""

    def __init__(self, atoms=None, definition=None, set_from_atoms=True):
        """Use Potential.from_atoms to reuse already created Potential object stored
        in atoms"""
        if definition is None:
            from .definitions.potential import potential_definition as definition
        self._atoms = atoms
        self._complete = False
        super().__init__(definition)
        if set_from_atoms and atoms:
            self.set_from_atoms()

    def _attach_atoms(self, atoms):
        """Make this potential the potential attached to ``atoms``.

        A Potential may be associated with atoms without being their attached
        potential (for example a temporary Potential created from atoms). This
        helper establishes the stronger, bidirectional attachment relationship.
        """
        atoms = SPRKKRAtoms.promote_ase_atoms(atoms)
        previous_atoms = self._atoms
        if (
            previous_atoms is not None
            and previous_atoms is not atoms
            and previous_atoms._potential is self
        ):
            previous_atoms._potential = None

        if atoms is not None:
            previous_potential = atoms._potential
            if (
                previous_potential is not None
                and previous_potential is not self
                and previous_potential._atoms is atoms
            ):
                previous_potential._atoms = None

        self._atoms = atoms
        if atoms is not None:
            atoms._potential = self
        return atoms

    def read_from_file(self, file, atoms=None, allow_dangerous=False):
        super().read_from_file(file, allow_dangerous=allow_dangerous)
        self.make_complete()
        if atoms is not False:
            self._attach_atoms(self.update_atoms(atoms or self._atoms))

    def make_complete(self):
        """Call this function, if you set manually all the properties necessary
        to create the atoms object"""
        self._complete = True

    @property
    def atoms(self):
        if self._atoms is None and self._complete:
            self._attach_atoms(self.update_atoms())
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        self._attach_atoms(atoms)

    def update_atoms(self, atoms=None):
        """Update the ASE object from the values contained in the sections
        of the potential"""
        atoms = SPRKKRAtoms.promote_ase_atoms(atoms or self._atoms)
        iodata = ReadIoData()
        readed = set()

        def do(section):
            nonlocal atoms
            if section._definition.name in readed:
                return
            readed.add(section._definition.name)
            for d in section._depends_on():
                if d in self:
                    do(self[d])
            atoms = section._update_atoms(atoms, iodata) or atoms

        for i in self:
            do(i)

        return atoms

    def set_from_atoms(self, atoms=None, io_data: WriteIoData = None):
        """Set the sections' values of the potential according to the given ASE atoms object.

        Parameters
        ----------
        atoms
          The atoms object, from which the data will be set. If it is None, the ``atoms`` property
          of the potential (``self.atoms``) is used.

        io_data
          The additional (in the time of the creation frozen) state of the atoms object,
          that contains e.g. numbering of the sites, atomic types etc.
          If is not set, it is created from the atoms.
        """
        atoms = atoms or self._atoms
        atoms = SPRKKRAtoms.promote_ase_atoms(atoms)
        atoms.update_sites_from_symbols()
        super().set_from_atoms(atoms, io_data)
        self.make_complete()

    def is_scf_converged(self, atoms=None):
        """Return whether this potential represents a converged SCF state."""
        status = self.SCF_INFO.SCFSTATUS()
        if status == "CONVERGED":
            return True
        if status != "START":
            return False

        atoms = atoms or self._atoms
        return atoms is not None and atoms.sprkkr_has_converged_data()

    @class_property
    @cache
    def potential_definition(cls):
        # import here to avoid circular import issues
        from .definitions import potential as definition

        return definition.potential_definition

    @staticmethod
    def from_file(
        filename,
        atoms=None,
        allow_dangerous=False,
        *,
        check_file_unchanged=False,
    ):
        """Create a potential from a given potential file.

        Parameters
        ----------
        check_file_unchanged
          If True, require ``filename`` to be a filesystem path and verify with
          cheap filesystem fingerprints that it was not replaced or rewritten
          while being parsed. The check is opt-in because normal potential
          loading does not require this provenance guarantee.
        """
        fingerprint = None
        path = None
        if check_file_unchanged:
            if not isinstance(filename, (str, bytes, PathLike)):
                raise TypeError(
                    "check_file_unchanged=True requires a filesystem path"
                )
            path = Path(filename)
            fingerprint = file_fingerprint(path)

        pd = Potential.potential_definition
        potential = pd.read_from_file(
            filename,
            atoms=atoms,
            allow_dangerous=allow_dangerous,
        )

        if check_file_unchanged:
            try:
                current = file_fingerprint(path)
            except OSError as exc:
                raise PotentialFileChanged(
                    f"Potential file changed while it was being read: {path}"
                ) from exc
            if current != fingerprint:
                raise PotentialFileChanged(
                    f"Potential file changed while it was being read: {path}"
                )
        return potential

    @classmethod
    def from_string(cls, string: str, atoms=None, allow_dangerous=False):
        """Create a potential from a string containing a content of a potential file."""
        return cls.from_file(StringIO(string), atoms, allow_dangerous)

    @classmethod
    def from_atoms(cls, atoms):
        """Create a potential, that describes the given atoms object."""
        SPRKKRAtoms.promote_ase_atoms(atoms)
        if atoms.has_potential():
            return atoms.potential
        pd = Potential.potential_definition
        return cls(atoms=atoms, definition=pd)

    def reset(self, update_atoms=True):
        copy = [i for i in self]
        for i in copy:
            if not i._definition.mandatory:
                i.reset()
        self.SCF_INFO.SCFSTATUS = "START"
        if update_atoms:
            self.update_atoms()

    def __repr__(self):
        return "SPRKKR POTENTIAL"

    def __str__(self):
        return "SPRKKR POTENTIAL"

    def save_to_file(
        self,
        file,
        atoms=None,
        *,
        validate: Union[str, bool] = "save",
        restart_scf: bool = False,
    ):
        """Save the potential, optionally restarting SCF only in serialized output.

        ``restart_scf`` is an execution/serialization intent, not potential state.
        All normal synchronization from atoms is performed first. The value
        written for ``SCFSTATUS`` is then overridden to ``START`` without
        changing the canonical status or discarding converged electronic data.
        """
        self.set_from_atoms()
        if atoms is not False:
            self.set_from_atoms(atoms)

        status = None
        if restart_scf:
            status = self.SCF_INFO.SCFSTATUS
            status.result = "START"
        try:
            return super().save_to_file(file, False, validate=validate)
        finally:
            if status is not None:
                status.clear_result()


# At last - to avoid circular import problem
