from ..common.conf_containers import Section
from ase.spacegroup import get_spacegroup

class AtomSection(Section):
  """ Section, that sets the ASE atoms properties """

  __setattr__ = object.__setattr__

  def __init__(self, definition, container):
     super().__init__(definition, container)
     self._hook_disabled = 0

  @property
  def _atoms(self):
     return self._container.atoms

  @property
  def _bravais_lattice(self):
     return self._container.bravais_lattice

  @property
  def _cell_spacegroup(self):
     return self._container.cell_spacegroup

  def _add(self, v):
      super()._add(v)
      v._hook = self._hook

  def clear(self, do_not_check_required=False):
      self._hook_disabled+=1
      try:
         super().clear(do_not_check_required)
         self._process()
      finally:
         self._hook_disabled-=1

  def set(self, values):
      self._hook_disabled+=1
      try:
         super().set(values)
         self._process()
      finally:
         self._hook_disabled-=1

  def _hook(self, obj):
       if not self._hook_disabled:
          self._process()

  def save_to_file(self, file):
      self._hook_disabled+=1
      try:
         self._set_from_atoms()
         super().save_to_file(file)
      finally:
         self._hook_disabled-=1
