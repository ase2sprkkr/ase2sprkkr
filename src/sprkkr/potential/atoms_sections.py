from ..common.conf_containers import Section

class AtomsSection(Section):
  """ Section, that sets/retrieve the ASE atoms properties.
      Its values should not be (generally) edited directly,
      since they will be overwrited during saving the POT file
  """

  __setattr__ = object.__setattr__

  def __init__(self, definition, container):
     super().__init__(definition, container)
     self._hook_disabled = 0

  @property
  def _atoms(self):
     """ Atoms object, from which the values are readed """
     return self._container._atoms_io_data.atoms

  @property
  def _atoms_io_data(self):
     """ Object for storing temporary datas during reading/writing """
     return self._container._atoms_io_data

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
