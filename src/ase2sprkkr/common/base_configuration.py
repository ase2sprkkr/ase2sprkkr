class BaseConfiguration:
  """ Common base class for configuration """

  def __init__(self, definition, container=None):
      self._definition = definition
      self._container = container

  def _get_path(self, include_root=False):
      name = self._definition.name
      if self._container and (include_root or self._container._container):
         return f'{self._container._get_path()}.{name}'
      return name

  def _get_root_container(self):
      return self._container._get_root_container() if self._container else self
