""" This module contains just a base class for both configuration values - :class:`Options<ase2sprkkr.common.options.Option>`
and configuration containers - :class:`Sections<ase2sprkkr.common.configuration_containers.Section>`.
"""

from typing import Union


class Configuration:
  """ The common base class for all configurations values and containers. I.e.
  for :class:`Options<ase2sprkkr.common.options.Option>` and :class:`Sections<ase2sprkkr.common.configuration_containers.Section>`.
  """

  def __init__(self, definition, container=None):
      """
      Create the object. Just sets the two properties from the parameters.

      Parameters
      ----------
      definition: ase2sprkkr.common.configuration_definitions.BaseDefinition
        Definition of this configuration object.

      container: ase2sprkkr.common.configuration_containers.ConfigurationContainer
        The container, that owns this configuration object.
      """

      self._definition = definition
      """
      The "definition" of the option or section. The definition determines
      the name(s), value type(s) etc... contained in the configuration object.
      Instance of :class:`ase2sprkkr.common.configuration_definitions.BaseDefinition`
      """

      self._container = container
      """
      The parent container. I.e. the container that holds this object (e.g.
      for a value it is the section that owns the value)
      Instance of :class:`ase2sprkkr.common.configuration_containers.ConfigurationContainer`
      """

  def _get_path(self, include_root=False):
      """ Return the dot-delimited path to the item in the configuration tree.

      E.g. the ``ENERGY`` option in the ``CONTROL`` section has the path
      ``CONTROL.ENERGY``
      """
      name = self.name
      if self._container and (include_root or self._container._container):
         return f'{self._container._get_path()}.{name}'
      return name

  def _get_root_container(self):
      """ Return the root object of the configuration tree.

      I.E. the object, that represents the whole configuration or problem-definition file
      """
      return self._container._get_root_container() if self._container else self

  @property
  def name(self):
      """ Return the name of the option/section. The name is defined by the definition
      of the object.

      Returns
      -------
      name: str
      The name of the object.
      """
      return self._definition.name

  def as_dict(self, only_changed:Union[bool,str]='basic'):
      """ Return the value of self, in the case of container as a dictionary. To be redefined in the descendants.

      Parameters
      ----------
      only_changed
        Return only changed values, or all of them?
        If True, return only the values, that differ from the defaults.
        If False, return all the values.
        The default value 'basic' means, return all non-expert values
        and all changed expert values.
      """
      raise NotImplementedError()

  def to_dict(self, only_changed:Union[bool,str]='basic'):
      """ Alias of the method :meth:`as_dict`. """
      return self.as_dict(only_changed)

  def show(self):
      """ Print the configuration, as it will be saved into the configuration/problem definition file. """
      print(self.to_string())

  @property
  def info(self):
      return self._definition.info()

  @property
  def doc(self):
      try:
         return self._definition.description()
      except AttributeError as e:
         raise Exception("Cannot retrieve documentation") from e

  def help(self, verbose=False, show_hidden=False):
      if verbose is True:
         verbose='all'
      elif verbose is False:
         verbose=True
      print(self._definition.description(verbose, show_hidden))

      global _help_warning_printed
      if not _help_warning_printed:
          import __main__ as main
          if verbose is True and not hasattr(main, '__file__'):  # I'm in repl
             print('\n You can use <Configuration>.help(True) for a more detailed description of the possible configuration options. Enjoy ASE2SPRKKR!\n')
          _help_warning_printed = True

  def __repr__(self):
      d = self._definition
      out = d.configuration_type_name
      out = out + ' ' + d.name.upper()
      return out

  def save_to_file(self, file, *, validate:Union[str, bool]='save'):
      """ Save the configuration to a file in a given format.

      This routine do some basic stuff and then call _save_to_file routine,
      that contains the implementation specific for the type of the
      configuration container/value.

      Parameters
      ----------
      file: str or file
        File to read the data from

      validate
        Validate the data in the container first and raise an exception,
        if there is an error (e.g. the the data are not complete).
        The string value can be used to select the type of validation
        ``save`` means the full check (same as the default value ``True``),
        use ``set`` to allow some missing values.
      """
      if validate:
         self.validate(validate)

      if not hasattr(file, 'write'):
         with open(file, "w") as file:
           out=self._save_to_file(file)
           file.flush()
      else:
         out=self._save_to_file(file)
         file.flush()
      return out

  def to_string(self):
      """
      Return the configuration (problem definition) in a string.

      Returns
      -------
      configuration:str
      The configuration, as it should be saved in a configuration/problem definition file.
      """
      from io import StringIO
      s = StringIO()
      self._save_to_file(s)
      return s.getvalue()


_help_warning_printed=False
