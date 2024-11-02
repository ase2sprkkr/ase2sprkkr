""" This module contains just a base class for both configuration values - :class:`Options<ase2sprkkr.common.options.Option>`
and configuration containers - :class:`Sections<ase2sprkkr.common.configuration_containers.Section>`.
"""

from typing import Union
from .warnings import DataValidityWarning, DataValidityError
import warnings


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

  def _as_dict(self, get):
      raise NotImplementedError()

  @staticmethod
  def as_dict_getter(only_changed:Union[bool,str]='basic', generated=False, copy=False):

      def get(self):
          d = self._definition
          if d.is_generated and not generated:
               return None
          if only_changed == 'explicit':
                if self._definition.is_generated and not generated:
                    return None
                v = self._unpack_value(self._value)
          elif only_changed and (only_changed!='basic' or d.is_expert) and not d.is_always_added:
               v,c = self.value_and_changed()
               if not c:
                    return None
          else:
               v = self(all_values=True)
          if v is not None:
              if copy:
                   v = self._definition.copy_value(v, all_values=True)
          return v

      return get

  def as_dict(self, only_changed:Union[bool,str]='basic', generated=False, copy=False, getter=None):
      """ Return the value of self, in the case of container as a dictionary. To be redefined in the descendants.

      Parameters
      ----------
      only_changed
        Return only changed values, or all of them?
        If True, return only the values, that differ from the defaults.
        If False, return all the values.
        The default value ``basic`` means, return all non-expert values
        and all the changed expert values.
        ``explicit`` means just the values, that were explicitly set (even
        if they are the same as the default value)
      """
      if not getter:
          getter = self.as_dict_getter(only_changed, generated, copy)
      return self._as_dict(getter)

  to_dict = as_dict

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

  def check_for_errors(self, validate='save', print=True):
      with warnings.catch_warnings(record = True) as lst:
          warnings.simplefilter("always", DataValidityWarning)
          self._validate(validate)
      if lst and print:
          for warning in lst:
              warnings.showwarning(
                  message=warning.message,
                  category=warning.category,
                  filename=warning.filename,
                  lineno=warning.lineno
              )

  def validate(self, why='save'):
      if why=='warning':
          return self._validate(why)
      with warnings.catch_warnings():
          warnings.simplefilter("error", DataValidityError)
          self._validate(why)

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
      if validate == 'warning':
          self.check_for_errors('save')
      else:
          self.validate(validate)

      if not hasattr(file, 'write'):
          with open(file, "w") as file:
              out=self._save_to_file(file)
              file.flush()
      else:
          out=self._save_to_file(file)
          file.flush()
      return out

  def to_string(self, *, validate:Union[str, bool]='warning'):
      """
      Return the configuration (problem definition) in a string.

      Parameters
      ----------
      Validate. How to validate before retrieving. See the method
      validate.
      Default 'warning' means the same as 'save', but only throw a warning.

      Returns
      -------
      configuration:str
      The configuration, as it should be saved in a configuration/problem definition file.
      """
      from io import StringIO
      s = StringIO()
      self.save_to_file(s, validate=validate)
      return s.getvalue()

  def _find_member(self, name, lower_case:bool=False, is_option=None):
      item = self._find_members(name, lower_case, is_option)
      try:
          return next(item)
      except StopIteration:
          return None


_help_warning_printed=False
