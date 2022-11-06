""" The classes for storing one configuration value. """
from typing import Union
from ..common.grammar_types import mixed
from .configuration import Configuration
from ..common.misc import as_integer

class Option(Configuration):
  """ Class for one option (a configuration value) of SPRKKR - either to
  be used as a part of InputParameters or Potential configuration.
  Usage:

  >>> from ase2sprkkr.sprkkr.calculator import SPRKKR
  >>> calculator = SPRKKR()
  >>> conf = calculator.input_parameters
  >>> conf.ENERGY.ImE = 5.
  >>> conf.ENERGY.ImE()
  5.0
  >>> conf.ENERGY.ImE.info
  'Configuration value ImE'
  >>> conf.ENERGY.ImE.help()                     # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
  Configuration value ImE
  <BLANKLINE>
  ImE : Energy (<Real> [Ry|eV]) ≝ 0.0  (optional)
  """
  def __init__(self, definition, container=None, value=None):
      """"
      Parameters
      ----------
      definition: ValueDefinition
          The value type of the option and its format (in potential and/or task file)

      container:
          The container, that owns the object

      value: mixed
          The value of the option.
      """
      super().__init__(definition, container)
      self._value = value
      self._hook = None
      self._definition.type.enrich(self)

  def __call__(self, all_values:bool=False):
      """
      Return the value of the option.

      Parameters
      ----------
      all_values: For numbered_array (see :class:`ConfigurationDefinition.is_numbered_array <ConfigurationDefinition>`,
      pass True as this argument to obtain array of all values. If False (the default) is given,
      only the 'wildcard' value (i.e. the one without array index, which is used for the all values
      not explicitly specified) is returned.
      """
      if self._value is not None:
          if self._definition.is_numbered_array and not all_values:
              return self._value.get('def', self.default_value)
          return self._value
      if self._definition.is_numbered_array and all_values and self.default_value is not None:
          return { 'def' : self.default_value }
      return self.default_value

  @property
  def default_value(self):
      """ Return default value for the option.
      The function is here, and not in the definition, since the default value can be given
      by callable, that accepts the Option as argument. This possibility is used in ase2sprkkr,
      when the default values of some options are generated from the underlined Atoms object
      """
      return self._definition.get_value(self)

  def set(self, value, *, unknown=None):
      """
      Set the value of the option.

      Parameters
      ----------
      value: mixed
        The new value of the option.

      unknown: str or None
        A dummy argument to make the method compatibile with
        ase2sprkkr.sprkkr.common.configuration_containers.ConfigurationContainer.set()
      """
      if value is None:
          return self.clear()
      if self._definition.is_numbered_array:
        if isinstance(value, dict):
           self.clear(do_not_check_required=value, call_hooks=False)
           for k,v in value.items():
               self._set_item(k, v)
        else:
           self._set_item('def', v)
      else:
        self._value = self._definition.convert_and_validate(value)
      self._post_set()

  def _post_set(self):
      """ Thus should be called after all modifications """
      if hasattr(self,'_result'):
        del self._result
      if self._hook:
        self._hook(self)

  def _check_array_access(self):
      """ Check, whether the option is numbered array and thus it can be accessed as array using [] """
      if not self._definition.is_numbered_array:
          raise TypeException('It is not allowed to access {self._get_path()} as array')

  def __setitem__(self, name, value):
      """ Set an item of a numbered array. If the Option is not a numbered array, throw an Exception. """
      self._check_array_access()

      if isinstance(name, (list, tuple)):
          for n in name:
            self._set_item(n, value)
      elif isinstance(name, slice):
         if slice.stop is None:
             raise KeyError("To get/set values in a numbered array using slice, you have to specify the end index of the slice")
         for n in range(name.stop)[name]:
             self._set_item(n, value)
      else:
         self._set_item(name, value)
      self._post_set()

  def _set_item(self, name, value):
      """ Set a single item of a numbered array. For internal use - so no sanity checks """
      if self._value is None:
         self._value = {}
      if name != 'def':
         try:
            name = as_integer(name)
         except TypeError as e:
            raise KeyError('Numbered array indexes can be only integers, lists or slices') from e

      if value is None:
         del self._value[name]
         if not self._value:
            self._value = None
      else:
         self._value[name] = self._definition.convert_and_validate(value)

  def __getitem__(self, name):
      """ Get an item of a numbered array. If the Option is not a numbered array, throw an Exception. """
      self._check_array_access()
      if isinstance(name, (list, tuple)):
          return [ self._getitem(n) for n in name ]
      elif isinstance(name, slice):
         if slice.stop is None:
             raise KeyError("To get/set values in a numbered array using slice, you have to specify the end index of the slice")
         return [ self._getitem(n) for n in range(name.stop)[name] ]
      return self._getitem(name)

  def _getitem(self, name):
      """ Get a single item from a numbered array. For internal use - so no sanity checks """
      if name != 'def':
         try:
            name = as_integer(name)
         except TypeError as e:
            raise KeyError('Numbered array indexes can be only integers, lists or slices') from e
      return None if self._value is None else self._value.get(name, None)


  def __hasitem__(self, name):
      self._check_array_access()
      return self._value is not None and name in self._value

  def get(self):
      """ Return the value of self """
      return self()

  @property
  def result(self):
      """ Return the result value.

      In some cases, the value of an option have to be translated for the output.
      E.g. the site can be given as site object, but the integer index is
      required in the output.

      In a such case, this property can be utilized: the value of the option is
      retained as is and the transformed value is stored in the result.
      """
      if hasattr(self, '_result'):
          return self._result
      return self(all_values=True)

  @result.setter
  def result(self, value):
      self._result = value

  def clear(self, do_not_check_required=False,call_hooks=True):
      """ Clear the value: set it to None """
      if not self._definition.type.has_value:
         return
      if self._definition.default_value is None and not do_not_check_required and self._definition.required:
         raise ValueError(f'Option {self._get_path()} must have a value')
      self._value = None
      if call_hooks:
        self._post_set()

  def is_changed(self) -> bool:
      """ True, if the value is set and the value differs from the default """
      return self.value_and_changed()[1]

  def _save_to_file(self, file):
      """ Write the name-value pair to the given file, if the value
      is set. """
      d = self._definition
      if not d.type.has_value:
          return d.write(file, None)
      value = self.result
      if value is None or (d.is_expert and self.is_it_the_default_value(value)):
          return
      return d.write(file, value)

  def validate(self, why='save'):
      d = self._definition
      if d.type.has_value:
        value = self()
        if value is None:
           if not d.is_optional:
               name = self._get_root_container()
               raise Exception(f'Value {self._get_path()} is None and it is not an optional value. Therefore, I cannot save the {name}')
        else:
           d.type.validate(value, why)

  @property
  def name(self):
      return self._definition.name

  def as_dict(self, only_changed: Union[bool,str]='basic'):
      d = self._definition
      if only_changed and (only_changed!='basic' or d.is_expert):
           v,c = self.value_and_changed()
           return v if c else None
      return self(all_values=True)

  def value_and_changed(self):
      """ Return value and whether the value was changed

          Returns
          -------
          value:mixed
            The value of the options (return all values for 'numbered array')

          changed:bool
            Whether the value is the same as the default value or not
      """
      if self._value is not None:
         return self._value, not self.is_it_the_default_value(self._value)
      if self._definition.is_numbered_array:
         return {'def' : self.default_value}, False
      else:
         return self.default_value, False

  def is_it_the_default_value(self, value):
      """ Return, whether the given value is the default value. For
      numbered array, only the wildcard value can be set and this value
      have to be the same as the default. """

      default = self.default_value
      d = self._definition
      if d.is_numbered_array:
          return 'def' in value and len(value) == 1 and \
                  d.type.is_the_same_value(value['def'], default)
      else:
          return d.type.is_the_same_value(value, default)


  def _find_value(self, name):
      if self._definition.name == name:
         return self

  def get_path(self):
      return self._get_path()

  def __repr__(self):
      v = self._value
      if v is None:
         v = self.default_value
         if v:
            o=' (default)'
            v=' = '+str(v)
         else:
            o='out'
            v=''
      else:
          v=' = ' + str(v)
          o=''

      return f"<Option {self._get_path()} of type {self._definition.type} with{o} value{v}>"

class CustomOption(Option):
  """ An user-added option (configuration value). It can be removed from the section. """

  def remove(self):
      """ Remove me from my "parent" section """
      self._container.remove_member(self._definition.name)

  @classmethod
  def factory(cls, value_definition, type = mixed):
      """ Returns factory function for the given value definition """

      def create(name, section):
          definition = value_definition(name, type)
          definition.removable = True
          return cls(definition, section)

      create.grammar_type = type
      return create
