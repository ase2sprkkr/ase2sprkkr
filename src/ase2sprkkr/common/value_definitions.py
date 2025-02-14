from .configuration_definitions import RealItemDefinition
from .options import Option, DangerousValue
from .grammar_types import GrammarType, type_from_type, type_from_value, Array, QString

import builtins
from typing import Union, Dict
import numpy as np
import pyparsing as pp
from .warnings import warnings, DataValidityError


class ValueModifier:
    """ If this class is given as a type of a Value, it will modify the definition
    of value somehow. It is responsibile to set the True type of the value """


class InheritingValueModifier(ValueModifier):
    """ The definition of the value will be inherited from this class as well. """

    _enriching_classes = {}

    def modify_definition(self, definition):
        dcls = definition.__class__
        cls = self._enriching_classes.get(dcls)
        if not cls:
            scls = self.__class__
            self._enriching_classes[dcls] = cls = \
                type(scls.__name__ + dcls.__name__, (scls, dcls), {})

        definition.__class__ = cls
        return self.type


class ValueDefinition(RealItemDefinition):

  result_class = Option

  name_in_grammar = None
  is_generated = False
  is_validated = None  # default is is_generated

  item_type = 'value'

  def __init__(self, name, type=None, default_value=None,
               default_value_from_container=None,
               written_name=None, alternative_names=None,
               fixed_value=None, init_by_default=False,
               result_is_visible=False, info=None, description=None,
               is_stored=None, is_hidden=False, is_optional=None, is_required=None,
               is_expert=False, is_repeated:Union[bool,str,RealItemDefinition.Repeated]=False,
               is_always_added:bool=None,
               name_in_grammar=None, name_format=None, expert=None,
               write_alternative_name:bool=False,
               write_condition=None, condition=None,
               result_class=None,
               ):
    """
    Definition of a configuration value.

    Parameters
    ----------
    name: str
      Name of the configuration value

    type: Optional[GrammarType|mixed]
      Configuration value data type.
      If it is set to anyting what is not derived from GrammarType, the given value is used as the default value
      and the data type is derived from it.
      If it is None, the default value have to be set using ``expert`` parameter.

    default: mixed
      Default value for the configuration option.
      Can accept callable (with option instance as an argument) - then the value will be determined
      at 'runtime' (possibly according to the other values of the section)

    written_name: str
      Name of the configuration value in the input file

    alternative_names: str or [str]
      Value can have an alternative name (that alternativelly denotes the value)

    fixed_value: mixed
      If it is given, this option have a fixed_value value (provided by this parameter),
      that can not be changed by an user.
      #TODO - currently, callback (as in default_value) is not supported

    init_by_default: bool
      If the value is not set, init it by default

    result_is_visible: bool
      If True, the result (see the result property of :class:`Option`) assigned to the
      option is visible if the value is get as its default value - in this mode the result
      can be used as a kind of default value, specific for given configuration object.

      If False, the result is visible to the user just using the result property
      and it should be some transformation of the value of the object (e.g.
      relative path of the given absolute path, id of assigned object etc.)

    is_stored: bool
      If True, the value is readed/writed from the output file.
      If None, set to False if the value is Generated.

    is_optional: bool or None
      If True, the value can be omited, if fixed order in the section is required
      None means True just if required is False (or it is determined to be False),
      see the ``required`` parameter.

    is_required: bool or callable or str
      Required option can not be set to None (however, a required one
      can be still be optional, if it has a default values).
      If required = None, it is set to True if both the conditions are met:

       * the value is not expert
       * the optional is not True and the option has not a default_value

      If required is str, it means the same as True, the string will be
      used as error message.

      If it is callable, it is evaluated on demand, with the Option as the
      argument of the function.

    is_hidden: bool
      The value is hidden from the user (no container.name access to the value).

    is_expert: Union[bool,mixed]
      Expert values are somewhat hidden (e.g. listed at end) from the user.
      Expert values are not exported to the result, if they are set to the
      default value.

    is_repeated
      Whether the value can apper more than once in the output. The result
      can be then (dense) array or (sparse) dict,
      see :class:`ValueDefinition.Repeated`.

    is_always_added
      If False, add the value, only if its value is not the default value.
      Default None means False for expert values, True for the others.

    name_in_grammar: bool or None
      The value in the conf file is prefixed by <name><name_value_delimiter>
      If None, the default type value (type.name_in_grammar) is used

    name_format: str or None
      The way how the name is written

    expert: Optional[mixed]
      If not None, set ``is_expert`` to True, ``default_value`` to the given value and
      ``required`` to False. Note, that also ``type`` can be determined from such given
      ``default_value``.

    write_alternative_name
       Wheter use the name or the (first) alternative name in the output.

    write_condition
       If defined, write the value, only if write_condition(the option) is True.

    condition
       If defined, the condition
        - the condition.parse_condition() is invoked, when given grammar element
          should be parsed. If it is False, the element is skipped
        - the condition() is invoked, when the elements of the container is listed
          to hide the inactive members

    result_class
       Redefine the class that holds data for this option/section
    """
    if default_value_from_container:
       default_value = lambda o: default_value_from_container(o._container)
    if expert is not None:
       if type is None:
          type = expert
       else:
          default_value=expert
       is_expert = True
       is_required = False
    elif type is None:
       raise TypeError("The data-type of the configuration value is required.")

    if is_always_added is None:
       self.is_always_added = not is_expert
    else:
       self.is_always_added = is_always_added

    self.init_by_default = init_by_default
    self.is_stored = not self.is_generated if is_stored is None else is_stored

    if fixed_value is None:
       self.is_fixed = False
    else:
       default_value = fixed_value
       self.is_fixed = True

    if isinstance(type, ValueModifier):
        type = type.modify_definition(self)

    if default_value is None and not isinstance(type, (GrammarType, builtins.type)):
       self.type = type_from_value(type, type_map = self.type_from_type_map)
       self.default_value = None if isinstance(type,dict) else self.type.convert(type)
    else:
       self.type = type_from_type(type, type_map = self.type_from_type_map)
       if default_value is not None:
          if not callable(default_value):
             default_value = self.type.convert(default_value)
       self.default_value = default_value

    assert isinstance(self.type, GrammarType), "grammar_type (sprkkr.common.grammar_types.GrammarType descendat) required as a value type"
    self.type.used_in_definition(self)

    self.is_repeated = self.Repeated.create(is_repeated)

    self.grammar_type = self.type
    if self.is_repeated.is_array:
        self.type = Array(self.type)

    if self.default_value is None and self.type.default_value is not None:
       self.default_value = self.type.default_value

    if result_is_visible:
       self.default_value = lambda o: o._result if hasattr(o, '_result') else default_value

    if is_required is None:
       is_required = not is_expert and (not is_optional and default_value is None)
    self.is_required = is_required

    if is_optional is None:
       is_optional = is_required is False

    super().__init__(
         name = name,
         written_name = written_name,
         alternative_names = alternative_names,
         is_optional = is_optional,
         is_hidden = is_hidden,
         is_expert = is_expert,
         name_in_grammar = name_in_grammar,
         info=info,
         description = description,
         write_alternative_name = write_alternative_name,
         name_format = name_format,
         write_condition = write_condition,
         condition = condition,
         result_class = result_class,
    )

    if self.name_in_grammar is None:
        self.name_in_grammar = self.type.name_in_grammar

    if self.is_repeated.is_numbered and not self.name_in_grammar:
       raise ValueError('Repeated numbered values have to have its name in the grammar')

  configuration_type_name = 'OPTION'

  type_from_type_map = {}
  """ Redefine this in descendants, if you need to create different types that the defaults to be
  'guessed' from the default values """

  def allow_duplication(self):
       """ Can be the item repeated in the output file """
       return self.is_repeated and not self.is_repeated.is_numbered

  @property
  def is_independent_on_the_predecessor(self):
      """ Some value have to be positioned just after their predecessor
      in the output.
      """
      return self.name_in_grammar or self.type.is_independent_on_the_predecessor

  def enrich(self, option):
      """ The Option can be enriched by the definition, e.g. the docsting can be extended. """
      self.type.enrich(option)

  def data_description(self, verbose:Union[bool,str]=False, show_hidden=False, prefix:str=''):
    """
    Return the description of the contained data type and their type.

    Parameters
    ----------
    verbose
      If ``False``, return only one-line string with a basic info.
      If ``True``, return more detailed informations.
      'verbose' means here the same thing as True

    show_hidden
      If `False``, do not show hidden members... which has no meaning for Values.

    prefix
      The string, with with each line will begin (commonly the spaces for the indentation).
    """
    out = f"{prefix}{self.name} : {self.type}"
    if callable(self.default_value):
        value = getattr(self.default_value,'__doc__', '<function>')
    else:
        value = self.get_value()
    if value is not None:
        out+=f" ≝ {value}"

    flags = []
    if self.is_optional:
       flags.append('optional')
    if self.is_hidden:
       flags.append('hidden')
    if self.is_expert:
       flags.append('expert')
    if self.is_expert == self.is_always_added:
       flags.append('always add' if self.is_expert else 'add non-default')
    if self.is_repeated:
       flags.append('array')
    if self.is_fixed:
       flags.append('read_only')
    if flags:
       flags = ', '.join(flags)
       out+= f"  ({flags})"

    if verbose:
       add = self.additional_data_description(prefix=prefix + self._description_indentation)
       if add:
          out+= '\n'
          out+= add
    return out

  def additional_data_description(self, verbose=False, show_hidden=False, prefix:str='') -> str:
    """ Return the additional runtime-documentation for the configuration value.
        E.g. return the possible choices for the value, etc...

        Parameters
        ----------
        verbose
          This parameter has no effect here. See :meth:`RealItemDefinition.data_description` for its explanation.

        show_hidden
          This parameter has no effect here. See :meth:`RealItemDefinition.data_description` for its explanation.

        prefix
          Prefix for the indentation of the description.

        Returns
        -------
        additional_data_description

          An additional description of the values accepted by this configuration option, retrieved from the documentation type.
    """
    return self.type.additional_description(prefix)

  def added_to_container(self, container):
       """ Hook called, when the object is assigned to the container (currently from the container
       constructor) """
       self.container = container
       self.type.added_to_container(container)

  def validate_type(self, item:bool):
      """ Return the DataType against which should be data validated.

      Parameters
      ----------
      item: if True, not the whole value is set, but only item of an array
            (in the case of repeated option, or e.g. the one with :class:`ase2sprkkr.common.grammar_types.Array`
            type)
      """
      if item:
          if self.grammar_type is self.type:
               return self.type.type
          return self.grammar_type
      return self.type

  def validate(self, opt, value, why='set', item=False):
      try:
          if value is None:
             req = opt.is_required
             if req:
                if req == 'save':
                    if why != 'save':
                        return True
                    else:
                        req = True
                if req is True:
                    raise ValueError(f"The value is required for {opt._get_path()}, it can't be None")
                else:
                    raise ValueError(req)
             return True
          if self.is_fixed and not np.array_equal(self.default_value, value):
              ValueError(f'The value of {opt._get_path()} is required to be {self.default_value}, cannot set it to {value}')
          self.validate_type(item).validate(value, self.get_path, why=why)
          self.validate_warning(value)
      except ValueError as e:
          DataValidityError.warn(str(e))

  def convert_and_validate(self, opt, value, why='set', item=False):
      with warnings.catch_warnings():
          warnings.simplefilter("error", DataValidityError)
          try:
              value = self.validate_type(item).convert(value)
              self.validate(opt, value, why, item)
              return value
          except ValueError as v:
              DataValidityError.warn(str(v))

  @property
  def value_name_format(self):
    return self.name_format

  @value_name_format.setter
  def value_name_format(self, value):
      self.name_format = value

  def __str__(self):
    out="<{}: {}>".format(self.name, str(self.type))
    try:
      val = self.get_value()
    except Exception:
      val = "<ERRORNEOUS VALUE>"
    if val is not None:
      out+= "={}".format(val)
    return out

  def __repr__(self):
    return str(self)

  def _grammar_of_value(self, delimiter, allow_dangerous=False):
    """ Return grammar for the (possible optional) value pair """
    type = self.grammar_type
    body = type.grammar(self.name)

    if self.is_fixed:
      def check_fixed(s, loc, x, body=body):
          if self.default_value.__class__ is np.ndarray:
             eq = np.array_equal(x[0], self.default_value)
          else:
             eq = x[0]==self.default_value
          if eq:
             return x
          message="The value of {} is {} and it should be {}".format(self.name, x[0], self.default_value)
          raise pp.ParseException(s,loc,message, body)
      body=body.copy().addParseAction(check_fixed)

    if allow_dangerous and hasattr(self, 'type_of_dangerous'):
        danger = pp.Forward()
        danger << self.type_of_dangerous.grammar(self.name + '_dangerous')
        danger.addParseAction(lambda x: DangerousValue(x[0], self.type_of_dangerous, False))
        body = body ^ danger

    if delimiter:
       body = delimiter + body

    optional, df, _ = type.missing_value()
    if optional:
      body = pp.Optional(body).setParseAction( lambda x: x or df )
    return body

  @property
  def _grammar(self):
     if not self.is_stored:
        return None
     return self._hooked_grammar

  def _create_grammar(self, allow_dangerous=False, name_in_grammar=None, name_value_delimiter=None, original=False):
    """ Return a grammar for the name-value pair """
    if self.output_definition is not self and not original:
        g = self.output_definition._grammar
        return g and g(allow_dangerous)

    name_in_grammar = name_in_grammar if name_in_grammar is not None else self.name_in_grammar
    if (name_value_delimiter is None and name_in_grammar) or name_value_delimiter is True:
       name_value_delimiter=self.grammar_of_delimiter

    body = self._grammar_of_value(name_value_delimiter, allow_dangerous)

    if name_in_grammar:
       nbody = self.formated_name.strip()
    else:
       nbody = ''

    if not self.grammar_type.missing_value()[0]:
        if nbody:
            nbody +=str(name_value_delimiter) or ' '
        nbody+=self.grammar_type.grammar_name()

    out = self._tuple_with_my_name(body, has_value=self.type.has_value,
                                         name_in_grammar=name_in_grammar)
    out.setName(nbody)
    return out

  def get_value(self, option=None):
     """ Return the default or fixed value of this option.

     The function can accept the Option (which of the definition is): if the default value is given by callable,
     this argument is passed to it. (E.g. to set the default value using some properties obtained from the
     configuration objects.
     """
     if self.default_value is not None:
        if callable(self.default_value):
           return self.type.convert(self.default_value(option))
        return self.default_value
     return None

  def _save_to_file(self, file, option, always=False, name_in_grammar=None, delimiter=''):
      value, write = option._written_value(always)
      if write:
          return self.write(file, value, name_in_grammar, delimiter=delimiter)
      else:
          return

  def write(self, file, value, name_in_grammar=None, delimiter=''):
     """
     Write the option to the open file

     Parameters
     ----------
     file
      The file to write to.

     value
      The value to write. It can be instance of DangerousValue, in such case
      It's own type is used to write the value.
     """
     if name_in_grammar is None:
         name_in_grammar = self.name_in_grammar

     def write(name, value):
         if name_in_grammar:
            if delimiter:
                file.write(delimiter)
            self.write_name(file, name)
            self.write_value(file, value, self.name_value_delimiter)
            return True
         else:
            if delimiter:
                deli = delimiter + self.prefix
            else:
                deli = self.prefix
            return self.write_value(file, value, deli)

     name = self.formated_name
     if self.is_repeated:
        nmb = self.is_repeated.is_numbered
        if self.is_repeated.type == self.Repeated.Type.DICT:
            if nmb == self.Repeated.Numbering.WITH_DEFAULT:
                written = ( (name + (str(i) if i!='def' else ''), v) for i, v in value.items() )
            else:  # Dict has to be numbered
                written = ( (name + str(i) , v) for i, v in value.items() )
        else:
            if nmb:
                written = ( (name + str(i + 1), v) for i, v in enumerate(value) )
            else:
                written = ( (name, v) for v in value )
        out = False
        for mname, val in written:
            if write(mname, val):
                delimiter = self.container.delimiter
                out = True
        return out
     else:
        return write(name, value)

  def write_value(self, file, value, delimiter=''):
     """ Write given value and a given delimiter before the value to a file """
     dangerous = isinstance(value, DangerousValue)
     if dangerous:
        type = value.value_type
        if not type:
           if hasattr(self, "type_of_dangerous"):
              type = getattr(self.type_of_dangerous, "string_type", QString)
           else:
              type = QString
        value = value()

     else:
        type = self.grammar_type

     if type.has_value:
         if value is None and not dangerous:
           value = self.get_value()
         if value is None:
           return False
         missing, df, _ = type.missing_value()
         write_value = not ( missing and df == value )
     else:
        value=None
        write_value=True

     if write_value:
         file.write(delimiter)
         type.write(file, value)
         return True
     else:
         return False

  def write_name(self, file, name, delimiter=''):
      file.write(delimiter + self.prefix)
      file.write(name)

  def remove(self, name):
     del self.section[name]
     return self

  def _generic_info(self):
      return f"Configuration value {self.name}"

  @property
  def can_be_repeated(self):
      return bool(self.is_repeated)

  def _get_copy_args(self)->Dict[str, str]:
       """
       Compute the dictionary that defines the attributes to create a copy of this object.

       Returns
       -------
       copy: Dict
          The returning dictionary has this structure:
          { name of the argument of the __init__ function : name of the object attribute }
       """
       out = super()._get_copy_args()
       if self.is_fixed:
           out['fixed_value'] = out['default_value']
       return out

  _copy_excluded_args = RealItemDefinition._copy_excluded_args + ['fixed_value', 'result_is_visible', 'default_value_from_container']

  def copy_value(self, value, all_values=False):
      """ Creates copy of the value

      Parameters
      ----------
      values
        The value to copy

      all_values
        Wheter, for a numbered array, a whole dict is supplied
      """
      if not all_values or not self.is_repeated.is_dict:
          return self.type.copy_value(value)
      return { k:self.type.copy_value(v) for k,v in value.items() }

  def check_array_access(self):
      """ Check, whether the option is array type (or repeated) and thus it can be accessed as array using [] """
      if not self.is_repeated and not self.type.array_access:
          raise TypeError('It is not allowed to access {self.get_path()} as array')
