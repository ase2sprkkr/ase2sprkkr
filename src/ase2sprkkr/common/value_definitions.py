from .configuration_definitions import RealItemDefinition, Validator
from .options import Option
from .dangerous_values import DangerousValue
from .grammar_types import GrammarType, type_from_type, type_from_value, Array, QString
from .warnings import warnings, DataValidityError, ValidationResult
from .repetition import Repeated
from .decorators import cached_property
from .grammar import as_delimiter

import builtins
import copy as copy_module
from typing import Iterable, List, Optional, Sequence, Tuple, Union, Dict, Any
import numpy as np
import pyparsing as pp



def _validate_required(
    option: Option, container: Any, why: str
) -> Optional[DataValidityError]:
    """Validate requiredness of ``option`` in ``container`` for phase ``why``."""
    requirement = option._definition.is_required
    if why == "set" and not callable(requirement):
        return None
    required = option.is_required()
    if not required or (required == "save" and why != "save"):
        return None
    value = option(unpack=False, all_values=True)
    if value is not None or isinstance(value, DangerousValue):
        return None
    if required is True or required == "save":
        message = f"The value is required for {option._get_path()}, it can't be None"
    else:
        message = required
    return DataValidityError(message)


def _validate_numbering(
    option: Option, container: Any, why: str
) -> Optional[DataValidityError]:
    """Validate conditional numbering of ``option`` for phase ``why``."""
    if option._definition.is_repeated.numbering_condition is None:
        return None
    value = option(unpack=False, all_values=True)
    return option._definition.numbering_issue(option, value)


class ValueModifier:
    """If this class is given as a type of a Value, it will modify the definition
    of value somehow. It is responsibile to set the True type of the value"""

    validators = ()


class InheritingValueModifier(ValueModifier):
    """The definition of the value will be inherited from this class as well."""

    _enriching_classes = {}

    def modify_definition(self, definition):
        modifiers = (*getattr(definition, "_modifiers", ()), self)
        base_classes = getattr(
            definition, "_base_classes", (definition.__class__,)
        )
        base_classes = (type(self), *base_classes)
        definition.__class__ = self._modified_class(base_classes)
        definition._modifiers = modifiers
        return self.type

    def copy(self, **kwargs: Any) -> Any:
        """Copy a modified definition together with its modifier objects."""
        out = super().copy(**kwargs)
        out._modifiers = tuple(
            copy_module.copy(modifier) for modifier in self._modifiers
        )
        modifier_validators = (
            validator
            for modifier in out._modifiers
            for validator in modifier.validators
        )
        out.validators += tuple(
            validator
            for validator in modifier_validators
            if validator not in out.validators
        )
        return out

    @classmethod
    def _modified_class(cls, base_classes: Tuple[type, ...]) -> type:
        """Return the cached dynamic class combining ``base_classes``."""
        modified = cls._enriching_classes.get(base_classes)
        if modified is None:
            modified = type(
                "".join(item.__name__ for item in base_classes),
                base_classes,
                {
                    "__init__": base_classes[-1].__init__,
                    "__reduce_ex__": _reduce_modified_definition,
                    "_base_classes": base_classes,
                },
            )
            cls._enriching_classes[base_classes] = modified
        return modified


def _reduce_modified_definition(
    definition: "ValueDefinition", protocol: int
) -> tuple:
    """Describe ``definition`` for pickle protocol ``protocol``."""
    state = definition.__dict__.copy()
    modifiers = state.pop("_modifiers")
    return (
        _rebuild_modified_definition,
        (definition._base_classes, modifiers),
        state,
    )


def _rebuild_modified_definition(
    base_classes: Tuple[type, ...], modifiers: Sequence[ValueModifier]
) -> "ValueDefinition":
    """Reapply ``base_classes`` and ``modifiers`` before restoring state."""
    definition_class = base_classes[-1]
    instance = definition_class.__new__(definition_class)
    instance.__class__ = InheritingValueModifier._modified_class(base_classes)
    instance._modifiers = modifiers
    return instance


class ValueDefinition(RealItemDefinition):
    result_class = Option
    delimiter = pp.Empty().set_name("")
    intrinsic_validators = RealItemDefinition.intrinsic_validators + (
        _validate_required,
        _validate_numbering,
    )

    name_in_grammar = None
    is_generated = False
    is_validated = None  # default is is_generated

    item_type = "value"

    def __init__(
        self,
        name,
        type=None,
        default_value=None,
        default_value_from_container=None,
        written_name=None,
        fixed_value=None,
        init_by_default=False,
        result_is_visible=False,
        info=None,
        description=None,
        is_stored=None,
        is_hidden=False,
        is_optional=None,
        is_required=None,
        is_expert=False,
        is_repeated: Union[bool, str, RealItemDefinition.Repeated] = False,
        repeated_delimiter=None,
        repeated_count=None,
        is_always_added: bool = None,
        name_in_grammar=None,
        name_format=None,
        expert=None,
        write_condition=None,
        condition=None,
        warning_condition=None,
        validators: Optional[Union[Validator, Iterable[Validator]]] = (),
        result_class=None,
        delimiter=None,
        delimiter_grammar=None,
        indent=None,
    ):
        """
        Definition of a configuration value.

        Parameters
        ----------
        name: str | Tuple(str)
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
          Required options must have a value when parsing or saving (however,
          a required option can still be optional in the grammar if it has a
          default value). Callable requiredness is also checked while setting,
          so changing a selector cannot leave its dependent option invalid.
          If required = None, it is set to True if both the conditions are met:

           * the value is not expert
           * the optional is not True and the option has not a default_value

          If required is str, it means the same as True, the string will be
          used as error message.

          If it is callable, it is evaluated on demand with the Option as its
          argument. The runtime section is available as ``option._container``.

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

        repeated_delimiter: string
            The delimiter used between repeated instances of self, if not repeated_with_name

        repeated_count: string or int or callable
            The number of repetitions. If string, the repetition is given from the parsed value at that exact
            path; ``..`` selects the parent section (for example ``..NPAN``). If callable, the function is called
            on the values of the local options named by its parameters.

        repeated_with_name: bool
            If True, the whole value-name pair is repeated


        is_always_added
          If False, add the value, only if its value is not the default value.
          Default None means False for expert values, True for the others.

        validators: callable or iterable of callables
          Semantic validators receiving the option, its runtime section and
          validation reason.

        name_in_grammar: bool or None
          The value in the conf file is prefixed by <name><name_value_delimiter>
          If None, the default type value (type.name_in_grammar) is used

        name_format: str or None
          The way how the name is written

        expert: Optional[mixed]
          If not None, set ``is_expert`` to True, ``default_value`` to the given value and
          ``required`` to False. Note, that also ``type`` can be determined from such given
          ``default_value``.

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

        delimiter
           If not None, use the specified name_value_delimiter instead of common one
        """
        if callable(default_value_from_container):
            default_value = default_value_from_container
            default_value_from_container = True
        else:
            default_value_from_container = default_value_from_container
        self.default_value_from_container = default_value_from_container
        self.result_is_visible = result_is_visible

        if expert is not None:
            if type is None:
                type = expert
            else:
                default_value = expert
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

        modifier = type if isinstance(type, ValueModifier) else None
        if modifier:
            type = modifier.modify_definition(self)

        if default_value is None and not isinstance(type, (GrammarType, builtins.type)):
            self.type = type_from_value(type, type_map=self.type_from_type_map)
            default_value = None if isinstance(type, dict) else self.type.convert(type)
        else:
            self.type = type_from_type(type, type_map=self.type_from_type_map)
            if default_value is not None:
                if not callable(default_value):
                    default_value = self.type.convert(default_value)
            default_value = default_value

        if default_value is None and self.type.default_value is not None:
            default_value = self.type.default_value

        self._default_value = default_value

        assert isinstance(self.type, GrammarType), (
            "grammar_type (sprkkr.common.grammar_types.GrammarType descendat) required as a value type"
        )
        self.type.used_in_definition(self)

        self.grammar_type = self.type
        if self.is_repeated.is_array:
            self.type = Array(self.type)

        if is_required is None:
            is_required = not is_expert and (not is_optional and default_value is None)
        self.is_required = is_required

        if is_optional is None:
            is_optional = is_required is False

        super().__init__(
            name=name,
            written_name=written_name,
            is_optional=is_optional,
            is_hidden=is_hidden,
            is_expert=is_expert,
            name_in_grammar=name_in_grammar,
            info=info,
            description=description,
            name_format=name_format,
            write_condition=write_condition,
            condition=condition,
            warning_condition=warning_condition,
            validators=validators,
            result_class=result_class,
            is_repeated=is_repeated,
            repeated_delimiter=repeated_delimiter,
            repeated_count=repeated_count,
            repeated_with_name=True,
            indent=indent
        )
        if modifier:
            self.validators += modifier.validators

        if delimiter is not None:
            self.delimiter = as_delimiter(
                delimiter if delimiter_grammar is None else delimiter_grammar,
                None if delimiter_grammar is None else delimiter,
            )

    configuration_type_name = "OPTION"

    default_repeated_type = Repeated.REPEATED

    type_from_type_map = {}
    """ Redefine this in descendants, if you need to create different types that the defaults to be
  'guessed' from the default values """

    @cached_property
    def name_in_grammar(self):
        return self.type.name_in_grammar

    @property
    def default_value(self):
        defval = self._default_value
        if self.default_value_from_container:
            defval = (lambda d: lambda o: d(o._container))(defval)
        if self.result_is_visible:
            return lambda o: o._result if hasattr(o, "_result") else defval
        return defval

    @default_value.setter
    def default_value(self, val):
        self._default_value = val

    def allow_duplication(self):
        """Can be the item repeated in the output file"""
        return self.is_repeated and not self.is_repeated.is_numbered

    @property
    def is_independent_on_the_predecessor(self):
        """Some value have to be positioned just after their predecessor
        in the output.
        """
        return self.name_in_grammar or self.type.is_independent_on_the_predecessor

    def enrich(self, option):
        """The Option can be enriched by the definition, e.g. the docsting can be extended."""
        self.type.enrich(option)

    def data_description(self, verbose: Union[bool, str] = False, show_hidden=False, prefix: str = ""):
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
            value = getattr(self.default_value, "__doc__", "<function>")
        else:
            value = self.get_value()
        if value is not None:
            out += f" ≝ {value}"

        flags = []
        if self.is_optional:
            flags.append("optional")
        if self.is_hidden:
            flags.append("hidden")
        if self.is_expert:
            flags.append("expert")
        if self.is_expert == self.is_always_added:
            flags.append("always add" if self.is_expert else "add non-default")
        if self.is_repeated:
            flags.append("array")
        if self.is_fixed:
            flags.append("read_only")
        if flags:
            flags = ", ".join(flags)
            out += f"  ({flags})"

        if verbose:
            add = self.additional_data_description(prefix=prefix + self._description_indentation)
            if add:
                out += "\n"
                out += add
        return out

    def additional_data_description(self, verbose=False, show_hidden=False, prefix: str = "") -> str:
        """Return the additional runtime-documentation for the configuration value.
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
        """Hook called, when the object is assigned to the container (currently from the container
        constructor)"""
        super().added_to_container(container)
        self.type.added_to_container(container)

    def validate_type(self, item: bool):
        """Return the DataType against which should be data validated.

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

    def validate(self, opt, value, why="set", item=False):
        try:
            if value is None:
                return True
            if self.is_fixed and not np.array_equal(self.default_value, value):
                raise ValueError(
                    f"The value of {opt._get_path()} is required to be {self.default_value}, cannot set it to {value}"
                )
            self.validate_type(item).validate(
                value, why=why, option=opt
            )
        except ValueError as e:
            DataValidityError.warn(str(e))

    def numbering_issue(
        self, opt: Option, value: Any
    ) -> Optional[DataValidityError]:
        """Return an error for invalid conditional numbering, if any."""
        repeated = self.is_repeated
        if repeated.numbering_condition is None or value is None or repeated.numbering_for(opt):
            return None
        try:
            length = len(value)
        except TypeError:
            return None
        if length <= 1:
            return None
        message = (
            f"{opt._get_path()} is unnumbered for the current configuration "
            f"and therefore can contain only one value; got {length}"
        )
        return DataValidityError(message)

    def convert_value(self, opt: Option, value: Any, item: bool = False) -> Any:
        """Convert ``value`` for ``opt``; ``item`` selects element conversion."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", DataValidityError)
            try:
                return self.validate_type(item).convert(value)
            except DataValidityError:
                raise
            except (TypeError, ValueError) as error:
                raise DataValidityError(str(error)) from error

    def validate_converted(
        self, opt: Option, value: Any, why: str = "set", item: bool = False
    ) -> List[ValidationResult]:
        """Validate converted ``value`` for ``opt`` and phase ``why``.

        ``item`` selects element validation rather than whole-value validation.
        """
        return ValidationResult.collect(self.validate, opt, value, why, item)

    def convert_and_validate(
        self, opt: Option, value: Any, why: str = "set", item: bool = False
    ) -> Any:
        """Strictly convert and validate ``value`` using the supplied arguments."""
        value = self.convert_value(opt, value, item)
        ValidationResult.report(
            self.validate_converted(opt, value, why, item), "raise"
        )
        return value

    @property
    def value_name_format(self):
        return self.name_format

    @value_name_format.setter
    def value_name_format(self, value):
        self.name_format = value

    def __str__(self):
        out = "<{}: {}>".format(self.name, str(self.type))
        try:
            val = self.get_value()
        except Exception:
            val = "<ERRORNEOUS VALUE>"
        if val is not None:
            out += "={}".format(val)
        return out

    def __repr__(self):
        return str(self)

    def _grammar_of_value(self, delimiter, allow_dangerous=False):
        """Return grammar for the (possible optional) value pair"""
        type = self.grammar_type
        body = type.grammar(self.name)

        if self.is_fixed:

            def check_fixed(s, loc, x, body=body):
                if self.default_value.__class__ is np.ndarray:
                    eq = np.array_equal(x[0], self.default_value)
                else:
                    eq = x[0] == self.default_value
                if eq:
                    return x
                message = "The value of {} is {} and it should be {}".format(self.name, x[0], self.default_value)
                raise pp.ParseException(s, loc, message, body)

            body = body.copy().add_parse_action(check_fixed)

        if allow_dangerous and hasattr(self, "type_of_dangerous"):
            danger = pp.Forward()
            danger << self.type_of_dangerous.grammar(self.name + "_dangerous")
            danger.add_parse_action(lambda x: DangerousValue(x[0], self.type_of_dangerous, False))
            body = body ^ danger

        if delimiter:
            body = delimiter + body

        optional, df, _ = type.missing_value()
        if optional:
            body = pp.Optional(body).set_parse_action(lambda x: x or df)
        return body

    @property
    def _grammar(self):
        if not self.is_stored:
            return None
        return self._hooked_grammar

    def _create_grammar(self, allow_dangerous=False, name_in_grammar=None, name_value_delimiter=None, original=False):
        """Return a grammar for the name-value pair"""
        if self.output_definition is not self and not original:
            g = self.output_definition._grammar
            return g and g(allow_dangerous)

        name_in_grammar = name_in_grammar if name_in_grammar is not None else self.name_in_grammar
        if (name_value_delimiter is None and name_in_grammar) or name_value_delimiter is True:
            name_value_delimiter = self.delimiter

        body = self._grammar_of_value(name_value_delimiter, allow_dangerous)

        if name_in_grammar:
            nbody = self.formated_name.strip()
        else:
            nbody = ""

        if not self.grammar_type.missing_value()[0]:
            if nbody:
                nbody += str(name_value_delimiter) or " "
            nbody += self.grammar_type.grammar_name()

        out = self._tuple_with_my_name(body, has_value=self.type.has_value, name_in_grammar=name_in_grammar)
        out.set_name(nbody)
        return out

    def get_value(self, option=None):
        """Return the default or fixed value of this option.

        The function can accept the Option (which of the definition is): if the default value is given by callable,
        this argument is passed to it. (E.g. to set the default value using some properties obtained from the
        configuration objects.
        """
        if self.default_value is not None:
            if callable(self.default_value):
                return self.type.convert(self.default_value(option))
            return self.default_value
        return None

    def _save_to_file(self, file, option, always=False, name_in_grammar=None, delimiter=""):
        value, write = option._written_value(always)
        if write:
            return self.write(file, value, name_in_grammar, delimiter=delimiter, option=option)
        else:
            return

    def write(self, file, value, name_in_grammar=None, delimiter="", option=None):
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

        if self.is_repeated.numbering_condition is not None:
            if option is None:
                raise ValueError("Writing a NUMBERED_IF value requires its runtime Option")
            issue = self.numbering_issue(option, value)
            if issue:
                raise issue

        def write(name, value):
            if name_in_grammar:
                if delimiter:
                    file.write(str(delimiter))
                self.write_name(file, name)
                self.write_value(file, value, str(self.delimiter))
                return True
            else:
                if delimiter:
                    deli = str(delimiter) + self.prefix
                else:
                    deli = self.prefix
                return self.write_value(file, value, deli)

        name = self.formated_name
        if self.is_repeated:
            nmb = self.is_repeated.numbering_for(option)
            if self.is_repeated.type == self.Repeated.Type.DICT:
                if nmb == self.Repeated.Numbering.WITH_DEFAULT:
                    written = ((name + (str(i) if i != "def" else ""), v) for i, v in value.items())
                else:  # Dict has to be numbered
                    written = ((name + str(i), v) for i, v in value.items())
            else:
                if nmb:
                    written = ((name + str(i + 1), v) for i, v in enumerate(value))
                else:
                    written = ((name, v) for v in value)
            out = False
            for mname, val in written:
                if write(mname, val):
                    delimiter = self.container.delimiter
                    out = True
            return out
        else:
            return write(name, value)

    def write_value(self, file, value, delimiter=""):
        """Write given value and a given delimiter before the value to a file"""
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
            try:
                write_value = not (missing and df == value)
            except TypeError:  # unyt fix
                write_value = True
        else:
            value = None
            write_value = True

        if write_value:
            file.write(str(delimiter))
            type.write(file, value)
            return True
        else:
            return False

    def write_name(self, file, name, delimiter=""):
        file.write(str(delimiter) + self.prefix)
        file.write(name)

    def remove(self, name):
        del self.section[name]
        return self

    def _generic_info(self):
        return f"Configuration value {self.name}"

    def _get_init_args_for_copy(self, **kwargs) -> Dict[str, Any]:
        """
        Compute the values for creating the copy.

        Returns
        -------
        copy: Dict
           The returning dictionary has this structure:
           { name of the argument of the __init__ function : name of the object attribute }
        """
        out = super()._get_init_args_for_copy(**kwargs)
        if self.is_repeated.is_array:
            # ``self.type`` is the storage Array wrapper added by __init__.
            # Passing it back to the constructor would wrap it for a second
            # time; the grammar type is the original per-item type.
            out["type"] = self.grammar_type
        if self.is_fixed:
            out["fixed_value"] = out["default_value"]
        if "delimiter" in self.__dict__:
            out["delimiter"] = self.delimiter
        return out

    _copy_excluded_args = RealItemDefinition._copy_excluded_args + [
        "fixed_value",
        "result_is_visible",
        "delimiter",
        "delimiter_grammar",
    ]

    def copy_value(self, value, all_values=False):
        """Creates the copy of the value

        Parameters
        ----------
        values
          The value to copy

        all_values
          Wheter, for a numbered array, a whole dict is supplied
        """
        if not all_values or not self.is_repeated.is_dict:
            return self.type.copy_value(value)
        return {k: self.type.copy_value(v) for k, v in value.items()}

    def check_array_access(self):
        """Check, whether the option is array type (or repeated) and thus it can be accessed as array using []"""
        if not self.is_repeated and not self.type.array_access:
            raise TypeError("It is not allowed to access {self.get_path()} as array")
