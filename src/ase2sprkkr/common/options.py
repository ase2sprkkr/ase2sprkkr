"""The classes for storing one configuration value."""

from __future__ import annotations

import copy
import numpy as np
from typing import Any, Optional, Union

from ..common.grammar_types import mixed
from .configuration import Configuration, UnknownMemberPolicy
from ..common.misc import as_integer
from .dangerous_values import DangerousValue
from .configuration_transaction import ConfigurationTransaction
from .warnings import (
    DataValidityError,
    InvalidValuePolicy,
    ReportInvalidPolicy,
    RetainInvalidPolicy,
    ValidationResult,
)


class _DiscardInvalidValue(Exception):
    """The proposed value was reported but must not be staged."""


class BaseOption(Configuration):
    """A base placeholder for a leaf element of a grammar file,
    both the a-value-holding ones (:class:`Option`) and
    dummy ones (Dummy)
    """

    def _save_to_file(self, file, always=False, name_in_grammar=None, delimiter=""):
        """Write the name-value pair to the given file, if the value
        is set."""
        return self._definition.output_definition._save_to_file(file, self, always, name_in_grammar, delimiter)

    def _find_members(self, name, is_option=None, lower_case=True):
        if self._definition.has_name(name, lower_case) and is_option is not False:
            yield self

    def get_path(self):
        return self._get_path()

    def _as_dict(self, get):
        return None

    def stage_clear(
        self,
        transaction: ConfigurationTransaction,
        *,
        check_required: bool = True,
    ) -> bool:
        """Return false because a placeholder has no state to stage."""
        return False


class Dummy(BaseOption):
    def _validate(self, why: str) -> bool:
        """Accept validation phase ``why`` because a dummy has no value."""
        return True

    def has_any_value(self):
        return False

    def __repr__(self):
        return f"<DUMMY {self._definition.name}>"


class DummyStub(Dummy):
    def _as_dict(self, get):
        if not self._definition.allowed(self._container):
            return None
        return get(self._container[self._definition.item])


class Option(BaseOption):
    """Class for one option (a configuration value) of SPRKKR - either to
    be used as a part of InputParameters or Potential configuration.
    Usage:

    >>> from ase2sprkkr.sprkkr.calculator import SPRKKR
    >>> calculator = SPRKKR()
    >>> conf = calculator.input_parameters
    >>> conf.ENERGY.ImE = 5.0
    >>> conf.ENERGY.ImE()
    unyt_quantity(5., 'Ry')
    >>> conf.ENERGY.ImE.info
    'Configuration value ImE'
    >>> conf.ENERGY.ImE.help()  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    Configuration value ImE
    <BLANKLINE>
    ImE : Energy (<Real> [Ry|eV]) ≝ 0.0 Ry (optional)
    >>> conf.ENERGY.ImE.set_dangerous("1J")
    >>> conf.ENERGY.ImE()
    '1J'
    """

    class Change:
        """Original state and finalization of a staged option value."""

        _MISSING_RESULT = object()

        def __init__(self, option: "Option") -> None:
            """Capture the current state and hook of ``option``."""
            self.option = option
            self.hook = option._hook
            self.old_value = option._value
            self.old_result = getattr(
                option, "_result", self._MISSING_RESULT
            )

        def commit(self) -> Any:
            if self.hook:
                return self.option, self.hook

        def rollback(self) -> None:
            self.option._value = self.old_value
            if self.old_result is self._MISSING_RESULT:
                self.option.__dict__.pop("_result", None)
            else:
                self.option._result = self.old_result

    class IndexedChange(Change):
        """Rollback state for an in-place indexed value change."""

        def __init__(self, option: "Option", key: Any) -> None:
            """Capture ``key`` and the surrounding state of ``option``."""
            super().__init__(option)
            self.key = key
            self.target = option(unpack=False, all_values=True)
            self.old_item = copy.deepcopy(self.target[key])

        def rollback(self) -> None:
            self.target[self.key] = self.old_item
            super().rollback()

    _NO_KEY = object()

    def __init__(self, definition, container=None, value=None):
        """ "
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
        self._hook = None
        self._definition.enrich(self)
        self._value = value

    def _value_or_default(self):
        d = self._definition
        if d.is_generated:
            return d.get_generated(self)
        if hasattr(self, "_result"):
            return self._result
        if self._value is not None:
            return self._value
        value = self.default_value
        if self._definition.is_repeated.is_dict and value is not None:
            return {"def": value}
        return value

    def __call__(self, all_values: bool = False, unpack=True):
        """
        Return the value of the option.

        Parameters
        ----------
        all_values: Control the behavior for the dict_like repeated values
        (see `is_repeated` attribute of :class:`ConfigurationDefinition`).
        Pass True as this argument to obtain dictionary
        of all values. If False (the default) is given, only the 'wildcard' value
        (i.e. the one without array index, which is used for the all values not explicitly specified)
        is returned.
        """
        d = self._definition
        value = self._value_or_default()
        if not d.is_generated and d.init_by_default and self._value is None:
            self._value = self._pack_value(value)
        if isinstance(value, DangerousValue) and unpack:
            value = value()
        if d.is_repeated.is_dict and not all_values:
            value = (
                value.get("def", self.default_value)
                if value is not None
                else self.default_value
            )
        if unpack:
            value = self._unpack_value(value)
        return value

    def is_dangerous(self):
        """Return, whether the option is set to a dangerous value, i.e. a value
        that bypass the validation."""
        return isinstance(self._value, DangerousValue)

    def set_dangerous(self, value, index=None):
        """Set the option to a dangerous value - i.e. to a value that bypass the
        type and value checks and enforcing.

        However, the type of such value is still checked by the proper mixed type.
        To completly bypass the check, set the value to an instance of DangerousValue
        class directly.
        """
        value = self._create_dangerous_value(value)
        if index is not None:
            self[index] = value
        else:
            self.set(value)

    def _create_dangerous_value(self, value):
        return DangerousValue(value, self._definition.type_of_dangerous)

    @property
    def default_value(self):
        """Return default value for the option.
        The function is here, and not in the definition, since the default value can be given
        by callable, that accepts the Option as argument. This possibility is used in ase2sprkkr,
        when the default values of some options are generated from the underlined Atoms object
        """
        return self._definition.get_value(self)

    def set(
        self,
        value: Any,
        *,
        unknown: UnknownMemberPolicy = None,
        retain_invalid: Optional[RetainInvalidPolicy] = None,
        report_invalid: Optional[ReportInvalidPolicy] = None,
    ) -> None:
        """Set a value using independent retention and reporting policies.

        ``retain_invalid`` accepts ``"none"``, ``"typed"`` or ``"all"``.
        The latter stores unconvertible values as :class:`DangerousValue`.
        ``report_invalid`` accepts ``"raise"``, ``"warn"`` or ``"ignore"``.
        ``unknown`` is accepted for the common container/option assignment
        interface; leaf options do not perform member lookup.
        """
        with self._mutation(
            "set", retain_invalid, report_invalid
        ) as (transaction, _policy):
            if value is None and not self._definition.is_generated:
                self.stage_clear(transaction)
            else:
                self.stage(
                    transaction, value, unknown=unknown
                )

    def stage(
        self,
        transaction: ConfigurationTransaction,
        value: Any,
        *,
        unknown: UnknownMemberPolicy = None,
        key: Any = _NO_KEY,
    ) -> bool:
        """Stage ``value`` or indexed ``key`` in ``transaction``.

        Conversion and validation issues are recorded in the policy owned by
        ``transaction``. ``unknown`` is accepted for compatibility with
        container staging.
        """
        policy = transaction.policy
        definition = self._definition
        if definition.is_generated:
            if key is self._NO_KEY:
                definition.stage_generated(self, transaction, value)
            else:
                definition.stage_generated(self, transaction, value, key)
            transaction._finish_stage()
            return True
        if key is not self._NO_KEY:
            return self._stage_item(transaction, key, value)

        def pack(item: Any) -> Any:
            return self._pack_proposed_value(item, policy)

        try:
            if value is None:
                converted = None
            elif definition.is_repeated.is_dict:
                if isinstance(value, dict):
                    converted = {}
                    rejected = False
                    for item_key, item in value.items():
                        try:
                            converted[item_key] = pack(item)
                        except _DiscardInvalidValue:
                            rejected = True
                    if rejected:
                        raise _DiscardInvalidValue
                else:
                    current = self(unpack=False, all_values=True)
                    converted = dict(current) if isinstance(current, dict) else {}
                    try:
                        converted["def"] = pack(value)
                    except (TypeError, ValueError):
                        converted.update(
                            (index, pack(item))
                            for index, item in enumerate(value, 1)
                        )
            elif definition.is_repeated:
                if isinstance(value, (str, bytes)):
                    value = (value,)
                else:
                    try:
                        iter(value)
                    except TypeError:
                        value = (value,)
                converted_items = [pack(item) for item in value]
                if (
                    isinstance(value, np.ndarray)
                    or definition.type.array_access
                    or len(converted_items) > 1
                ):
                    converted = np.asarray(converted_items)
                else:
                    converted = converted_items
            else:
                converted = pack(value)
        except _DiscardInvalidValue:
            transaction._finish_stage()
            return False

        change = self.Change(self)
        self._value = converted
        self.__dict__.pop("_result", None)
        transaction.push(change)
        transaction._finish_stage()
        return True

    def _stage_item(
        self, transaction: ConfigurationTransaction,
        key: Any,
        value: Any,
    ) -> bool:
        """Stage converted ``value`` at ``key`` in ``transaction``."""
        definition = self._definition
        policy = transaction.policy
        try:
            converted = self._pack_proposed_value(
                value,
                policy,
                item=not definition.is_repeated and definition.type.array_access,
            )
        except _DiscardInvalidValue:
            transaction._finish_stage()
            return False

        if self._value is None or hasattr(self, "_result"):
            current = copy.deepcopy(self._value_or_default())
            current = self._unpack_value(current)
            change = self.Change(self)
            self._value = self._pack_value(current)
            self.__dict__.pop("_result", None)
            target = self._value
        else:
            change = self.IndexedChange(self, key)
            target = self._value

        try:
            target[key] = converted
        except Exception:
            change.rollback()
            raise

        transaction.push(change)
        self.__dict__.pop("_result", None)
        transaction._finish_stage()
        return True

    def add_hook(self, hook):
        self._hook = hook

    def _check_array_access(self):
        """Check, whether the option is array type (or repeated) and thus it can be accessed as array using []"""
        return self._definition.check_array_acces(self)

    def __setitem__(self, name, value):
        """Set an item of a numbered array. If the Option is not a numbered array, throw an Exception."""
        d = self._definition
        if d.is_generated:
            with self._mutation("set") as (transaction, _policy):
                self.stage(transaction, value, key=name)
            return

        d.check_array_access()
        if not d.is_repeated.is_dict:
            with self._mutation("set") as (transaction, _policy):
                self.stage(transaction, value, key=name)
            return

        current = self(unpack=False, all_values=True)
        proposed = {} if current is None else d.copy_value(current, all_values=True)

        def set_item(key: Any, item: Any) -> None:
            if not (d.is_repeated.is_numbered.has_default and key == "def"):
                try:
                    key = as_integer(key)
                except TypeError as exc:
                    raise KeyError("Numbered array indexes can be only integers, lists or slices") from exc
                if key < 1:
                    raise KeyError("Numbered array indexes has to be greater than zero")
            if item is None:
                proposed.pop(key, None)
            else:
                proposed[key] = item

        if isinstance(name, (list, tuple)):
            for key in name:
                set_item(key, value)
        elif isinstance(name, slice):
            try:
                count = len(value)
                step = name.step or 1
                start = name.start or 1
                stop = name.stop or start + step * count
                for key, item in zip(range(start, stop, step), value):
                    set_item(key, item)
            except (TypeError, ValueError):
                if name.stop is None:
                    raise KeyError(
                        "To set a numbered array slice to one value, the slice end has to be specified"
                    )
                for key in range(name.start or 1, name.stop, name.step or 1):
                    set_item(key, value)
        else:
            set_item(name, value)
        self.set(proposed or None)

    def __getitem__(self, name):
        """Get an item of a numbered array. If the Option is not a numbered array, throw an Exception."""
        d = self._definition
        if d.is_generated:
            return d.get_generated(self, name)
        d.check_array_access()
        if not d.is_repeated.is_dict:
            return self()[name]

        if isinstance(name, (list, tuple)):
            return [self._getitem(n) for n in name]
        elif isinstance(name, slice):
            if name.stop is None:
                if self._value is None:
                    stop = 2
                else:
                    try:
                        stop = max(i for i in self._value if i != "def") + 1
                    except ValueError:
                        stop = 2
                name = slice(max(1, name.start or 1), stop, name.step)
            return [self._getitem(n) for n in range(name.start, name.stop, name.step or 1)]
        return self._getitem(name)

    def _getitem(self, name):
        """Get a single item from a numbered array. For internal use - so no sanity checks"""
        if name != "def":
            try:
                name = as_integer(name)
            except TypeError as e:
                raise KeyError("Numbered array indexes can be only integers, lists or slices") from e
        if self._value is None:
            return self.default_value
        if name in self._value:
            out = self._value[name]
        elif "def" in self._value:
            out = self._value["def"]
        else:
            return self.default_value
        return self._unpack_value(out)

    def _unpack_value(self, value):
        """Unpack potentionally dangerous values."""
        if isinstance(value, DangerousValue):
            value = value()
        r = self._definition.is_repeated
        if r:
            if r.is_dict and isinstance(value, dict):
                value = {i: v() if isinstance(v, DangerousValue) else v for i, v in value.items()}
            elif isinstance(value, list):
                value = [v() if isinstance(v, DangerousValue) else v for v in value]
        return value

    def _pack_value(self, value):
        """Validate the value, if it's to be."""
        if isinstance(value, DangerousValue):
            """ The dangerous value is immutable, checked during its creation """
            return value
        return self._definition.convert_and_validate(self, value)

    def _pack_proposed_value(
        self,
        value: Any,
        policy: InvalidValuePolicy,
        *,
        item: bool = False,
    ) -> Any:
        """Pack ``value`` according to the transaction's policy."""
        if isinstance(value, DangerousValue):
            return value
        definition = self._definition
        try:
            converted = definition.convert_value(self, value, item=item)
        except DataValidityError as issue:
            if policy.retain_all:
                policy.add(issue)
                return DangerousValue(
                    value,
                    getattr(definition, "type_of_dangerous", None),
                    validate=False,
                )
            policy.discard(issue)
            raise _DiscardInvalidValue from issue

        issues = definition.validate_converted(self, converted, item=item)
        retain = policy.retain_typed or not any(
            isinstance(issue, DataValidityError) for issue in issues
        )
        record = policy.add if retain else policy.discard
        for issue in issues:
            record(issue)
        if not retain:
            raise _DiscardInvalidValue
        return converted

    def __hasitem__(self, name):
        d = self._definition
        d.check_array_access()
        if not d.is_repeated.is_dict:
            return name in self()

        if self._value is None:
            return False
        value = self._unpack_value(self._value)
        return name in value

    def get(self):
        """Return the value of self"""
        return self()

    @property
    def result(self):
        """Return the result value.

        In some cases, the value of an option have to be translated for the output.
        E.g. the site can be given as site object, but the integer index is
        required in the output.

        In a such case, this property can be utilized: the value of the option is
        retained as is and the transformed value is stored in the result.
        """
        if hasattr(self, "_result"):
            if isinstance(self._result, DangerousValue):
                return self._result.value
            return self._result
        return self(all_values=True)

    @result.setter
    def result(self, value):
        self._result = value

    def clear_result(self):
        if hasattr(self, "_result"):
            del self._result

    def stage_clear(
        self,
        transaction: ConfigurationTransaction,
        *,
        check_required: bool = True,
    ) -> bool:
        """Stage clearing in ``transaction`` without global validation.

        ``check_required`` rejects clearing a required option without a default.
        """
        policy = transaction.policy
        if self._definition.is_generated:
            return self.stage(transaction, None)
        if not self._definition.type.has_value:
            return False
        if (
            self._definition.default_value is None
            and check_required
            and self.is_required()
        ):
            issue = DataValidityError(
                f"Option {self._get_path()} must have a value."
            )
            if policy.retain_typed:
                policy.add(issue)
                return self.stage(transaction, None)
            policy.discard(issue)
            transaction._finish_stage()
            return False
        return self.stage(transaction, None)

    def is_changed(self) -> bool:
        """True, if the value is set and the value differs from the default"""
        return self.value_and_changed()[1]

    def is_set(self) -> bool:
        """True, if the value is set (even equal to the default value)"""
        return self._value is not None

    def _written_value(self, always=False):
        """
        Parameters
        ----------
        always:
          Skip all condition checking

        Returns
        -------
        write value: Any
          The value to be written
        write: bool,
          Whether to write the value or not
        """
        d = self._definition
        if not d.is_stored:
            return None, False

        if not always:
            if not d.allowed(self._container):
                return None, False
        if not d.write_condition(self):
            return None, False

        if not d.type.has_value:
            return None, True

        if self.is_dangerous():
            return self._value, self._value() is not None

        value = self.result
        missing, _, np = d.type.missing_value()
        if np.__class__ is value.__class__ and np == value:
            return value, False
        if value is None or (not d.is_always_added and self.is_it_the_default_value(value)):
            return value, False
        return value, True

    def is_required(self) -> Union[bool, str]:
        """Return requiredness for this option."""
        r = self._definition.is_required
        if not r:
            return False
        if callable(r):
            return r(self)
        return r

    def _validate(self, why: str) -> None:
        """Validate this option for phase ``why``."""
        d = self._definition

        container = self._container
        if why == "parse" and d.validate_parsed:
            d.validate_parsed(self)
        if not d.allowed(container):
            return

        has_stored_value = (
            not d.is_generated
            and d.type.has_value
        )

        validate_value = (
            has_stored_value and
            why == 'save' and (
                d.is_validated
                if d.is_validated is not None
                else True )
        )

        if has_stored_value and (validate_value or d.is_repeated):
            value = self(unpack=False, all_values=True)
        else:
            value = None

        if validate_value:
            if d.is_repeated:
                if value is None:
                    validation_items = (value,)
                elif isinstance(value, DangerousValue):
                    validation_items = ()
                elif d.is_repeated.is_dict:
                    validation_items = tuple(value.values())
                else:
                    validation_items = value
            else:
                validation_items = (value,)

            for item in validation_items:
                if not isinstance(item, DangerousValue):
                    d.validate(self, item, why)

        else:
             value = None

        if d.is_repeated and value is not None:
            ValidationResult.emit(
                d.repeated_count.validate(self, len(value), why)
            )

        d.run_validators(self, container, why)

    @property
    def name(self):
        return self._definition.name

    def _as_dict(self, get):
        if not self._definition.allowed(self._container):
            return None
        return get(self)

    def value_and_changed(self):
        """Return value and whether the value was changed

        Returns
        -------
        value:mixed
          The value of the options (return all values for 'numbered array')

        changed:bool
          Whether the value is the same as the default value or not
        """
        d = self._definition
        if d.is_generated:
            return self(), False

        value = self._unpack_value(self._value)
        if value is not None:
            return value, not self.is_it_the_default_value(value)
        if d.is_repeated.is_numbered.has_default and self.default_value is not None:
            return {"def": self.default_value}, False
        else:
            return self.default_value, False

    def is_it_the_default_value(self, value):
        """Return, whether the given value is the default value. For
        numbered array, only the wildcard value can be set and this value
        have to be the same as the default."""
        d = self._definition
        if d.is_generated:
            return True

        default = self.default_value
        if d.is_repeated.is_numbered.has_default:
            return "def" in value and len(value) == 1 and d.type.is_the_same_value(value["def"], default)
        else:
            return d.type.is_the_same_value(value, default)

    def has_any_value(self):
        return self.result is not None

    def __len__(self):
        return len(self())

    def __iter__(self):
        return iter(self())

    def __bool__(self):
        return True

    def __repr__(self):
        if self._definition.is_generated:
            return f"<Generated value {self._get_path()}>"
        else:
            v = self._value
            o = None

        if o is None and v is None:
            v = self.default_value
            if callable(v):
                v = "fn()"
            if v is not None:
                o = " (default)"

        if v is None:
            if o:
                o = "out" + o
            else:
                o = "out"
            v = ""
        else:
            o = ""
            v = " = " + str(v)
            if len(v) > 20:
                v = f"{v[:10]}...{v[-10:]}"

        type = "(generated)" if self._definition.is_generated else f"of type {self._definition.type}"
        return f"<Option {self._get_path()} {type} with{o} value{v}>"


class CustomOption(Option):
    """An user-added option (configuration value). It can be removed from the section."""

    def remove(self):
        """Remove me from my "parent" section"""
        self._container.remove_member(self._definition.name)

    @classmethod
    def factory(cls, value_definition, type=mixed):
        """Returns factory function for the given value definition"""

        def create(name, section):
            definition = value_definition(name, type)
            definition.removable = True
            return cls(definition, section)

        create.grammar_type = type
        return create
