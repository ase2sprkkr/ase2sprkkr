from __future__ import annotations

from .configuration_containers import BaseConfigurationContainer
from .configuration import UnknownMemberPolicy
from typing import Any, Iterable, Mapping, Optional, Union
from .warnings import (
    DataValidityError,
    ReportInvalidPolicy,
    RetainInvalidPolicy,
    ValidationReason,
)
from .configuration_definitions import BaseDefinition
from .configuration_transaction import ConfigurationTransaction
import numpy as np


class RepeatedConfigurationContainer(BaseConfigurationContainer):
    """A container for configuration (problem-definition) options and/or sections.

    Options in the configuration (problem-definition) files are grouped to
    sections, sections are then grouped in a configuration file object.
    This is a base class for these containers.
    """

    class ContentsChange:
        """Original structure of a staged repeated section."""

        def __init__(
            self,
            section: "RepeatedConfigurationContainer",
            proposed: "RepeatedConfigurationContainer",
        ) -> None:
            """Capture ``section`` before installing ``proposed`` children."""
            self.section = section
            self.old_values = section._values
            self.old_parents = tuple(
                (child, child._container) for child in proposed.values()
            )

        def commit(self) -> None:
            pass

        def rollback(self) -> None:
            section = self.section
            section._values = self.old_values
            for child, parent in self.old_parents:
                child._container = parent
            for child in section.values():
                child._container = section

    def __init__(self, definition, container=None):
        """Create the container and its members, according to the definition"""
        super().__init__(definition, container)
        """
        The members of the container, in a form of ``{obj.name : obj}``
        """
        self._values = [] if self._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION else {}

    def __getitem__(self, name):
        """
        The members of the container are accesible using ``container["member name"]`` notation.
        """
        return self._values[name]

    def __len__(self):
        return len(self._values)

    def __bool__(self):
        return True

    def add(self, id=None):
        out = self._definition.create_object(self, repeated=False)
        if (id is None) != (self._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION):
            raise ValueError("Non-none Id for repeated list is allowedonly for dict-like repeated containers")
        if self._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION:
            self._values.append(out)
        else:
            self._values[id] = out
        return out

    def __contains__(self, name):
        """The check for existence of a member with the given name."""
        if self._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION:
            ln = len(self._values)
            return isinstance(name, int) and abs(name) < len(self) or name == -ln
        return name in self._values

    def stage_clear(
        self,
        transaction: ConfigurationTransaction,
        *,
        check_required: bool = True,
    ) -> bool:
        """Stage an empty collection in ``transaction``.

        ``check_required`` is accepted for the common staged-clear protocol.
        """
        empty = (
            []
            if self._definition.is_repeated
            == BaseDefinition.Repeated.LIST_SECTION
            else {}
        )
        return self.stage(transaction, empty)

    def get(self, name=None, unknown="find"):
        """
        Get the value, either of self or of a child of a given name.

        Parameters
        ----------
        name: None or str
          If None, return contained values as a dictionary.
          Otherwise, return the value of the member with the given name.

        unknown: str or None
          If unknown == 'find' and there is no member with a given name,
          try to find the first such in descendant conainers.

        Return
        ------
        value: mixed
        """

        if name is None:
            return self.as_dict()
        if "." in name:
            section, name = name.split(".", 1)
            return self._values[section].get(name)
        try:
            val = self._values[name]
        except (IndexError, KeyError, TypeError):
            raise KeyError(f"No {name} member of {self}") from None
        return val.get()

    def set(
        self,
        values: Optional[Union[Mapping[Any, Any], Iterable[Any], str]] = None,
        value: Any = None,
        *,
        unknown: UnknownMemberPolicy = "find",
        validation_reason: ValidationReason = "set",
        retain_invalid: Optional[RetainInvalidPolicy] = None,
        report_invalid: Optional[ReportInvalidPolicy] = None,
        merge: bool = False,
        **kwargs: Any,
    ) -> None:
        """Set repeated values and validate the complete proposal.

        By default, replace the complete collection. With ``merge=True``,
        append values to a list-like collection or update a dictionary-like
        collection by key.

        Invalid-value retention and reporting are controlled independently by
        ``retain_invalid`` (``"none"``, ``"typed"`` or ``"all"``) and
        ``report_invalid`` (``"raise"``, ``"warn"`` or ``"ignore"``).
        ``values`` supplies the children, or one child name paired with
        ``value``. ``unknown`` controls child assignment lookup and
        ``validation_reason`` selects ``"set"`` or ``"parse"`` validation.
        """
        with self._mutation(
            validation_reason, retain_invalid, report_invalid
        ) as (transaction, _policy):
            self.stage(
                transaction,
                values,
                value,
                unknown=unknown,
                merge=merge,
                **kwargs,
            )

    def stage(
        self,
        transaction: ConfigurationTransaction,
        values: Optional[Union[Mapping[Any, Any], Iterable[Any], str]] = None,
        value: Any = None,
        *,
        unknown: UnknownMemberPolicy = "find",
        merge: bool = False,
        **kwargs: Any,
    ) -> bool:
        """Stage repeated ``values`` in ``transaction``.

        ``unknown`` controls descendant assignment. ``merge``
        updates the current collection instead of replacing it; ``kwargs`` are
        forwarded to descendant staging.
        """
        proposed = self._definition.create_object(self._container)
        if merge:
            proposed._values = self._values.copy()
        if isinstance(values, str):
            if self._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION:
                values = [value]
            else:
                values = {values: value}
        elif value is not None:
            raise ValueError(
                "If value argument of Container.set method is given, the values have to be string name of the value"
            )
        elif values is None:
            values = ()

        try:
            items = values.items()
        except AttributeError:
            items = enumerate(values)

        for key, child_values in items:
            if (
                merge
                and proposed._definition.is_repeated
                == BaseDefinition.Repeated.DICT_SECTION
                and key in proposed._values
            ):
                child = proposed._values[key]
            else:
                child = proposed.add(
                    None
                    if proposed._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION
                    else key
                )
            child.stage(
                transaction,
                child_values,
                unknown=unknown,
                **kwargs,
            )

        change = self.ContentsChange(self, proposed)
        self._values = proposed._values
        for child in self.values():
            child._container = self
        transaction.push(change)
        return True

    def __iter__(self):
        """Iterate over all members of the container"""
        if self._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION:
            yield from range(len(self))
        else:
            yield from self._values.keys()

    def items(self):
        if self._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION:
            return enumerate(self._values)
        return self._values.items()

    def values(self):
        if self._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION:
            return self._values
        return self._values.values()

    def _as_dict(self, only_changed: Union[bool, str] = "basic", generated: bool = False, copy=False):
        """
        Return the content of the container as a dictionary.
        Nested containers will be transformed to dictionaries as well.

        Parameters
        ----------
        only_changed
          Return only changed values, or all of them?
          If True, return only the values, that differ from the defaults.
          If False, return all the values.
          The default value 'basic' means, return all non-expert values
          and all changed expert values.

        generated: bool
          Add generated values
        """
        if self._definition.is_repeated == BaseDefinition.Repeated.LIST_SECTION:
            out = [v.as_dict(only_changed, generated, copy) for v in self]
            if out and out[-1] is not None:
                for i in range(len(out) - 1, 0, -1):
                    if out[i] is not None:
                        out = out[: i + 1]
                    else:
                        out = None
        else:
            out = {}
            for k, v in self.items():
                value = v.as_dict(only_changed, generated, copy)
                if value is not None:
                    out[k] = value
        return out or None

    def is_changed(self):
        for i in self.values():
            if i.is_changed():
                return True
        return False

    def __repr__(self):
        return super().__repr__() + "[]"

    def _save_to_file(self, file, always=False, name_in_grammar=None, delimiter="") -> bool:
        """Save the content of the container to the file (according to the definition)

        Parameters
        ----------
        file: file
          File object (open for writing), where the data should be written

        always:
          Do not consider conditions

        Returns
        -------
        something_have_been_written
          If any value have been written return True, otherwise return False.
        """
        out = False
        for i in self.values():
            d = self._definition
            if d._save_to_file(file, i, always, name_in_grammar, delimiter=delimiter):
                name_in_grammar = False
                delimiter = d.repeated_delimiter
                out = True
        return out

    def _validate(self, why: str):
        """Validate the configuration data. Raise an exception, if the validation fail.

        Parameters
        ----------
        why
          Type of the validation. Possible values
          ``save`` - Full validation, during save.
          ``set`` - Validation on user input. Allow required values not to be set.
          ``parse`` - Validation during parsing - some check, that are enforced by the parser, can be skipped.
        """
        if why == "save" and not self._definition.is_optional and not self.has_any_value():
            DataValidityError.warn(f"Non-optional section {self._definition.name} has no value to save")

        values = self.values()
        self._definition.repeated_count.validate(self, len(values), why)
        for item in values:
            item._validate(why)


    def values_of(self, name):
        ln = len(self)
        if not ln:
            return np.empty((0,))
        return np.fromiter(
            (self[i][name]() for i in self), count=ln, dtype=self[0][name]._definition.type.numpy_dtype()[0]
        )
