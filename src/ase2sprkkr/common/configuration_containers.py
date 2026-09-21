"""In this file the common containers of configuration values are,
either for task or potential.

Configuration containers are classes, that holds configuration values
and other containers, and are able to write them to a configuration file,
and that are results of parsing of a configuration file.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from .configuration import Configuration, UnknownMemberPolicy
from .file_utils import filename_from_file
import copy
import itertools
from typing import Callable, Iterator, Mapping, Optional, Union, Any
from .configuration_transaction import ConfigurationTransaction
from .warnings import (
    DataValidityError,
    ReportInvalidPolicy,
    RetainInvalidPolicy,
    ValidationReason,
)


class DisabledAttributeError(AttributeError):
    """This exception is raised, if the attribute of a container exists,
    but it is disabled. E.g. because it has no sense for the current data.
    """


class BaseConfigurationContainer(Configuration, ABC):
    """Configuration container, that holds members, either in classical way
    (see :class:ConfigurationContainer) or treat them in a special way
    """

    def copy(self, copy_values: bool = False):
        """Create a copy of the container

        Parameters
        ----------
        copy_values
          If true, the copy of values is done, so their modifications do not affects the container.
          (e.g. for numpy arrays)
        """
        d = self._definition
        values = self.as_dict(copy=copy_values, only_changed=True)
        out = d.result_class(definition=d)
        out.set(
            values,
            unknown="add",
            retain_invalid="all",
            report_invalid="ignore",
        )
        return out

    @abstractmethod
    def values(self):
        """Return the configuration objects contained by this container."""

    def has_any_value(self) -> bool:
        """
        Return True if any member of the section has value.

        Return
        ------
          has_any_value: bool
              True, if no value in the container is set, False otherwise
        """
        for i in self.values():
            if i.has_any_value():
                return True
        return False

    def get_members(
        self,
        name: Any = None,
        unknown: UnknownMemberPolicy = "find",
        is_option: Optional[bool] = None,
        lower_case: bool = True,
        *,
        accept: Optional[Callable[[Configuration], bool]] = None,
    ) -> Iterator[Configuration]:
        """Return members selected by a path.

        ``None`` and leading parent operators are common to all containers;
        each concrete container implements lookup of its own members in
        :meth:`_get_members`.
        """
        if name is None:
            yield self
            return
        if isinstance(name, str) and name.startswith(".."):
            parent = self._container
            # A concrete item of a repeated section is held in a technical
            # collection carrying the same definition.  It must not add a
            # level to explicit paths, because no such level exists in the
            # definition tree.
            if (
                parent is not None
                and parent._definition is self._definition
            ):
                parent = parent._container
            if parent is not None:
                yield from parent.get_members(
                    name[2:],
                    unknown,
                    is_option,
                    lower_case,
                    accept=accept,
                )
            return
        yield from self._get_members(
            name,
            unknown,
            is_option,
            lower_case,
            accept=accept,
        )

    @abstractmethod
    def _get_members(
        self,
        name: Any,
        unknown: UnknownMemberPolicy,
        is_option: Optional[bool],
        lower_case: bool,
        *,
        accept: Optional[Callable[[Configuration], bool]],
    ) -> Iterator[Configuration]:
        """Return members owned directly or indirectly by this container."""
        raise NotImplementedError

    def get_member(
        self,
        name: Any,
        *,
        unknown: UnknownMemberPolicy = "find",
        is_option: Optional[bool] = None,
        lower_case: bool = True,
        accept: Optional[Callable[[Configuration], bool]] = None,
    ) -> Configuration:
        """Return the first matching member, raising ``KeyError`` if absent.

        ``unknown``, ``is_option``, ``lower_case`` and ``accept`` have the same
        meaning as in :meth:`get_members`.
        """
        members = self.get_members(
            name,
            unknown,
            is_option,
            lower_case,
            accept=accept,
        )
        try:
            return next(members)
        except StopIteration:
            raise KeyError(f"No member with name {name} in {self}") from None

    @property
    def definition(self):
        """The definition of the section.

        Returns
        -------
        ase2sprkkr.common.configuration_definitions.ContainerDefinition
        The definition of the section. I.e. the object that defines, which configuration values
        are in the section, their default values etc.
        """
        return self._definition


class ConfigurationContainer(BaseConfigurationContainer):
    """A container for configuration (problem-definition) options and/or sections.

    Options in the configuration (problem-definition) files are grouped to
    sections, sections are then grouped in a configuration file object.
    This is a base class for these containers.

    """

    class MemberChange:
        """Attachment of a custom member staged on a container."""

        def __init__(
            self, container: "ConfigurationContainer", member: Configuration
        ) -> None:
            """Record the custom ``member`` attached to ``container``."""
            self.container = container
            self.member = member

        def commit(self) -> None:
            pass

        def rollback(self) -> None:
            self.container.remove_member(self.member.name)

    class ParsedValuesChange:
        """Temporary record of the values explicitly present in parsed input."""

        _MISSING = object()

        def __init__(
            self, container: "ConfigurationContainer", values: Mapping[str, Any]
        ) -> None:
            """Expose parsed ``values`` temporarily on ``container``."""
            self.container = container
            self.old_values = getattr(
                container, "_parsed_values", self._MISSING
            )
            container._parsed_values = values

        def commit(self) -> None:
            self._restore()

        def rollback(self) -> None:
            self._restore()

        def _restore(self) -> None:
            if self.old_values is self._MISSING:
                self.container.__dict__.pop("_parsed_values", None)
            else:
                self.container._parsed_values = self.old_values

    def __init__(self, definition, container=None):
        """Create the container and its members, according to the definition"""
        super().__init__(definition, container)
        """
      The members of the container, in a form of ``{obj.name : obj}``
      """
        self._init_members_from_the_definition()

    def _init_members_from_the_definition(self):
        self._members = {}
        self._lowercase_members = {}
        """
      Non-hidden members of the containers, accesible via sanitized names.
      I.e. via names with whitespaces and other special characters replaced by underscore.
      These sanitized names are then used as names for "attributes" of this container, to
      make the members accesible via ``<container>.<member>`` notation.
      """
        self._interactive_members = {}
        for v in self._definition.members():
            if v.create_object:
                self._add(v.create_object(self))

    def items(self):
        """Members of the container. I.e. the options of the section, or sections
        of the configuration file e.t.c.

        Returns
        -------
        members: dict
        A dictionary of the shape ``{ name : member }``

        """
        return self._members

    def _get_attribute_member(self, name: str) -> Configuration:
        """
        Return a directly accessible member for attribute access.

        Search both the original and sanitized interactive names, and reject
        members which are hidden or disabled by their condition.
        """
        members = self.__dict__.get("_members")
        if members is None:
            raise AttributeError(name)
        if name in members:
            out = members[name]
        elif name in self._interactive_members:
            out = self._interactive_members[name]
        else:
            raise AttributeError(f"No {name} member of {self._definition}")
        d = out._definition
        if d.is_hidden:
            raise DisabledAttributeError(
                f"member {name} of {self} is not directly accessible. "
                "Probably it"
                "s a hidden attribute used for some kind of logic, "
                "for which a direct access has no sense. If you really need "
                'an access to the attribute, you can use the "container['
                "name"
                ']" notation.'
            )
        allowed = d.allowed(self)
        if not allowed:
            if allowed is False:
                raise DisabledAttributeError(
                    f"member {name} of {self} is not accessible for "
                    "the current data. It is probably not available or has no sense "
                    "in this particular case (e.g. data file does not contain needed "
                    "data for it). If you eally need an access to the attribute, you can use the "
                    '"container['
                    "name"
                    ']" notation.'
                )
            else:
                raise DisabledAttributeError(str(allowed))
        return out

    def __getattr__(self, name):
        """
        The members of the container are accesible as attributes of the container, too.
        Either using their normal, or ``sanitized`` names.
        """
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            out = self._get_attribute_member(name)
        except AttributeError as e:
            if isinstance(e, DisabledAttributeError):
                msg = str(e)
                cls = DisabledAttributeError
            else:
                msg = f"There is no value with name {name} in {self}.\nMaybe, you want to add a custom value using the add method?"
                cls = AttributeError
            raise cls(msg) from e
        return out

    def __getitem__(self, name):
        """
        The members of the container are accesible using ``container["member name"]`` notation.
        """
        if isinstance(name, tuple):
            if not name:
                raise KeyError("An empty tuple not allowed as a key.")
            out = self._members[name[0]]
            ll = len(name)
            if ll == 1:
                return out
            if ll == 2:
                return out[name[1]]
            return out[name[1:]]

        return self._members[name]

    def _get(self, name, default=None):
        return self._members.get(name, default)

    def __dir__(self):
        """
        Expose the interactive_members in the container attribute listing.
        Interactive_members are the non-hidden members identified by their sanitized names.
        """
        # def ok(member):
        #    d = member._defintion
        #    return not d.condtion or d.condition(self)

        members = (i for i in self._interactive_members.keys())
        if self._definition.dir_common_attributes:
            members = itertools.chain(members, super().__dir__())
        return members

    def __contains__(self, name):
        """The check for existence of a member with the given name."""
        if isinstance(name, tuple):
            ll = len(name)
            if ll == 0:
                return False
            n = name[0]
            if n not in self._members:
                return False
            if ll == 1:
                return True
            member = self._members[n]
            if ll == 2:
                return name[1] in member
            else:
                return name[1:] in member
        return name in self._members

    def stage_clear(
        self,
        transaction: ConfigurationTransaction,
        *,
        check_required: bool = True,
    ) -> bool:
        """Stage descendant clearing in ``transaction`` without validation.

        ``check_required`` is forwarded to descendant options.
        """
        changed = False
        for member in self._members.values():
            if member._definition.is_generated:
                continue
            changed = member.stage_clear(
                transaction,
                check_required=check_required,
            ) or changed
        return changed

    def _get_members(
        self,
        name: str,
        unknown: UnknownMemberPolicy,
        is_option: Optional[bool],
        lower_case: bool,
        *,
        accept: Optional[Callable[[Configuration], bool]],
    ) -> Iterator[Configuration]:
        """
        Get all the members of given name. According to ``unknown`` parameter,
        either only from self, or from any child containers, too.

        Parameters
        ----------
        name: str
          Name or dotted path of the requested member.

        unknown: str or None
          If unknown == 'find' and there is no member with a given name,
          try to find the first such-named item (case insensitive)
          in the descendant conainers.

        is_option: bool
          If set, limit to either Option or non-option items

        lower_case: bool
          If true, try to search for the lower-cased name

        accept: callable
          If given, return only members for which this predicate is true.
          A rejected direct member does not prevent searching descendants.

        Return
        ------
        value: mixed
        """
        if "." in name:
            section_name, child_name = name.split(".", 1)
            members = self.get_members(
                section_name,
                unknown,
                lower_case=lower_case,
            )
            for member in members:
                if member._definition.is_option:
                    continue
                yield from member.get_members(
                    child_name,
                    unknown,
                    is_option,
                    lower_case,
                    accept=accept,
                )
            return

        member = self._members.get(name)
        if member is None and lower_case:
            member = self._lowercase_members.get(name.lower())

        if (
            member is not None
            and (is_option is None or member._definition.is_option is is_option)
            and (accept is None or accept(member))
        ):
            yield member
            return
        elif unknown == "find":
            search_name = name.lower() if lower_case else name
            for child in self:
                for found_member in child._find_members(
                    search_name, is_option, lower_case
                ):
                    if accept is None or accept(found_member):
                        yield found_member

    def get(self, name=None, unknown="find", is_option=True):
        """
        Get the value, either of self or of a child of a given name.

        Parameters
        ----------
        name: None or str
          If None, return contained values as a dictionary.
          Otherwise, return the value of the member with the given name.

        unknown: str or None
          If unknown == 'find' and there is no member with a given name,
          try to find the first such-named item (case insensitive)
          in the descendant conainers.
          unknown == 'find_exact' do the same, case sensitive.

        Return
        ------
        value: mixed
        """
        if name is None:
            return self.as_dict(only_changed=False, generated=True)
        item = self.get_member(name, unknown=unknown, is_option=is_option)
        return item.as_dict(only_changed=False, generated=True)

    def set(
        self,
        values: Optional[Union[Mapping[str, Any], str]] = None,
        value: Any = None,
        *,
        unknown: UnknownMemberPolicy = "find",
        validation_reason: ValidationReason = "set",
        retain_invalid: Optional[RetainInvalidPolicy] = None,
        report_invalid: Optional[ReportInvalidPolicy] = None,
        **kwargs: Any,
    ) -> None:
        """Set values using independent retention and reporting policies.

        ``retain_invalid`` accepts ``"none"``, ``"typed"`` or ``"all"``;
        ``report_invalid`` accepts ``"raise"``, ``"warn"`` or ``"ignore"``.
        Parsing defaults to retaining typed invalid values and warning. Other
        assignments default to discarding invalid values and raising.

        ``values`` is a mapping or one member name paired with ``value``;
        ``unknown`` selects lookup, addition, ignoring or failure behavior.
        ``validation_reason`` is ``"set"`` or ``"parse"``. Additional keyword
        arguments are treated as member assignments.
        """
        with self._mutation(
            validation_reason, retain_invalid, report_invalid
        ) as (transaction, _policy):
            self.stage(
                transaction,
                values,
                value,
                unknown=unknown,
                **kwargs,
            )

    def stage(
        self,
        transaction: ConfigurationTransaction,
        values: Optional[Union[Mapping[str, Any], str]] = None,
        value: Any = None,
        *,
        unknown: UnknownMemberPolicy = "find",
        **kwargs: Any,
    ) -> bool:
        """
        Stage values in a transaction. This is the sole assignment-routing
        implementation for configuration containers.

        Usage:

        > input_parameters.set({'NITER': 5, 'NE': [10]})
        or
        > input_parameters.set(NITER=5, NE=[10])

        Parameters
        ----------

        values:
          Dictionary of values to be set, or the name of the value, if the value is given.

        value:
          Value to be set. Setting this argument require to pass string name to the values argument.

        unknown: 'add', 'find' or None
          How to handle unknown (not known by the definition) parameters.
          If 'find', try to find the values in descendant containers.
          If 'add', add unknown values as custom values.
          If None, throw an exception.
          Keyword only argument.

        transaction: ConfigurationTransaction
          Transaction receiving the staged changes.

        **kwargs: dict
          The values to be set (an alternative syntax as syntactical sugar)
        """
        if isinstance(values, str):
            assignments = {values: value}
        elif value is not None:
            raise ValueError(
                "If value argument of Container.set method is given, the values have to be string name of the value"
            )
        elif values is None:
            assignments = {}
        else:
            try:
                values.items()
            except AttributeError:
                raise ValueError(
                    "Only a mapping can be assigned to a section."
                ) from None
            assignments = values

        if kwargs:
            try:
                assignments = copy.copy(assignments)
                assignments.update(kwargs)
            except (AttributeError, TypeError):
                assignments = dict(assignments)
                assignments.update(kwargs)

        policy = transaction.policy
        if policy.why == "parse":
            transaction.push(self.ParsedValuesChange(self, assignments))

        changed = False
        for name, assignment in assignments.items():
            changed = self.stage_value(
                transaction,
                name,
                assignment,
                unknown=unknown,
            ) or changed
        return changed

    def stage_value(
        self,
        transaction: ConfigurationTransaction,
        name: str,
        value: Any,
        *,
        unknown: UnknownMemberPolicy = "find",
    ) -> bool:
        """Resolve ``name`` and stage ``value`` using the supplied policies."""

        def accepts_value(member: Configuration) -> bool:
            """Return whether ``member`` accepts the proposed value."""
            return member._definition.accept_value(value)

        try:
            member = self.get_member(
                name,
                unknown=unknown,
                accept=accepts_value,
            )
        except KeyError:
            if unknown == "ignore":
                return False
            if unknown != "add":
                raise
            if "." in name:
                raise KeyError(
                    f"Cannot add dotted path {name} to {self}; "
                    "set the value through its containing section"
                ) from None
        else:
            return member.stage(
                transaction, value, unknown=unknown
            )

        with transaction.savepoint() as savepoint:
            change = transaction.push(self.stage_member(name))
            changed = change.member.stage(
                transaction, value, unknown="add"
            )
            if not changed:
                savepoint.rollback()
        return changed

    def add(self, name: str, value=None):
        """
        Add custom value to the container

        Parameters
        ----------
        name: str
          Name of the added value

        value: value
          Value of the added value
        """
        with self._mutation("set") as (transaction, _policy):
            change = transaction.push(self.stage_member(name))
            if value is not None:
                change.member.stage(
                    transaction, value, unknown="add"
                )

    def stage_member(self, name: str) -> "ConfigurationContainer.MemberChange":
        """Attach custom member ``name`` and return its rollback change."""
        custom_class = getattr(self._definition, "custom_class", None)
        if not custom_class:
            raise TypeError(
                f"Can not add custom members to a configuration class "
                f"{self._definition}"
            )
        if name in self._members:
            raise TypeError(
                f"Section member {name} is already in the section "
                f"{self._definition}"
            )
        member = custom_class(name, self)
        change = self.MemberChange(self, member)
        self._add(member)
        return change

    def remove_member(self, name: str):
        """
        Remove a (previously added) custom value from the container
        """
        cclass = getattr(self._definition, "custom_class", False)
        if not cclass:
            raise TypeError("Can not remove items of {}".format(name))
        if not getattr(self._members[name], "remove"):
            raise KeyError("No custom member with name {} to remove".format(name))
        member = self._members[name]
        del self._members[name]

        iname = name.lower()
        if self._lowercase_members.get(iname) is member:
            del self._lowercase_members[iname]

        df = member._definition
        if not df.is_hidden:
            for iname in df.interactive_names():
                if self._interactive_members.get(iname) is member:
                    del self._interactive_members[iname]
                    iname = iname.lower()
                    if self._lowercase_members.get(iname) is member:
                        del self._lowercase_members[iname]

    def __iter__(self):
        """Iterate over all members of the container"""
        yield from self._members.values()

    def values(self):
        """Return the contained configuration values."""
        return self._members.values()

    def _values(self):
        """Iterate over all members of the container"""
        yield from self.values()

    def _as_dict(self, get):
        """
        Return the content of the container as a dictionary.
        Nested containers will be transformed to dictionaries as well.

        Parameters
        ----------
        get
          This function will be applied to the options to (possible) obtain
          the values
        """
        out = {}
        for i in self._values():
            value = i._as_dict(get)
            if value is not None:
                out[i._definition.real_name] = value
        return out or None

    def _find_members(self, name: str, is_option=None, lower: bool = False):
        """
        Iterates over a value of a given name in self or in any
        of owned subcontainers.

        Parameters
        ----------
        name: str
        A name of the sought options

        is_option:bool
        If True, find only options, if False, find only the others.

        lower:bool
        If True, find an option with given lowercased name (case insensitive)

        Returns
        -------
        value:typing.Optional[ase2sprkkr.common.options.Option]
        The first option of the given name, if such exists. ``None`` otherwise.
        """
        if is_option is not True and name == (self.name.lower() if lower else self.name):
            yield self
        for i in self._values():
            if i._definition.is_hidden:
                continue
            yield from i._find_members(name, is_option, lower)

    def _add(self, member):
        name = member.name
        self._members[name] = member

        lname = name.lower()
        if lname not in self._lowercase_members:
            self._lowercase_members[lname] = member

        df = member._definition

        if not df.is_hidden:
            for iname in df.interactive_names():
                if iname not in self._interactive_members:
                    self._interactive_members[iname] = member
                    iname = iname.lower()
                    if iname not in self._lowercase_members:
                        self._lowercase_members[iname] = member

    def is_changed(self):
        for i in self._values():
            if i.is_changed():
                return True
        return False

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
        return self._definition._save_to_file(file, self, always, name_in_grammar, delimiter)

    def __setattr__(self, name, value):
        """Setting the (unknown) attribute of a section sets the value of the member
        with a given name"""
        if name[0] == "_" or name in self.__dict__ or hasattr(getattr(self.__class__, name, None), "__set__"):
            super().__setattr__(name, value)
        else:
            val = self._get_attribute_member(name)
            val.set(value)

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
        for item in self._values():
            if (
                not item._definition.is_option
                and not item._definition.allowed(self)
            ):
                continue
            item._validate(why)
        self._definition.run_validators(self, self, why)

class BaseSection(ConfigurationContainer):
    """A section of SPRKKR configuration - i.e. part of the configuration file."""


class Section(BaseSection):
    """A standard section of a task or potential (whose content is predefinded by SectionDefinition)"""

    @property
    def definition(self):
        """The definition of the section.

        Returns
        -------
        ase2sprkkr.common.configuration_definitions.ContainerDefinition
        The definition of the section. I.e. the object that defines, which configuration values
        are in the section, their default values etc.
        """
        return self._definition


class CustomSection(BaseSection):
    """Custom task section. Section created by user with no definition"""

    def remove(self):
        """Remove the custom section from the parent container"""
        self._container.remove(self.name)

    @classmethod
    def factory(cls, definition_type):
        """Create a factory for custom values.

        Parameters
        ----------
        definition_type: ase2sprkkr.common.configuration_definitions.BaseDefinition
          Type (definitions) of the custom values created by
          the resulting function

        Return
        ------
        factory: callable
          Factory function of the signature (name: str, container: ase2sprkkr.common.configuration_containers.ConfigurationContainer)
          that created a custom value or section of the given definition

        """

        def create(name, container):
            definition = definition_type(name)
            definition.removable = True
            return cls(definition, container)

        return create


class RootConfigurationContainer(ConfigurationContainer):
    """Base class for data of configuration/problem-definition files

    In addition to container capabilities, it can read its data from/to file.
    """

    name_in_grammar = False

    _filename = None

    def read_from_file(self, file, clear_first: bool = True, allow_dangerous: bool = False):
        """Read data from a file

        Parameters
        ----------
        file: str or file
          File to read the data from

        clear_first
          Clear the container first.
          Otherwise, the data in the sections that are not present in the
          file are preserved.
        allow_dangerous
          Allow to load dangerous_values, i.e. the values that do not pass the requirements for the input values (e.g. of a different type or constraint-violating)
        """
        values = self._definition.parse_file(file, allow_dangerous=allow_dangerous)
        with self._mutation("parse") as (transaction, _policy):
            if clear_first:
                self.stage_clear(transaction, check_required=False)
            self.stage(
                transaction,
                values,
                unknown="add",
            )
        self._filename = filename_from_file(file, None)

    @property
    def original_filename(self):
        return self._filename

    def find(self, name, unknown="find", is_option=True, lower_case=True, first=True):
        """Find a configuration value of a given name in the owned sections"""
        if first:
            return self.get_member(
                name,
                unknown=unknown,
                is_option=is_option,
                lower_case=lower_case,
            )
        return self.get_members(name, unknown, is_option, lower_case)
