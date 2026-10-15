from .configuration_definitions import RealItemDefinition, BaseDefinition, Validator
from .grammar import as_delimiter, delimitedList
from .misc import dict_first_item
from .repeated_configuration_containers import RepeatedConfigurationContainer
from .configuration_containers import Section
from .repetition import Repeated
from .decorators import add_to_signature, cache
from .parsing_results import dict_from_parsed


import pyparsing as pp
from typing import Iterable, Optional, Union
import re
import inspect
from io import StringIO


class ContainerDefinition(RealItemDefinition):
    """Base class for a definition (of contained data, format, etc)
    of either a whole configuration file
    (e.g. :class:`InputParameters<ase2sprkkr.input_parameters.input_parameters.InputParameters>` or
    e.g. :class:`Potential<ase2sprkkr.potentials.potentials.Potential>`) or
    its :class:`Section<ase2sprkkr.common.configuration_containers.Section>`.
    """

    force_order = False
    """ Force order of its members """

    delimiter = pp.Empty().set_name("")

    is_option = False
    """Container definitions create sections or configuration roots, not options."""

    value_name_format = None
    """ The (print) format, how the name is written """

    write_last_delimiter = True

    item_type = "section"

    @staticmethod
    def _dict_from_named_values(args, items=None):
        """auxiliary method that creates dictionary from the arguments"""
        items = items or {}
        for value in args:
            items[value.name] = value
        return items

    dir_common_attributes = True
    """ In dir listing, show the common 'object' attributes """

    def __init__(
        self,
        name,
        members=[],
        written_name=None,
        info=None,
        description=None,
        is_optional=False,
        is_hidden=False,
        is_expert=False,
        has_hidden_members=False,
        name_in_grammar=None,
        force_order=None,
        name_regex=False,
        result_class=None,
        is_repeated=False,
        repeated_delimiter=None,
        repeated_count=None,
        repeated_with_name=False,
        write_condition=None,
        warning_condition=None,
        validators: Optional[Union[Validator, Iterable[Validator]]] = (),
    ):
        """
        Definition of container (e.g. section of an input file).
        For the rest of the parameters see the :class:`RealItemDefinition`.

        Parameters
        ----------
        has_hidden_members: bool
          If true, this section is not intended for a direct editing

        force_order: bool
          If True, the items has to retain the order, if False, the items can be in the input file in any order.

        validators: callable or iterable of callables
          Semantic validators passed to :class:`RealItemDefinition`.
        """

        super().__init__(
            name=name,
            written_name=written_name,
            is_optional=is_optional,
            is_hidden=is_hidden,
            is_expert=is_expert,
            name_in_grammar=name_in_grammar,
            info=info,
            description=description,
            name_regex=name_regex,
            result_class=result_class,
            is_repeated=is_repeated,
            repeated_delimiter=repeated_delimiter,
            repeated_count=repeated_count,
            repeated_with_name=repeated_with_name,
            write_condition=write_condition,
            warning_condition=warning_condition,
            validators=validators,
        )

        if not isinstance(members, dict):
            members = self._dict_from_named_values(members)

        if self.value_name_format:
            for i in members.values():
                i.value_name_format = self.value_name_format
        self._members = members

        self.has_hidden_members = has_hidden_members
        if force_order is not None:
            self.force_order = force_order
        if isinstance(self.delimiter, str):
            self.delimiter = as_delimiter(self.delimiter)

    default_repeated_type = Repeated.LIST_SECTION

    configuration_type_name = "SECTION"
    """ Name of the container type in the runtime documentation """

    def added_to_container(self, container):
        """Attach this container and propagate the callback to its children."""
        super().added_to_container(container)
        for member in self._members.values():
            member.added_to_container(self)

    def allow_duplication(self):
        return self.is_repeated

    def __repr__(self):
        return f"<{self.configuration_type_name} {self.name}>"

    def data_description(self, verbose: Union[bool, str, int] = False, show_hidden: bool = False, prefix: str = ""):
        """
        Return the runtime documentation for the configuration described by this object.

        Parameters
        ----------
        verbose
          If ``False``, only one line with the section name and basic info is returned.
          If ``True``, the items contained in the section are listed.
          If ``'all'``, add also detailed info about all descendants.
          If an ``int`` is given, print detailed informations about n levels. I.e. ``1`` is the same as ``True``

        show_hidden
          If False, do not show hidden members.

        prefix
          The string, with with each line will begin (commonly the spaces for the indentation).
        """

        def container_name():
            out = self.configuration_type_name
            if out:
                out += " "
            out += self.name
            return out

        out = f"{prefix}{container_name()}"
        flags = []
        if self.force_order:
            flags.append("fixed-order")
        if self.is_hidden:
            flags.append("hidden")
        if self.is_optional:
            flags.append("optional")
        if self.is_expert:
            flags.append("expert")
        if self.is_repeated:
            flags.append("repeated")
        if flags:
            flags = ", ".join(flags)
            out += f" ({flags})"

        if verbose:
            if isinstance(verbose, int):
                verbose -= 1
            else:
                verbose = verbose if verbose == "all" else False

            add = self.additional_data_description(verbose, show_hidden, prefix)
            if self.info_in_data_description:
                info = self.info(False)
                if info:
                    add = prefix + info + "\n\n" + add + "\n"
            if add:
                out += " contains:"
                under = prefix + "-" * len(out) + "\n"
                out = f"{prefix}{out}\n{under}{add}"
        return out

    def additional_data_description(self, verbose: Union[bool, str, int] = False, show_hidden=False, prefix: str = ""):
        """
        Return the description (documentation for runtime) of the items in the container.

        Parameters
        ----------
        verbose
          If ``True``, include detailed description of the children.
          If ``'all'``, include even detailed description.
          If ``int`` is given, print detailed informations up to n levels.

        show_hidden
          If False, do not show hidden members.

        prefix
          The string, with with each line will begin (commonly the spaces for the indentation).
        """
        cprefix = prefix + self._description_indentation
        out = []

        def write(i):
            s = i.data_description(verbose, show_hidden, cprefix)
            if not i.info_in_data_description:
                if "\n" not in s:
                    info = i.info(False)
                    if info:
                        s = s + (" " * (max(40 - len(s), 0) + 2)) + info
                else:
                    ccprefix = cprefix + i._description_indentation
                    s += "\n\n"
                    s += ccprefix + i.info(False).replace("\n", "\n" + ccprefix)
                    s += "\n"
            out.append(s)

        expert = False
        for i in self:
            if i.is_hidden and not show_hidden:
                continue
            if not i.is_expert:
                write(i)
            else:
                expert = True

        if expert:
            out.append(f"{cprefix}\n{cprefix}Expert options:")
            out.append(cprefix + "--------------")
            cprefix += self._description_indentation

            for i in self:
                if i.is_expert:
                    if i.is_hidden and not show_hidden:
                        continue
                    write(i)

        return "\n".join(out)

    def __iter__(self):
        return iter(self._members.values())

    def members(self):
        return self._members.values()

    def names(self):
        return self._members.keys()

    def __getitem__(self, key):
        return self._members[key]

    def __setitem__(self, key, value):
        self._members[key] = value

    def __contains__(self, key):
        return key in self._members

    def get_member(self, name):
        """Return the definition at the given explicit dotted path.

        A leading ``..`` selects the parent container, so ``..N`` resolves
        ``N`` in the direct parent and ``....N`` resolves it two levels up.
        """
        if name.startswith(".."):
            if self.container is None:
                raise KeyError(f"No parent of {self}")
            return self.container.get_member(name[2:])
        if "." in name:
            section_name, child_name = name.split(".", 1)
            return self._members[section_name].get_member(child_name)
        return self._members[name]

    def remove(self, name):
        del self._members[name]
        return self

    def copy(self, add=[], remove=[], defaults={}, **kwargs):
        """Get the args to copy the container.
        The method is just repeated to add named argumets.

        Parameters
        ----------
        """
        arguments = self._get_init_args_for_copy(add, remove, defaults, **kwargs)
        return self.__class__(**arguments)

    def _get_init_args_for_copy(self, add=[], remove=[], defaults={}, **kwargs):
        """Get the args to copy the container."""
        members = dict(((k, i.copy()) for k, i in self._members.items()))
        for i in remove:
            del members[i]
        members.update(self._dict_from_named_values(add))
        for i, v in defaults.items():
            members[i].default_value = members[i].type.convert(v)

        out = super()._get_init_args_for_copy(**kwargs)
        out["members"] = members
        return out

    def copy_member(self, name, **kwargs) -> BaseDefinition:
        """Copy a member, allowing to redefine its properties.

        Returns
        -------
        new_member: BaseDefinition
          The newly created member
        """
        out = self._members[name].copy(**kwargs)
        self._members[name] = out
        return out

    def _grammar_of_values(self, allow_dangerous: bool = False, delimiter=None):
        if self.custom_class:
            custom_value = self.custom_member_grammar(self.excluded_names_condition())
        else:
            custom_value = None
        delimiter = delimiter or self.delimiter

        def repeated_grammars():
            """If the item can be repeated, do it here - we don't know, whether there is a fixed order in any way
            (e.g. the item is followed by the items without name in grammar)
            """
            for i in self._members.values():
                g = i._grammar and i._grammar(allow_dangerous)
                if not g:
                    continue
                if i.is_repeated and i.repeated_with_name:
                    g = i.repeated_count.grammar(g, i.repeated_delimiter if i.repeated_delimiter is not None else delimiter)
                yield i, g

        def grammars():
            """This function iterates over the items of the container, joining all the without name_in_grammar with the previous ones."""
            it = iter(repeated_grammars())
            head_item, grammar_chain = next(it)

            for item, grammar in it:
                if item.is_independent_on_the_predecessor:
                    yield head_item, grammar_chain
                    head_item, grammar_chain = item, grammar
                else:
                    add = delimiter + grammar
                    if item.condition and self.force_order:
                        add = item.condition.prepare_grammar(item, add)
                    if item.is_optional:
                        add = pp.Optional(add)
                    grammar_chain = grammar_chain + add
            yield head_item, grammar_chain

        if self.force_order:
            init = pp.Empty().set_name("<INIT>")

            def set_loc(loc, toks):
                init.location = loc

            init.set_parse_action(set_loc)

            first = pp.Empty().add_condition(lambda loc, toks: loc == init.location).set_name("<FIRST>")
            if custom_value:
                cvs = pp.ZeroOrMore(custom_value + delimiter).set_name("<custom...>")
                after = delimiter + cvs
            else:
                after = pp.Forward() << delimiter
            after.add_condition(lambda loc, toks: loc != init.location)
            inter_cvs = (first | after).set_name("<?DELIM>")
            inter = first | delimiter.copy().add_condition(lambda loc, toks: loc != init.location)

            def sequence():
                for head, g in repeated_grammars():
                    if head.is_independent_on_the_predecessor:
                        delim = inter_cvs
                    else:
                        delim = inter
                    g = delim + g
                    if head.condition:
                        g = head.condition.prepare_grammar(head, g)
                    if head.is_optional:
                        g = pp.Optional(g)
                    yield g

            values = pp.And([i for i in sequence()])

            if custom_value:
                if not self._first_section_has_to_be_first():
                    values = cvs + values
                values += pp.ZeroOrMore(delimiter + custom_value)
            values = init + values

        else:
            it = grammars()
            # store the first fixed "chain of sections"
            first = self._first_section_has_to_be_first() and next(it)[1]
            # the rest has any order
            values = pp.MatchFirst([i for head, i in it])
            if custom_value:
                values |= custom_value
            values = delimitedList(values, delimiter)
            if first:
                values = first + pp.Optional(delimiter + values)

        values.set_parse_action(lambda x: dict_from_parsed(x.asList()))
        if self.write_last_delimiter:
            values = values + pp.Optional(delimiter)
        if self.is_repeated and not self.repeated_with_name:
            values = self.repeated_count.grammar(values, self.repeated_delimiter)
            #values = pp.Group(values, aslist=True).set_name(f"<{self.name}[]>")
            values = pp.Group(values, aslist=True).set_name(f"<{self.name}[]>")
        return values

    def _allow_duplicates_of(self, name):
        """Can a given element (identified by name) have more values in the parsed results?
        (However, not all definitions have to specify allow_duplicates, just the ones
        that have a value). For the others, this function raises an error.
        """
        return self[name].allow_duplication()

    def _create_grammar(self, allow_dangerous=False):
        delimiter = self.delimiter
        values = self._grammar_of_values(allow_dangerous, delimiter)
        out = self._tuple_with_my_name(values, delimiter)
        out.set_name(self.name)
        return out

    @classmethod
    @cache
    def delimited_custom_value_grammar(cls):
        """Return the grammar for the custom child with delimiter.
        The delimiter can delimite it either from the previous child or from the section name."""

        return cls.child_class.delimiter + cls.custom_value_grammar()

    custom_name_characters = pp.alphanums + "_-()"
    """ Which characters can appears in an unknown child (value/section) name """

    @classmethod
    def custom_member_grammar(cls, name_condition=None):
        """Grammar for the custom - unknown - child"""
        name = pp.Word(cls.custom_name_characters).set_parse_action(lambda x: x[0].strip())
        if name_condition:
            name.add_condition(name_condition)
        out = (name + cls.delimited_custom_value_grammar()).set_parse_action(lambda x: tuple(x))
        out.set_name(cls.custom_value_name)
        return out

    def all_member_names(self):
        for i in self:
            yield from i.all_names_in_grammar()

    def excluded_names_condition(self):
        """Add the condition to the element, that
        its value is not any of given names"""
        names = set((_ending_numbers.sub("", i).upper() for i in self.all_member_names()))

        if not names:
            return

        def cond(x):
            striped = _ending_numbers.sub("", x[0]).upper()
            return striped not in names

        return cond

    def _first_section_has_to_be_first(self):
        """Has/ve the first child(s) in an unordered sequence fixed position?"""
        return not dict_first_item(self._members).is_independent_on_the_predecessor

    def parse_file(self, file, return_value_only=True, allow_dangerous=False):
        """Parse the file, return the parsed data as dictionary"""
        grammar = self.grammar(allow_dangerous)
        out = grammar.parse_file(file, parse_all=True)
        return self.parse_return(out, return_value_only)

    def parse(self, string, whole_string=True, return_value_only=True, allow_dangerous=False):
        """Parse the string, return the parsed data as dictionary"""
        grammar = self.grammar(allow_dangerous)
        out = grammar.parse_string(string, parse_all=True)
        return self.parse_return(out, return_value_only)

    def parse_return(self, val, return_value_only: bool = True):
        """Clean up the parsed values (unpack then from unnecessary containers)

        Parameters
        ----------
        return_value_only
          Return only value, not name - value tuple

        """
        val = val[0]
        if return_value_only:
            val = val[1]
        return val

    async def parse_from_stream(
        self, stream, up_to, start=None, whole_string=True, return_value_only=True, allow_dangerous=False
    ):
        """
        Parse string readed from asyncio stream.
        The stream is readed up to the given delimiter
        """

        result = await stream.readuntil(up_to)
        result = result[: -len(up_to)].decode("utf8")
        if start:
            result = start + result
        return self.parse(result, whole_string)

    def read_from_file(self, file, allow_dangerous=False, **kwargs):
        """Read a configuration file and return the parsed Configuration object"""
        out = self.result_class(definition=self, **kwargs)
        out.read_from_file(file, allow_dangerous=allow_dangerous)
        return out

    def read_from_dict(self, values, **kwargs):
        out = self.result_class(definition=self, **kwargs)
        out.set(values, unknown="add")
        return out

    def read_from_string(self, string, allow_dangerous=False, **kwargs):
        return self.read_from_file(StringIO(string), allow_dangerous, **kwargs)

    repeated_class = RepeatedConfigurationContainer
    """ Class for the repeated sections """

    def create_object(self, container=None, repeated: bool = True):
        """Create an instance (section)

        container: BaseConfigurationContainer
            To which container the created object will belong

        repeated:
            Has meaning only for a is_repeated section. Then, if it is True,
            a Container for repeated values of the section is returned.
            Otherwise, the container just for one instance of a section is
            returned.
        """
        self._finalize()
        if repeated and self.is_repeated:
            return self.repeated_class(self, container)
        return super().create_object(container)

    def _save_to_file(self, file, value, always=False, name_in_grammar=None, delimiter="") -> bool:
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
        if not always:
            if not self.write_condition(value._container) or not self.allowed(value._container):
                return

        if self.is_expert:
            if not value.is_changed():
                return False
        else:
            if not value.has_any_value():
                return False
        if name_in_grammar is None:
            name_in_grammar = self.name_in_grammar

        if delimiter:
            file.write(str(delimiter))
        if name_in_grammar:
            file.write(self.formated_name)
            file.write(str(self.delimiter))

        members = iter(value)
        if self.write_last_delimiter:
            for o in members:
                if o._save_to_file(file, always):
                    file.write(str(self.delimiter))
        else:
            delimiter = ""
            for o in members:
                if o._save_to_file(file, always, delimiter=delimiter):
                    delimiter = self.delimiter

        return True


_ending_numbers = re.compile("[0-9]*$")


class SectionDefinition(ContainerDefinition):
    """Base class for definition of the sections in Pot or InputParameters files.

    It just redefine a few properties/methods to values/behavior typical for the sections
    """

    result_class = Section

    @property
    def values(self):
        return self._members

    custom_value_name = "CUSTOM_VALUE"
    """ Just the name that appears in the grammar, when it is printed."""

    @classmethod
    @cache
    def delimited_custom_value_grammar(cls):
        gt = cls.custom_class.grammar_type
        # here the child (Value) class delimiter should be used
        out = cls.child_class.delimiter + gt.grammar()
        optional, df, _ = gt.missing_value()
        if optional:
            out = out | pp.Empty().set_parse_action(lambda x: df)
        return out

    def _generic_info(self):
        return f"Configuration section {self.name}"

    def accept_value(self, value) -> bool:
        if isinstance(value, dict):
            return True
        return self.is_repeated and isinstance(value, Iterable)


class ConfigurationRootDefinition(ContainerDefinition):
    """From this class, the definition of the format of a whole configuration file should be derived."""

    write_last_delimiter = False
    """ Do not print additional newline after the last section """

    name_in_grammar = False
    """ The configuration files has commonly no "name" of its content, they
   just contains their content.

   However, in some cases the name_in_grammar could be used for some kind of
   prefix in the file, however, it is better to use a fixed value for this purpose.
   """
    item_type = "configuration"

    @classmethod
    def definition_from_dict(cls, name, defs=None):
        """
        Create instance of the definition from a dictionary, creating
        the sections (and values) definitions recursively.
        """

        def gen(i):
            section = defs[i]
            if not isinstance(defs, SectionDefinition):
                section = cls.child_class(i, section)
            return section

        if defs is None:
            defs = name
            name = cls.__name__

        return cls((gen(i) for i in defs))

    @add_to_signature(ContainerDefinition.__init__, prepend=True)
    def __init__(self, name, members=[], **kwargs):
        if not members and not isinstance(name, str):
            members = name
            name = self.__class__.__name__
        super().__init__(name, members, **kwargs)
        self._finalize()

    @property
    def sections(self):
        return self._members

    custom_value_name = "CUSTOM_SECTION"
    """ Just the name that appears in the grammar, when it is printed."""

    def _tuple_with_my_name(self, expr, delimiter=None):
        """Do not create tuple (name, value) for the root class."""
        return expr

    def parse_return(self, val, return_value_only=True):
        """Clean up the parsed values (unpack then from unnecessary containers)

        There is no name in the parsed results (see how
        ConfigurationRootDefinition._tuple_with_my_name is redefined).
        """
        val = val[0]
        return val

    def _create_grammar(self, allow_dangerous=False):
        """Returns the grammar to parse the configuration file.

        This method just tweaks the grammar (generated by the common container implementation) to ignore comments,
        so the comments would be ignored just once.
        """
        out = super()._create_grammar(allow_dangerous)
        out = self.add_ignored(out)
        return out

    def add_ignored(self, grammar):
        grammar = pp.Suppress(pp.Regex(r"(\s*\n)*")) + grammar
        grammar = grammar.ignore("#" + pp.restOfLine + pp.LineEnd())
        return grammar

    def _generic_info(self):
        return "Configuration"
