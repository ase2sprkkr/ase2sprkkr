from .configuration_definitions import RealItemDefinition, BaseDefinition
from .grammar import delimitedList
from .misc import dict_first_item
from .repeated_configuration_containers import RepeatedConfigurationContainer
from .configuration_containers import Section
from .decorators import cache

import pyparsing as pp
import numpy as np
from typing import Union
import itertools
from collections.abc import Iterable
import inspect
from io import StringIO

# This serves just for dealing with various pyparsing versions
_parse_all_name = 'parse_all' if \
  'parse_all' in inspect.getfullargspec(pp.Or.parseString).args \
  else 'parseAll'


def dict_from_parsed(values, allowed_duplicates):
    """ Create a dictionary from the arguments.
    From duplicate arguments create numpy arrays.
    Moreover, if there is key of type (a,b), it will be transformed to subdictionary.
    Such a keys do not allow duplicates.

    >>> dict_from_parsed( [ ('x', 'y'), (('a','b'), 1 ), (('a', 'c'), 2) ], [] )
    {'x': 'y', 'a': {'b': 1, 'c': 2}}
    >>> dict_from_parsed( [ ('x', 1), ('x', '2') ], [] ) # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    pyparsing.exceptions.ParseException: There are non-unique keys: x
    >>> dict_from_parsed( [ ('x', 1), ('x', 2) ], ['x'] )
    {'x': array([1, 2])}
    """
    out = {}
    duplicates = []
    errors = []

    if isinstance(allowed_duplicates, list):
       _allowed_duplicates = lambda x: x in allowed_duplicates
    else:
       _allowed_duplicates = allowed_duplicates

    def add(out, k, key, check):
         if k in out:
            if key in duplicates:
               out[k].append(v)
            else:
               if not _allowed_duplicates(check):
                   return errors.append(check)
               out[k] = [ out[k], v ]
               duplicates.append(k)
         else:
            out[k] = v

    for k, v in values:
        if isinstance(k, tuple):
           if not k[0] in out:
              o = out[k[0]] = {}
           elif not isinstance(out[k[0]], dict):
              raise pp.ParseException(f"There is duplicate key {k[0]}, however, I should created "
                                "both arrays and dict from the data, which is not possible")
           else:
              o = out[k[0]]
           add(o, k[1], k, k[0])
        else:
           add(out, k, k, k)
    for k in duplicates:
        if isinstance(k, tuple):
            o=out[k[0]]
            k=k[1]
        else:
            o=out
        o[k] = np.array(out[k])
    if errors:
       errors = ','.join(errors)
       raise pp.ParseException(f"There are non-unique keys: {errors}")
    return out


def add_excluded_names_condition(element, names):
    """ Add the condition to the element, that
    its value is not any of given names """
    if not names:
       return
    names = set((i.upper() for i in names))
    element.addCondition(lambda x: x[0].upper() not in names)


class ContainerDefinition(RealItemDefinition):
    """ Base class for a definition (of contained data, format, etc)
    of either a whole configuration file
    (e.g. :class:`InputParameters<ase2sprkkr.input_parameters.input_parameters.InputParameters>` or
    e.g. :class:`Potential<ase2sprkkr.potentials.potentials.Potential>`) or
    its :class:`Section<ase2sprkkr.common.configuration_containers.Section>`.
    """

    force_order = False
    """ Force order of its members """

    value_name_format = None
    """ The (print) format, how the name is written """

    write_last_delimiter = True

    @staticmethod
    def _dict_from_named_values(args, items=None):
        """auxiliary method that creates dictionary from the arguments"""
        items = items or {}
        for value in args:
           items[value.name] = value
        return items

    def __init__(self, name, members=[], alternative_names=[], info=None, description=None,
                 is_optional=False, is_hidden=False, is_expert=False,
                 has_hidden_members=False, name_in_grammar=None, force_order=None,
                 write_alternative_name:bool=False, result_class=None,
                 is_repeated=False
                 ):
       """
       Definition of container (e.g. section of an input file).
       For the rest of the parameters see the :class:`RealItemDefinition`.

       Parameters
       ----------
       has_hidden_members: bool
         If true, this section is not intended for a direct editing
       is_repeated: bool or string
         The section can be repeated. The name of the section appears only once on the beginning (this differs from ValueDefinition.is_repeated #TODO - merge the meaning of the swtich).
         If a non-empty string is given, the values are divided by the string.

       force_order: bool
         If True, the items has to retain the order, if False, the items can be in the input file in any order.
       """

       super().__init__(
           name = name,
           alternative_names = alternative_names,
           is_optional = is_optional,
           is_hidden = is_hidden,
           is_expert = is_expert,
           name_in_grammar = name_in_grammar,
           info = info,
           description = description,
           write_alternative_name = write_alternative_name,
           result_class = result_class
       )

       if not isinstance(members, dict):
          members = self._dict_from_named_values(members)

       if self.value_name_format:
          for i in members.values():
              i.value_name_format = self.value_name_format
       self._members = members
       for i in self._members.values():
           i.added_to_container(self)

       self.has_hidden_members = has_hidden_members
       if force_order is not None:
          self.force_order = force_order
       self.is_repeated = is_repeated

    configuration_type_name = 'SECTION'
    """ Name of the container type in the runtime documentation """

    def __repr__(self):
        return f"<{self.configuration_type_name} {self.name}>"

    def data_description(self, verbose:Union[bool,str,int]=False, show_hidden:bool=False, prefix:str=''):
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
               out+=' '
            out+=self.name
            return out

        out = f"{prefix}{container_name()}"
        flags = []
        if self.force_order:
           flags.append('fixed-order')
        if self.is_hidden:
           flags.append('hidden')
        if self.is_optional:
           flags.append('optional')
        if self.is_expert:
           flags.append('expert')
        if self.is_repeated:
           flags.append('repeated')
        if flags:
           flags = ', '.join(flags)
           out+=f" ({flags})"

        if verbose:
           if isinstance(verbose, int):
              verbose-=1
           else:
              verbose = verbose if verbose=='all' else False
           add = self.additional_data_description(verbose, show_hidden, prefix)
           if add:
              out+=' contains:'
              under=prefix + "-" * len(out) + '\n'
              out=f"{prefix}{out}\n{under}{add}"
        return out

    def allow_duplication(self):
        return self.is_repeated

    def additional_data_description(self, verbose:Union[bool,str,int]=False, show_hidden=False, prefix:str=''):
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
        cprefix=prefix + self._description_indentation
        out = []

        def write(i):
           s = i.data_description(verbose, show_hidden, cprefix)
           if not '\n' in s:
              info = i.info(False)
              if info:
                 s = s + (' ' * (max(40 - len(s), 0) + 2)) + info
           else:
              ccprefix = cprefix + i._description_indentation
              s+='\n\n'
              s+= ccprefix + i.info(False).replace('\n', '\n' + ccprefix)
              s+='\n'
           out.append(s)

        expert = False
        for i in self:
            if i.is_hidden and not show_hidden:
               continue
            if not i.is_expert:
               write(i)
            else:
               expert=True

        if expert:
          out.append(f'\n{cprefix}Expert options:')
          out.append(cprefix + '--------------')
          cprefix+=self._description_indentation

          for i in self:
              if i.is_expert:
                 if i.is_hidden and not show_hidden:
                    continue
                 write(i)

        return '\n'.join(out)

    def __iter__(self):
        return iter(self._members.values())

    def members(self):
        return self._members.values()

    def names(self):
        return self._members.keys()

    def __getitem__(self, key):
        return self._members[key]

    def __setitem__(self, key, value):
        self._members[key]=value

    def __contains__(self, key):
        return key in self._members

    def remove(self, name):
        del self._members[name]
        return self

    def copy(self, args=[], items=[], remove=[], defaults={}, **kwargs):
        """ Copy the section with the contained values modified by the arguments."""
        members = dict( ( (k,i.copy()) for k,i in self._members.items() ) )
        for i in remove:
            del members[i]
        members.update(self._dict_from_named_values(args, items))
        for i,v in defaults.items():
            members[i].default_value = members[i].type.convert(v)

        default = { k: getattr(self, v) for k,v in self._get_copy_args().items() }
        default.update(kwargs)
        default['members'] = members
        return self.__class__(**default)

    def copy_member(self, name) -> BaseDefinition:
        """ Copy a member, allowing to redefine its properties.

            Returns
            -------
            new_member: BaseDefinition
              The newly created member
        """
        out = self._members[name].copy()
        self._members[name] = out
        return out

    def all_member_names(self):
        return itertools.chain.from_iterable(
            i.all_names_in_grammar() for i in self
        )

    def _grammar_of_values(self, allow_dangerous:bool=False, delimiter=None):
       if self.custom_class:
          custom_value = self.custom_member_grammar(self.all_member_names())
       else:
          custom_value = None
       delimiter = delimiter or self.grammar_of_delimiter

       def repeated_grammars():
           """ If the item can be repeated, do it here - we don't know, whether there is a fixed order in any way
           (e.g. the item is followed by the items without name in grammar)
           """
           for i in self._members.values():
               g = i._grammar and i._grammar(allow_dangerous)
               if not g:
                   continue
               if i.can_be_repeated:
                   dlmtr = delimiter if i.can_be_repeated is True else i.can_be_repeated
                   g = delimitedList(g, dlmtr)
               yield i,g

       def grammars():
         """ This function iterates over the items of the container, joining all the without name_in_grammar with the previous ones. """
         it = iter(repeated_grammars())
         head_item, grammar_chain = next(it)

         for item, grammar in it:
             if item.is_independent_on_the_predecessor:
               yield head_item, grammar_chain
               head_item, grammar_chain = item, grammar
             else:
               add = delimiter + grammar
               if item.is_optional:
                  add = pp.Optional(add)
               if item.condition and self.force_order:
                  add = item.condition.prepare_grammar(item, add)
               grammar_chain = grammar_chain + add
         yield head_item, grammar_chain

       if self.force_order:

           init = pp.Empty()

           def set_loc(loc, toks):
               init.location = loc
           init.setParseAction(set_loc)

           first = pp.Empty().addCondition(lambda loc, toks: loc == init.location)
           if custom_value:
              cvs = pp.ZeroOrMore(custom_value + delimiter).setName('<custom...>')
              after = delimiter + cvs
           else:
              after = pp.Forward() << delimiter
           after.addCondition(lambda loc, toks: loc != init.location)
           inter_cvs = (first | after).setName('<?DELIM>')
           inter = (first | delimiter.copy().addCondition(lambda loc, toks: loc != init.location))

           def sequence():
               for head,g in repeated_grammars():
                   if head.is_independent_on_the_predecessor:
                      delim = inter_cvs
                   else:
                      delim = inter
                   g = delim + g
                   if head.is_optional:
                      g = pp.Optional(g)
                   if head.condition:
                      yield head.condition.prepare_grammar(head, g)
                   else:
                      yield g

           values = pp.And([ i for i in sequence()])

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
           values = pp.MatchFirst([i for head,i in it])
           if custom_value:
               values |= custom_value
           values = delimitedList(values, delimiter)
           if first:
               values = first + pp.Optional(delimiter + values)

       values.setParseAction(lambda x: dict_from_parsed(x.asList(), self._allow_duplicates_of))

       if self.validate:
          def _validate(s, loc, value):
              # just pass the dict to the validate function
              is_ok = self.validate(MergeDictAdaptor(value[0], self), 'parse')
              if is_ok is not True:
                if is_ok is None:
                   is_ok = f'Validation of parsed data of {self.name} section failed'
                raise pp.ParseException(s, loc, is_ok)
              return value
          values.addParseAction(_validate)
       if self.is_repeated:
          rdelim = delimiter
          if self.is_repeated is not True:
              rdelim = rdelim + pp.Literal(self.is_repeated)
          values = pp.delimitedList(values, rdelim)
          values.addParseAction(lambda x: [x.asList()])
       return values

    def _allow_duplicates_of(self, name):
        """ Can a given element (identified by name) have more values in the parsed results?
            (However, not all definitions have to specify allow_duplicates, just the ones
            that have a value). For the others, this function raises an error.
        """
        return self[name].allow_duplication()

    def _create_grammar(self, allow_dangerous=False):
       delimiter = self.grammar_of_delimiter
       values = self._grammar_of_values(allow_dangerous, delimiter)
       out = self._tuple_with_my_name(values, delimiter)
       out.setName(self.name)
       return out

    validate = None
    """ A function for validation of just the parsed result (not the user input) """

    @classmethod
    @cache
    def delimited_custom_value_grammar(cls):
        """ Return the grammar for the custom child with delimiter.
        The delimiter can delimite it either from the previous child or from the section name."""

        return cls.child_class.grammar_of_delimiter + cls.custom_value_grammar()

    custom_name_characters = pp.alphanums + '_-()'
    """ Which characters can appears in an unknown child (value/section) name """

    @classmethod
    def custom_member_grammar(cls, value_names = []):
       """ Grammar for the custom - unknown - child """
       name = pp.Word(cls.custom_name_characters).setParseAction(lambda x: x[0].strip())
       add_excluded_names_condition(name, value_names)
       out = (name + cls.delimited_custom_value_grammar()).setParseAction(lambda x: tuple(x))
       out.setName(cls.custom_value_name)
       return out

    def _first_section_has_to_be_first(self):
       """ Has/ve the first child(s) in an unordered sequence fixed position? """
       return not dict_first_item(self._members).is_independent_on_the_predecessor

    def parse_file(self, file, return_value_only=True, allow_dangerous=False):
       """ Parse the file, return the parsed data as dictionary """
       grammar = self.grammar(allow_dangerous)
       out = grammar.parseFile(file, **{ _parse_all_name: True } )
       return self.parse_return(out, return_value_only)

    def parse(self, string, whole_string=True, return_value_only=True, allow_dangerous=False):
       """ Parse the string, return the parsed data as dictionary """
       grammar = self.grammar(allow_dangerous)
       out = grammar.parseString(string, **{ _parse_all_name: whole_string } )
       return self.parse_return(out, return_value_only)

    def parse_return(self, val, return_value_only:bool=True):
        """ Clean up the parsed values (unpack then from unnecessary containers)

        Parameters
        ----------
        return_value_only
          Return only value, not name - value tuple

        """
        val = val[0]
        if return_value_only:
           val = val[1]
        return val

    async def parse_from_stream(self, stream, up_to, start=None, whole_string=True, return_value_only=True, allow_dangerous=False):
        """
        Parse string readed from asyncio stream.
        The stream is readed up to the given delimiter
        """

        result = await stream.readuntil(up_to)
        result = result[:-len(up_to)].decode('utf8')
        if start:
           result = start + result
        return self.parse(result, whole_string)

    def read_from_file(self, file, allow_dangerous=False, **kwargs):
        """ Read a configuration file and return the parsed Configuration object """
        out = self.result_class(definition = self, **kwargs)
        out.read_from_file(file, allow_dangerous=allow_dangerous)
        return out

    def read_from_dict(self, values, **kwargs):
        out = self.result_class(definition = self, **kwargs)
        out.set(values, unknown='add')
        return out

    def read_from_string(self, string, allow_dangerous=False, **kwargs):
        return self.read_from_file(StringIO(string), allow_dangerous, **kwargs)

    def validate(self, container, why:str='save'):
        self.validate_warning(container)
        return True

    repeated_class = RepeatedConfigurationContainer
    """ Class for the repeated sections """

    def create_object(self, container=None, repeated:bool=True):
        """ Create an instance (section)

        container: BaseConfigurationContainer
            To which container the created object will belong

        repeated:
            Has meaning only for a is_repeated section. Then, if it is True,
            a Container for repeated values of the section is returned.
            Otherwise, the container just for one instance of a section is
            returned.
        """
        if repeated and self.is_repeated:
            return self.repeated_class(self, container)
        return super().create_object(container)

    def _save_to_file(self, file, value, always=False, name_in_grammar=None, delimiter='')->bool:
        """ Save the content of the container to the file (according to the definition)

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
            if not self.write_condition(self) or (self.condition and not self.condition(self)):
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
           file.write(delimiter)
        if name_in_grammar:
           file.write(self.name)
           file.write('\n')

        members = iter(value)
        if self.write_last_delimiter:
           for o in members:
               if o._save_to_file(file, always):
                   file.write(self.delimiter)
        else:
           delimiter = ''
           for o in members:
               if o._save_to_file(file, always, delimiter=delimiter):
                   delimiter=self.delimiter

        return True


class SectionDefinition(ContainerDefinition):
   """ Base class for definition of the sections in Pot or InputParameters files.

       It just redefine a few properties/methods to values/behavior typical for the sections
   """

   result_class = Section

   @property
   def values(self):
       return self._members

   custom_value_name = 'CUSTOM_VALUE'
   """ Just the name that appears in the grammar, when it is printed."""

   @classmethod
   @cache
   def delimited_custom_value_grammar(cls):
        gt = cls.custom_class.grammar_type
        # here the child (Value) class delimiter should be used
        out = cls.child_class.grammar_of_delimiter + gt.grammar()
        optional, df, _ = gt.missing_value()
        if optional:
           out = out | pp.Empty().setParseAction(lambda x: df)
        return out

   def _generic_info(self):
      return f"Configuration section {self.name}"

   def accept_value(self, value) -> bool:
       if isinstance(value, dict):
           return True
       return self.is_repeated and isinstance(value, Iterable)


class ConfigurationRootDefinition(ContainerDefinition):
   """ From this class, the definition of the format of a whole configuration file should be derived.

   """
   write_last_delimiter = False
   """ Do not print additional newline after the last section """

   name_in_grammar = False
   """ The configuration files has commonly no "name" of its content, they
   just contains their content.

   However, in some cases the name_in_grammar could be used for some kind of
   prefix in the file, however, it is better to use a fixed value for this purpose.
   """

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

       return cls(( gen(i) for i in defs))

   def __init__(self, name, members=[], **kwargs):
       if not members and not isinstance(name, str):
          members = name
          name = self.__class__.__name__
       super().__init__(name, members, **kwargs)

   @property
   def sections(self):
       return self._members

   custom_value_name = 'CUSTOM_SECTION'
   """ Just the name that appears in the grammar, when it is printed."""

   def _tuple_with_my_name(self, expr, delimiter=None):
       """ Do not create tuple (name, value) for the root class. """
       return expr

   def parse_return(self, val, return_value_only=True):
        """ Clean up the parsed values (unpack then from unnecessary containers)

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
       out=super()._create_grammar(allow_dangerous)
       out=self.add_ignored(out)
       return out

   def add_ignored(self, grammar):
       grammar = pp.Suppress(pp.Regex(r'(\s*\n)*')) + grammar
       grammar = grammar.ignore("#" + pp.restOfLine + pp.LineEnd())
       return grammar


class MergeDictAdaptor:
    """ This class returns a read-only dict-like class
    that merge values from a container and from the
    definition of a section """

    def __init__(self, values, definition):
        self.values = values
        self.definition = definition

    def __hasitem__(self, name):
        return self.values.__hasitem__(name) or \
               self.definition.__hasitem__(name)

    def __getitem__(self, name):
        try:
            return self.values[name]
        except KeyError:
            return self.definition[name].get_value()

    def __repr__(self):
        return f'Section {self.definition.name} with values {self.values}'
