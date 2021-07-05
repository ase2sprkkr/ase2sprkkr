""" In this file the common containers of configuration values are,
either for task or potential """

from ..common.misc import OrderedDict
from ..common.grammar_types import mixed
from .options import Option
import pyparsing as pp
from .conf_common import ConfCommon
import itertools
class ConfContainer(ConfCommon):
  """ Custom task section. Section created by user with no definition """

  def __init__(self, definition, container=None):
      super().__init__(definition, container)
      self._members = OrderedDict()
      for v in definition.members():
          self._add(v.create_object(self))

  def members(self):
      return self._members

  def _get_member(self, name):
      if not name in self._members:
          raise AttributeError(f'No {name} member of {self._definition}')
      out = self._members[name]
      if out._definition.is_hidden:
          raise AttributeError(f'Member {name} of {self._definition} is not directly accessible')
      return out

  def __getattr__(self, name):
      out = self._get_member(name)
      return out

  def __getitem__(self, name):
    return self._members[name]

  def __iter__(self):
      yield from self._members.values()

  def __dir__(self):
      return itertools.chain(self._members.keys(), super().__dir__())

  def __contains__(self, name):
      return name in self._members

  def clear(self, do_not_check_required=False):
      if not do_not_check_required and self._definition.required:
         raise ValueError(f'Section {self._definition.name} cannot be cleared')
      for i in self._members.values():
          i.clear(True)

  @property
  def name(self):
      return self._definition.name

  def set(self, values, unknown='add'):
      """
      Set the values of the container

      Parameters
      ----------

      values: dict
        Values to be set

      unkwnown: 'add', 'find' or None
        If add, add unkwnown values as Custom values.
        If find, try to find the values in descendant containers.
      """

      for i,v in values.items():
          if not i in self._members:
             if unknown == 'find':
                value = self.find(i)
                if value:
                   value.set(v)
                   continue
             if not unknown == 'add':
                raise KeyError("No option with name {} in {}".format(i, str(self)))
                continue
             self.add(i, v)
          else:
             self._members[i].set(v)

  def add(self, name, value=None):
      if not getattr(self._definition, 'custom_class', False):
         raise TypeError(f'Can not add custom members to a configuration class {self._definition}')
      if name in self._members:
         raise TypeError(f'Section member {name} is already in the section {self._definition}')
      cc = self._definition.custom_class
      self._add(cc(name, self))
      if value is not None:
          self._members[name].set(value)

  def remove(self, name):
      cclass = getattr('custom_class', self._definition, False)
      if not cclass:
         raise TypeError("Can not remove items of {}".format(name))
      if not getattr(self._members[name], 'remove'):
         raise KeyError("No custom member with name {} to remove".format(name))
      del self._members[name]

  def __iter__(self):
      yield from self._members.values()

  def to_dict(self, dct=None):
      out = OrderedDict()
      for i in self:
          i.to_dict(out)
      if dct is not None and out:
          dct[self.name] = out
      return out

  def _find_value(self, name):
      for i in self:
          if i._definition.is_hidden:
             continue
          out = i._find_value(name)
          if out:
             return out

  def _add(self, member):
      self._members[member.name] = member


class BaseSection(ConfContainer):
  """ A section of SPRKKR configuration  """

  def __setattr__(self, name, value):
      if name[0]=='_':
        super().__setattr__(name, value)
      else:
        val = self._get_member(name)
        val.set(value)

  def has_any_value(self):
      for i in self:
        if i() is not None:
           return True
      return False

  def save_to_file(self, file):
      if not self.has_any_value():
         if not self._definition.is_optional:
            raise ValueError(f"Non-optional section {self._definition.name} has no value to save")
         return
      if self._definition.name_in_grammar:
         file.write(self._definition.name)
         file.write('\n')
      for o in self:
          if o.save_to_file(file):
             file.write(self._definition.delimiter)
      file.flush()

  @property
  def seciton_name(self):
      return self._definition.name


class Section(BaseSection):
  """ A standard section of a task or potential (whose content is predefinded by SectionDefinition) """

  @property
  def definition(self):
      return self._definition


class CustomSection(BaseSection):
  """ Custom task section. Section created by user with no definition """

  def remove(self):
      self._container.remove(self.name)

  @classmethod
  def factory(cls, definition_type):
      def create(name, container):
          definition = definition_type(name)
          definition.removable = True
          return cls(definition, container)
      return create


class RootConfContainer(ConfContainer):

  def save_to_file(self, file):
      """ Save the configuration to a file in a format readable by SPR-KKR """
      if not hasattr(file, 'write'):
         with open(file, "w") as file:
           return self.save_to_file(file)

      it = iter(self)
      i = next(it)
      if i:
        i.save_to_file(file)
        for i in it:
          file.write(self._definition.delimiter)
          i.save_to_file(file)
      file.flush()

  def read_from_file(self, file, clear_first=True):
      values = self._definition.parse_file(file)
      #except Exception as e:
      #   print(e)
      #   breakpoint()
      #   print(e)
      if clear_first:
         self.clear(True)
      self.set(values)

  def find(self, name):
      """ Find a configuration value of a given name in the owned sections """
      return self._find_value(name)
