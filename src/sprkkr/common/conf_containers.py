from ..common.misc import OrderedDict
from ..common.grammar_types import mixed
from .options import Option
import pyparsing as pp

class ConfContainer:
  """ Custom task section. Section created by user with no definition """

  def __init__(self, definition):
      self._definition = definition
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

  def set(self, values, allow_add=True, ):
      for i in values:
          if not i in self._members:
             if not allow_add:
                raise KeyError("No option with name {} in {}".format(i, str(self)))
             self.add(i, values[i])
          else:
             self._members[i].set(values[i])

  def add(self, name, value=None):
      if not getattr(self._definition, 'custom_class', False):
         raise TypeError(f'Can not add custom members to a configuration class {self._definition}')
      if name in self._members:
         raise TypeError(f'Section member {name} is already in the section {self._definition}')
      cc = self._definition.custom_class
      self._add(cc(self, name))
      if value is not None:
          self._members[name].set(value)

  def _add(self, member):
      self._members[member.name] = member

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

class BaseSection(ConfContainer):
  """ A section of SPRKKR configuration  """

  def __init__(self, definition, container=None):
      super().__init__(definition)
      self._container = container

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
         return
      file.write(self._definition.name)
      file.write('\n')
      for o in self:
          if o.save_to_file(file):
             file.write(self._definition.delimiter)

  @property
  def seciton_name(self):
      return self._definition.name


class Section(BaseSection):

  @property
  def definition(self):
      return self._definition


class CustomSection(BaseSection):
  """ Custom task section. Section created by user with no definition """

  def remove(self):
      self._container.remove(self.name)

  @classmethod
  def factory(cls, definition_type):
      def create(container, name):
          definition = definition_type(name)
          definition.removable = True
          return cls(definition, container)
      return create


class RootConfContainer(ConfContainer):

  def save_to_file(self, file):
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

  def read_from_file(self, file, clear_first=True):
      values = self._definition.grammar().parseFile(file, True)
      #except Exception as e:
      #   print(e)
      #   breakpoint()
      #   print(e)
      assert len(values) == 1
      values = values[0]
      if clear_first:
         self.clear(True)
      self.set(values)
