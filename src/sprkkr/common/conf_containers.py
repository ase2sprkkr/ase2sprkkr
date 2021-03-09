from collections import OrderedDict
from ..common.grammar_types import mixed
from .options import Option

class ConfContainer:

  _item_class = Option

  def __init__(self, definition):
      self._definition = definition
      self._members = OrderedDict()
      for v in definition.members():
          self._members[v.name] = self._item_class(v)

  def __getattr__(self, name):
      return self._members[name]

  def __getitem__(self, name):
    return self._members[name]

  def __iter__(self):
      yield from self._members.values()

  def __contains__(self, name):
      return name in self._members

  def clear(self):
      for i in self._members.values():
          i.clear()

  @property
  def name(self):
      return self._definition.name

  def set(self, values, allow_add=True):
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
      cc = self._definition.custom_class
      self._members[name] = cc(self, name)
      if value:
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

class BaseSection(ConfContainer):
  """ A section of SPRKKR configuration  """

  def __setattr__(self, name, value):
      if name[0]=='_':
        super().__setattr__(name, value)
      else:
        self._options[name].set(value)

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
  def __init__(self, container, definition):
      super().__init__(definition)
      self._container = container

  def remove(self):
      self._container.remove(self.name)

  @classmethod
  def factory(cls, definition_type):
      def create(container, name):
          definition = definition_type(name)
          definition.removable = True
          return cls(container, definition)
      return create



class RootConfContainer(ConfContainer):

  _item_class = Section

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
      assert len(values) == 1
      values = values[0]
      if clear_first:
         self.clear()
      self.set(values)


