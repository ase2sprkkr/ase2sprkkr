from .configuration_containers import BaseConfigurationContainer, DictAdaptor
from typing import Union, Any, Dict


class RepeatedConfigurationContainer(BaseConfigurationContainer):
    """ A container for configuration (problem-definition) options and/or sections.

    Options in the configuration (problem-definition) files are grouped to
    sections, sections are then grouped in a configuration file object.
    This is a base class for these containers.
    """

    def __init__(self, definition, container=None):
         """ Create the container and its members, according to the definition """
         super().__init__(definition, container)
         """
         The members of the container, in a form of ``{obj.name : obj}``
         """
         self._values = {}

    def __getitem__(self, name):
         """
         The members of the container are accesible using ``container["member name"]`` notation.
         """
         return self._values[name]

    def __len__(self):
         return len(self._values)

    def __bool__(self):
         return True

    def add(self, id):
         out=self._definition.create_object(self, repeated=False)
         self._values[id]=out
         return out

    def __contains__(self, name):
        """ The check for existence of a member with the given name."""
        return name in self._values

    def clear(self, do_not_check_required=False, call_hooks=True, generated=None):
        """
        Erase all values (or reset to default) from all options in the container
        (ad subcontainers)

        Parameters
        ----------
        do_not_check_required: bool
          Do not check validity of the values after clearing. If ``False`` (default)
          is passed as this argument, the required option without a default value
          (or a section containing such value) throw an exception, which prevents the
          clearing (neverthenless, previous values in the section will be cleared anyway).

        call_hooks: bool
          If False, the cleared values do not raise theirs hooks

        generated: bool
          If True
        """
        self._values = {}

    def get(self, name=None, unknown='find'):
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
        if '.' in name:
           section, name = name.split('.')
           return self._values[section].get(name)
        if name in self._values:
           val = self._values[name]
        else:
           val = None
        if not val:
           raise ValueError(f"No {name} member of {self}")
        return val.get()

    def set(self, values:Union[Dict[str,Any],str,None]={}, value=None, *, unknown='find', error=None, **kwargs):
        """
        Set the value(s) of parameter(s). Usage:

        > input_parameters.set({'NITER': 5, 'NE': [10]})
        or
        > input_parameters.set(NITER=5, NE=[10])

        Parameters
        ----------

        values:
          Dictionary of values to be set, or the name of the value, if the value is given.

        value:
          Value to be set. Setting this argument require to pass string name to the values argument.

        unkwnown: 'add', 'find' or None
          How to handle unknown (not known by the definition) parameters.
          If 'find', try to find the values in descendant containers.
          If 'add', add unknown values as custom values.
          If None, throw an exception.
          Keyword only argument.

        **kwargs: dict
          The values to be set (an alternative syntax as syntactical sugar)
        """
        if values.__class__ is str:
           values = { values : value }
        elif value is not None:
           raise ValueError("If value argument of Container.set method is given,"
           " the values have to be string name of the value")

        try:
            items = values.items()
        except AttributeError:
            items = enumerate(values)

        for k,v in items:
           if not k in self._values:
               self.add(k).set(v)
           else:
               self._values[k].set(v)

    def __iter__(self):
        """ Iterate over all members of the container """
        yield from self._values.keys()

    def items(self):
        return self._values.items()

    def values(self):
        return self._values.values()

    def as_dict(self, only_changed:Union[bool,str]='basic', generated:bool=False, copy=False):
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
        out = {}
        for k,v in self.items():
            value = k.as_dict(only_changed, generated, copy)
            if value is not None:
                out[k] = value
        return out or None

    def is_changed(self):
        for i in self.values():
            if i.is_changed():
                return True
        return False

    def _save_to_file(self, file, always=False, name_in_grammar=None, delimiter='')->bool:
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
        delim = '' if self._definition.is_repeated is True else self._definition.is_repeated
        out = False
        for i in self.values():
            if self._definition._save_to_file(file, i, always, name_in_grammar, delimiter=delimiter):
                name_in_grammar=False
                delimiter=delim
                out=True
        return out

    def validate(self, why:str='save'):
        """ Validate the configuration data. Raise an exception, if the validation fail.

        Parameters
        ----------
        why
          Type of the validation. Possible values
          ``save`` - Full validation, during save.
          ``set`` - Validation on user input. Allow required values not to be set.
          ``parse`` - Validation during parsing - some check, that are enforced by the parser, can be skipped.
        """
        for i in self.values():
            self._definition.validate(DictAdaptor(i), why)
            if why == 'save' and not self._definition.is_optional and not self.has_any_value():
                raise ValueError(f"Non-optional section {self._definition.name} has no value to save")
        for o in self.values():
            o.validate(why)
