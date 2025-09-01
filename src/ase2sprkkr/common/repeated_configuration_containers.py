from .configuration_containers import BaseConfigurationContainer
from typing import Union, Any, List
from .warnings import DataValidityError


class RepeatedConfigurationContainer(BaseConfigurationContainer):
    """ A group, that can be repeated """

    def __init__(self, definition, container=None):
         """ Create the container and its members, according to the definition """
         super().__init__(definition, container)
         """
         The members of the container, in a form of ``{obj.name : obj}``
         """
         self._values = []

    def __getitem__(self, name):
         """
         The members of the container are accesible using ``container["member name"]`` notation.
         """
         return self._values[name]

    def __len__(self):
         return len(self._values)

    def __bool__(self):
         return True

    def add(self):
         out=self._definition.create_object(self, repeated=False)
         self._values.append(out)
         return out

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
        self._values = []

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
           return self._values[int(section)].get(name)
        name = int(name)
        if name in self._values:
           val = self._values[name]
        else:
           val = None
        if not val:
           raise ValueError(f"No {name} member of {self}")
        return val.get()

    def set(self, values:Union[List,None]={}, value=None, *, unknown='find', error=None, **kwargs):
        self._set(values, value, unknown=unknown, error=error, **kwargs)

    def _set(self, values:Union[List,None]={}, value=None, *, unknown='find', error=None, **kwargs):
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
        self._values = []

        if values:
            for v in values:
               self.add().set(v)

    def __iter__(self):
        """ Iterate over all members of the container """
        yield from self._values

    def items(self):
        yield from enumerate(self._values)

    def values(self):
        yield from self._values

    def _as_dict(self, only_changed:Union[bool,str]='basic', generated:bool=False, copy=False):
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
        out = [
            v.as_dict(only_changed, generated, copy)
            for v in self._values
        ]
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
        out = False
        for i in self.values():
            d=self._definition
            if d._save_to_file(file, i, always, name_in_grammar, delimiter=delimiter):
                name_in_grammar=False
                delimiter=d.repeated_delimiter
                out=True
        return out

    def _validate(self, why:str='save'):
        """ Validate the configuration data. Raise an exception, if the validation fail.

        Parameters
        ----------
        why
          Type of the validation. Possible values
          ``save`` - Full validation, during save.
          ``set`` - Validation on user input. Allow required values not to be set.
          ``parse`` - Validation during parsing - some check, that are enforced by the parser, can be skipped.
        """
        if why == 'save' and not self._definition.is_optional and not self.has_any_value():
            DataValidityError.warn(f"Non-optional section {self._definition.name} has no value to save")
        for o in self.values():
            o._validate(why)
