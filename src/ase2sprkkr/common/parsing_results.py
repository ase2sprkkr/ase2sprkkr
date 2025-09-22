import pyparsing as pp


class Result:
    """ When data are parsed, some (repeated) options can yield
    special keys, results in agregating """

    def result(self):
        return self.value


class Key:
    """ A base class for items, that have to be treated in a special way """

    NONE = lambda x:x
    """ The regular items has no special key """

    def __init__(self, key):
        self.key = key

    def get(self, too):
        if self.key in too:
            out = too[self.key]
            if out.__class__ is not self.ResultClass:
                raise pp.ParseException(f"Conflicting values for item {self.key}")
            return out
        out = self.ResultClass()
        too[self.key] = out
        too.process.append(self.key)
        return out


class IgnoredKey(Key):
    """ This key is totaly and silently ignored. Use it
    for the keys, that should be written, but not read.
    """

    def add(self, too, val):
        return


class ValidateKey(Key):
    """ This key is totaly and silently ignored. Use it
    for the keys, that should be written, but not read.
    """
    def add(self, too, val):
        too.checks.append(val)


class SubKey(Key):
    """ A base class for items that have subkeys """

    def __init__(self, key, sub):
        self.key = key
        self.sub = self.convert(sub)

    def convert(self, sub):
        return int(sub)


class DictKey(SubKey):

    class ResultClass(Result):

        def __init__(self):
            self.value = {}

    def add(self, too, value):
        out=self.get(too)
        if self.sub in out.value:
             raise pp.ParseException(f"Duplicate key {self.sub} in {self.key}")
        out.value[self.sub] = value


class DefDictKey(DictKey):

    def convert(self, sub):
        return sub if sub == 'def' else int(sub)


class ArrayKey(SubKey):

    class ResultClass(Result):

        NOT_SET = object()

        def __init__(self):
            self.value = []

        def result(self):
            return self.value

    def add(self, too, val):
        out = self.get(too)
        ln = len(out.value)
        i = self.sub - 1
        if ln < i:
            out.value+=[out.NOT_SET] * (i - ln + 1)
        if ln == i:
            out.value.append(val)
        else:
            if out.value[i] is not out.NOT_SET:
                raise pp.ParseException(f"Duplicate key {self.sub} in {self.key}")
            out.value[i]=val


class RepeatedKey(Key):

    class ResultClass(Result):

        def __init__(self):
            self.value = []

    def add(self, too, val):
        self.get(too).value.append(val)


class Values(dict):
    """ Result of dict_from_parsed: dictionary with list of checks on the parsed values."""
    def __init__(self):
        super().__init__()
        self.checks = []
        self.process = []

    def to_dict(self):
        return { i: j.to_dict() if isinstance(j, Values) else j
                    for i,j in self.items() }


def dict_from_parsed(values):
    """ Create a dictionary from the arguments.
    From duplicate arguments create numpy arrays.
    Moreover, if there is key of type (a,b), it will be transformed to subdictionary.
    Such a keys do not allow duplicates.

    >>> dict_from_parsed( [ ('x', 'y'),  ((DictKey('a', 1)), 1 ), ((DictKey('a', 3)), 2) ] )
    {'x': 'y', 'a': {1: 1, 3: 2}}
    >>> dict_from_parsed( [ ('x', 1), ('x', '2') ] ) # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    pyparsing.exceptions.ParseException: There are non-unique keys: x
    >>> dict_from_parsed( [ (RepeatedKey('x'), 1), (RepeatedKey('x'), 2) ] )
    {'x': [1, 2]}
    """
    out = Values()
    duplicates = set()
    errors = []

    def add(key, value):
        if isinstance(key, Key):
            key.add(out, value)
        elif key in out:
            duplicates.add(k)
        else:
            out[key] = value

    for k, v in values:
        try:
           add(k, v)
        except Exception as e:
           errors.append(e)

    for key in out.process:
        out[key] = out[key].result()

    for i in out.checks:
         try:
           i(out)
         except Exception as e:
           errors.append(e)

    if duplicates:
        duplicates = ", ".join((i.upper() for i in duplicates))
        errors.append(pp.ParseException(f"There are duplicate items named {duplicates}"))

    if errors:
        if len(errors) == 1:
            raise errors[0]
        errors = '\n'.join((str(e) for e in errors))
        raise pp.ParseException(f"There are errors in data: {errors}")
    return out
