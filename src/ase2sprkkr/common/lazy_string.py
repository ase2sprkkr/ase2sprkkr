from .decorators import cached_property


class LazyString:
    """ Lazy string is a string, evaluated only by demand """

    def __init__(self, value):
        self._value = value

    @cached_property
    def value(self):
        return self._value()

    def __str__(self):
        return self.value

    def __add__(self, other):
        return LazyString(lambda: self.value + other)

    def __iadd__(self, other):
        if 'value' in self.__dict__:
            del self.__dict__['value']
        self.value = lambda: self.value + other

    def __radd__(self, other):
        return LazyString(lambda: other + self.value)

    def __iter__(self):
        return iter(self.value)

    def __getattr__(self, name):
        return getattr(self.value, name)
