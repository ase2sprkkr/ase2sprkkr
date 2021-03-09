import functools
_x = object()

def lazy_value(fce):
    """ The decorator for only once computed value """
    x = []

    """ Hack for staticmethod decorator, which is in fact binded by the descriptor protocol """
    if isinstance(fce, staticmethod):
        fce = fce.__func__

    @functools.wraps(fce)
    def cached_fce():
        if not x:
           x.append(fce())
        return x[0]
    return cached_fce
