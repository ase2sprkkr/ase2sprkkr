def check_symmetry(pbc, symmetry, error, err_class = ValueError):
    """ Check the symmetry of the PCB, raise an exception, if they not satisfy the requirements """
    if not (pbc == symmetry).all():
       if callable(error):
           error = error()
       raise err_class(error)
