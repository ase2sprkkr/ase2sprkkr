"""
Functionality for enriching docstrings.
"""
import importlib


def verbatim(s):
    return '.. code-block:: text\n\n ' + s.replace('\n', '\n ')


def process_input_parameters_definition(name):
    """
    Add runtime-generated information about the task input parameters to the
    docstring of the module.
    """
    m=importlib.import_module(name)
    ip = m.input_parameters
    m.__doc__ = m.__doc__ + """\n

**Description of the sections and parameters**

""" + verbatim(ip.description(verbose='all'))
