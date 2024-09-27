"""
Functionality for enriching docstrings.
"""


def verbatim(s):
    return '.. code-block:: text\n\n ' + s.replace('\n', '\n ')


def process_input_parameters_definition(module, definition):
    """
    Add runtime-generated information about the task input parameters to the
    docstring of the module.
    """
    module.__doc__ = (getattr(module, '__doc__', None) or '') + """\n

**Description of the sections and parameters**

""" + verbatim(definition.description(verbose='all'))
