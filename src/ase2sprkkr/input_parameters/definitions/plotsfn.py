"""Native PLOTSFN selector."""
from ._plot_selectors import plot_parameters

def input_parameters_definition():
    return plot_parameters("PLOTSFN", 'Plot the radial spherical-harmonic expansion coefficients of shape functions from an existing FULLPOT potential, for each mesh. ASA has no such output.')
