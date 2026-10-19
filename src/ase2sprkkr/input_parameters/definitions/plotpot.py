"""Native PLOTPOT selector."""
from ._plot_selectors import plot_parameters

def input_parameters_definition():
    return plot_parameters("PLOTPOT", 'Plot the stored radial potential and exchange field; FULLPOT includes nonspherical potential components.')
