"""Native PLOTRHO selector."""
from ._plot_selectors import plot_parameters

def input_parameters():
    return plot_parameters("PLOTRHO", 'Plot the stored charge and spin densities; FULLPOT includes nonspherical density components.')
