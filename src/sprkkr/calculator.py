# coding: utf-8
# vim: set et sw=4 ts=4 sts=4 nu ai cc=80 fdm=indent mouse=a:

"""
Module calculator
=================

This module contains different classes used to define a new calculator for
specific spectroscopies understood by MsSpec.

For more informations about SPR-KKR, follow this
`link <http://olymp.cup.uni-muenchen.de/index.php?option=com_content&
view=article&id=8%3Asprkkr&catid=4%3Asoftware&Itemid=7&lang=en>`__

"""

import os
import numpy as np
from ase.calculators.calculator import FileIOCalculator

from .misc import LOGGER
from .calcio import InputFile, PotFile



class SPRKKR(FileIOCalculator):
    command = 'kkrscf < kkr.inp > kkr.out'
    implemented_properties = []
    
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label="kkr", atoms=None, **kwargs):    
        """
        Construct SPRKKR calculator.

        """
        
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        LOGGER.debug(f'Output directory is: {self.directory}')
        LOGGER.debug(f'Output prefix is   : {self.prefix}')

        filename = os.path.join(self.directory, self.prefix + ".inp")
        self.input = InputFile(filename)
        self.atoms = None


    def set_atoms(self, atoms):
        self.atoms = atoms

    def write_input(self, atoms, properties=None, system_changes=None):
        # this will create directories
        FileIOCalculator.write_input(self, self.atoms)
        # create the input file
        self.input.write()

        # create the start pot file
        potfilename = self.input.filename.replace(".inp", ".pot")
        potfile = PotFile(self.atoms, filename=potfilename)
        potfile.write()

    def read_results(self):
        raise NotImplementedError

    def scf(self):
        input_filename = os.path.abspath(self.input.filename)
        output_filename = input_filename.replace(".inp", ".out")
        pot_filename = os.path.abspath(os.path.join(self.directory, self.prefix + ".pot"))
        self.command = "kkrscf < " + input_filename + " > " + output_filename
        self.calculate(self.atoms, None, None)


        

