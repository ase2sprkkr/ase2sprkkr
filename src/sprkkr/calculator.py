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
from ase.units import Rydberg

from .misc import LOGGER
from .calcio import InputFile, PotFile

class SPRKKR(FileIOCalculator):
    command = 'kkrscf < kkr.inp > kkr.out'
    implemented_properties = ['energy']
    
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label="kkr", task='scf',atoms=None, **kwargs):    
        """
        Construct SPRKKR calculator.

        """
        
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        LOGGER.debug(f'Output directory is: {self.directory}')
        LOGGER.debug(f'Output prefix is   : {self.prefix}')
        

        self.atoms = None
        self.task = task

        self.inpfile=os.path.join(self.directory, self.prefix + ".inp")
        self.output= self.inpfile.replace(".inp", ".out")
        self.potfile=self.inpfile.replace(".inp", ".pot")
        self.sysfile=self.inpfile.replace(".inp", ".sys")
        
        LOGGER.debug(f'INP FILE:{self.inpfile}')
        LOGGER.debug(f'POT FILE:{self.potfile}')        
        LOGGER.debug(f'OUT FILE:{self.output}')                
        self.input = InputFile(self.inpfile,self.task)
        if (self.task.upper()=='SCF'):
            input_filename = os.path.abspath(self.input.filename)
            output_filename = input_filename.replace(".inp", ".out")
            pot_filename = os.path.abspath(os.path.join(self.directory, self.prefix + ".pot"))
            self.command = "mpirun.openmpi -np 4 kkrscfMPI  " + input_filename + " > " + output_filename

        else:
            print("TASK {} not implemeted in ASE",format(self.task))
            raise NotImplementedError    

    def set_potfile(self,filename):
        self.potfile=filename
        self.pot=PotFile(self.atoms, filename=self.potfile)
        
    def set_inpfile(self,filename):
        self.inpfile=filename
        self.input= InputFile(self.inpfile,self.task)

    def set_outfile(self,filename):
        self.outfile=filename
            
    def set_atoms(self, atoms):
        self.atoms = atoms
 
    def write_pot(self):
        self.pot=PotFile(self.atoms, filename=self.potfile,sysfilename=self.sysfile)
        self.pot.write()
        self.pot.write_sys()
        
    def write_input(self, atoms, properties=None, system_changes=None):
        # this will create directories
        FileIOCalculator.write_input(self, self.atoms)
        potline=os.path.join(self.prefix + ".pot")
        self.input.control_section.set(POTFIL=potline)
        # create the input file
        self.input.write()
        # create the start pot file
        self.write_pot()

    @staticmethod
    def _skip_lines(fd, num):
        for ii in range(num):
            line = next(fd)
        return line

    @staticmethod
    def _skip_lines_to(fd, key):
        while 1:
            try:
                line = next(fd)

            except StopIteration:
                return ''

            if key in line:
                return line

    def read_output(self,filename):
        out = {
            'it' : [],
            'ERR' : [],
            'ETOT' : [],
            'EF' : [],
            'M' : [],
            'converged' : [],
            'atom_confs' : [],
        }
        atom_confs = {}
        with open(filename, 'r') as fd:
            for line in fd:
                if 'ERR' in line and 'EF' in line:
                    items = line.split()
                    out['it'].append(int(items[0]))
                    out['ERR'].append(float(items[2]))
                    out['EF'].append(float(items[5]))
                    out['M'].append((float(items[10]), float(items[11])))
                    items = self._skip_lines(fd, 1).split()
                    out['ETOT'].append(float(items[1]))
                    flag = items[5] == 'converged'
                    out['converged'].append(flag)

                    out['atom_confs'].append(atom_confs)
                    atom_confs = {}

                    self._skip_lines(fd, 1)

                elif 'SPRKKR-run for:' in line:
                    run = line.replace('SPRKKR-run for:', '').strip()
                    out['run'] = run

                elif ' E= ' in line:
                    atom = line.split()[-1]
                    akeys = self._skip_lines(fd, 1).split()
                    line = self._skip_lines_to(fd, 'sum').split()
                    avals = list(map(float, line[1:8])) + [float(line[9])]
                    line = self._skip_lines(fd, 1).split()
                    akeys.append(line[0])
                    avals.append(line[1])
                    atom_conf = dict(zip(akeys, avals))
                    atom_confs[atom] = atom_conf

        return out

    def read_results(self):
       outstrg=self.read_output(self.output)
       lastiter=len(outstrg['it'])
       self.niter=lastiter
       converged = outstrg['converged'][lastiter-1]
       if not converged:
           raise RuntimeError('SPRKKR did not converge! Check ' + self.output)

       self.results['raw_output'] = outstrg
       self.results['energy']=outstrg['ETOT'][lastiter-1]*Rydberg
#        if not converged:
#            raise RuntimeError('ELK did not converge! Check ' + self.out)
#        self.read_energy()
#        if self.parameters.get('tforce'):
#            self.read_forces()
#        self.width = self.read_electronic_temperature()
#        self.nbands = self.read_number_of_bands()
#        self.nelect = self.read_number_of_electrons()
#        self.niter = self.read_number_of_iterations()
#        self.magnetic_moment = self.read_magnetic_moment()

   

    def read(self):
        raise NotImplementedError
        
        
        

#    def scf(self):
#        input_filename = os.path.abspath(self.input.filename)
#        output_filename = input_filename.replace(".inp", ".out")
#        pot_filename = os.path.abspath(os.path.join(self.directory, self.prefix + ".pot"))
#        self.command = "mpirun.openmpi -np 4 kkrscfMPI  " + input_filename + " > " + output_filename
#        self.calculate(self.atoms, None, None)


        

