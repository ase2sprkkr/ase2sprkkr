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
from .calcio import _skip_lines, _skip_lines_to, InputFile, PotFile

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
        self.atoms = None
        self.task = task
        if restart is None:
            if input_parameters.upper() == 'SCF':
                self.restart = False
            else:
                self.restart = True

        self.inpfile=os.path.join(self.prefix + ".inp")
        self.outfile= self.inpfile.replace(".inp", ".out")
        self.potfile=self.inpfile.replace(".inp", ".pot")
        self.sysfile=self.inpfile.replace(".inp", ".sys")

        LOGGER.debug(f'Calculation will run in  directory: {self.directory}')
        LOGGER.debug(f'Prefix is   : {self.prefix}')
        LOGGER.debug(f'INP FILE:{self.inpfile}')
        LOGGER.debug(f'POT FILE:{self.potfile}')
        LOGGER.debug(f'OUT FILE:{self.outfile}')
        LOGGER.debug(f'OUT FILE:{self.sysfile}')

        self.input = InputFile(filename=self.inpfile, task=self.task,directory=self.directory)

        if (self.input_parameters.upper()=='SCF'):
            self.command = "mpirun.openmpi -np 4 kkrscfMPI  " + self.inpfile + " > " + self.outfile
            self.command = "echo"
        elif (self.input_parameters.upper()=='PHAGEN'):
            self.command = "mpirun.openmpi -np 4 kkrscfMPI  " + self.inpfile + " > " + self.outfile
        elif (self.input_parameters.upper()=='GEN'):
            self.command = "mpirun.openmpi -np 4 kkrgenMPI  " + self.inpfile + " > " + self.outfile
        else:
            print("TASK {} not implemeted in ASE",format(self.task))
            raise NotImplementedError

    def set_potfile(self,filename):
        self.potfile=filename
        if not self.restart:
            self.pot=PotFile(self.atoms, filename=self.potfile,sysfilename=self.sysfile,directory=self.directory)
            self.pot.write()
            self.pot.write_sys()
        else:
            LOGGER.debug(f'POTENTIAL FILE WILL BE USED:{self.potfile}')


        LOGGER.debug(f'POT FILE:{self.potfile}')
    def set_inpfile(self,filename):
        self.inpfile=filename
        self.input= InputFile(filename=self.inpfile, task=self.task,directory=self.directory)
        LOGGER.debug(f'INP FILE:{self.inpfile}')
    def set_outfile(self,filename):
        self.outfile=filename
        LOGGER.debug(f'OUT FILE:{self.outfile}')
    def set_atoms(self, atoms):
        self.atoms = atoms

    def set_command(self,task):
        self.task=input_parameters.upper()
        print("TASK:",self.task)
        if (self.input_parameters.upper()=='SCF'):
            self.command = "mpirun.openmpi -np 4 kkrscfMPI  " + self.inpfile + " > " + self.outfile
        elif (self.input_parameters.upper()=='PHAGEN'):
            self.restart=True
            self.command = "kkrscf <  " + self.inpfile + " > " + self.outfile
        elif (self.input_parameters.upper()=='PHAGEN'):
            self.restart=True
            self.command = "kkrgen <  " + self.inpfile + " > " + self.outfile
        else:
            print("TASK {} not implemeted in ASE",format(self.task))
            raise NotImplementedError

    def write_input(self, atoms, properties=None, system_changes=None):
        # this will create directories
        FileIOCalculator.write_input(self, self.atoms)
        self.input.control_section.set(POTFIL=self.potfile)
        # create the input file
        self.input.write()
        # create the start pot file
        self.set_command(self.task)
        self.set_potfile(self.potfile)

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
                    items = _skip_lines(fd, 1).split()
                    out['ETOT'].append(float(items[1]))
                    flag = items[5] == 'converged'
                    out['converged'].append(flag)

                    out['atom_confs'].append(atom_confs)
                    atom_confs = {}

                    _skip_lines(fd, 1)

                elif 'SPRKKR-run for:' in line:
                    run = line.replace('SPRKKR-run for:', '').strip()
                    out['run'] = run

                elif ' E= ' in line:
                    atom = line.split()[-1]
                    akeys = _skip_lines(fd, 1).split()
                    line = _skip_lines_to(fd, 'sum').split()
                    avals = list(map(float, line[1:8])) + [float(line[9])]
                    line = _skip_lines_to(fd, 'E_band').split()
                    akeys.append(line[0])
                    if len(line) >= 2:
                       avals.append(line[1])
                    atom_conf = dict(zip(akeys, avals))
                    atom_confs[atom] = atom_conf

        return out

    def read_results(self):
       outstrg=self.read_output(os.path.join(self.directory,self.outfile))
       lastiter=len(outstrg['it'])
       self.niter=lastiter
       self.converged = outstrg['converged'][lastiter-1]
       if not self.converged:
           raise RuntimeError('SPRKKR did not converge! Check ' + self.outfile)

       self.results['raw_outfile'] = outstrg
       self.results['energy']=outstrg['ETOT'][lastiter-1]*Rydberg

    def phagen(self):
         self.calculate(self.atoms, None, None)

    def kkrgen(self):
         self.calculate(self.atoms, None, None)

    def kkrspec(self):
         self.calculate(self.atoms, None, None)

    def kkrchi(self):
         self.calculate(self.atoms, None, None)

    def read(self):
        raise NotImplementedError
