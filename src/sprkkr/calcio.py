# coding: utf-8
# vim: set et sw=4 ts=4 sts=4 nu ai cc=80 fdm=indent mouse=a:

"""
Module calcio
=============

"""

import numpy as np
from datetime import datetime
from collections import OrderedDict

from ase import Atom
from ase.units import Bohr
from ase.spacegroup import get_spacegroup

from .misc import LOGGER, get_occupancy


class Section(object):
    def __init__(self, name):
        self.name = name
        self.variables = []
        self.switches = []

    def add_variables(self, *variables):
        # dynamically add an attribute
        for variable in variables:
            self.variables.append(variable)

    def set(self, **kwargs):
        names = [_.name.upper() for _ in self.variables]
        for kw, val in kwargs.items():
            if kw.upper() not in names:
                self.set_custom(**{kw:val})
            else:
                for v in self.variables:
                    if kw.upper() == v.name.upper():
                        v.value = val

    def set_custom(self, **kwargs):
        LOGGER.info(f"Adding custom variables to the \"{self.name}\" section")
        for kw, val in kwargs.items():
            v = Variable(kw)
            v.value = val
            self.add_variables(v)

    def __str__(self):
        # for pretty printing, get the indentation
        indent = len(self.name)
        s = self.name.upper()
        for i, variable in enumerate(self.variables):
            if i > 0:
                s += " " * indent
            s += f" {variable}\n"
        return s
        

class Variable(object):
    def __init__(self, name, default=None, fmt="{:s}", dim=0, 
                 always_show=False):
        self.name = name
        self.fmt = fmt
        self.dim = dim
        self.default = default
        self.value = self.default
        self.always_show = always_show

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, val):
        self._value = val

    def __str__(self):
        if self.value is None:
            raise NameError(f"You must set a value for \"{self.name}\"")
        if self.value == self.default and not(self.always_show):
            return ""

        s   = val = ""
        key = self.name.upper()
        if self.dim > 0:
            try:
                val += ",".join([self.fmt.format(_) for _ in self.value])
            except:
                val += self.fmt.format(self.value)
            finally:
                s = key + "={" + val + "}"
        else:
            s = key + "=" + self.fmt.format(self.value)
        return s


class InputFile(object):
    def __init__(self, filename="test.inp"):
        self.filename = filename
        control_section = Section("CONTROL")
        control_section.add_variables(
            Variable("DATASET"),
            Variable("ADSI"),
            Variable("POTFIL"),
            Variable("PRINT", fmt="{:d}", default=0))

        sites_section = Section("SITES")
        sites_section.add_variables(
            Variable("NL", fmt="{:d}", default=3, dim=1, always_show=True))

        tau_section = Section("TAU")
        tau_section.add_variables(
            Variable("BZINT", default="POINTS", always_show=True),
            Variable("NKTAB", fmt="{:d}", default=300, always_show=True)
            )


        self.sections = OrderedDict({
                        "control": control_section, 
                         "tau"    : tau_section,
                         "sites"  : sites_section})
        for key, section in self.sections.items():
            setattr(self, key + "_section", section)

    def set(self, section_name, **variables):
        if not(hasattr(self, section_name + "_section")):
            LOGGER.info(f"Adding custom section \"{section_name}\"")
            setattr(self, section_name + "_section", Section(section_name))

        section = getattr(self, section_name + "_section")
        section.set(**variables)
        self.sections["section_name"] = section
        
    def write(self):
        with open(self.filename, "w") as fd:
            fd.write(str(self))

    def __str__(self):
        now = datetime.now()
        s  = "#"*80 + "\n"
        s += "#  SPR-KKR input file   {}\n".format(self.filename)
        s += "#  created by sprkkr Python package on {}\n".format(now.ctime())
        s += "#"*80 + "\n"
        s += "\n"

        for key, section in self.sections.items():
            s += str(section)
            s += "\n"

        return s





class LatticeData(object):
    sym = OrderedDict({
            'aP': 'C_i',
            'mP': 'C_2h', 'mS': 'C_2h',
            'oP': 'D_2h', 'oS': 'D_2h', 'oI': 'D_2h', 'oF': 'D_2h',
            'tP': 'D_4h', 'tI': 'D_4h',
            'hR': 'D_3d',
            'hP': 'D_6h',
            'cP': 'O_h', 'cI': 'O_h', 'cF': 'O_h'})

    def __init__(self, atoms):
        cell = atoms.get_cell()
        bl = cell.get_bravais_lattice()
        sg = get_spacegroup(atoms)
        ps = bl.pearson_symbol

        self.bravais_number = list(self.sym.keys()).index(ps) + 1
        self.schoenflies_symbol = self.sym[ps]
        self.basis = 0

        self.bravais = cell.get_bravais_lattice()
        self.alat = bl.a / Bohr
        self.boa = cell.cellpar()[1] / cell.cellpar()[0]
        self.coa = cell.cellpar()[2] / cell.cellpar()[0]
        self.rbas = sg.scaled_primitive_cell

    def __repr__(self):
        desc = self.bravais.longname
        s = "<{}({})>".format(self.__class__.__name__, desc)
        return s


class SiteData(object):
    SITES_LIST = []

    def __init__(self, coords, occupancy=[]):
        self.x, self.y, self.z = coords
        self.occupancy = occupancy
        self.noq = len(occupancy)
        self.qmtet = 0.
        self.qmphi = 0.
        self.__class__.SITES_LIST.append(self)
        self.id = len(self.__class__.SITES_LIST)

    @classmethod
    def get_sites(cls):
        return cls.SITES_LIST

    def __repr__(self):
        s = "<SiteData(id={:d}, noq={:d})>".format(self.id, self.noq)
        return s


class TypeData(object):
    TYPES_LIST = []

    def __init__(self, symbol, concentration=1., site=None, mesh=None):
        self.symbol = symbol
        self.atomic_number = Atom(symbol).number
        self.nat = 1
        self.concentration = concentration
        self.mesh = None
        self.site = site
        self.__class__.TYPES_LIST.append(self)
        self.id = len(self.__class__.TYPES_LIST)

    @classmethod
    def get_types(cls):
        return cls.TYPES_LIST

    def __repr__(self):
        s = "<TypeData(\'{}\', id={:d}, concentration={:.4f}, site={:d})>".format(
            self.symbol, self.id, self.concentration, self.site.id)
        return s


class MeshData(object):
    MESH_ID = 0

    def __init__(self):
        self.__class__.MESH_ID += 1
        self.id = self.MESH_ID
        self.type = ""

    def __repr__(self):
        s = "<{}(id={:d}, type=\"{}\")>".format(self.__class__.__name, self.id,
             self.type)
        return s


class ExponentialMeshData(MeshData):
    def __init__(self):
        super().__init__()
        self.type = "exponential"
        self.rws = 2.6
        self.rmt = 2.2
        self.jws = 721
        self.jmt = 0
        self.r_1 = 1e-6
        self.dx  = 2e-2


class PotFile(object):
    VERSION_STRING = "VERSION  3  (13.11.2000)"

    def __init__(self, atoms, filename=None, title=None, system=None):
        self.filename = filename
        self.atoms = atoms
        self.title = None

        self.system = system
        if system is None:
            self.system = atoms.get_chemical_formula()
        
        self.nm = 1
        self.nq = 1
        self.nt = 2
        self.ns = 2
        self.ne = 30
        self.irel = 2
        self.ibzint = 2
        self.nktab = 0

        self.info = None
        self.xc_pot = "VWM"
        self.scf_alg = "BROYDEN2"
        self.scf_iter = 0
        self.scf_mix = 2e-1
        self.scf_tol = 1e-5
        self.breitint = False
        self.orbpol = None
        self.extfield = False
        self.blcoupl = False
        self.bext = 0.
        self.kmrot = 0
        self.qmvec = np.array([0,0,0])
        self.ef = 0

        self.ld = LatticeData(atoms)
        self.md = ExponentialMeshData()
        self.all_at = self.get_atom_types()
        self.all_sd = SiteData.get_sites()

        # set the mesh to the atom types
        for at in self.all_at:
            at.mesh = self.md

    def get_atom_types(self):
        atoms = self.atoms
        sg = get_spacegroup(atoms)
        sites = sg.unique_sites(atoms.get_scaled_positions())
        occupancy = get_occupancy(atoms)
        all_types = []
        for i, sitexyz in enumerate(sites):
            # create a SiteData object
            sd = SiteData(sitexyz) 
            for symbol, concentration in occupancy[i].items():
                # create a TypeData object
                td = TypeData(symbol, concentration, site=sd)
                sd.occupancy.append(td)
                all_types.append(td)
        return all_types
        
    def write(self):
        with open(self.filename, "w") as fd:
            fd.write(str(self))


    def __str__(self):
        separator = "*"*80 + "\n"
        s  = 'SPRKKR    SCF-start data created by Python sprkkr package '
        s +=  datetime.now().ctime() + "\n"
        s += "{:<10}{}\n".format("FORMAT", self.VERSION_STRING)
        s += "{:<10s}{}\n".format("TITLE", self.title)
        s += "{:<10s}{}\n".format("SYSTEM", self.system)

        s += "{:<18}{:>2d}\n".format("NM", self.nm)
        s += "{:<18}{:>2d}\n".format("NQ", self.nq)
        s += "{:<18}{:>2d}\n".format("NT", self.nt)
        s += "{:<18}{:>2d}\n".format("NS", self.ns)
        s += "{:<18}{:>2d}\n".format("NE", self.ne)
        s += "{:<18}{:>2d}{:>8}\n".format("IREL", self.irel, "FREL")
        s += "{:<18}{:>2d}\n".format("IBZINT", self.ibzint)
        s += "{:<18}{:>2d}\n".format("NKTAB", self.nktab)

        s += "{:<10}{}\n".format("INFO", "NONE" if self.info is None else self.info)
        s += "{:<10}{}\n".format("XC-POT", self.xc_pot)
        s += "{:<10}{}\n".format("SCF-ALG", self.scf_alg)
        s += "{:<18}{:d}\n".format("SCF-ITER", self.scf_iter)
        s += "{:<18}{:1.10f}\n".format("SCF-MIX", self.scf_mix)
        s += "{:<18}{:1.10f}\n".format("SCF-TOL", self.scf_tol)

        s += "{:<12}{}\n".format("BREITINT", "F" if not self.breitint else "T")
        s += "{:<10}{}\n".format("ORBPOL", "NONE" if self.orbpol is None else
                                 self.orbpol)
        s += "{:<12}{}\n".format("EXTFIELD", "F" if not self.breitint else "T")
        s += "{:<12}{}\n".format("BLCOUPL", "F" if not self.breitint else "T")
        s += "{:<18}{:1.10f}\n".format("BEXT", self.bext)
        s += "{:<18}{:>2d}\n".format("KMROT", self.kmrot)
        s += ("{:<10}" + "{:>20.10f}"*3 + "\n").format("QMVEC", *self.qmvec)
        s += "{:<10}{:>20.10f}\n".format("EF", self.ef)


        desc = self.ld.bravais.longname
        pg = ""
        s += "{:<18}{:>2d}  {:>31}  {}  {}\n".format("BRAVAIS", 
             self.ld.bravais_number, desc, pg, self.ld.schoenflies_symbol)
        s += "{:<18}{:>2d}\n".format("BASIS", self.ld.basis)
        s += "{:<18}{:1.10f}\n".format("ALAT", self.ld.alat)
        s += "{:<18}{:1.10f}\n".format("BOA", self.ld.boa)
        s += "{:<18}{:1.10f}\n".format("COA", self.ld.coa)
        for i, row in enumerate(self.ld.rbas):
            s += ("R-BAS A({:d})" + "{:>20.10f}" * 3 + "\n").format(i+1, *row)

        s += separator

        s += "MESH INFORMATION\n"
        s += "MESH-TYPE {}\n".format(self.md.type.upper())
        s += "MESH{:>16d}\n".format(self.md.id)
        s += "{:<18s}{:1.10f}\n".format("RWS", self.md.rws)
        s += "{:<18s}{:1.10f}\n".format("RMT", self.md.rmt)
        s += "{}{:>17d}\n".format("JWS", self.md.jws)
        s += "{}{:>17d}\n".format("JMT", self.md.jmt)
        s += "{:<18s}{:1.10f}\n".format("R(1)", self.md.r_1)
        s += "{:<18s}{:1.10f}\n".format("DX", self.md.dx)

        s += separator
    
        s += "SITES{:>17d}\n".format(len(self.all_sd))
        for sd in self.all_sd:
            s += ("{:>4d} Q= ({:>3.10f},{:>3.10f},{:>3.10f})"
                  "   NOQ={:>3d}").format(sd.id, sd.x, sd.y, sd.z, sd.noq)
            s += "  IT,CONC="
            for at in sd.occupancy:
                s += "{:>3d} {:1.3f}".format(at.id, at.concentration)
            s += "\n"
        for sd in self.all_sd:
            s += "{:>4d} QMTET = {:>5.10f} QMPHI = {:5.10f}\n".format(sd.id, 
                  sd.qmtet, sd.qmphi)

        s += separator

        s += "MADELUNG MATRIX SMAD\n\n"

        s += separator

        for iat, at in enumerate(self.all_at):
            s += "TYPE{:>16d}\n".format(at.id)
            s += "{}\n".format(at.symbol)
            s += "Z{:>19d}\n".format(at.atomic_number)
            s += "NAT{:>17d}\n".format(at.nat)
            s += "CONC{:>26.10f}\n".format(at.concentration)
            s += "MESH{:>16d}\n".format(at.mesh.id)
            s += "SITES Q{:>13d}\n".format(at.site.id)
            if iat < len(self.all_at) - 1:
                s += separator

        return s
