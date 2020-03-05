# coding: utf-8
# vim: set et sw=4 ts=4 sts=4 nu ai cc=80 fdm=indent mouse=a:

"""
Module calcio
=============

"""

import numpy as np
from datetime import datetime
from collections import OrderedDict
from ast import literal_eval

from ase import Atom
from ase.units import Bohr
from ase.spacegroup import get_spacegroup

from .misc import LOGGER, get_occupancy
from .data import get_sprkkr_input

_tr = str.maketrans('{}', '[]')

def _parse_value(val):
    aux = val.strip().translate(_tr)
    try:
        out = literal_eval(aux)

    except (ValueError, SyntaxError):
        aux2 = aux.split()
        if len(aux2) == 2: # Treat value with a unit.
            try:
                out = literal_eval(aux2[0])

            except (ValueError, SyntaxError):
                out = aux

        else:
            out = aux

    return out

def load_parameters(lines):
    """
    Load parameters.

    Notes
    -----
    The following assumptions should hold:
    - The section name is the first word on a line.
    - The sections are separated by empty lines or comments.
    - Each entry is on a separate line, with the exception of the first entry
      following the section name.

    Returns
    -------
    pars : dict
        The parameters as {section_name : section_dict}. When a parameter takes
        no value, the returned value is '@defined'.
    """
#    with open(filename) as fd:
#        lines = fd.readlines()

    # Replace comments.
    lines = [line.strip() if not line.startswith('#') else ''
             for line in lines]
    # Split sections.
    sections = '\n'.join(lines).split('\n\n')

    # Parse sections.
    pars = OrderedDict()
    for section in sections:
        section = section.strip()
        if not len(section): continue

        # Get section name.
        items = section.split('\n')
        name = items[0].split()[0]
        items[0] = items[0].replace(name, '')
        name = name.lower()

        # Parse section items.
        vals = OrderedDict()
        for item in items:
            aux = item.split('=')
            key = aux[0].strip()
            if len(aux) == 2: # item is key = value
                val = _parse_value(aux[1])

            else: # item is flag
                val = '@defined'

            vals[key] = val

        pars[name] = vals

    return pars

def make_sections(pars):
    sections = OrderedDict()
    for name, vals in pars.items():
        section = Section(name)

        for key, item in vals.items():
            dim = 0
            if isinstance(item, int):
                fmt = '{:d}'

            elif isinstance(item, float):
                fmt = '{:f}'

            elif isinstance(item, list):
                fmt = '{}'
                dim = len(item)

            else:
                fmt = '{:s}'

            var = Variable(key, default=item, fmt=fmt, dim=dim,
                           always_show=True)
            section.add_variables(var)

        sections[name] = section

    return sections

def compare_sections(s1, s2):
    is_eq = {name : (section == s2[name]) for name, section in s1.items()}
    return is_eq

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

    def __eq__(self, other):
        if len(self.variables) != len(other.variables):
            return False

        if len(self.switches) != len(other.switches):
            return False

        for ii, var in enumerate(self.variables):
            if var != other.variables[ii]:
                return False

        return True

class Variable(object):
    def __init__(self, name, default=None, fmt="{:s}", dim=0,
                 always_show=True):
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

    def __eq__(self, other):
        is_eq = ((self.name == other.name)
                 and (self.value == other.value)
                 and (self.dim == other.dim))
        return is_eq

class InputFile(object):
    def __init__(self, filename="test.inp", defaults_filename=None,task='scf'):
        self.filename = filename
        self.set_defaults(filename=defaults_filename)

    def set_defaults(self, sections=None, filename=None):
        if sections is not None:
            self.sections = sections

        elif filename is not None:
            pars = load_parameters(filename)
            self.sections = make_sections(pars)

        else: # Empty sections - to remove?
            control_section = Section("CONTROL")
            sites_section = Section("SITES")
            tau_section = Section("TAU")
            scf_section=Section("SCF")

            self.sections = OrderedDict({
                "control": control_section,
                "tau"    : tau_section,
                "scf"    : scf_section,
                "sites"  : sites_section,
            })

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
        self.irefq=1
        self.imq=1
        self.qmtet = 0.
        self.qmphi = 0.
        self.__class__.SITES_LIST.append(self)
        self.id = len(self.__class__.SITES_LIST)
        self.itoq=0


    @classmethod
    def get_sites(cls):
        return cls.SITES_LIST

    def __repr__(self):
        s = "<SiteData(id={:d}, noq={:d})>".format(self.id, self.noq)
        return s


class TypeData(object):
    TYPES_LIST = []

    def __init__(self, symbol, concentration=1., iq=None,site=None, mesh=None,iqat=None):
        self.symbol = symbol
        self.atomic_number = Atom(symbol).number
        self.nat = 1
        self.concentration = concentration
        self.mesh = None
        self.site = site
        self.iqat = iqat
        self.__class__.TYPES_LIST.append(self)
        self.id = iq 
        #len(self.__class__.TYPES_LIST)

    @classmethod
    def get_types(cls):
        return cls.TYPES_LIST

    def __repr__(self):
        s = "<TypeData(\'{}\', id={:d}, concentration={:.4f}, site={:d},iqat={})>".format(
            self.symbol, self.id, self.concentration, self.site.id,self.iqat)
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
    VERSION_STRING = "6  (21.05.2007)"

    def __init__(self, atoms, filename=None, title=None, system=None):
        self.filename = filename
        self.atoms = atoms
        self.title = None
        self.system = system
        if system is None:
            self.system = atoms.get_chemical_formula()
        self.package='SPRKKR'
        
        self.nm = 1
        self.nq = 1
        self.nt = 2
        self.ns = 2
        self.ne = 30
        self.irel = 3
        self.ibzint = 2
        self.nktab = 0
        self.scf_status='START'
        self.info = None
        self.xc_pot = "VWM"
        self.scf_alg = "BROYDEN2"
        self.fullpot=False
        self.nonmag=False
        self.semicore=False
        self.lloyd='F'
        self.scf_iter = 0
        self.scf_mix = 2e-1
        self.scf_tol = 1e-5
        self.breitint = False
        self.orbpol = None
        self.extfield = False
        self.blcoupl = False
        self.bext = 0.
        self.kmrot = 0
        self.rmsavv=999999.
        self.rmsavb=999999.
        
        self.qmvec = np.array([0,0,0])
        self.ef = 999999.
        self.vmtz=0.7
        self.cartesian=True
        self.basscale = np.array([1.0,1.0,1.0])
        
        self.vref=4.0000
        self.rmtref=0.00000

        self.tabncore=[0,0,0,2,2,2,2,2,2,2,2,10,10,10,10,10,10,10,10,18,18,
                       18,18,18,18,18,18,18,18,18,28,28,28,28,28,28,28,36,
                       36,36,36,36,36,36,36,36,36,36,46,46,46,46,46,46,46,
                       54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,
                       68,68,68,68,68,68,68,68,78,78,78,78,78,78,78,86,86,
                       86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,
                       86,86,86,86,86]

        # create lattice, mesh, site and types from ase.atoms
        self.ld = LatticeData(atoms)
        self.all_at = self.get_atom_types()
        self.all_sd = SiteData.get_sites()


        # set the mesh to the atom types
        self.md = ExponentialMeshData()
        for at in self.all_at:
            at.mesh = self.md

    def get_atom_types(self):
        atoms = self.atoms
#Use symmetry in order to find equivalet sites 
        sg = get_spacegroup(atoms)
        unique_sites=sg.unique_sites(atoms.get_scaled_positions())
        all_sites = atoms.get_scaled_positions()
        occupancy = get_occupancy(atoms)
        type_table=np.zeros(shape=(len(unique_sites),len(all_sites)),dtype=int)
        tol=1e-5
        iqat={}
        iqref=[]
        itoq={}
        for i in range(len(all_sites)):
            iqref.append(-1)
        if len(unique_sites)!=len(all_sites):
            LOGGER.info(f"==========================================")
            LOGGER.info(f"Symmetry adaped list of types will be used")
            LOGGER.info(f"==========================================")

            for ui in range(len(unique_sites)):
                su,k=sg.equivalent_sites(unique_sites[ui])
                lst=[]
                for i in range(len(su)):    
                    for j in range(len(all_sites)):
                        diff=np.abs(su[i]-all_sites[j])<=[tol,tol,tol]
                        if diff.all():
                            type_table[ui][j]=1
                            lst.append(j)
                            iqref[j]=ui
                            if occupancy[j]!=occupancy[ui]:
                                LOGGER.warning(f"\n OCCUPANCY NOT CONFORM WITH THE CRYSTAL SYMMETRY \n")
                                LOGGER.warning(f"Occupation of SITE {j:d} will be taken from SITE {ui:d} \n")
                                occupancy[j]=occupancy[ui]           
                iqat[ui]=lst
            LOGGER.info(f"Table of equivalent types \n \"{type_table}\"")
        else:
            for j in range(len(all_sites)):
                type_table[j][j]=1
                iqat[j]=[j]
                iqref[j]=j
#       
        all_types = []
        dsf={}
        it=0

        for ui in range(len(unique_sites)):
            lsf=[]
            for symbol, concentration in occupancy[ui].items():
                lsf.append(it)
                it=it+1
            dsf[ui]=lsf
        for i, sitexyz in enumerate(all_sites):
            ui=iqref[i]
            sd = SiteData(sitexyz)
            sd.noq=len(occupancy[ui])
            sd.itoq=ui
            ll=0
            lsf=dsf[ui]
            for symbol, concentration in occupancy[ui].items():
                td = TypeData(symbol, concentration,iq=i, site=lsf[ll],iqat=iqat[ui])
                sd.occupancy.append(td)
                all_types.append(td)
                ll=ll+1
        key_max = max(dsf.keys(), key=(lambda k: dsf[k]))
        self.nq=len(all_sites)
        self.nt=max(dsf[key_max])+1
 
        return all_types
        
    def write(self):
        with open(self.filename, "w") as fd:
            fd.write(str(self))


    def __str__(self):
        separator = "*"*80 + "\n"
        #**************************************************
        s  = separator
        s += 'HEADER    SCF-start data created by ase2sprkkr'
        s +=  datetime.now().ctime() + "\n"
        #**************************************************
        s += separator
        #**************************************************
        s += "{:<10s}{}\n".format("TITLE", self.title)
        s += "{:<10s}{}\n".format("SYSTEM", self.system)
        s += "{:<10s}{}\n".format("PACKAGE", self.package)
        s += "{:<10}{}\n".format("FORMAT", self.VERSION_STRING)
        #**************************************************
        s += separator
        #**************************************************
        s += "{}\n".format("GLOBAL SYSTEM PARAMETER")
        s += "{:<18}{:>2d}\n".format("NQ", self.nq)
        s += "{:<18}{:>2d}\n".format("NT", self.nt)
        s += "{:<18}{:>2d}\n".format("NM", self.nm)
        s += "{:<18}{:>2d}\n".format("IREL", self.irel)
        #**************************************************
        s += separator
        #**************************************************
        s += "{}\n".format("SCF-INFO")
        s += "{:<10}{}\n".format("INFO", "NONE" if self.info is None else self.info)
        s += "{:<10}{}\n".format("SCFSTATUS",self.scf_status)  
        s += "{:<10}{}\n".format("FULLPOT", "F" if not self.fullpot else "T")
        s += "{:<12}{}\n".format("BREITINT", "F" if not self.breitint else "T")  
        s += "{:<10}{}\n".format("NONMAG","F" if not self.nonmag else "T")
        s += "{:<10}{}\n".format("ORBPOL", "NONE" if self.orbpol is None else
                                self.orbpol)
        s += "{:<12}{}\n".format("EXTFIELD", "F" if not self.extfield else "T")
        s += "{:<12}{}\n".format("BLCOUPL", "F" if not self.blcoupl else "T")
        s += "{:<18}{:1.10f}\n".format("BEXT", self.bext)  
        s += "{:<10}{}\n".format("SEMICORE",self.semicore)
        s += "{:<10}{}\n".format("LLOYD",self.lloyd)
        s += "{:<18}{:>2d}\n".format("NE", self.ne)  
        s += "{:<18}{:>2d}\n".format("IBZINT", self.ibzint)
        s += "{:<18}{:>2d}\n".format("NKTAB", self.nktab)
        s += "{:<10}{}\n".format("XC-POT", self.xc_pot)
        s += "{:<10}{}\n".format("SCF-ALG", self.scf_alg)
        s += "{:<18}{:d}\n".format("SCF-ITER", self.scf_iter)
        s += "{:<18}{:1.10f}\n".format("SCF-MIX", self.scf_mix)
        s += "{:<18}{:1.10f}\n".format("SCF-TOL", self.scf_tol)
        s += "{:<18}{:1.10f}\n".format("RMSAVV", self.rmsavv)
        s += "{:<18}{:1.10f}\n".format("RMSAVB", self.rmsavb)
        s += "{:<10}{:>20.10f}\n".format("EF", self.ef)
        s += "{:<10}{:1.10f}\n".format("VMTZ", self.vmtz)
        s += separator
        #**************************************************
        s += "{}\n".format("LATTICE")    
        s += "{}\n".format("SYSDIM    3D")
        s += "{}\n".format("SYSTYPE   BULK")
        desc = self.ld.bravais.longname
        pg = ""
        s += "{:<18}{:>2d}  {:>31}  {}  {}\n".format("BRAVAIS", 
             self.ld.bravais_number, desc, pg, self.ld.schoenflies_symbol)
        s += "{:<12}{:1.10f}\n".format("ALAT", self.ld.alat)
        for i, row in enumerate(self.ld.rbas):
            s += ("A({:d})" + "{:>20.10f}" * 3 + "\n").format(i+1, *row)       
        s += separator
        #**************************************************
        s += "{}\n".format("SITES")       
        s += "{:<12}{}\n".format("CARTESIAN", "F" if not self.cartesian else "T")
        s += ("{:<10}" + "{:>15.10f}"*3 + "\n").format("BASSCALE", *self.basscale)
        s += "{}\n".format("        IQ      QX              QY              QZ")
        for sd in self.all_sd:
            s += ("{:>10d}"+"{:>15.10f}"*3+"\n").format(sd.id, sd.x, sd.y, sd.z)
        s += separator
        #**************************************************
        s += "{}\n".format("OCCUPATION")
        s += "{}\n".format("        IQ     IREFQ       IMQ       NOQ  ITOQ  CONC")
        for sd in self.all_sd:
            s += ("{:>10d}"+
                  "{:>10d}"*3).format(sd.id, sd.irefq,sd.imq,sd.noq)
            for at in sd.occupancy:
                if at.id+1 == sd.id:
                    s += "{:>3d} {:1.3f}".format(at.site+1,at.concentration)
            s += "\n"
        s += separator
        #**************************************************
        s += "{}\n".format("REFERENCE SYSTEM")    
        s += "{:<18}{:>2d}\n".format("NREF", self.nm)
        s += "{}\n".format("      IREF      VREF            RMTREF")
        for iref in range(self.nm):
            s += ("{:>10d}"+"{:>15.10f}"*2+"\n").format(iref+1, self.vref, self.rmtref)
        s += separator
        #**************************************************
        s += "{}\n".format("MAGNETISATION DIRECTION")       
        s += "{:<18}{:>2d}\n".format("KMROT", self.kmrot)
        s += ("{:<10}" + "{:>15.10f}"*3 + "\n").format("QMVEC", *self.qmvec)
        s += "{}\n".format("         IQ      QMTET           QMPHI")
        for sd in self.all_sd:
            s += "{:>10d} {:>15.10f} {:15.10f}\n".format(sd.id, 
                  sd.qmtet, sd.qmphi)
        s += separator
        #**************************************************
        s += "{}\n".format("MESH INFORMATION")       
        s += "MESH-TYPE {}\n".format(self.md.type.upper())
        s += "{}\n".format("   IM      R(1)            DX         JRMT      RMT        JRWS      RWS")
        s += "{:>5d}".format(self.md.id)
        s += "{:>18.10f}".format(self.md.r_1)
        s += "{:>18.10f}".format(self.md.dx)
        s += "{:>5d}".format(self.md.jmt)
        s += "{:>18.10f}".format(self.md.rmt)
        s += "{:>5d}".format(self.md.jws)
        s += "{:>18.10f}\n".format(self.md.rws)
        s += separator
        #**************************************************
        s += "{}\n".format("TYPES")
        s += "{}\n".format("   IT     TXTT        ZT     NCORT     NVALT    NSEMCORSHLT")
       
        itdone= [False for i in range(self.nt)]
        
        
        for iat, at in enumerate(self.all_at):
            if (not itdone[at.site]):
                s += "{:>5d}".format(at.site+1)
                s += "{:>9}".format(at.symbol)
                s += "{:>10d}".format(at.atomic_number)
                s += "{:>10d}".format(self.tabncore[at.atomic_number])
                s += "{:>10d}".format(at.atomic_number-self.tabncore[at.atomic_number])
                s += "{:>10d}\n".format(0)
                itdone[at.site]=True
                
        s += separator

        return s
