""" Thys module contain functions to write xband sysfiles.

TODO: Just a-bit-cleaned-up old legacy code. Please, use it with a propper caution.
"""

from .sprkkr_atoms import SPRKKRAtoms
from ..physics.lattice_data import LatticeData
from ..common.unique_values import UniqueValuesMapping


def sysfile_content(atoms, filename='<unknown>'):
    """ Return the content of the xband sysfile for the given atoms object.

    TODO
    ----
    Not thoroughly tested nor guaranted to be fully correct
    """
    atoms = SPRKKRAtoms.promote_ase_atoms(atoms)

    class SiteData(object):

        def __init__(self, iq, icl, site, coors):
            self.iq =iq     # id of the site
            self.icl = icl  # id of the unique site
            self.site = site
            self.coors = coors
            self.occupancy = []

        def __repr__(self):
            s = "<SiteData(iq={:d}, noq={:d})>".format(self.iq, len(self.site.occupation))
            return s

    class TypeData(object):

        def __init__(self, symbol, atomic_number, concentration=1., it=None, mesh=None, iqs=None):
            self.symbol = symbol
            self.atomic_number = atomic_number
            self.concentration = concentration
            self.it = it
            self.iqs = iqs

        def __repr__(self):
            s = "<TypeData(\'{}\', id={:d}, concentration={:.4f}, site={:d},equivalent_sites={})>".format(
                self.symbol, self.id, self.concentration, self.site.id,self.equivalent_sites)
            return s

    # create lattice site and types from ase.atoms
    ld = LatticeData(atoms)
    atomic_types = []     # types
    sites = []            # sites

    # Use symmetry in order to find equivalet sites
    uv = UniqueValuesMapping.from_values([ i.site_type for i in atoms.sites])
    equivalent_sites=uv.indexes(start_from = 1)   # mapping unique_sites => sites
    iqref=uv.mapping  # mapping sites => unique_sites

    for i, (site, coors) in enumerate(zip(atoms.sites, atoms.positions)):
        sd = SiteData(i + 1, iqref[i], site, coors)
        sites.append(sd)

    it = 0
    for site, icl in uv.value_to_class_id.items():
        for symbol, concentration in site.occupation.items():
            it+=1
            iqs = equivalent_sites[icl]
            td = TypeData(str(symbol), symbol.atomic_number, concentration,it=it, iqs=iqs)
            atomic_types.append(td)
            for iq in equivalent_sites[icl]:
              sites[iq - 1].occupancy.append(it)

    filestring = "system data-file created by python ase2sprkkr \n"
    filestring += filename + "\n"
    filestring += "xband-version\n"
    filestring += "5.0\n"
    # It would be really cool to support lower dimensions...one day.
    filestring += "dimension\n"
    filestring += "3D\n"
    filestring += "Bravais lattice\n"
    filestring += " ".join(map(str, ld.pearson.xband_data() )) + "\n"
    filestring += "space group number (ITXC and AP)\n"
    filestring += "%5i%5i" % (ld.sgno,ld.apno) + "\n"
    filestring += "structure type\n"
    filestring += "UNKNOWN\n"
    filestring += "lattice parameter A  [a.u.]\n"
    filestring += "%18.12f\n" % ld.alat
    filestring += "ratio of lattice parameters  b/a  c/a\n"
    filestring += "%18.12f%18.12f\n" % (ld.boa,ld.coa)
    filestring += "lattice parameters  a b c  [a.u.]\n"
    filestring += "%18.12f%18.12f%18.12f\n" % (ld.alat,ld.blat,ld.clat)
    filestring += "lattice angles  alpha beta gamma  [deg]\n"
    filestring += "%18.12f%18.12f%18.12f\n" % (ld.alpha,ld.beta,ld.gamma)
    filestring += "primitive vectors     (cart. coord.) [A]\n"
    for vec in ld.rbas:
        for p in vec:
            filestring += "%18.12f" % p
        filestring += "\n"
    # Get number of sites and fill out with empty spheres if the sites are not fully filled
    filestring += "number of sites NQ\n"
    filestring += "%3i\n" % (len(sites))
    filestring += " IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ\n"

    rws=0.0
    angmom=4
    for sd in sites:
        filestring += "%3i%4i%18.12f%18.12f%18.12f  %18.12f%4i%5i " % (sd.iq,sd.icl,sd.coors[0],sd.coors[1],sd.coors[2],rws,angmom,len(sd.occupancy))
        for at in sd.occupancy:
           filestring += "%3i" % (at)
        filestring+="\n"
    filestring+="number of sites classes NCL \n"
    filestring += "%3i\n" % ( len(equivalent_sites) )
    filestring+="ICL WYCK NQCL IQECL (equivalent sites)\n"
    for icl, iqs in equivalent_sites.items():
        filestring += "%3i   %1s%5i" % (icl,'-',len(iqs))
        for key in iqs:
            filestring += "%3i" % (key)
        filestring+="\n"

    filestring += "number of atom types NT\n"
    filestring += "%3i\n" % len(atomic_types)

    filestring += " IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"
    for iat, at in enumerate(atomic_types):
        filestring += " %2i%4i  %8s%5i%6.3f" % (at.it,at.atomic_number,at.symbol,len(at.iqs),at.concentration)
        for key in at.iqs:
            filestring += "%3i" % (key)
        filestring+="\n"

# Average Wigner-Seitz radi
# =============================================================================
#         rws = 0.
#         iq = 0
#         icl = 0
#         itoq = 0
#         for a in cell.atomdata:
#             icl += 1
#             itoqs = []
#             for sp in a[0].species:
#                 itoq += 1
#                 itoqs.append(itoq)
#             for b in a:
#                 iq += 1
#                 if minangmom:
#                     angmom = max(max([ed.angularmomentum[ed.elementblock[spcs]] for spcs in b.species])+1,minangmom)
#                 else:
#                     angmom = max([ed.angularmomentum[ed.elementblock[spcs]] for spcs in b.species])+1
#                 v = mvmult3(cell.latticevectors,b.position)
#                 filestring += "%3i%4i%18.12f%18.12f%18.12f  %18.12f%4i%5i "%(iq,icl,v[0],v[1],v[2],rws,angmom,len(a[0].species))
#                 for i in itoqs:
#                     filestring += "%3i"%i
#                 filestring += "\n"
#         filestring += "number of sites classes NCL\n"
#         filestring += "%3i\n"%len(len(atomic_types))
#         filestring += "ICL WYCK NQCL IQECL (equivalent sites)\n"
#         iq = 0
#         icl = 0
#         for a in cell.atomdata:
#             icl += 1
#             filestring += "%3i   %1s%5i"%(icl,'-',len(a))
#             for b in a:
#                 iq += 1
#                 filestring += "%3i"%iq
#             filestring += "\n"
#         filestring += "number of atom types NT\n"
#         len(atomic_types) = 0
#         for a in cell.atomdata:
#             len(atomic_types) += len(a[0].species)
#         filestring += "%3i\n"%len(atomic_types)
#         filestring += " IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"
#         iq = 0
#         it = 0
#         for a in cell.atomdata:
#             corr = 0
#             for sp,conc in a[0].species.iteritems():
#                 it += 1
#                 filestring += " %2i%4i  %8s%5i%6.3f"%(it,ed.elementnr[sp],sp,len(a),conc)
#                 iq -= corr*len(a)
#                 for b in a:
#                     iq += 1
#                     filestring += "%3i"%iq
#                 corr = 1
#                 filestring += "\n"
# =============================================================================
    return filestring


def write_sysfile(atoms, filename:str):
    with open(filename, "w") as fd:
       fd.write(str(sysfile_content(atoms, filename)))
