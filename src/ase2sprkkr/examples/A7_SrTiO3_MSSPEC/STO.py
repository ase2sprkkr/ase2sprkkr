"""
SCF calculations for  SrTiO3 in combination with msspec code to generate
photoelectron diffraction.

Call with ``msspec`` and/or ``sprkkr`` commandline argument.
"""

from ase.lattice.tetragonal import SimpleTetragonalFactory
from ase.visualize import view
from ase2sprkkr import SPRKKR
import numpy as np
import sys
import glob


def main(args):

    # Define a Perovskite Factory class
    class PerovskiteFactory(SimpleTetragonalFactory):
        bravais_basis = [[0, 0, 0.0], [0.5, 0.5, 0.5],[0., 0.5, 0.5], [0.5, 0.5, 0],
                         [0.5, 0.0, 0.5]]
        element_basis = (0, 1, 2, 2, 2)

    Perovskite = PerovskiteFactory()

    # Generate the base STO cell
    a0 = 3.905
    STO = Perovskite(('Sr', 'Ti', 'O'),
                     latticeconstant={'a': a0, 'c/a': 1.},
                     size=(1,1,1))
    if 'view' in args:
        view(STO)

    # ########## SPRKKR part
    if 'sprkkr' in args:
        calc = SPRKKR(atoms=STO,mpi=['mpirun','-np','16'])
        calc.input_parameters.set(NL=3)
        calc.input_parameters.SCF.MIX=0.01
        calc.input_parameters.SCF.NITER=1000
        calc.input_parameters.CONTROL.add('NOSYM', True)

        out_scf=calc.calculate()
        calc.input_parameters='PHAGEN'
        calc.calculate(potential=out_scf.potential_filename)

    # ######### MsSpec part
    if 'msspec' in args:

        from msspec.calculator import MSSPEC
        from msspec.utils import SPRKKRPotential, hemispherical_cluster, \
                                 get_atom_index

        ###########################################################################
        # Build the SrTiO3 cluster                                                #
        ###########################################################################
        # We start by loading the potential since it updates info in the STO cell
        pot = SPRKKRPotential(STO, "SrTiO3_scf.pot_new",
                              *glob.glob("*PHAGEN.pot"))

        # tag each atom in the STO object to easily set the emitter in the
        # hemispherical_cluster function
        [atom.set('tag', ('Sr', 'Ti', 'O').index(atom.symbol)) for atom in STO]

        # Create a hemispherical cluster centered on a Ti emitter from the
        # STO primitive cell used in the SPRKKR step.
        # For this example we use 4 planes and the Ti emitter (tag #1) is
        # located in the 2nd plane (numbering starts at 0)
        cluster = hemispherical_cluster(STO,
                                        planes=4,
                                        emitter_plane=1, emitter_tag=1)

        # The created cluster is centered on the emitter atom, so defining
        # the absorber attribute is straightforward:
        cluster.absorber = get_atom_index(cluster, 0, 0, 0)

        ###########################################################################
        # Set up the PhotoElectron Diffraction calculator                         #
        ###########################################################################
        # Create the calculator
        calc = MSSPEC(spectroscopy='PED', algorithm='expansion', folder='PED')

        # minimalistic set of parameter for XPD scan
        calc.set_atoms(cluster)
        calc.calculation_parameters.scattering_order = 2
        # We need to impose a maximum number of tl to use because there is still
        # a memory bug that I need to investigate. Anyway, usually 25 is not that
        # bad for the kind of atoms and the photon energy of lab sources...
        calc.tmatrix_parameters.lmax_mode = 'imposed'
        calc.tmatrix_parameters.lmaxt = 25

        data = None
        for source_energy in np.arange(500., 1501., 100.):
            calc.source_parameters.energy = source_energy
            # Run a polar scan with the default MuffinTin potential for Ti(2p3/2)
            calc.tmatrix_parameters.potential = 'muffin_tin'
            data = calc.get_theta_scan(level='2p3/2', data=data)

            # Now we use the previously generated SPRKKR potential and run the same
            # calculation
            calc.tmatrix_parameters.potential = pot
            data = calc.get_theta_scan(level='2p3/2', data=data)

            # To better see the differences, plot the normalized signal on the
            # same graph
            # add a new dataset for storing normalized values
            dset = data.add_dset(f"comparison at {source_energy} eV")
            # make a copy of previous scans values
            theta, muffintin_cs, sprkkr_cs = (data[-3].theta.copy(),
                                            data[-3].cross_section.copy(),
                                            data[-2].cross_section.copy())
            # divide by the max
            sprkkr_cs /= sprkkr_cs.max()
            muffintin_cs /= muffintin_cs.max()
            # add a new dataset with those values
            dset.add_columns(theta=theta, sprkkr=sprkkr_cs, muffintin=muffintin_cs)
            # add a view for this dataset
            dview = dset.add_view('Comparison',
                                title=('Comparison of XPD polar scans for '
                                       r'Ti(2p3/2) at $h\nu$ = {:.0f} eV'.format(
                                           source_energy)),
                                xlabel=r'$\Theta (\degree$)',
                                ylabel='Normalized Signal (a.u.)')
            dview.select('theta', 'sprkkr', legend='With SPRKKR potential')
            dview.select('theta', 'muffintin', legend='With internal MT potential')
            dview.set_plot_options(autoscale=True)

        # Pop up the final result
        data.view()

    if len(args) <= 1:
        print("Please specify either 'sprkkr', 'msspec' keywords or both "
              "of them on the command line")


# Just run the script only when directly called from command line

if __name__ == "__main__":
    main(sys.argv)
