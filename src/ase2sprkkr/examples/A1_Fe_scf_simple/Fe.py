"""
Simple SCF calculation for bcc Fe.
"""


def main():
    from ase.build import bulk
    from ase2sprkkr.sprkkr.calculator import SPRKKR

    # Create structure
    atoms = bulk('Fe')

    # choose sprkkr calculator
    calculator = SPRKKR(atoms=atoms,mpi=['mpirun','-np','4'])

    # perform scf calculations
    out=calculator.calculate()
    # out object includes results
    print(out.energy)
    print(len(out.iterations))
    print(out.iterations[-1]['error'])
    print(out.last_iteration['moment'])


# Just run the script only when directly called from command line

if __name__ == "__main__":
    main()
