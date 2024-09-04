"""
Demonstration how to chage input SPRKKR file.
"""


def main():
    from ase.build import bulk
    from ase2sprkkr.sprkkr.calculator import SPRKKR

    atoms = bulk('Fe')

    calculator = SPRKKR(atoms=atoms,mpi=['mpirun','-np','4'])
    # print out input parameters. This will be used when calling
    # calculator to create corresponding input file. Please note
    # that potential filename will be set by method calculate
    print("ORIGINAL INPUT PARAMETERS==================")
    print(calculator.input_parameters.to_dict())
    print("===========================================")
    # There are several way to modify parameters
    # Lets modify l-expanstion (this will nodify first occurence of NL)
    calculator.input_parameters.set(NL=3)
    calculator.input_parameters.set(NE=32)
    # Lest also change mixing
    calculator.input_parameters.SCF.MIX=0.20
    calculator.input_parameters.ENERGY.ImE=0.0
    # calculator.input_parameters.ENERGY.GRID=[5,3]
    print("NEW INPUT PARAMETERS=======================")
    print(calculator.input_parameters.to_dict())
    print("===========================================")
    out=calculator.calculate()
    print(out.energy)
    print(len(out.iterations))
    print(out.iterations[-1]['error'])
    print(out.last_iteration['moment'])


# Just run the script only when directly called from command line

if __name__ == "__main__":
    main()
