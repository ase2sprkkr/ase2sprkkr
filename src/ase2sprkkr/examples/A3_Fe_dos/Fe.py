"""
SCF calculation for Fe using kkrscf and consequent DOS by kkrgen.
"""


def main():
    from ase.build import bulk
    from ase2sprkkr.sprkkr.calculator import SPRKKR

    atoms = bulk('Fe')

    print("FIST STEP: SELF CONSISTENT CALCULATIONS====")
    calculator = SPRKKR(atoms=atoms,mpi=True)
    calculator.input_parameters.set(NL=3)
    calculator.input_parameters.SCF.MIX=0.20
    calculator.input_parameters.ENERGY.ImE=0.0
    calculator.input_parameters.ENERGY.GRID=[5,3]
    calculator.input_parameters.set(NE=32)
    out=calculator.calculate()

    print(out.energy)
    print(len(out.iterations))
    print(out.iterations[-1]['error'])
    print(out.last_iteration['moment'])
    print("SECOND STEP: CALCULATION OF DOS============")
    # Lets now calculate DOS
    # First we need to change task (there are several input data tabulated for
    # various tasks to help user.

    calculator = out.calculator
    calculator.input_parameters='bsf'
    calculator.input_parameters.ENERGY.EMAX=0.2
    calculator.calculate()
    calculator.input_parameters='DOS'

    print("NEW INPUT PARAMETERS=======================")
    print(calculator.input_parameters.to_dict())
    print("===========================================")

    calculator.input_parameters.set(NE=300)
    calculator.input_parameters.set(NL=3)

    # Pass a newly converged potential to the DOS calculation
    calculator.calculate(potential=out.potential_filename)

    # For the processing of the results of the DOS task use xband
    calculator.run_xband()


# Just run the script only when directly called from command line

if __name__ == "__main__":
    main()
