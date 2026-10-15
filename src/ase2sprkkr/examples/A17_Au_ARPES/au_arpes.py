"""Calculate Au bulk DOS, angle-resolved photoemission and a Fermi map."""


def main():
    from ase.build import bulk
    from ase2sprkkr import SPRKKR

    # The SCF step produces the converged potential required by all subsequent
    # spectroscopy tasks in this self-contained workflow.
    gold = bulk("Au", "fcc", a=4.08)
    scf_options = {
        "CONTROL.KRMT": 4,
        "TAU.NKTAB3D": 1000,
        "SCF.VXC": "VWN",
        "SCF.FULLPOT": False,
        "SCF.NITER": 2000,
        "SCF.MIX": 0.1,
        "SCF.TOL": 1e-5,
    }
    calculator = SPRKKR(atoms=gold, options=scf_options)
    calculator.calculate()

    # A DOS calculation is a useful check of the converged bulk potential.
    dos_result = calculator.calculate(
        input_parameters="dos",
        options={
            "ENERGY.ImE": 0.001,
            "ENERGY.EMINEV": -12.0,
            "ENERGY.EMAXEV": 7.0,
            "ENERGY.NE": 1000,
            "TAU.NKTAB": 2500,
        },
    )
    dos_result.dos.plot()

    # One-step ARPES includes photon geometry, the surface orientation and the
    # detector scan in a single calculation.
    arpes_options = {
        "TAU.NKTAB": 500,
        "ENERGY.EMINEV": -4.0,
        "ENERGY.EMAXEV": 2.0,
        "ENERGY.EWORK_EV": 5.0,
        "ENERGY.NE": 300,
        "ENERGY.IMV_INI_EV": 0.02,
        "ENERGY.IMV_FIN_EV": 2.0,
        "TASK.MILLER_HKL": [1, 1, 1],
        "SPEC_PH.THETA": 45,
        "SPEC_PH.PHI": 90,
        "SPEC_PH.EPHOT": 50,
        "SPEC_STR.N_LAYDBL": [12, 12],
        "SPEC_STR.N_LAYER": 30,
        "SPEC_EL.THETA": [-20, 20],
        "SPEC_EL.NT": 200,
        "SPEC_EL.SPOL": 4,
    }
    arpes_result = calculator.calculate(
        input_parameters="arpes", options=arpes_options
    )
    arpes_result.spc.plot()

    # Reuse the same experimental geometry on a two-dimensional k grid at EF.
    fermi_options = {
        **arpes_options,
        "ENERGY.EMINEV": 0.0,
        "ENERGY.EMAXEV": 0.0,
        "ENERGY.NE": 1,
        "SPEC_EL.KA": [-0.5, -0.5],
        "SPEC_EL.K1": [1, 0],
        "SPEC_EL.K2": [0, 1],
        "SPEC_EL.NK1": 120,
        "SPEC_EL.NK2": 120,
        "SPEC_EL.PSPIN": [0, 0, 1],
    }
    fermi_result = calculator.calculate(
        input_parameters="arpes", options=fermi_options
    )
    fermi_result.spc.plot()


if __name__ == "__main__":
    main()
