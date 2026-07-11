"""Calculate and plot Ni L2,3 XAS and XMCD spectra of Ni2FeGa."""


def make_ni2fega():
    """Return the ordered primitive Ni2FeGa Heusler cell."""
    import numpy as np
    from ase2sprkkr import SPRKKRAtoms

    cell = 5.748 * np.array(
        [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
    )
    return SPRKKRAtoms(
        symbols=["Ni", "Ni", "Fe", "Ga"],
        scaled_positions=[
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.75],
            [0.50, 0.50, 0.50],
            [0.00, 0.00, 0.00],
        ],
        cell=cell,
        pbc=True,
    )


def main():
    from glob import glob
    from ase2sprkkr import OutputFile, SPRKKR

    # XAS is a post-SCF task, so first obtain a converged potential.  Keeping
    # this step here makes the example independent of every other example.
    calculator = SPRKKR(atoms=make_ni2fega())
    scf_result = calculator.calculate()

    # IT=1 selects Ni and CL='2p' selects the spin-orbit split L2,3 edge.
    xas_options = {
        "CONTROL.DATASET": "Ni2FeGa",
        "TAU.NKTAB": 500,
        "TASK.IT": 1,
        "TASK.CL": "2p",
        "TASK.FRAMETET": 0,
        "TASK.FRAMEPHI": 0,
        "ENERGY.NE": 180,
        "ENERGY.GRID": 6,
        "ENERGY.EMAX": 40.0,
    }
    calculator.calculate(
        potential=scf_result.potential_filename,
        input_parameters="xas",
        options=xas_options,
    )

    # SPR-KKR writes one RAT file for each selected absorbing edge.
    rat_files = glob("Ni2FeGa_XAS_*.rat")
    if len(rat_files) != 1:
        raise RuntimeError(f"Expected one XAS RAT file, found {rat_files}")
    OutputFile.from_file(rat_files[0]).plot()


if __name__ == "__main__":
    main()
