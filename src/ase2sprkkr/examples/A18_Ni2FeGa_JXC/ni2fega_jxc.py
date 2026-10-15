"""Calculate Ni2FeGa exchange couplings and export them for UppASD."""


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
    from ase2sprkkr import SPRKKR

    # Exchange parameters require a converged magnetic potential.  Performing
    # SCF here keeps the example independent of the CPA and XAS examples.
    calculator = SPRKKR(atoms=make_ni2fega())
    scf_result = calculator.calculate()

    jxc_options = {
        "CONTROL.DATASET": "Ni2FeGa",
        "TAU.NKTAB": 1000,
        "MODE.LLOYD": False,
        "TASK.CLURAD": 4.5,
        "SITES.NL": 4,
        "MODE": "SP-SREL",
        "ENERGY.NE": 32,
    }
    result = calculator.calculate(
        potential=scf_result.potential_filename,
        input_parameters="jxc",
        options=jxc_options,
    )

    # Inspect all three exchange-tensor components, report SPR-KKR's mean-field
    # Curie temperature and create jfile/posfile/momfile/inpsd for UppASD.
    result.jxc.plot(layout=(1, 3))
    print("Mean-field Curie temperature:", result.mean_field_curie_temperature)
    result.write_uppasd_files()


if __name__ == "__main__":
    main()
