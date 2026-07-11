"""Scan Fe/Co substitution in Ni2FeGa with the coherent potential approximation."""


def make_ni2fega():
    """Build the ordered primitive cell used as the scan template."""
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
    import numpy as np
    from ase2sprkkr import SPRKKR

    template = make_ni2fega()
    results = []

    # CPA represents substitutional disorder in the primitive cell, avoiding
    # a separate large supercell for every concentration.
    for fe_fraction in np.linspace(0.0, 1.0, 21):
        atoms = template.copy()
        atoms.sites[2].site_type.occupation.set(
            {"Fe": fe_fraction, "Co": 1.0 - fe_fraction}
        )
        result = SPRKKR(atoms=atoms).calculate()
        results.append((fe_fraction, result.energy))

    for fe_fraction, energy in results:
        print(f"Fe fraction={fe_fraction:5.2f}, total energy={energy}")


if __name__ == "__main__":
    main()
