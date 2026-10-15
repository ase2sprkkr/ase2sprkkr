"""Add empty spheres to open diamond Si for an atomic-sphere calculation."""


def main():
    from ase.build import bulk
    from ase2sprkkr import SPRKKR

    # Diamond Si has large interstitial voids.  It is therefore a natural
    # example of a structure whose ASA basis benefits from empty spheres.
    silicon = bulk("Si", "diamond", a=5.43)
    calculator = SPRKKR(atoms=silicon)

    # The simplest form asks ASE2SPRKKR to choose all settings automatically:
    # calculator.save_input(potential_file="Si.pot", empty_spheres=True)

    # A two-item tuple limits the permitted empty-sphere radii:
    # calculator.save_input(
    #     potential_file="Si.pot",
    #     empty_spheres=(0.2, 2.0),
    # )

    # A dictionary exposes the complete set of useful tuning parameters.  This
    # example only writes the input, so no SPR-KKR executable is launched.
    calculator.save_input(
        potential_file="Si.pot",
        empty_spheres={
            "min_radius": 0.2,
            "max_radius": 2.0,
            "max_spheres": 300,
            "verbose": True,
        },
    )


if __name__ == "__main__":
    main()
