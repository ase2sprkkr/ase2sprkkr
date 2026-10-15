"""Relax MoSe2 with GPAW and reuse the ASE structure for SPR-KKR ARPES."""


def main():
    from ase.filters import FrechetCellFilter
    from ase.optimize import BFGS
    from ase.spacegroup import crystal
    from ase2sprkkr import SPRKKR
    from gpaw import FermiDirac, GPAW, PW

    # GPAW is an optional dependency.  The relaxed object remains an ordinary
    # ASE Atoms instance and can therefore be passed directly to ASE2SPRKKR.
    atoms = crystal(
        symbols=["Mo", "Se"],
        basis=[(1 / 3, 2 / 3, 1 / 4), (1 / 3, 2 / 3, 0.621)],
        spacegroup=194,
        cellpar=[3.288, 3.288, 12.900, 90, 90, 120],
    )
    atoms.calc = GPAW(
        mode=PW(650),
        xc="PBE",
        kpts={"size": (8, 8, 4), "gamma": True},
        occupations=FermiDirac(0.05),
        txt="mose2_bulk_relax.txt",
    )
    optimizer = BFGS(
        FrechetCellFilter(atoms),
        trajectory="mose2.traj",
        logfile="mose2.log",
    )
    optimizer.run(fmax=0.02)
    atoms.write("mose2_relaxed.cif")

    # First converge the SPR-KKR potential, then reuse it by changing the task
    # on the same calculator object.
    calculator = SPRKKR(atoms=atoms)
    calculator.calculate()
    calculator.change_task("arpes")
    result = calculator.calculate()
    result.spc.plot()


if __name__ == "__main__":
    main()
