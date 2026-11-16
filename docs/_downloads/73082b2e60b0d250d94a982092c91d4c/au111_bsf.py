"""Calculate the Bloch spectral function of a semi-infinite Au(111) surface."""


def main():
    from ase.build import fcc111
    from ase2sprkkr import SPRKKR
    from ase2sprkkr.sprkkr.build import semiinfinite_system

    # Build a central interaction zone between bulk Au and vacuum.
    bulk_slab = fcc111("Au", size=(1, 1, 3), periodic=True)
    surface = semiinfinite_system(
        atoms=bulk_slab, repeat=[3, 3], atoms2=None, axis=2
    )

    bulk_options = {
        "CONTROL.KRMT": 2,
        "SITES.NL": 4,
        "TAU.KKRMODE": "TB",
        "TAU.NKTAB3D": 1000,
        "TAU.NKTAB2D": 100,
        "TAU.CLURAD": 2.7,
        "SCF.VXC": "VWN",
        "SCF.NITER": 500,
        "SCF.MIX": 0.1,
        "SCF.TOL": 1e-5,
    }
    calculator = SPRKKR(atoms=surface, options=bulk_options)
    bulk_result = calculator.calculate()

    # The surface SCF starts from the converged bulk potential and uses gentler
    # mixing because convergence of the interaction zone is more demanding.
    surface_options = {
        **bulk_options,
        "TAU.NKTAB2D": 71,
        "SCF.NITER": 2000,
        "SCF.MIX": 0.01,
    }
    surface_calculator = SPRKKR(
        potential=bulk_result.potential_filename, options=surface_options
    )
    surface_result = surface_calculator.calculate()

    # Evaluate the energy-dependent BSF along a path through the surface BZ.
    bsf_options = {
        "TAU.NKTAB3D": 1000,
        "TAU.NKTAB2D": 71,
        "TAU.CLURAD": 2.7,
        "ENERGY.ImE": 0.002,
        "ENERGY.EMINEV": -4.0,
        "ENERGY.EMAXEV": 2.0,
        "ENERGY.NE": 400,
        "TASK.NK": 300,
        "TASK.KA": [[-0.66667, 0.0, 0.0], [0.0, 0.0, 0.0]],
        "TASK.KE": [[0.0, 0.0, 0.0], [0.66667, 0.0, 0.0]],
    }
    result = surface_calculator.calculate(
        potential=surface_result.potential_filename,
        input_parameters="bsfek",
        options=bsf_options,
    )
    result.bsf.plot(IQ=[0, 1, 2, 10, 11, 12])


if __name__ == "__main__":
    main()
