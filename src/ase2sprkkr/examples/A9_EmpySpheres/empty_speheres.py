def main():

    from ase.build import bulk
    from ase2sprkkr import SPRKKR

    atoms=bulk('Cu','fcc', a=3.6, orthorhombic=True)
    atoms = atoms[:1]

    calc=SPRKKR()

    #default values
    #out = calc.calculate(atoms=atoms, options={'NITER': 1 }, empty_spheres=True)

    #min and max radius
    #out = calc.calculate(atoms=atoms, options={'NITER': 1 }, empty_spheres=(0.5, 1.5))

    #custom arguments of ase2sprkkr.bindings.empty_spheres, resulting either in calling
    #es_finder, if it is available, or xband empty_spheres routines, called from the function
    #empty_spheres in ase2sprkkr.bindings.xband.spheres (compiled from spheres.pyx)
    out = calc.calculate(atoms=atoms, options={'NITER': 1 }, empty_spheres={'max_radius': 2., 'verbose' : True, 'max_spheres': 300})
    print(atoms.positions)

if __name__ == '__main__':
    main()
