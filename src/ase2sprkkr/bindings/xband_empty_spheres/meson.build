subdir('symmetry')
subdir('spheres')

py.extension_module(
           'symmetry',
           'symmetry.pyx',
           dependencies :  python_dep,
           link_with: symmetry_lib,
           install: true,
           subdir: 'ase2sprkkr/bindings/xband_empty_spheres'
)

py.extension_module(
           'spheres',
           'spheres.pyx',
           dependencies :  python_dep,
           link_with: spheres_lib,
           install: true,
           subdir: 'ase2sprkkr/bindings/xband_empty_spheres'
)
