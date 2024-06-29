subdir('symmetry')
subdir('spheres')

py.extension_module(
           'symmetry',
           'symmetry.pyx',
           dependencies :  python_dep,
           link_with: symmetry_lib,
           install: true,
           include_directories: numpy_include,
           subdir: 'ase2sprkkr/bindings/xband'
)

py.extension_module(
           'spheres',
           'spheres.pyx',
           dependencies :  python_dep,
           link_with: spheres_lib,
           install: true,
           include_directories: numpy_include,
           subdir: 'ase2sprkkr/bindings/xband'
)
