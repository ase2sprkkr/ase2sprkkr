python = import('python')

py = python.find_installation(pure: false)
python_dep = py.dependency()
python_exe = py.full_path()
glob = run_command(python_exe, 'sources.py', '[.](py|in)', check: true)

sources = glob.stdout().strip().split('\n\n')
common_fortran_args = [ '--std=legacy', '-cpp' ]
numpy_include = run_command(python_exe,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()


dir = ''
foreach current : sources
    if dir != ''
       foreach j: current.split('\n')
         message('>'/ j / '<')
       endforeach
       py.install_sources(
         current.split('\n'),
         subdir : dir
       )
       dir = ''
     else
       dir = current
     endif
endforeach

subdir('ase2sprkkr/bindings/xband')
