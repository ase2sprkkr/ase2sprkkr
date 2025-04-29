""" Example of an ARPES calculation and plotting its results. """

import numpy as np
from ase.atoms import Atoms
from ase2sprkkr import SPRKKR


def main():
    a=Atoms(symbols='H',cell=np.array([(1.,0,0),(0,1,0),(0,0,1)]))
    a.pbc=True
    xx=SPRKKR(atoms=a)
    out = xx.calculate(options={'NITER':2})
    calc = out.calculator

    # running ARPES TASK
    calc.change_task('arpes')
    calc.set('NE', 10)
    calc.set('NKTAB', 30)
    calc.set('NKTAB2D', 30)
    calc.set('NKTAB3D', 30)
    out = out.calculator.calculate()
    print(out.spc.ENERGY())
    # plotting
    out.spc.plot()
    out.spc.TOTAL.plot(filename='out.png')

    # Direct reading of earlier calculated output files
    from ase2sprkkr import OutputFile
    of2=OutputFile.from_file('H_arpes_ARPES_data.spc')
    of2.ENERGY()

    # One can do calculations with SPC files
    print((out.spc).ENERGY())
    print((of2 + out.spc).ENERGY())
    print((of2 - out.spc).TOTAL())
    (of2 - out.spc).plot()


if __name__ == '__main__':
    main()
