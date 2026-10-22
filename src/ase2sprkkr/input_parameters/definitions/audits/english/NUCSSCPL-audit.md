# NUCSSCPL — SPR-KKR 9.7 and ASE2SPRKKR audit

NUCSSCPL has an input definition **prepared in `tasks-workflow`**, absent from `develop`. The native SPR-KKR 9.7 implementation has demonstrated defects; adding an input definition does not fix its calculation.

## Purpose and entry

NUCSSCPL counts tensor nuclear spin-spin bonds A_ij and especially Ruderman's Kittel's contribution R_ij. INITVAR:491–496 sets TASK=NUCSSCPL, RUNELOOP=true and ICORE=0. This is not a renaming to E-INTEGRAL.

GEN:358–386 prepares off-diagonal TAUIJ over INIT_MOD_TAUIJ_STAR and requires IBZINT=2 (BZINT=POINTS). TASK therefore contains CLURAD, NSHLCLU, SELECT and EXCLUDE from a common reader. If NSHLCLU is entered, it takes precedence over CLURAD; both selectrs simultaneously lead to STOP. Lists have a capacity of100 and contain chemical symbols.

TASK.QREPONLY has not been added: INITVAR reads it inside the branch CALC_FMAG/JXC/DMI/JALF/JORB/BSDOSIJ, where it does not enter the branch NUCSSCPL itself. Copying the same offer as for BSDOSIJ would be wrong.

NUCSSCPL:77–83 requires NONMAG, IREL=3, KMROT=0, one energy path, point core and GRID=5. GEN:1043–1047 stops for FULLPOT. The definition is valid for NONMAG, MODE, CONTROL.FINNUC, energy path, BZINT, inexplicable NE and exclusiveness of selectors. FULLPOT from potential and KMROT must be verified while running; Python definition does not know them.

Isotopes are not set by input task. TABMUQIN (table.f) selects for Z table isotope with non-zero magnetic moment; the table covers Z=0..110 and contains also zero items. The definition therefore does not invent options ISOTOPE, MNUC or INUC.

## Defaults and their writing

| Option |SPR-KKR in absence|Definitions|
| --- | --- | --- |
| CONTROL.NONMAG |Commonly false, this task then STOP|true, explicit|
| CONTROL.FINNUC | false | false, omitted |
| MODE |Commonly FREL, IREL=3| None, omitted |
| TAU.BZINT |common general choice|POINTS, explicit|
| TASK.CLURAD |without NSHLCLU does not have a usable default, STOP|2,2 ALAT, explicit|
| TASK.NSHLCLU |Reroute absence to CLURAD| omitted |
| TASK.SELECT/EXCLUDE |blank lists| omitted |
| ENERGY.GRID |NUCSSCPL does not have its own initialization branch; the routine requires5|explicitly5|
| ENERGY.NE | MOD_ENERGY inicializuje0; NEMAX=0 vede k STOP |explicitly30|
| ENERGY.EMIN | -0,2 Ry |explicitly -0,2 Ry|
| ENERGY.EMAX |General Sentinel replacement for EF, INIT_MOD_ENERGY:182| omitted, None |
| ENERGY.ImE | 0,01 Ry |explicitly0,01 Ry|

Explicit differences NONMAG, GRID, NE and CLURAD prepare the required input; omitted FINNUC, MODE and EMAX correspond to native behaviour. The radius2,2 is not shown to be converged and should be checked against larger clusters. NE=30 is a common integration option, numerical convergence not verified. For GRID=5 EPATH, Gauss's quadrature is on a complex arc and ImE is not used in this branch; therefore, this is not a physical difference against a common explicit ImE=0 in ASE2SPRKKR.

## Mistakes in Fortran

1. **Deallocation non-allocated fields:** for each energy except the last energy, DEALLOCATE(AITJT,RITJT,DITJT,NITJT,DVECITJT,XVEC,YVEC) shall be performed in line 341. Their first ALLOCATE is in the finishing branch 366–368. NE=30 thus encounters unprepared field detachment after the first energy. Change NE to 1 is not a substantive repair of integration workflow.
2. ** Loss of matrix of previous types:** loop DO IT=1,NT in line 115 calls CINIT(M*M*3*NT,ME) to 133. It zeroes the whole ME for all types at each IT; only the last type remains filled after the loop. Pair contributions of other types are lost.
3. **Overridden radial mesh:** IM=IMT(IT) is replaced immediately IM=1 (118–120). Reading and integration of all types thus use JRWS(1), independently of their assignment of mesh.
4. **Dangerous diagnostic division:** CHECK is fixed true. RAT divides MEZZ/ME for non-zero XA+XB, which does not guarantee non-zero denominator ME (198–200). Output filters AIJ/A11 and RIJ/R11 (440 and454) divides the zero reference diagonal member.
5. **Incorrect plot scaling:*** by loop over pairs is used again for graphs FTP last set for a specific pair. Multiple isotope moments require a factor according to the given IT,JT. In addition, YMIN=YMAX*FTP and subsequently YMAX=YMIN*FTP (547–548) overwrite the lower limit of the upper and upper limit multiplied twice.
6. **MPI is not completed:** AIJINT/RIJINT are allocated only at IECURR=1, without MPI reduction and without process resolution in output. Energy divided between processes does not have a common initialization or a complete integral. The definition uses mpi=False.

The main file is named DATASET_A_ij.dat; other outputs are fixed RK.dat and tab_aij.dat and Aij/Rij graphs. Python parser results are not implemented. Due to the above errors support input does not justify the claim that native numerical task works.

## Verification scope

The prepared Python definition has test coverage for serialization and selected validation rules. Counts from old standalone patch runs are not the test results of the current `tasks-workflow` branch. No numerical native SPR-KKR reference calculation is claimed here.
