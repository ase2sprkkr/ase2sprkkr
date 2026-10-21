# BSDOSIJ — current definition audit

`bsdosij.py` defines experimental site-resolved back-scattering DOS. The native and input-parameter constraints below refer to SPR-KKR 9.7; numerical accuracy has not been established.
## Input and defaults

TASK BSDOSIJ starts site-resolved back-scattering DOS in kkrgen; initvar.f:540–554 changes internal task to E-INTEGRAL. GEN:358–386 prepares off-diagonal TAUIJ and requires IBZINT=2, i.e. BZINT=POINTS. INIT_MOD_TAUIJ_STAR reads TASK.CLURAD or NSHLCLU, lists of chemical symbols SELECT/EXCLUDE and uses QREPONLY from INITVAR.

| Option |Native missed behaviour|Definitions and writing|
| --- | --- | --- |
| GRID |E-INTEGRAL sets 5, BSDOSIJ then requires 3|explicitly 3; necessary repair of non-functional native default|
| NE | 30 |explicitly 300 for energy resolution; convergence not verified|
| ImE | 0,01 Ry |explicitly 0,01 Ry|
| EMIN | -0,2 Ry |explicitly -0,2 Ry|
| EMAX | EF |omitted, None delegate to EF|
| CLURAD |no applicable default; without NSHLCLU the program stops|explicitly 2,2 ALAT as well as existing JXC definition|
| NSHLCLU |absence selects CLURAD|omitted; it takes precedence over CLURAD|
| SELECT / EXCLUDE |blank lists|omitted; mutually exclusive, not more than 100 symbols|
| QREPONLY |false in absence of a keyword|false; does not write|
| MODE |Commonly FREL|omitted; other MODE refused|

Differences GRID, NE and CLURAD are explicit and not erroneous omitted defaults. Native GRID=5 is not better for this task: leads to STOP. NE=30 may be sufficient for gross control, but 300 provides a softer energy mesh; suitability must verify the calculation. The radius 2,2 is not universally converged. ImE remains native.

The definition rejects NONMAG, other MODE, other GRID, two energy routes, SPLITSS, non-positive NE and other BZINT. In addition, the source requires KMROT=0, NKM=NKMMAX and representative IQAT(1,IT); these conditions depend on potential and the module is does not verify. QREPONLY changes central sites but does not guarantee the fulfilment of all these checks.

## Serious findings SPR-KKR

1. **The matrix component:** bsdosij.f:319–349 iterates IC=1..3, but the first CMATMUL receives BAR_M_JT_JQ without offset IC. It transmits the beginning of the field, i.e. the first component, while the second local matrix uses IC. This seems to confuse the components of the contribution to JQ.
2. **Surfaced symmetrical rotation:** condition JQ.EQ.IQREPQ(JQ) .OR.1.EQ.1 in line 325 is always true. Rotation branch ME_ROT_MREG is unattainable.
3. **Chained output indexing:** J_ITJT is allocated as (NIJ,NT,NT,3,3), but lines434–436 store J_ITJT(IE,I,IT,JT,IC) as if it had an energy dimension. IE may exceed NIJ, I exceeds NT and JT exceeds 3. Subsequent reading uses indices (I,IT,JT,component,component). I does not read correctly stored values if the limits are observed accidentally; the blanks are not initialized.
4. **The output is copied by exchange coupling:*** file is named DATASET_XCPLTEN_Jij.dat and other outputs/planes are designated as exchange links. Final processing does not go through all energy BSDOS_IJ. The name and intended DOS behaviour thus do not correspond to the completed output workflow.

These errors does not fix Python domain file. Support input is experimental; numerical results SPR-KKR 9.7 are not confirmed by this audit.

## MPI
mpi=True corresponds to the implemented reduction BSDOS_IJ over MPI_SUM_C (bsdosij.f:112–126) and output at the process0. Numerical accuracy alone MPI has not been tested.
