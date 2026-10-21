# JXC — current task audit

`jxc.py` is the shared JXC/DMI input definition. For JXC, native `initvar.f:535–545` selects magnetic `E-INTEGRAL`; `gen.f:930–996,1193–1255` branches between SP-SREL scalar exchange and FREL tensor exchange. The restrictions below concern the supported periodic workflow, not a numerical certification.

## Scope and corrections

JXC is already supported. INITVAR reads the keyword JXC, sets CALC_JXC and internal task E-INTEGRAL. GEN:930–990 distinguishes SP-SREL (IREL=2) with XCPLJ0/XCPLJIJ and FREL (IREL=3) with XCPLTENSOR. NREL/SREL are not the same as SP-SREL; the initial warning_condition also prohibited only the first two modes. New tests explicitly verify that SP-SREL remains accepted.

The integrated schema exposes the native TASK.NSHLCLU, SELECT, EXCLUDE and QREPONLY. INIT_MOD_TAUIJ_STAR:125–161 reads lists of chemical symbols, prohibits their current use and gives shell count preference over CLURAD when present NSHLCLU. QREPONLY reads INITVAR inside E-INTEGRAL branches; affects the list of central sites in INIT_MOD_TAUIJ_STAR. CLURAD is in units ALAT, not angstroms.

The current validation rejects NREL/SREL/NONMAG, not positive NE, more energy routes, unaudited integration and current SELECT/EXCLUDE. Mode and NONMAG warnings accompany the present validation. For the standard periodic calculation, BZINT=POINTS is explicit. GEN:358–386 for off-diagonal TAUIJ requires IBZINT=2.

The audited range is normal periodic JXC on GRID=5. In the source there are other branches embedded-cluster, occupation-weighted GRID3 (WGT_OCC), MULTISPIN, SLC_JXC and SLC_TRQ. They do not have a complete audit; therefore, the validator is not a claim that SPR-KKR generally never knows GRID3 or cluster. Those expanded native branches are outside the present input contract.

## Defaults and their writing

| Option | SPR-KKR |Definitions|
| --- | --- | --- |
| TASK.CLURAD |without NSHLCLU a value is required; a generic default is not|explicit2,2 ALAT, left|
| NSHLCLU |selects CLURAD in absence| None, omitted |
| SELECT/EXCLUDE |blank lists| None, omitted |
| QREPONLY |false in absence| false, omitted |
| ENERGY.GRID |5 for E-INTEGRAL|explicit5, left|
| ENERGY.NE |30 for E-INTEGRAL|explicit 30|
| ENERGY.ImE |generally 0.01 Ry|explicit 0.01 Ry|
| ENERGY.EMIN | -0,2 Ry |explicit -0,2 Ry, left|
| ENERGY.EMAX |General Sentinel replacement for EF|omitted, taken over EF|
| MODE |Commonly FREL|None, omitted; SP-SREL with writes when entering|
| TAU.BZINT |general choice depends on workflow|explicitly POINTS for documented off-diagonal TAU flow|

The integrated definition sets NE=30 and ImE=0.01 Ry, matching these source defaults. GRID=5 still requires physical convergence checks; no numerical comparison has been performed. Radio2,2 is a convenient start, not a guarantee of the converted sum J_ij. Convergence NE, the lower limits, k-mesh and cluster remain on numerical verification.

## Serious ALE in SPR-KKR

1. **Out rotation according to symmetry at FREL:** XCPLTENSOR:312 has JQ.EQ.IQREPQ(JQ) .OR.1.EQ.1. Thus, the second branch with rotation ME_ROT_MREG will never run. For non-representative neighboring site can be used matrix in the wrong framework. QREPONLY restricts central site, does not automatically solve all neighbors. This is not what I call breaking up every single-atomic collinear run; it is an important mistake in the general tenor workflow.
2. **FREL has potential conditions:** NONMAG=false, KMROT=0, NKM=NKMMAX and representative IQAT(1,IT) are directly controlled XCPLTENSOR:74–80,155. Python knows explicit NONMAG, other conditions require real potential. Cannot enter QREPONLY to certify all its properties.
3. **More routes are not processed in general:** routines use ETAB/WETAB first routes and have NEPATH=1 guard. Later enforcement E-INTEGRAL does not change this condition to support split contour.
4. **Cluster SP-SREL is not a normal periodic flow:** CLUXCPLJIJ reads scratch unit19 but TAU_DRIVE normally calls TAU_REAL_SPACE with KTAUIJ=0 which shuts down its export. Without an alternative unit preparation audit19, it is not right to promise this workflow only by switching BZINT=CLUSTER. FREL distinguishes embedded-cluster separately.
5. **Higher spin and spin-lattice branches have their own interface:** GEN reads MULTISPIN/SLC_JXC/SLC_TRQ and selects other routines. The current schema does not expose these branches as fully supported tasks.

Unlike BSDOSIJ uses XCPLTENSOR for JC correct offset BAR_M_JT_JQ(1,1,JC) and output J_ITJT has consistent pair indexing. Errors from BSDOSIJ cannot be mechanically attributed JXC only due to a similar source.

## MPI and output

mpi=True is maintained and supported: J_IJINT is zeroed at the first call, MPI_SUM_C performs reduction and GEN:1193–1241 has an additional finishing call. SP-SREL has its own XCPLJIJ reduction. Numerical serial/MPI match has not been measured.

The FREL tensor output is in meV (DIMAG(J_IJINT)/pi*RY_EV*1000); SP-SREL uses scalar Jij and a separate result family. No numerical verification of the output parser or downstream Hamiltonian conventions is claimed.

## Current output contract

The shared `jxc.py` manifest selects the scalar `{dataset}_J_ij.dat` and `{dataset}_Jij_*.agr` for SP-SREL JXC, or the `{dataset}_XCPLTEN_Jij.dat`, `{dataset}_XCPLTEN_Dij.dat` and matching graph families for FREL JXC. When `TASK.DMI` is also enabled, it adds `DMIVEC_Dij` files. Native `initvar.f:535–545` selects the shared E-INTEGRAL computation and `gen.f:951–996` dispatches the modes; `dmivector.f:416–418` confirms the `TRIM(DATSET)` prefix for DMI. As `kkrmain.f:239–245` appends nonempty `ADSI` to native `DATSET`, the current manifest, which uses `CONTROL.DATASET()` alone, can misname the outputs with nonempty ADSI. This is a filename-manifest discrepancy, not proof that JXC numerics fail.

## Verification scope

Input serialization and validation have dedicated regression tests; no fresh test-suite run or numerical SPR-KKR calculation is claimed here.
