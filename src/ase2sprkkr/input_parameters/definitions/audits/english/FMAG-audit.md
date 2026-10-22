# FMAG — SPR-KKR 9.7 and ASE2SPRKKR audit

FMAG is **prepared in `tasks-workflow`**, not present in the published `develop` branch. The native SPR-KKR 9.7 behavior is analyzed independently of this Python definition.

## What actually starts the input FMAG

TASK FMAG activates CALC_FMAG and initvar.f:534–550 renamed internal TASK to E-INTEGRAL. Sets NOWRDOS=true, RDTAUDIR=false, WRTAUDIR=false, ICORE=0 and IHYPER=0. GEN for E-INTEGRAL turns on writing of both wave functions, CALCINT and GETIRRSOL, sets ICORE=1 and NEPATH=1 (gen.f:296–303). The ICORE setting in INITVAR is therefore not final.

GEN calls DENSITY for each energy and after integration (gen.f:1078–1084,1180–1186). DENSITY accumulates charge, spin and orbital density, creates dens and fmag outputs and ends the program using STOP_REGULAR. FMAG does not have its own momentum-mesh parameter reader or specific TASK options for form factor. The new definition does not invent QMAX or NQ. The general QREPONLY switch is read in the shared E-INTEGRAL branch, but DENSITY itself does not use it; this definition therefore does not expose it as an FMAG parameter.

## Defaults and writing

| Parameter |Native behaviour|Definitions| Serialization |
| --- | --- | --- | --- |
| ENERGY.GRID |5 for internal E-INTEGRAL| 5 |explicit|
| ENERGY.NE |30 for internal E-INTEGRAL| 30 |explicit|
| ENERGY.ImE |Common default 0,01 Ry, without FMAG override| 0 Ry |explicit|
| ENERGY.EMIN | -0,2 Ry | -0,2 Ry |explicit|
| ENERGY.EMAX |EF in absence|None, delegation per source| omitted |
| MODE |normally FREL when omitted| None, FREL | omitted |
| TASK |FMAG flag, no other form-factor options| FMAG |explicit|

GRID=5 and NE=30 correspond to init_mod_energy.f:89–92. ImE=0 Ry is an explicit option ASE2SPRKKR for complex arc integration, the same as common definition ENERGY used JXC. It differs from the general native EIMAG=0,01 Ry; because it is registered, it is not incorrect omitted default. Native 0,01 is not automatically better: the construction and convergence of the integration contour itself is important. Numerical comparison of both values was not performed. On the contrary, omitted EMAX must leave a true EF, not a fixed value.

The definition has kkrgen and mpi=True. DENSITY contains MPI barrier and MPI_SUM_R for density and OBS (density.f:254–264), subsequent files writes MPI_ID=0. This branch completes the output before general postprocessing MPI by checking in GEN. The parallel run was not numerically verified.

## Source findings

1. **FULLPOT is actually prohibited even if the original check can miss task.** check_input.f:119–123 controls TASK='FMAG' and GRID=5. The internal renaming to E-INTEGRAL can skip this branch. However, GEN stops directly in CALC_FMAG branches FULLPOT (1078–1081 and1180–1183). The definition documents ASA; FULLPOT is a potential property that the general validator InputParameters does not know.
2. **The integration contour should be occupied by the energy interval.** GRID=5 is native the requirement of the original FMAG check. The definition preserves it and warns for another GRID instead of claiming that the source always rejects such input. Convergence depends on the lower between the number of points.
3. **Only the first energy path.** GEN for E-INTEGRAL sets NEPATH=1 and DENSITY takes WETAB(IECURR,1). Definition rejects SPLITSS and two items GRID/NE. It is a protection against ignored second path; no specialized split-density workflow is introduced here.
4. **Normalisation by dividing by zero.** density.f:397–400 divides FAVR, FSPN and FORB by their values at q=0 without control. Zero total, spin or orbital moment can cause NaN or Inf. The input definition does not have calculated moments, so in general it cannot detect the problem. For explicit NREL/SREL warns against normalization nonmagnetic result. Not even FREL guarantees non-zero moments, of course.
5. **Fixed momentum mesh.** NQNEU=min(200,100)=100, QNEUMAX=4*pi*A0_ANG and output axis x=q/(4*pi*A0_ANG) runs from 0 to 1 (density.f:380–403,430–431). NQNEUMAX=200 is not the number of points actually written or input default. The range is not read from TASK.
6. **Type size and MPI reduction are not equally wide.** Densities are reset to ITBOT..ITTOP (density.f:83–91), but MPI_SUM_R has a length NRMAX*NTMAX. For calculation limited to a part of types, it can reduce even uninitialized unused parts of fields. The specific impact on the results requires running verification.
7. **Relativist reader is unconditional.** DENSITY calls WAVFUN_READ_REL, but also contains its own branches IREL<=1 and a warning of dimensions (density.f:107,121–156). Therefore, I am not introducing an unwarranted total ban on modes other than full relativity. Their numerical accuracy is not verified.
8. **The output is normalized curves.** FAVR is the sum of spin and orbital contributions before separate normalization. After normalization generally does not apply FAVR=FSPN+FORB. Parser or further processing shall not add up these normalized curves without the original moments.

The combination FMAG with DMFT/LDA+U is not rejected in the prepared definition without direct native ban; their physical consistency with input potential and the kernel calculation has not been numerically verified. Patch does not create specialized parser XMGR files and does not fix Fortran.

## Verification scope

The prepared Python definition has test coverage for serialization and selected validation rules. Counts from old standalone patch runs are not the test results of the current `tasks-workflow` branch. No numerical native SPR-KKR reference calculation is claimed here.
