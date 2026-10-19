# JORB audit and new definition

Basic: supplied SPR-KKR 9.7 and currently fetched ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. JORB does not have its own definition. Single patch 0018-jorb.patch adds jorb.py and test_jorb.py against this base; does not include earlier tasks.

## Meaning and Input

JORB activates CALC_JORB. initvar.f:539–550 switches internal TASK to E-INTEGRAL and sets NOWRDOS=true, RDTAUDIR=false and WRTAUDIR=false. This is not about reading existing DOS or TAU directory. GEN for E-INTEGRAL sets writing regular and irregular wave functions, CALCINT=true, GETIRRSOL=true, ICORE=1 and NEPATH=1.

V every energy calls GEN ORBCURDNS s WETAB(IE,IEPATH) (gen.f:1098–1099); integration followed by final ORBCURDNS s EFCORRECT=true (1586–1588). The routine describes the relativistic paramagnetic current density of valence electrons in Atomic units, in full-potential formalism. The result is neither the Stoner parameter nor the magnetic form factor.

ORBCURDNS does not have its own TASK reader. The new definition therefore does not add non-existent options for network or components. General QREPONLY reads the common branch E-INTEGRAL but is not documented as a separate parameterization JORB output; I do not display it here as JORB option.

## Defaults and their writing

| Option |Native behaviour|Definitions| Serialization |
| --- | --- | --- | --- |
| ENERGY.GRID |5 for E-INTEGRAL| 5 |explicit|
| ENERGY.NE |30 for E-INTEGRAL| 30 |explicit|
| ENERGY.ImE |Common default 0,01 Ry, without JORB override| 0 Ry |explicit|
| ENERGY.EMIN | -0,2 Ry | -0,2 Ry |explicit|
| ENERGY.EMAX |EF in absence| None | omitted |
| MODE |normally full relativity in absence| None, FREL | omitted |
| TASK |JORB flag| JORB |explicit|

GRID=5 and NE=30 correspond directly to init_mod_energy.f:89–92. ImE=0 is an explicit option corresponding to a common integration definition ENERGY, different from a general native EIMAG=0,01. Therefore it is not wrong omitted default. Native 0,01 is not clearly better without convergence comparison; for a complex arc depends on the contura itself. Deleted EMAX correctly takes EF from potential.

Another GRID will trigger an alert to the necessary occupied integration interval, not a hard ban on all alternative contours. NE is to be positive. SPLITSS or two items GRID/NE are rejected because E-INTEGRAL sets one path. ORBCURDNS contains branches SPLITSS, but do not themselves provide completed two-konture workflow in this branch GEN.

## Resources limitations and findings

1. **MPI does not complete the final JORB output.** After the energy loop GEN:1355–1368 continues in MPI only for the listed tasks or CALC_JXC/CALC_DMI. Internal E-INTEGRAL with only CALC_JORB condition fails and stops before the final call in line 1588. In addition, ORBCURDNS does not reduce itself JORB_LMSCT between processes. The definition has mpi=False. The presence of a common MPI energy loop therefore does not mean functional parallel JORB.
2. **Relativistic assumptions.** ORBCURDNS uses WAVFUN_READ_REL and relativistic indices/couplings; has no own IREL branch for scalar-relativistic solution. Conservatively rejects the definition NREL/SREL/SP-SREL. This is not a documented separate native check JORB but a restriction on workflow based on the implementation used.
3. **FULLPOT is intended formalism, not documented ban ASA.** The header ORBCURDNS states full-potential formalism. CHECK_INPUT does not have JORB ban ASA or FULLPOT. The definition therefore does not invent a ban by mere comment. Numerical accuracy for ASA or FULLPOT has not been verified; the potential must provide the necessary angular representation.
4. **Angular cutoff has an internal condition.** AME_NAB_RLMC stops NKM>UBOUND(CGC,1) and NLAMAX<2*NLMAX (orbcurdns.f:462–466). These are derived internal dimensions; alone SITES.NL is not sufficient to verify them. I don't create non-existent TASK.NLAMAX.
5. **Automatic diagnostic output fort.*.*** For odd l and abs(m)=1 routine writes per unit 1000+LMR without explicit OPEN with name (lines 374–386). Fortran runtime typically creates fort.<unit>. These files are not differentiated by DATASET or type; multiple types add data to the same unit. They describe diagnostic harmonics, they are not a clear main output.
6. **Final correction is not input constant JORB.** WE at EFCORRECT=true is SHFTEF; at normal loop WEINP. RENORMALIZE if applicable, multiplys the stored JORB_LMSCT value SCLNOS. GEN initializes SHFTEF=0 and RENORMALIZE=false or adjusts them to other branches (e.g. Lloyd). The definition does not exhibit SHFTEF/SCLNOS as undocumented TASK options.
7. **Exit NE is not the number of density rows.** WRHEAD gets NETAB(1), but the data are written by types, LM and radial points, with three cartesian components. Parser cannot consider header NE as a number of data rows.

Main files have names DATASET0_JORB.dat, DATASET0_JORB_UP.dat and DATASET0_JORB_DN.dat. The total output is the sum of two spin components. For each type there are IRTOP, NLAMAX-1 and NLMAMAX, followed by LM blocks of radial values. These are not individual atomic currents without harmonious development.

DMFT/LDA+U are not rejected in the definition without a direct source ban; ORBCURDNS distinguishes different left/right solutions. The consistency of these combinations has not been numerically verified. Patch does not fix Fortran does not even add specialized parser.

## Tests

9 new cases: registration, serial declaration, defaults, omitted EMAX, file roundtrip and pickle; three ways of entering split contours; alerts for another GRID; three incompatible relativistic modes; zero number of energy.

Full input_parameter set: **117 passed**, Python 3.12, 5,30 s.

```bash
PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 /workspace/scratch/e79ca15d9801/test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q
```

Passed git diff --check and checked the patch reverse application. Numerical SPR-KKR run not started.
