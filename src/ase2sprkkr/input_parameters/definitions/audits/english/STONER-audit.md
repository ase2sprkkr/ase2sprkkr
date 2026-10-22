# STONER audit and new definition

Compared to supplied sources SPR-KKR 9.7 and currently loaded branches ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. V this branch STONER does not have its own definition. Patch 0012-stoner.patch is independent against the given commit, does not include previous tasks.

## Entrance and supported workflow

STONER starts kkrgen; TASK does not read any of its own parameters. initvar.f:384–390 sets ICORE=0, NOWRDOS=true and RDTAUDIR=true. gen.f:1592–1610 compiles DOS weights and calls SUSC_GAMMA.

The definition contains CONTROL, TAU, MODE, STRCONST, ENERGY, TASK and SITES. It exhibits conservative workflow with one real energy at EF and relativistic wave functions. Alternative energy scans generally permit the source, but STONER then uses current DOBS_LTX and wave functions; the definition rejects them so that the result does not indicate the last scan energy as EF. This is a restriction of a new definition, not a claim of native ban of scans.

## Defaulty

| Option | Native STONER default |Definitions| Serialization |
| --- | --- | --- | --- |
| ENERGY.GRID | 0 | 0 |explicit|
| ENERGY.NE | 1 | 1 |explicit|
| ENERGY.ImE | 0 Ry | 0 Ry |explicit|
| ENERGY.EMIN | EF |None, delegation to SPR-KKR| omitted |
| ENERGY.EMAX | EF |None, delegation to SPR-KKR| omitted |
| ENERGY.EMINEV/EMAXEV |n.a.| None | omitted |
| ENERGY.SPLITSS | false | false | omitted |

Energy defaults are directly in init_mod_energy.f:108–113. From the common definition ENERGY, the explicit EMIN=-0,2 Ry: overwritten native EF shall not be taken over. Absolute EF is unknown when creating input, therefore both absolute bounds remain absent. Zero relative offsets are also allowed; native the reader applies relative limits together in the presence of both. Non-zero offsets and explicit absolute limits reject this definition.

There is no difference between physical native default and default ASE2SPRKKR. General SITES.NL may be undefined and left to potential; the validator only controls explicit values. The actual NLQ potential points must be sufficient for this workflow.

## Finding and restricting

1. **MPI does not complete STONER.**gen.f:1355–1368 has a list tasks allowed after an energy loop. STONER is not there and MPI run ends STOP_REGULAR before call SUSC_GAMMA. The definition has mpi=False.
2. **Relativistic format is used without branching according to IREL.** SUSC_GAMMA reads WAVFUN_READ_REL (susc_gamma.f:108) and works with relativistic indices kappa/mj. The definition rejects NREL, SREL and SP-SREL. This is a safe restriction according to the format used; there is no separate ban on these modes for STONER.
3. **D channel is used without Brooks OP.** NLOP=3 is set at start SUSC_GAMMA and GAMMAL(IR,NLOP,IS) is used without condition (row 207). SUSC_ORBPOL also writes PRE(4), with PRE in size 2*(NLMAX-1). NL<3 can therefore mean out-of-bounds array access or read the unfilled channel. Definition rejects explicit NL<3. Native branch BROOKS-F uses NLOP=4 and needs f channel; through shared MODE.OP this special value is not exposed and this branch has not been tested.
4. **Kernel is firmly VWN.** SUSC_GAMMA:236 calls directly SUSC_VWN. The result therefore does not automatically use the kernel of the same function as SCF potential if it was created by another function. The definition does not invent TASK.VXC.
5. **Little DOS is not physically neutral fallback.** SUSC_GAMMA:87–90 at NEFT<=1e-8 sets all NEFTL=1 and NEFT=NL. The source himself points out that Stoner's enforcement will be unrealistic. This cannot be detected in advance only from the input.
6. **Normalization has unprotected division.** SUSC_ORBPOL:100 divides BRACAH/ORBSQINT(NLOP,IS) without control of the zero denominator; diagnostic output SUSC_GAMMA:290–291 also divides RHORINT. These are numerical resource risks, not repairable defaults input.
7. **Unbalanced cutoffs deserves a numerical audit.** gen.f meets NEFTL according to global NL, SUSC_GAMMA processes them according to local NLQ. Without verifying initialization and specific potential, I do not claim that such configurations are incorrect.

FULLPOT is not explicitly prohibited for STONER in check_input. I do not create a ban only according to the spherical structure of this routine; however, the physical completeness of its FULLPOT result is not verified. The same applies to combinations with correlations. Numerical run SPR-KKR has not been performed nor has there been a specialized parser output.

## Tests

12 new cases: registration, native energy and omitted limits, file roundtrip and pickle; six failed energy settings; zero relative limits; three incompatible wave modes; too small a cutoff and its correction.

Full input_parameter set: **120 passed** (Python 3.12, 4,98 s).

```bash
PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 /workspace/scratch/e79ca15d9801/test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q
```

Passed git diff --check and checked the patch reverse application. The tests verify the definition and serialization, not the numerical accuracy of the Stoner parameter.
