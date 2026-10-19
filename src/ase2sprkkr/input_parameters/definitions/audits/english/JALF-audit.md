# JALF audit and new definition

Compared to supplied sources SPR-KKR 9.7 and currently fetched branches ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. JALF does not have its own definition. A separate patch 0019-jalf.patch adds jalf.py and test_jalf.py against that base, without previous tasks.

## Calculation and input

JALF uses kkrgen and routine CURDNS: valence current density from Dirac matrix alpha, developed into real spherical harmonics with cartesian components. It's not a paramagnetic gradient contribution JORB. CURDNS simultaneously evaluates magnetic moment from current.

TASK JALF activates CALC_JALF and renames internal task to E-INTEGRAL (initvar.f:538–550). GEN calls AMEBI and allocates/zeroes JALF_LMCT (gen.f:331–336). For each energy it calls CURDNS with an integration weights WETAB, after a loop final CURDNS with EFCORRECT=true (1091–1092 and1580–1582). CURDNS does not have its own reader TASK; definition does not add non-existent networking, vector potential or rotational symmetry parameters.

## Defaulty

| Option |Native behaviour|Definitions| Serialization |
| --- | --- | --- | --- |
| ENERGY.GRID |5 for E-INTEGRAL| 5 |explicit|
| ENERGY.NE |30 for E-INTEGRAL| 30 |explicit|
| ENERGY.ImE |Common default 0,01 Ry| 0 Ry |explicit|
| ENERGY.EMIN | -0,2 Ry | -0,2 Ry |explicit|
| ENERGY.EMAX |EF in absence| None, delegace | omitted |
| MODE |full relativity in normal absence| None, FREL | omitted |
| TASK |JALF flag| JALF |explicit|

GRID/NE correspond to init_mod_energy.f:89–92. ImE=0 Ry is an explicit option ASE2SPRKKR, the same as in the common integration definition ENERGY. It differs from native 0,01 Ry; therefore it is not incorrect omitted default. Without a numerical convergence comparison, it cannot be said that native value is better. Deleted EMAX correctly takes EF.

The definition requires positive NE; warns for GRID different from the default arc 5. Native source generally allows other contours, so there is no blanket ban. SPLITSS and two items GRID/NE are rejected: GEN for E-INTEGRAL sets NEPATH=1. Local branches SPLITSS in CURDNS do not themselves show functional split-contour workflow.

## Finding and restricting

1. **MPI misses completion.** GEN:1355–1368 for theCALC_JALF itself stops before the final branch 1580–1582. CURDNS alone does not have MPI reduction JALF_LMCT. The definition therefore uses mpi=False. A mere parallel energy loop JALF does not complete the main file.
2. **Diagnostics accesses an array out of bounds at small cutoff.** CURDNS allocates JDNST(NRMAX,NLABIMAX,-1:+1,NTMAX) but fixed loop LR=1,3 writes JDNST(IR,LR,LAM,1) (curdns.f:368–375). KKRMAIN sets NLABIMAX=NLMAX (471). For NLMAX<3 therefore the diagnosis approaches out of bounds. Definition rejects explicit SITES.NL with max(NL)<3; cutoff taken from potential must additionally check running workflow.
3. **Diagnostic curdns files mix types/meshes.** The file uses R(IR,1) and JDNST(...,1), but IRTOP stays looped over type values of the last type. In addition, Mesh 1 may not be mesh type1 because IMT(1) may be different. This error is in auxiliary curdns1..3; the main DATASET0_JALF.dat selects IM=IMT(IT) and IRTOP for each type.
4. **Diagnostic names are not according to DATASET.** curdns1, curdns2 and curdns3 are fixed names, so different data are colliding in one directory. In addition, unfilled type can be read into the type1 diagnostics if ITBOT>1.
5. **Rotal symmetry is a compatible-time option.** FORCE_ROTATIONAL_SYMMETRY is PARAMETER=false (curdns.f:45–46). It is not input default that the user can change by entering TASK. Active output retains general harmonic; reduced JDNST serves magnetic torque and diagnostics.
6. **Vector potential is not calculated by default.** CALCVECPOT is set at each call false (curdns.f:93). After extracting the current, the routine returns to the line 379; the following extensive branch for A is not available without changing the source. The definition does not pretend to be support TASK.CALCVECPOT.
7. **Relativistic workflow.** CURDNS uses WAVFUN_READ_REL and AME_RLMSIGC without own scalar-relativistic resolution. The definition conservatively rejects NREL/SREL/SP-SREL; it is not documented separately native check JALF. Numerical correctness of inexplicably relativistic inputs has not been demonstrated.
8. **Corections EF and renormalizations are program status.** CURDNS for EFCORRECT will use SHFTEF as weight and multiply stored current SCLNOS at RENORMALIZE. These values are transmitted GEN, not loaded as TASK options JALF. The definition does not expose them as fictitious parameters.

FULLPOT or ASA are not explicitly prohibited for JALF in the source; it would be incorrect to impose an analogy ban with FMAG. Numerical validation of both modes did not take place here. The same applies to combinations with DMFT/LDA+U; the routine distinguishes the left and right solutions, but this patch their consistency does not verify and does not introduce a ban without a background.

## Output

The main file DATASET0_JALF.dat has the header J_alpha and subsequent TYPE IT blocks with IRTOP, NLAMAX-1 and NLMAMAX. For each harmonica and radial point writes three Cartesian components. Header NETAB(1) is the number of integration energies, not the number of density lines. The current at accumulation contains a factor C_AU (curdns.f:255–257). The magnetic moment is declared in the commentary in Bohr magnetons and uses a prefactor PRE and UNITS=0,5 (342–364); numerical units and the sign of the moment have not been independently validated.

## Tests

10 new cases: registration, serial MPI declaration, defaults, omitted EMAX, file roundtrip and pickle; three forms of split contours; other alerts GRID; three incompatible relativistic modes; zero NE; small cutoff and repair.

Full input_parameter set: **118 passed**, Python 3.12, 5,57 s.

```bash
PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 /workspace/scratch/e79ca15d9801/test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q
```

Passed git diff --check and checked the patch reverse application. Numerical SPR-KKR run not started. Patch does not fix Fortran does not even add specialized parser current density.
