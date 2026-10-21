# FSCAT — audit of the current implementation

The integrated `fscat.py` input definition exposes the native FSCAT calculation through `kkrgen`. Source-level indexing limitations below remain relevant; the definition does not repair native Fortran.

## Real interface

`initvar.f:176–177` recognizes TASK FSCAT; `gen.f` calls FSCAT in KKRGEN. Routine `fscat.f` does not read SECTION_SET/SECTION_FIND, i.e. does not have task-specific parameters. I do not create ITSEL, NGEO, DIRINI, DIRFIN or output switches that do not control the input.

All atomic types count. NGEO=3 is fixed; the initial directions are [001], [100], [111]/sqrt(3), the final direction is always the opposite. These are three back-diversity geometrys, not user-adjustable general angles. The output for each type and geometry contains four spin constituents of amplitudes and their phases.

The definition uses CONTROL, TAU, MODE, STRCONST, ENERGY, TASK and SITES; `executable=kkrgen`, `mpi=False`. Custom loops of routines do not divide between wounds and create the same output names. MPI or numerical run has not been verified.

| Option |Native behaviour|Definitions|
|---|---|---|
| ENERGY.GRID |first path ETAB|explicit 1, direct real energy path|
| ENERGY.NE |NETAB(1), excluding separate task-specific default|explicit 300, positive number|
| ENERGY.ImE |init_mod_energy converts to 0|default 0; at non-zero notification selection|
| ENERGY.EMIN | -0.2 Ry |explicit -0.2 Ry, identical s native default|
| ENERGY.EMAX |at missed actual EF|None, is omitted from the input and delegate EF SPR-KKR|
| ENERGY.EMINEV/EMAXEV |relative limits to EF in eV|Common pair interface ENERGY|
| ENERGY.SPLITSS |reset to false|false; when warning is switched on|
| SITES.NL |context cutoff|None, is omitted from the input; explicit values are checked|
| MODE.MODE |default full relativity|other alert schemes, not deemed native ban|

GRID=1 and NE=300 are written to the input: they are practical options ASE2SPRKKR, not claims native defaults. NETAB starts at {0,0} (`mod_energy.f:43–44`) and init_mod_energy for FSCAT does not set a specific number of points; it is not appropriate to keep it on this initial state. For quick checking, NE can be reduced. EMAX and NL are omitted defaults and maintain context SPR-KKR. Native EMAX=EF is left here; for unfilled states, the user must increase the upper limit, for example relative EMINEV/EMAXEV.

`init_mod_energy.f:315–320` overwrites ImE and turns off SPLITSS. FSCAT uses only the first ETAB path. The transcription and the other way of the definition are highlighted. Unlike PSHIFT, it does not convert the energy itself before SSITE to its fair part: it transmits the whole ETAB (`fscat.f:162–166`). Reset ImE can therefore not be issued as a guarantee that any other contour will have all the energy real. Default GRID=1 gives a direct real path.

The output energy axis is **E-EF in eV** (`fscat.f:163`), unlike the absolute axisPSHIFT. Graphs are created using XMGRHEAD s GEO1/GEO2/GEO3; a specialized output parser is not established by this audit.

## Restrictions and validation

FULLPOT explicitly reject `check_input.f:70–72` and corresponding FSCAT branch GEN. It is recognised from the potential that the definition of input itself does not load; the prohibition is documented and enforces it SPR-KKR. DMFT and LDAU forbids `gen.f:211–222`, the definition also rejects them.

FSCAT compiles relativistic kappa and spin channels but check_input does not contain a direct requirement IREL=3. To NREL/SREL/SP-SREL the definition draws attention to unverified combinations. For IREL<2, NK/NKM init_mod_angmom is set differently than for relativistic basis, so their use in relativistic construction is suspicious; numerical impact is not verified here.

FSCAT overwrites NLT and NLQ to global NL and NKMQ to NKM (`fscat.f:79–85`). Thus, the explicit non-uniform cutoff does not remain effective; the definition draws attention. If global NL=1, CALC_CHPLM allocates T(0:0,0:0) but writes T(1,1) and T(1,0); LEGENDRE additionally writes P(2) to field length 1. The explicit s-only cutoff therefore rejects the definition. At omitted NL without potential, this check cannot be performed.

## Serious errors in source SPR-KKR

1. **Incorrect index of spherical harmonics and reading out of bounds.** FSCAT uses `LM=L*(L+1)+L+ML+1` (`fscat.f:131`). CALC_CHPLM saves m=0 to `L*(L+1)+1` and m components with a ±m (`calc_chplm.f:155–162`). In FSCAT is additional +L. For example, if NL=3, l=2, m=2 occurs LM=11, but the harmonics have only 9 elements. This is not just a cosmetic error: it can devalue amplitudes or trigger bounds error. In addition, zero Clebsch-Gordan coefficients must be omitted before accessing harmonics, because unacceptable m=l+1 may also occur in their zero components.
2. **The pre-treatment phase calculation of zero.** After `AA=ABS(BB)` is immediately calculated `LOG(BB/AA)` (`fscat.f:207–209`); the AA>1e-8 check is followed afterwards. At BB=0 it is divided 0/0. Later allocation of zeros to tables does not prevent the previous invalid operation and PHI is still used for graph limits. The correct order requires to calculate PHI until inside the nonzero branch.
3. **Restriction NL=1 is not controlled in the source.** CALC_CHPLM and LEGENDRE record elements out of bounds for l_max=0 as above. The definition limits explicit input but does not remove the library routine error or harmonic index error for NL>=2.
4. **The output text is misleading.** STOP_REGULAR states `phase shifts calculated` (`fscat.f:284`) although the calculated outputs are the matrix of the dispersion amplitude. Local auxiliary FT is calculated but not used; there is no reason to advertise it as an available output.

The definition does not issue a syntax valid input for a numerically reliable calculation. Especially the first two problems require correction SPR-KKR; they are not solved by adding task-specific parameters ASE2SPRKKR.

## Output contract checked against native SPR-KKR

`fscat.f:87–99,237–245` uses three hard-coded geometries (`GEO1`–`GEO3`) and `XMGRHEAD(DATSET,...)` to produce plot files. The current Python manifest declares `{dataset}_GEO?_*.agr`; the glob is intended to match the three native families. Because `kkrmain.f:239–245` appends nonempty `ADSI` to `DATSET`, the actual native prefix includes ADSI. The manifest uses `CONTROL.DATASET()` alone, so it may not match the output for nonempty ADSI. No new claim is made about numerical correctness.

## Verification scope

The input definition and its focused tests are present in `src/ase2sprkkr/input_parameters/`. This audit does not claim a new test-suite run or numerical SPR-KKR validation.
