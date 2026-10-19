# POSANI — audit SPR-KKR 9.7 / ASE2SPRKKR

Order 31. SPR-KKR: PUB9.7_260420_18h13.tar.gz.
ASE2SPRKKR refactor: 267512c4b8ccbbfd37af7ede898febc066b7e3ee.
Single patch 0031-posani.patch against this commit; does not require prior
COMPTON patch and does not add another task.

## Stav

The definition POSANI has been missing so far. Added definitions for standard periodic
FREL calculation of momentum density of positrons, over kkrgen, with one
Energy contour GRID=5 and electron potential.

Preliminary classification: **works but has significant caveats** for limited workflow
with NQ=NT; **broken** for general systems NQ≠NT at work with positron TAUQQ
(see below). Native task exists and is not disabled. Physical values were not
verified numerically during.

## Real workflow and inputs

`gen.f:399–426` distinguishes POSANI, POSANIPREP, POSANITAU and MODE.POSITRON.
The priority is POSANI > POSANIPREP > POSANITAU > POSITRON. They're different steps:

1. Electron potential and electron TAU supporting documents shall correspond to the system,
the energy and the wave function convention used.
2. POSANIPREP prepares DATASET_pos_dat.pan with positron status. This task
and POSANITAU will be assessed separately; this definition does not trigger them.
3. POSANI reads **electronic potential**, TAU directory and this prepared file.
`gen.f:431–449` for POSANI does not call POSANIPOT. POSVC is therefore not the input of this
step; it belongs to the preparation of positronal potential.

`compton.f:424–426` calls POSANIREAD. Positron file name is fixed
DATASET_pos_dat.pan (`posaniread.f:38`), not path derived from TAUDIR.
The change DATASET thus changes the expected name of the positron file; TAUPATH does not redirect it.
The definition deliberately does not invent a non-existent keyword for its any way.

Common momentum grid inputs in TASK:
MEMODE=1/2, NPMAX, PNMAX, PPMAX, NPN, NPP, USE_PEQUIV,
PNVEC, PPVEC1, PPVEC2, MD2D, EMD2D, SHFT, SHFT_PP, SHFT_PN,
INSIDE_SPHERE and SQQ_RESOLVED. It reads them COMPTON and INIT_MOD_COMPTON_PMESH.
TAU.TAUDIR and TAU.TAUPATH set electron TAU supporting documents.

The definition allows NONMAG, requires FREL and POINTS. GRID=5 and one contour are
reduction of the audited workflow, not the claim that native routine prohibits all
other contours. ONLY_CORE and USE_ZJOVZZ are not exposed; they require further audit.

## Important difference versus COMPTON

INITVAR detects COMPTON before INIT_MOD_ENERGY but POSANI is searched for up to
v GEN (`gen.f:417`); TASK is switched to POSANI in `gen.f:522`.
INIT_MOD_ENERGY is already called in `kkrgen.f:1039`.
Therefore, only TASK POSANI **does not mean automatic take-over** COMPTON
GRID=5/5, NE=30/50 and SPLITSS=true.

The definition sets CONTROL.RDTAU=true and always writes keyword; INITVAR then turns off
normal energy loop (`initvar.f:687`). Step Using Existing Electron
TAU podklady. Without RDTAU this definition shall not pass at validation.

COMPTON subsequently rereads SPLITSS from TASK (`compton.f:210`). Definitions
POSANI exhibits ENERGY.SPLITSS and TASK.SPLITSS=false and rejects them,
that the preparation of energy and later calculation are consistent. Unlike the previous one
COMPTON task there is no SPLITSS keyword does not write.

## Defaults and writing

| Parameter |Definitions|Native effective behaviour without parameter|Assessment|
|---|---|---|---|
| RDTAU |true, always writes|INITVAR default false; POSANI does not set it itself|explicit workflow option, not match native default|
| GRID, NE |5, 30, entered|is not POSANI-specific default before initialization of energy; it depends on general state and TAU backgrounds|explicit option, must match the electron cache|
| ImE |0.01 Ry, writes| 0.01 Ry |the same|
| EMIN |-0.2 Ry, writes| -0.2 Ry |the same|
| EMAX | nevystaven, is omitted from the input |Sentinel is converted to EF in INIT_MOD_ENERGY:182|the same effective behaviour|
| SPLITSS v ENERGY and TASK | false, is omitted from the input | false v tomto postupu |same; switch rejected as other workflow|
| MEMODE |1, writes| 1 |the same|
| NPMAX |250, writes|250; internal limit 10000|Same default; definition requires 2..10000|
| PNMAX, PPMAX |0, 10 a.u., entered|the same|PNMAX=0 gives one longitudinal point and causes a warning|
| NPN, NPP |50, 50, entered|the same; at PNMAX=0 the effective NPN is changed to 1|Same input default|
| USE_PEQUIV |0, writes| 0 |Same, validation limits 0..2|
|directional vectors|singles, omitted from the input|Normal z, transverse x/y|the same|
|Other flags momentum grids|false not to be published| false |the same|
| TAUDIR / TAUPATH |singles, omitted from the input|DATASET.tau.dir / current directory|the same|

Defaults that are not written do not mask different native selection.
NE=30 is an explicit starting number, not an automatic native POSANI default.
More basic than its isolated value is a match to an electron TAU directory.
Native PNMAX=0 is not the whole profile; left with notice, the range must the user
zvolit. Convergence NE, momentum grids and NPMAX is not confirmed by tests.

## Error documented: number of seats versus types

`posaniprep.f:65–69` allocates TAUQQ with dimensions NQMAX,NQMAX. Ready File
but writes TAUQQ in cycles **IQA/IQB=1..NT** (`:268–272`).
`posaniread.f:88–94` similarly reads TAUQQPOS in cycles 1..NT although dummies
the field has dimensions NQMAX,NQMAX.

Following annihilation routine runs through physical sites IQA/IQB=1..NQ and uses
TAUQQPOS(IQA,IQB) (`compton.f:1043–1079`). This field before reading is not
completely initialized.

- NQ>NT: part of the elements used is neither stored nor loaded at all and remains uninitialized.
- NT>NQ: preparation/reading treat other indices as locations; at NT>NQMAX
  jde o out-of-bounds array access. Enlargement NQMAX does not itself correct the meaning of these elements.
- NQ=NT removes this particular discrepancy of numbers; does not guarantee all other
assumptions or numerical correctness.

Common mixed occupancy or more equivalent sites of one type may have NQ≠NT.
Python's entry without potential doesn't determine this relationship. The restriction is therefore explicitly
in the definition and audit description; it cannot be promised that validation InputParameters system
with incorrect numbers automatically stops. Native Fortran has not been corrected.

## Other significant ALE

- POSANIPREP uses only s_1/2 positron partial wave
(`posaniprep.f:91–93`); POSANIREAD suppresses the imaginary part ZGPOS
  (`posaniread.f:67–70`). It is native approximation, not defaults input.
- Reader reads rows in order of types/places and does not verify structure conformity, radical
meshes or order of types with current potential. Exchange of prepared documents
no definition tests are detected.
- For directional vectors native the code normalizes but does not ortogonalize the specified pair.
The definition requires the final nonzero three-component vectors and perpendicularity.
Ignored PPVEC2 at PNVEC rejects; without PNVEC requires both transverse vectors.
- Shift flags shall be evaluated with priority; the definition shall not exceed one.
Slided even grids use ±INT((N-1)/2) and lose edges; in shifted
direction required odd size.
- EMD2D Multiplys the longitudinal vector twice (`init_mod_compton_pmesh.f:183–185`
and `:273–279`). Allowed only PNMAX=0, other cuts rejected.
- SQQ_RESOLVED not implemented for MD2D/EMD2D
(`compton.f:354–360`); the definition of the combination rejects the silent switch-off.
- At INSIDE_SPHERE and PNMAX=0 is PFEMAX virtually zero (1e-10).
- MEMODE=1 requires at least two interpolating points for division NP-1.
- The routine switches GF_CONV_RH to ZJ and SOLVER_FP to BS and to FULLPOT
  (`compton.f:185–194`). The supporting documents shall correspond to the convention used.

## MPI and output

POSANI is allowed in the final MPI control GEN. Common COMPTON routine has
Energy distribution and reduction (`compton.f:581–619`, `:1273–1296`).
The last energy gets the rank 0 in MPI_DISTRIBUTE corresponding to the core calculation
contribution to the rate 0. For each active contura, it shall be NE>=NPROCS.
mpi=True has a source justification, the actual MPI run has not started.

Outputs are DATASET_BND.mpa, DATASET_COR.mpa, DATASET_TOT.mpa;
2D outputs use extension _pan2D. The files distinguish the sum and the spin difference.
The interpretation of approximations and numerical standardisation shall verify the subsequent calculation.

## Tests

24 new tests: defaults, writing RDTAU and absence SPLITSS, automatic loading
POSANI, pickle, TAU paths, vectors, NONMAG, modes, both SPLITSS, wrong directions,
valid transverse pair, shifted grids, EMD2D, parameter ranges, NE, contour/BZINT,
RDTAU and not supported SQQ_RESOLVED s 2D.

Whole set InputParameters: **132 passed**, Python 3.12. Four warning:
two expected alerts PNMAX=0 in the general test defaults and two existing ones
NumPy 2.5 DeprecationWarning. Diff check and patch back app went through.
Numerical SPR-KKR has not been started. Tests do not verify the preparation of the positron file,
Match cache/potential, NQ/NT, physical value or numerical MPI result.

## Files

- src/ase2sprkkr/input_parameters/definitions/posani.py
- src/ase2sprkkr/input_parameters/test/test_posani.py
- 0031-posani.patch
