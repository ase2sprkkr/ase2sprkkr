# POSANITAU — audit SPR-KKR 9.7 / ASE2SPRKKR

Order 33. SPR-KKR: PUB9.7_260420_18h13.tar.gz.
ASE2SPRKKR refactor: 267512c4b8ccbbfd37af7ede898febc066b7e3ee.
Single patch 0033-posanitau.patch against this commit.

## Status and relevance

Missing definitions supplemented for standard periodic 3D ASA/FREL calculation
positronic anihilation speeds and service life over kkrgen. Common categories
audit: **works but has significant caveats**. Native task is active, numerical
the correctness is not verified by this audit.

POSANITAU is not the same as POSANI or POSANIPREP. Doesn't read DATASET_pos_dat.pan,
does not create its TAUQQ export and therefore does not accept an error NT instead of NQ from this export.
Works with electron potential/hustota, builds positronic potential,
integrates positronic density and evaluates its overlap with electron density.

## Connection, inputs and outputs

`gen.f:414–426` recognizes keywords after initial initialization of energy;
POSANI and POSANIPREP have priority before POSANITAU. The definition only writes POSANITAU.
`gen.f:431–449` calls POSANIPOT and prohibits updating EF.
`gen.f:458–476` requires GRID=5 and stores a copy of the electron density before its
positronic density replacement later.

The definition shall show:
- TASK.POSVC: BM / DR / IPM, volba ASA electron-positron correlation potential.
- ENERGY: one GRID=5 contour, NE=30, ImE=0.01 Ry, absolute EMIN/EMAX.
- TAU: standard POINTS and common k-strip parameters.
- MODE: FREL; NONMAG is allowed.
- SCF.FULLPOT: explicit true is rejected as this audit covers ASA.
Native FULLPOT a branch exists and is not marked off. Derived Mode
directly from the potential is not verifiable by validating itself InputParameters.

POSANILAMBDA does not have its own keywords for selecting the amplifier factor, angular
Mesh or normalization. The definition doesn't invent such parameters.

`gen.f:1508–1564` completes density correction, writes DATASET_positron.rho,
force graphs PLOTPRS and for FULLPOT also PLOT2DPR, then calls POSANILAMBDA.
Speed and service life are printed on standard output: rate [1/s], lifetime [s]
and overall lifetime also [ps]. Fields RATET have a distinction by type.

## Defaults and writing

| Parameter |Definitions and writing|Native effective behaviour|Assessment|
|---|---|---|---|
| GRID | writes 5 |GEN requires 5 but keyword POSANITAU does not add its own early initialization|explicit required option|
| NE | writes 30 |is not independent POSANITAU default before GEN|explicit start integration, convergence verify|
| ImE | writes 0.01 Ry |Common default 0.01 Ry|the same|
| EMIN | writes -0.2 Ry |Common default -0.2 Ry|same number, suitability for positron band not guaranteed|
| EMAX |not entered, is omitted from the input; can be changed explicitly|Sentinel converted to EF in INIT_MOD_ENERGY:182|the same effective behaviour; EF from the electron input may not be appropriate positron band limit|
| POSVC |BM, writes| BM |Taken native default|
| SPLITSS | false, is omitted from the input | false |same; the switch is outside defined workflow|
| MODE |without keyword, FREL| FREL |the same|
| SCF.FULLPOT | false, is omitted from the input |may be based on potential|input default native does not overwrite the mode; use ASA is subject to definition|
| BZINT / NKTAB |POINTS /common ASE default 250, entered|Common native initialization|explicit choice of standard workflow, not guarantee of convergence|

Deleted EMAX does not hide another native default. The user can set it explicitly;
the physical suitability of integration limits is highlighted by the validator. Relative EMINEV
or EMAXEV the definition rejects as an unaudited selection related to electron EF;
it is not a general native ban on relative energy parameters.

## Significant ALE and documented errors

### POSVC does not change the durability enhancer

POSVC controls correlation potential in POSANIPOT (`posanipot.f:240–252`).
However, thenihilation speed always uses POSANI_GFACTOR
(`posanilambda.f:92`, for FULLPOT `:250`). This routine has a single fixed polynomial
v rs (`posani_gfactor.f:38–39`) and does not read at all POSVC.
**POSVC=IPM therefore does not mean life calculation with amplifier factor 1.**
Changing potential while lifetime enhancement remains firmly established.
This is a significant limitation of physical interpretation, not a keyword mapping error.

### Normalization can mask an inappropriate integration area

ASA branch POSANILAMBDA first integrates the positron charge with NAT and CONC
(`posanilambda.f:60–71`) and normalizes density to total charge 1 (`:74–78`).
Zero/small charge leads to STOP but not zero does not verify correct belt area.
Normalizace sama nepotvrzuje vhodnost EMIN/EMAX ani konvergenci NE.

### Cache density suppresses a new energy calculation

`gen.f:500–504` reads DATASET_positron.rho and turns off RUNELOOP at success.
Change EMIN, EMAX, NE or POSVC thus does not have to mean new density: can be repeated
evaluate the old file. POSANIRHOIO does not verify physical origin of cache or conformity
with current parameters; u FULLPOT also skips headers without verifying their content.
The definition does not delete this file and does not automatically address its disability.

### Inversion of zero speed of individual type

Only total RATE>=1e-8 is tested at the end POSANILAMBDA. Then for everyone
the active type prints 1/RATET(IT), without separate control RATET(IT).
The total positive speed does not prevent zero speeds of one type; individual
the output can then contain infinity or cause a floating-point exception.
Total lifetime has control RATE, problem is in type listing.

### Conversion range and potential

POSANIPOT requires 3D and contains branches dependent on SCFVXC; unsupported
The functional may end up STOP. This definition does not overwrite electron VXC; and
Audit does not verify all its variants. Using compatible electron
potential and density remain a condition.
FULLPOT The calculation uses a fixed angular quadrature NGAU=20×20; its total
the calculation of shape functions and convergence are not audited in this task.

## MPI

Definition mpi=False: the audited path is serial. GEN performs final correction
a writing new density subject to MPI_ID=0 but subsequent fences and POSANILAMBDA
are not protected by the same condition in this block. Keyword POSANITAU additionally can
leave the internal TASK=NONE that passes the general final MPI control.
Without verification of density reduction and subsequent lifecycle can not be promised safe
parallel calculation. This is not explicit native ban MPI.

## Tests

15 new tests: defaults and omitted EMAX, automatic recognition task,
pickle, BM/DR/IPM, explicit EMAX, NREL/SREL/SP-SREL, FULLPOT, other GRID,
more contours NE, relative energy limits, CLUSTER and allowed NONMAG.
They control alerts on the integration area and the meaning POSVC.

Whole set InputParameters: **123 passed**, Python 3.12. Four warning:
two expected alerts POSANITAU in the general test defaults and two existing ones
NumPy 2.5 DeprecationWarning. Diff check and patch back app went through.
Numerical SPR-KKR has not been started. Tests do not verify lifetime, density, convergence,
cache, electron-positron models or numerical MPI results.

## Files

- src/ase2sprkkr/input_parameters/definitions/posanitau.py
- src/ase2sprkkr/input_parameters/test/test_posanitau.py
- 0033-posanitau.patch
