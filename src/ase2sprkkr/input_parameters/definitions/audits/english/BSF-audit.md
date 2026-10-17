# 41. BSF audit of common definition BSFEK/BSFKK

Already supported task; revised common definition. Preliminary: **has fully meaningful***, with limitations defaults and native geometry. The audit focuses on the option E-K/K-K, numbers and geometry; it is not a numerical verification or a complete check of all tables KPATH for all Bravais grids.

## Opravy

- NK1/NK2 are now Integer(min=2). Native BLOCHSF when parameterization K-K the map divides NK1−1 and NK2−1 (blochsf.f:939–940); the previous plain int allowed division by zero. NK for E-K path has minimum 2, NE minimum 1.
- GRID default is dynamic: for E-K remains 3, for K-K is 0. Earlier GRID=3/NE=1 at native initialization results in a single point integration weight not different. The change is actually being written out. After explicit setting GRID, the user can of course deviate from default.
- K-K validation rejects endless/NAN directions and zero or collinear plane. The check takes place against KA, not only by K1/K2 alone.
- Fixed description K1/K2: end points; native subtracts KA from them. If the relative directions are not perpendicular, the other orthogonalizes. Thus, such an assignment does not have to be drawn by the originally intended parallelogram (blochsf.f:438–455).

## Defaults: written versus omitted

E-K: preserved NE=200, GRID=3, NK=300, ImE=0.001 Ry, EMIN=−0.2 Ry and EMAX=−1 Ry. **EMAX=−1 is real energy, not sentinel for EF.**Native init_mod_energy uses substitution EF to EMAX<−99. The default range −0.2→−1 is therefore decreasing and surprising for the normal spectrum. This is an explicit warning ASE default, not a correction of an ignored parameter native reader. It is appropriate to set its own range, possibly in the future default, to change to a physically understandable range to EF; it is not possible to determine the appropriate absolute top energy without potential. The default is now maintained according to the audit assignment defaults.

Native NK at reading E-K is 800, ASE writes 300. This is a valid explicit deviation; 800 may be better for fine paths, but the final number is divided between segments and may be adjusted according to k-table capacity. Neither 300 nor 800 is a guarantee of convergence.

K-K: generated energy EMIN=0.7 Ry is explicit ASE option, EMAX copy EMIN; it is not automatically EF. Written NE=1 and newly GRID=0 correspond to one point. The correct energy or relative energy EMINEV should be set to find the Fermi area. Unwritten KA native completes with zero start; Python validation uses the same beginning. NK1/NK2/K1/K2 do not have an imaginary defaults and must be supplied by the user.

Existing mpi=True left. BLOCHSF has a parallel energy process and GEN allows it; however does not verify has a numerical conformity MPI result. MODE/CPA/STRCONST/SITES remain common sections unchanged.

## Other restrictions

E-K is selected by the actual energy count NE>1, K-K at NE=1; names BSFEK/BSFKK are ASE aliases and hints, not two different native dispatch keys. TASK with writes as BSF, CONTROL ADSI as BLOCHSF. Existing validation warns in conflict with the required alias and NE.

In E-K mode native divides points by segment length, with a minimum of two points per segment. Degenerate explicit segments and zero total length can lead to problematic distribution. Their complete control, together with the geometry of each KPATH/Bravais combination, is not solved in this patch.

For 2D native in K-K the branch zeros the z-folder K1/K2 before deduction KA. Input-only check does not know the potential dimension, so the plane valid according to specified 3D vectors can degenerate after native transcription. For 2D, explicit consistent plane input is required.

At native reading TAU cache can take BLOCHSF energy from cache; the consistency of external files is not verified here. This audit does not meet Advanced BSFSPN/BSFORB/BSFETD/BSFALF orEPHOT.

## Verification

11 new tests, whole set InputParameters **119 passed**, two existing NumPy DeprecationWarning messages. Tests control K-K GRID=0, equality EMIN/EMAX, round-trip/pickle, conservation E-K defaults, non-zero KA with endpoints, minimum numbers, degeneration, colinerite and NN. Existing tests continue to verify the sharing of aliases, regime change and basic grammar.

Command: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

Sources: blochsf.f:220–268,310–357,419–455,939–940; init_mod_energy.f:182; epath.f (GRID=3 single point branch). ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Single patch only for BSF. No numerical SPR-KKR run.

## Current integration status (2026-10-17)

The source audit above remains a historical record of the original definition baseline. On branch `tasks`, the shared `bsf` definition now models the TAU directory as an optional reusable cache and explicitly manifests the native persistent results `{dataset}_BLOCHSF.bsf` and `{dataset}_BLOCHSF_spol.bsf` for `BSF`, `BSFEK`, and `BSFKK`. The aliases do not get separate filenames because all three select native `ADSI=BLOCHSF`.

This file-contract integration is present at commit `032d932b2e05aa3a540ce232c98522b1fb7065c6`. The historical test count above has not been rewritten; the later manifest regression test covers all three public names.
