# OPTABS — SPR-KKR 9.7 and ASE2SPRKKR audit

OPTABS is **prepared in `tasks-workflow`**, not present in `develop`. Its native SPR-KKR 9.7 optical-absorption path is subject to the documented limits.

## Definitions

CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES; minimum SCF for refusal FULLPOT. RSP forbids this branch in FULLPOT. The audited range is FREL and BZINT=POINTS; native OPT_ABS also allows IBZINT=4 and 6 which are not certified here. NONMAG is permissible. Many-body optics is not covered by this audit.

TASK offers NEB, NOMEGA, OMEGA1/2, TEMPEL, NKT and Euler FRAMETET/FRAMEPHI/FRAMEGAM. ENERGY has GRID=0, NE=1 only for previous initialization and ImE for actual calculation. OPT_ABS itself creates two energy meshes, sets EMIN=EF−OMEGA2 and does not use an energy range from ENERGY. NEA and TASK.IGRID are read, but the subsequent mesh overwrites their meaning; the definition therefore does not offer them. OMEGA2 is **v eV**, temperatures in K, ImE in Ry, angles in degrees.

## Specific errors and restrictions

1. **Non-initialised KBT:** at present TASK and TEMPEL≤1 is not local KBT assigned. It then controls the thermal extensity and enters the calculation NKT_ACTUAL. So zero temperature in this source is not a safe input. The validation requires a final TEMPEL>1 K; the actual T=0 procedure would need a correction of Fortran.
2. **NOMEGA is not the actual frequency count:** OMEGATAB is first allocated according to the entered NOMEGA, then the NOMEGA is converted to NINT(OMEGA2/DELTA_E), i.e. the original NEB. Entering NOMEGA<NEB causes writing per allocated field. The definition requires NOMEGA≥NEB and describes the parameter as capacity. NEB≥3 prevents silent native clamp.
3. **Pairing energy:*** final and initial networks have step DELTA_E=OMEGA2/NEB. Transitions shall be counted only in conformity ERYDA -ERYDB with OMEGATAB in tolerance 1e−5 Ry. A general shift OMEGA1 can eliminate a large part or all transitions. The definition limits the normal workflow to OMEGA1=0 which native replaces DELTA_E. The output network thus does not start with zero but OMEGA2/NEB. Special comunsurable shifts are not offered.
4. Native replaces almost zero ImE value 0.001 Ry; definition requires positive final ImE and explicit writes 0.01 Ry. The width must be selected according to the desired spectrum.
5. Native contains MPI distribution of energy points and reduction of spectral fields. They were not numerically certified throughout the string, therefore mpi=False; this is not a claim that native MPI of the branch does not necessarily work.
6. The thermal extension increases both energy meshes, but the number of frequencies remains original NEB. Higher temperature, too soft net or large NKT can significantly increase memory and time of calculation. Converging NEB, k-points and ImE. The tensor framework and physical interpretation of the results require control for a specific structure.

## Defaults: real versus written

| Parameter |Effective native|ASE definition and writing|
|---|---|---|
| NEB | 30, minimum 3 |30, entered|
| NOMEGA |initial capacity 100, actual number then NEB|100, written up; capacity validation|
| OMEGA1 | 0 eV → DELTA_E |0, entered|
| OMEGA2 | 8 eV |8, entered|
| TEMPEL |Reader default 1e6 K at TASK|300 K, written down|
| NKT |Locally 6, reader at TASK overwrites to 99|6, entered|
| ENERGY GRID/NE |own optical network created again|0/1, written only for initialisation|
| ImE |Almost zero replaced 0.001 Ry|0.01 Ry, written|
| EMIN/EMAX |EMIN converts OPT_ABS, its network determines TASK|not in the definition|
|Euler angles|0, frame changed only with the key present|Single, omitted|

TEMPEL=300 and NKT=6 are explicit deviations, not hidden defaults. They are more suitable for a normal start than native 1e6/99; the experimental temperature must be set by the user. Meaning of the Reader verified in scan_input.f: default different from the Sentinel 9999 overwrites the value of the missing key. Native initial DATA value of Euler angles agrees with their unwritten status.

## Verification scope

The prepared Python definition has test coverage for serialization and selected validation rules. Counts from old standalone patch runs are not the test results of the current `tasks-workflow` branch. No numerical native SPR-KKR reference calculation is claimed here.
