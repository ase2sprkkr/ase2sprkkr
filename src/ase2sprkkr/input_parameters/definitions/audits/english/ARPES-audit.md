# 43. ARPES audit SPR-KKR 9.7 and ASE2SPRKKR

Already supported task in kkrspec. Fixed existing definitions; preliminary **have fully meaningful***, with important geometry limitations and other draw modes. The audit focuses on modern SPEC_INPUT reader, energy units, angular/to network and their defaults. It is not a complete numerical certification of one-step photoemission, all surface barriers or aliases AIPES/SPLEED/BAND.

## Opravy

- ENERGY.EWORK_EV correctly described as **output work in eV**, not internal potential. SPEC_INPUT reads and converts /RY_EV; sets EWORK in the appropriate branch when entering.
- NT, NP and NK1..NK4 are positive integers; NK1=1 remains permissible. Reader has NT -1/NP -1 protects the condition >1 for one for the explicit branch. Required minimum 2 would unduly eliminate the valid native procedure.
- K1 requires NK1, at the same time K2 also NK2. Negative and zero numbers previously could lead to empty or invalid-sized networks. Angular range compliance check and NT/NP remained in place.
- For ARPES rejected K3/K4/NK3/NK4: in supplied SPEC_INPUT reader these keys are not read anywhere. In the public Python definition, they remain for reverse parsing, but validation does not pretend their effect.
- Final non-negative EWORK_EV and both IMV_*_EV, positive final EPHOT; final electron angles and k-vectors. Zero control is not arbitrary.
- N_LAYDBL must have exactly two non-negative integers. N_LAYER and NLAT_G_VEC are positive integers. INPVER/STRVER limited to 0/1. At INPVER=0, the new control of the modern reader is omitted because the data comes from an old external input.inp.

## Defaults: written versus omitted

| Parameter | ASE default |Comparison and writing|
|---|---|---|
| EWORK_EV | 4.2 eV |written; overwrites workfunction from SCF/Potential|
| IMV_INI_EV / IMV_FIN_EV | 0.05 / 2 eV |written; native photoemission start VIL/VIH is 0.1/1 eV|
| NE | 300 |written; native photoemission start EMESH=200|
| EMINEV / EMAXEV |−8 /+5 eV versus EF|written; explicit spectral range, not general physical limit|
| NL | 4 |written; explicit angular cutoff|
| NT / NP | 1 / 1 |written, consistent with zero angle range|
| Electron THETA/PHI |0..0 for ARPES|task-specific default; SPLEED has 45/270 degrees|
| EPHOT | 25 eV |written; explicit experimental option|
| N_LAYER | 50 |written; native photoemission start 30|
| NLAT_G_VEC | 57 |written; native start 30|
| N_LAYDBL | 10 / 10 |written; native has a typo LANZ1=10 twice instead of initialization LANZ2|
| INPVER / STRVER | 1 / 1 |written down; modern input and structure generation|

These listed deviations have been left. N_LAYER=50 or more reciprocal vectors can improve convergence at the cost of time; it is not possible to recommend native smaller numbers as better without the system. Output work and control are sensitive material/experimental parameters: the value of potential can be more suitable than generic 4.2 eV, but its use requires this explicit key to be omitted. Default is not just formatting settings.

Expert photon polarization and spin-amplitude defaults are generally not recorded on entry. This audit does not change their logic and does not constitute a complete certification of all combinations POL_P/NPOL/ICIRC/IDREH or POL_E/POL0/POL0L. Reader rejects some combinations based on the presence of alternative keys; numerical matching default value alone is not sufficient to replace writing and omission.

## Specific native errors and restrictions

1. **IMV_INI v Ry has a typo:** spec_inputs.f:1257 calls SECTION_SET_REAL('IMV_INI',VIH,...), i.e. writes to the final state variable instead of VIL. ASE offers IMV_INI_EV, the branch of which writes correctly to VIL; the wrong Ry alias is not added. VIL/VIH reader subsequently converts to internal spectroscopic input in eV.
2. **Double LANZ1 initialization:** reader sets LANZ1=10 twice. Explicit N_LAYDBL={10,10} in ASE sets both variables and is therefore better than blank missing fields with expectations of both native defaults.
3. K1/K2 in this reader **sun vectors are added to KA**, not endpoints as in BLOCHSF. The number of points determines NT/NP and the given vectors are used for k-mesh. Names cannot be interpreted uniformly between tasks.
4. Native photon convention gives normal incidence THETA=180 degrees. It is not possible to automatically transfer geometric meaning THETA from other tasks. The specifics of the crystal frame and sample tilt must be verified for a specific surface.
5. Validation does not verify existence and content of external input.inp/struc.inp at INPVER/STRVER=0, consistency of surface structure with potential and convergence of Layer Doubling and reciprocal Cutoff.
6. Executive kkrspec and existing parallel settings left; numerical MPI/serial calculation not performed. An audit performed here ARPES is not a confirmation of all other TASK options sharing the same definition.

## Tests

14 new tests; whole set InputParameters **122 passed**, two existing NumPy DeprecationWarning messages. Verified defaults, executable, round-trip/pickle, valid single point scan, missing NK1, invalid count/length, negative workfunction, NNdamping, zero photon energy, invalid INPVER and ignored third/fourth directions. Existing angular scan tests remain successful.

Command: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

The sourcee: spec_inputs.f:1198–1330,1420–1509,1660–1783; spec_rmbarpot.f (SURF_BAR). ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Single patch only for ARPES. Without numerical SPR-KKR run.
