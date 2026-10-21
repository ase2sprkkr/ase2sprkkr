# CHIXAS — current definition audit

The integrated `chixas.py` definition describes field-induced XAS dichroism in the serial ASA/FREL path. External physical susceptibility data remain mandatory, and no numerical SPR-KKR validation is claimed.
## Mapping and definitions

INITVAR recognizes CHIXAS; RSP calls CHIXAS (susc_xas.f). This is not an alias of normal XAS: SUSC_XASCALC explicitly has a purpose FIELD-INDUCED MCXD. The absorption type is entered TASK.ITXRAY, not XAS TASK.IT. Core shell is entered CL.

The definition contains CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES and minimum SCF section for refusal FULLPOT. RSP CHIXAS in FULLPOT explicitly prohibits. FREL, POINTS and one GRID=3 are the audited spectroscopic range chosen, not the explicit prohibition of each other native network. Native CHIXAS re-establishes the energy path and sets NEPATH=1; SPLITSS therefore we do not offer. NE>=2 prevents a single point GRID=3 with an undefined weight.

TASK also contains NEF (numbered dictionary → NEF1, NEF2, ...), CHISPN, CHIVV, ABSUNIT and PRINT. Not created imaginary physical defaults. Without NEF and both susceptibility validation fails; parser/skeleton itself can be created and then supplemented.

## Mandatory external data and validation

- NEF for each potential type is l-resolved DOS to EF. Native if missing NEF for any type ends STOP. The number of ingredients shall not exceed native NLMAX and shall correspond to both l-channels and standardisation used SUSC_GAMMA. This input-only validation does not detect from potential.
- CHISPN and total CHIVV are per-type susceptibility in 10^-6 cm^3/mol. SUSC_XASINIT converts them into atomic units. Total CHIVV divides by two and divides equally between spin channels.
- The definition requires the same length CHISPN/CHIVV and continuous NEF1..NEFn, controls the final values, non-negative DOS and ITXRAY to that extent. Actual NT is not read here in potential; the conformity of the input fields is not proof of each other's compliance with potential. The numbers in the tests are only grammar fixes, not material data.
- Spin-resolved CHIVVUP/CHIVVDN are not made available to avoid conflicting combinations with total CHIVV. Native parser temporarily overwrites global NT by the number of values found. Bad lengths can give zero fallbacks and other condition; it is necessary to deliver exactly all potential types.

## Defaults: co is written to the input

| Parameter | Native |ASE definition|
|---|---|---|
| ITXRAY |explicitly required|1, written down; need to choose the correct type|
| CL | 2P |2P, written|
| GRID / NE |for prefix CHI: 5 / 30|3 / 180, written for energy differentiated spectrum|
| EMIN | EMINDEFAULT |−0.2 Ry, written|
| EMAX | EF |4 Ry, written; explicit spectroscopic option|
| ImE | 0.01 Ry |0.01 Ry, written|
| ABSUNIT | false |false, not written|
| NEF, CHISPN, CHIVV |NEF mandatory, susceptibility without data can end up zero|without defaults, required by validation|
| TASK.PRINT | input IPRINT |not entered, not written|

Listed GRID/NE/EMAX are deliberate deviations from prefix CHI default. They are more suitable as the start of spectral mesh than integration only to EF, but the range, resolution and width must be set for a specific edge. EMAX=4 is not proof that it covers the required spectrum, and −0.2..4 Ry may unnecessarily include a large interval. Native default cannot be blindly recommended as better only because it is CHI prefix.

TASK.OP is omitted from the input: CHIXAS then sets local CHIORBPOL=BROOKS. This variable is not automatically the same as the general MODE.OP determining the correction of potential. Unwritten ABSUNIT accurately preserves native behavior. RDTAUFIL begins with sources false; restart and read old text TAU/TAUTAU caches are not offered by this definition.

## Next ALE

Native CHIXAS is not a self-consistent derivation of susceptiblite from material: first it needs NEF, then it calculates GAMMA and then uses the entered CHI in induced spectroscopy. Connection with real external data and core waves has not been numerically verified.

RDABXC/LASTITRXC, TTTFMT, surface/weak, cache and spin-resolved orbital procedures are not audited. At POINTS native in the local k-loop branch switches IBZINT from 2 to 6; this cannot be confused with the user normal TB workflow only by number. The definition sets mpi=False; the parallel string has not been certified; it is not a proven MPI error.
