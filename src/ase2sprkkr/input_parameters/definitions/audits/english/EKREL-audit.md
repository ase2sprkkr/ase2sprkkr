# EKREL — audit of the current implementation

The integrated `ekrel.py` input definition exposes the native `kkrgen` E(k) and Fermi-surface modes. The numerical and geometry restrictions below are derived from SPR-KKR 9.7 (`PUB9.7_260420_18h13`), not from the presence of a Python validator.

## Parameters TASK and defaults

EKREL calculates E(k) for arranged systems, or seeks Fermi's surface. INITVAR recognizes keyword EKREL and turns off normal RUNELOOP (`initvar.f:459–467`); GEN then calls its own EKREL. The routine has its own energy and loop.

| Parameter | Native default / interpretation |Definitions|
|---|---|---|
| NK |1000, number of points required| >=2, default is omitted from the input |
| EMIN |0.01 absolute Ry| default is omitted from the input |
| EMAX |1.0 absolute Ry| default is omitted from the input |
| NE |501 search energy| >=2, default is omitted from the input |
| F-SURF | false | flag |
| ESLICE |EF, at F-SURF only|None and not writing delegate EF; explicit value in Ry|
| KPATH | 1 |shared travel descriptions BSF/KDIRTAB; default is omitted from the input|
| NKDIR |0, use KDIRTAB|0..9; non-zero means manual path|
| KA1..KA9 / KE1..KE9 |during the manual journey both ends mandatory|three components in cartesian units 2*pi/a, to be completed explicitly|

Source: `ekrel.f:92–154`. EMIN, EMAX and NE are read from **TASK**, not ENERGY. The limits are absolute Ry. EMINEV/EMAXEV in general ENERGY sections do not replace these TASK limits.

The manual path is entered for example TASK.NKDIR=1, TASK.KA1=[0.,0.,0.], TASK.KE1=[1.,0.,0.]. Specific numbered options are created by a short declarative loop for native range 1..9. The validation checks completeness, finiteness, non-zero segments and unused indices. There is no hidden shortening NKDIR>9: definition rejects it while native routine cuts it to 9.

KPATH uses the same kdirs table as BSF, therefore the definition takes over the existing description dictionary. The availability of a specific number depends on the Bravais grid; the actual validation of input without its structure does not verify. The range of supported grids is not only extended by adding this definition.

## Fermiho povrch

F-SURF creates a Celobrilouin k grid for BXSF output. TASK.NK is default density requirement, not guaranteeing the exact result of points: the source recounts it, including double border points (`ekrel.f:270–305`). Manual k-path does not control the resultant sampling here. KDIRTAB is called before branches F-SURF, and therefore also for Fermi's surface can apply its grid limitation.

ESLICE remains EF, which may lie outside default search interval 0.01..1 Ry. Without potential, this relationship cannot be verified. Explicit ESLICE validator checks against TASK limits. The source actually permits equality with limits, although its error meshage claims strict inequality (`ekrel.f:132–134`); the definition corresponds to a condition `EMIN <= ESLICE <= EMAX`.

## Joint sections and modes

The definition has CONTROL, TAU, MODE, STRCONST, ENERGY, TASK and SITES; executable=kkrgen. `mpi=True` corresponds to the actual distribution of own to the loop using MPI_SET_WORKERS_ARRAY and IPROCK (`ekrel.f:435–448`). MPI or numerical startup was not verified here.

ENERGY has explicit GRID=1, NE=1, ImE=0, EMIN=-0.2 Ry; EMAX=None is omitted from the input and delegate EF. ENERGY serves a common initialization and does not constitute a search interval EKREL. INIT_MOD_ENERGY requires a non-zero initial number of points before later task-specific adjustments. The multiple ENERGY.NE definition draws attention to; the number of own search energy controls TASK.NE. SITES.NL=None is omitted from the input and preserves potential context.

Not recorded TASK defaults correspond to native routine. Positive EMIN=0.01 Ry is native choice, not forced choice ASE2SPRKKR; due to the poles and behaviour of E=0 I leave it. For a specific energy window it is necessary to set both TASK limits. Unwritten ESLICE must remain contextual, cannot be replaced by constant numeric EF.

**CPA alloys are natively rejected** (`ekrel.f:81–82`); the source recommends Bloch spectral function. The input definition without potential cannot determine NCPA, therefore the restriction documents and relies on native control. Does not create unnecessary CPA section. FULLPOT does not have a general ban in check_input for EKREL; nor does GEN prohibit DMFT/LDAU for this task. EKREL has its own reading DMFT self-energy (`:384–399`). The definition therefore does not invent their prohibition, but numerical compatibility has not been verified here.

## Mistakes and risks in the source

1. **NE>=2 is checked only for F-SURF, although division NE-1 is general.** `ETOL=(EMAX-EMIN)/(NEK-1)` (`ekrel.f:179`) is also performed for ordinary E(k). The definition therefore requires TASK.NE>=2 always.
2. **Source generally does not control the current energy window.** At EMIN=EMAX, ETOL=0 is created and the subsequent check of passage E=0 may divide by zero; moreover, reverse limits violate the search assumption. Definition requires final EMIN<EMAX.
3. **Incorrect meshage ESLICE does not correspond to the test.** Text states EMIN<ESLICE<EMAX, but the condition rejects only values strictly outside the interval. Equality is not natively forbidden and neither does the definition forbid it.
4. **NK is not the exact size of the resulting path.** Segments are divided by lengths and each receives at least two points (`:204–227`); F-SURF has its own grid. The definition describes NK as a requirement, not a guaranteed number.
5. **Handless zero segments or total zero length are not properly treated.** Dividing by SUMD in the source can divide by zero. The definition controls specific manual segments; availability and geometry of preset paths verify SPR-KKR according to grid.

The transition over E=0 source tries to solve by reducing NEK by 1 (`ekrel.f:181–187`). The actual amount of energy does not therefore have to correspond to the declared NE. Definition of this behaviour does not fix.

## Output contract checked against native SPR-KKR

For ordinary E(k), `kkrmain.f:1202–1208` opens `TRIM(DATSET)//'.bnd'` using `WRHEAD(...,'DISPERSION',...)`; this is the data written by `ekrel.f` on unit 10. With `F-SURF` enabled, `ekrel.f:697–703,785–790` **also** writes `TRIM(SYSTEM)//'.bxsf'`. Its basename is derived from `SYSTEM`, not from `DATASET`. Native `kkrmain.f:239–245` appends nonempty `ADSI` to `DATSET`, so the `.bnd` name contains ADSI when set, unlike the BXSF filename. The current Python task manifest declares only `{dataset}.bnd`; thus the optional BXSF is absent and the `.bnd` declaration disagrees with native naming for nonempty ADSI. These are manifest issues, not a finding that the native computation fails.

## Verification scope

The input definition and its focused tests are present in `src/ase2sprkkr/input_parameters/`. This audit does not claim a new test-suite run or numerical SPR-KKR validation.
