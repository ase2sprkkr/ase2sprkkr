# 36. CHILANDAU — SPR-KKR 9.7 / ASE2SPRKKR

So far missing task for Landau's orbital magnetic susceptibility. Added definition for **serial periodic ASA FREL procedure**. Preliminaryly: the serial branch works with essential limitations; **MPI the branch is broken**. Numerical calculation not carried out. It's not a joint final audit.

## Definitions and validation

Executive kkrchi; native INITVAR recognizes CHILANDAU, RSP calls SUSC_LANDAU. Added CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK and SITES; minimum optional SCF section allows to refuse FULLPOT, which RSP explicitly prohibits.

Required BZINT=POINTS (native IBZINT=2), one GRID=5 contour without SPLITSS and positive energy. SUSC_LANCALC otherwise at NEPATH!=1 or GRID!=5 silently returns without calculation. It returns at KMROT!=0; KMROT is a property native status/potential and does not detect it. FREL is a limitation of the audited procedure: angular and wavefunction routines use relativistic CGC and spinor channels, not found explicit STOP for all other modes. NONMAG is permissible © it is not a spin susceptibility requiring a magnetic base state.

TASK contains CLURAD, NSHLCLU, SELECT, EXCLUDE, QREPONLY and PRINT. SELECT/EXCLUDE are mutually exclusive; radius is final non-negative. NSHLCLU including zero takes precedence over radio. The definition does not offer all possible native calculation modifications.

## Defaults

**Written out:** GRID=5, NE=30, ImE=0.01, EMIN=−0.2 and CLURAD=2.2. These are explicit practical choices ASE; number of points and radius require convergence. Native INIT_MOD_TAUIJ_STAR does not have a usable default radius: without CLURAD and NSHLCLU ends STOP Therefore, the published radius is better than its omission, but the value 2.2 is not a physical guarantee of convergence.

**Unwritten:** EMAX is actually left to native EF; QREPONLY=False agrees to MOD_SITES.TREAT_QREP_ONLY=False. NSHLCLU is not forced to zero because its mere presence changes the choice of cluster. PRINT is single: SUSC_LANDAU after initialization sets IPRINT=0 and derives CHIPRINT from it. FULLPOT=False is not written; if the potential itself activates FULLPOT, native RSP it still rejects it. Python default does not switch full-potential potential to ASA itself.

## Specific errors and limitations native sources

1. **MPI: uninitialized field.** SUSC_LANCALC allocates SAVE angular field only at IECURR=1. SUSC_LANDAU divides energy between wounds. The other energy ratio does not have these fields allocated (susc_lancalc.f:66–88; susc_landau.f:193–214).
2. **MPI: reduction after the end of the field.** CHILANLT has a shape (0:NLMAX,NTMAX,2), but MPI_SUM_R starts at CHILANLT(0,1,2) and transmits the length of both blocks. Since the beginning of the second block, only half of the required length (susc_landau.f, final MPI total) remains. This is an independent memory error. mpi=False is therefore a deliberate default definition; explicit enforcement MPI does not correct this error via calculator.
3. ** Continuous outputs.** Condition IECURR.LT.-NETAB(1) in SUSC_LANCALC is always false for positive energy indices. The output of the susceptiblite is printed at each energy, not at the end. V serial run is required to use the last block; earlier blocks are only partial integrals.
4. **Reading TAU branch.** RDTAUFIL opens and reads the text file header, but the energy loop lacks the corresponding reading TAU and the calculation TAU is skipped at RDTAUFIL. This reuse procedure is not offered or certified by this definition.
5. The calculation uses the final sum of the neighbours and the specific Cartesian component. Convergence, geometric range, multiplication at QREPONLY and interpretation for general anisotropics need to be verified. Lower dimensionality and non-standard rotational states are not covered by this audit.


## Verification

12 new tests: discovery/defaults/serial executable, round-trip, pickle, NSHLCLU=0 and selectors, NONMAG, non-audited mode rejection, FULLPOT, BZINT, GRID, SPLITSS, invalid radius and NE, reciprocal exclusiveity SELECT/EXCLUDE. Full set InputParameters: **120 passed**, two existing NumPy DeprecationWarning messages.

Command: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

The sourcee: initvar.f:453–454; rsp.f:134–147; susc_landau.f; susc_lancalc.f; susc_laname.f; susc_lanme.f; init_mod_tauij_star.f:150–180; mod_sites.f:73. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Single patch only for CHILANDAU.
