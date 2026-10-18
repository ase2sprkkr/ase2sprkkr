# 42. SCF — audit ASE2SPRKKR proti SPR-KKR 9.7

Already supported task. Fixed Task Definition SCF; preliminary **has fully meaningful**. The audit controls the main SCF management, algorithms and defaults. It is not a certification of all FULLPOT, many-body, cluster or non-linear branches. The Numerical Solver hasn't been launched.

## Corrections to actual errors

1. ASE offered ALG=PULAY but scf.f:427–437 accepted only BROYDEN2, ANDERSON and TCHEBY. PULAY leads to SCFALG not found. Task copy ALG now offers only these three names; PULAY is not automatically mapped to ANDERSON because such an alias native reader does not have.
2. TCHEBY s FULLPOT native explicitly rejects; validation added. FULLPOT can also be activated by potential, which input-only validation itself does not detect.
3. NITER and ISTBRY are positive integers. ITDEPT limited to 3..40: SCFBROYDEN converts values below 3 to 40 and limits values above 40 to capacity. The definition now does not pretend that these requirements will be implemented accurately.
4. MIX, MIXOP and TOL must be final positive values; EFGUESS final energy. The upper limit MIX was not arbitraryly set at 1 because overrelaxation is not merely a type error.
5. NE is positive; single point GRID=2/3 is rejected due to native integration defects. Other advanced grid combinations are not fully controlled by this patch.
6. Description EFGUESS specified: This is the initial estimate EF in Ry and the omission of the initial search when creating starting potential. It does not prohibit all subsequent updates EF in SCF.

Corrections shall be made through the drawing copies SCF and ENERGY to avoid changing other definitions using the common sections. The common section ALG may still contain PULAY in other draw definitions; this is not a confirmation of its native support.

## Defaults and their writing

Normal values NITER=200, MIX=0.2, TOL=1e−5, ISTBRY=1 and ITDEPT=40 are written to the input and correspond to initialisation in scf.f:360–371. BROYDEN2 corresponds to MOD_SCF default and is default of the draw definition. VXC options, their default and metadata have not been modified; MOD_SCF initializes VWN, but specific potential and reader can continue to influence the status input. Energy default GRID=5/NE=30 corresponds to SCF to INIT_MOD_ENERGY; EMIN=−0.2 is written ASE option EMAX is omitted from the input and leaves native EF.

**MIXOP is an important unwritten case:** scf.f first sets SCFMIX=0.2 and SCFMIXOP=SCFMIX, then reads MIX and MIXOP. At MIX=0.5 and single MIXOP it is therefore not possible to claim that MIXOP will be 0.5; for normal non-export procedure, 0.2 remains. ASE does not impose another hidden default, MIXOP is single and is omitted from the input. The description of this dependence is now explicitly highlighted. To complete explicit MIXOP is better if the user wants the same value as changed MIX; automatic change default would change native behavior.

EFGUESS is not recorded without input; native is determined by the presence of the key, not by any replacement energy. FULLPOT=False is omitted from the input; it is not a promise that the full potential will change to ASA. TAU/CPA/STRCONST/SITES/MODE are common sections, here without interference.

## Other restrictions and branches

- Native EXPORT_HARTREE_POTENTIAL sets NITER=1, MIX=1 and skips reading VXC/NITER/MIX. The table of normal defaults is not transferable to this branch.
- ANDERSON mapuje IMIX=7, BROYDEN2 IMIX=4. TCHEBY uses a separate procedure; the algorithm name is not synonymous with the same numerical mix.
- Exhaustion NITER is not evidence of convergence. It is necessary to control the real convergence state, not just the existence of new potential or successful input parsing. This patch does not adjust output parser or calculator lifecycle.
- Native offers further management as SCFRESTART, MIXRHO, SIM, MAGCOUPL, SEMICORE and 2D Layer charge parameters. They are not completed here and their entire workflow is not covered by this audit.
- mpi=True preserved. SCF has parallel solutions and initial search; numerical match MPI/serial in this audit has not been tested.

## Tests

15 new tests; whole set InputParameters **123 passed**, two existing NumPy DeprecationWarning messages. Verified by three supported algorithms in round-trip and pickle, executable/MPI option, maintenance of main defaults, rejection PULAY, invalid numbers and history, NN mix, inexplicable tolerance, infinite EFGUESS, zero NE, single point GRID=3 and TCHEBY+FULLPOT.

Command: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

The sourcee: scf.f:360–437; scfbroyden.f:54–55,105,121–123; scfinitpot.f:443–455; mod_scf.f:33; init_mod_energy.f:76–79. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Single patch only for SCF.
