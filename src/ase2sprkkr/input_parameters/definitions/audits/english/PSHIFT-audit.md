# PSHIFT audit and new definition InputParameters

Background: attached SPR-KKR 9.7 (PUB9.7_260420_18h13) and current ASE2SPRKKR `refactor`, verified fetch, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`. PSHIFT does not have a default module in this revision. Single patch adds `definitions/pshift.py` and `test/test_pshift.py`; does not need previous patch.

## Entry

`initvar.f:173–174` recognizes TASK PSHIFT, `gen.f:1426–1431` calls the routine PSHIFT in KKRGEN. `pshift.f` does not read any task-specific parameters; it is not ITSEL, L, PRINT or the output switch. The routine counts all atomic types, uses CONTROL.PRINT and cutoff from potential/SITES.NL. The definition adds common CONTROL, TAU, MODE, STRCONST, ENERGY, TASK and SITES.

Uses `kkrgen`, `mpi=False`: loops PSHIFT do not contain the distribution of work between the wounds and output names are shared. MPI or numerical run was not verified.

| Option |Behavior / default SPR-KKR|Definitions|
|---|---|---|
| TASK | pouze keyword PSHIFT |no imaginary own parameters|
| MODE.MODE |missed means full relativity|FREL/None, NREL, SREL and SP-SREL; all have source branches|
| ENERGY.GRID |takes the fair part of the first ETAB path|explicit 1, direct equidistant real path|
| ENERGY.NE |NETAB(1), Global Initial Value 0 without task-specific default|explicit 300, positive number|
| ENERGY.ImE |reset to 0|default 0; non-zero notification|
| ENERGY.EMIN | -0.2 Ry (`init_mod_energy.f:49,74`) |explicit -0.2 Ry, identical default|
| ENERGY.EMAX |at missed actual EF|explicit 1.0 Ry, see below|
| ENERGY.EMINEV/EMAXEV |limit relative to EF in eV|common pair interface ENERGY preserved|
| ENERGY.SPLITSS |reset to false|false, when warning is switched on|
| SITES.NL |cutoff from potential|None, is omitted from the input; context remains SPR-KKR|

ENERGY.ImE and SPLITSS overwrites `init_mod_energy.f:315–320`. PSHIFT only works with the first way and the real part of energy (`pshift.f:66,79,123,368`). The second way of definition is pointed out. Other GRID is not natively forbidden but can generate repeated or reversed real energy; this is inappropriate for phase expansion. The definition therefore recommends GRID=1 and only warns other ways.

### Differences defaults

The upper limit **EMAX=1.0 Ry** is intentionally different from native omissions that are used EF. PSHIFT looks for resonance under absolute energy 1 Ry; such window is more useful for diagnosis than stopping at EF, which can miss unfilled resonances. Therefore, here native default I do not automatically consider to be better. This is an explicit option ASE2SPRKKR, not an incorrectly modeled omitted default. For another window, the user sets absolute or relative limits.

NE=300 and GRID=1 are also explicitly written practical defaults. Native NETAB starts at zero (`mod_energy.f:43–44`), init_mod_energy has no specific PSHIFT default number of points; this is not appropriate to take over. Deleted SITES.NL remains contextual and therefore preserves native behavior.

## Modes and outputs

- FREL: direct DIRBS calculation of phase shifts and logarithmic derivative; output panels include nonmagnetic channels and selected magnetic channels.
- NREL/SREL: one spin channel, phase shifts and the real/imaginal part of t-matic.
- SP-SREL: two spin channels plus an estimate of exchange resonance separation.
- The output energy axis is **absolute eV**, `Re(ETAB)*RY_EV`, not `(E-EF)*RY_EV` (`pshift.f:79`). The input EMINEV/EMAXEV does not change this meaning output of the axis.
- Relativist branch calls DIRBS directly without SOCTL (`:396–400`); general MODE.SOC its SOC does not flatten. Definition at explicit selection in FREL warns. It does not replace this branch with SOCPAR, which is actively rewriting SOCTL.
- Output graphs remain without specialized parser; patch adds InputParameters and uses general output reader.

FULLPOT explicitly rejects `check_input.f:80–82` and `gen.f:1428–1429`. It is detected from the potential; the input alone cannot reliably validate it without loading the potential. The restriction is documented and enforced by SPR-KKR. DMFT and LDAU rejects GEN (`gen.f:211–222`); definitions also. There is no non-existent ban NREL/SREL/SP-SREL.

## Mistakes and risks in the source

1. **Resonance split can use missing resonance.** V nonrelativistic branch is RESONANCE common symptom for both spins (`:213–232`). Just find one resonance; the other ERES can remain -9999.90 Ry and still be deducted. The same problem is in the relativistic branch (`:519–545`): a common symptom for different m_j and DESO additionally uses ERES of the previous kappa channel without checking for it. Numerical fracturing may not be a physical result.
2. **Resonance detection does not control actual pi/2.** Test only current phase >pi/2 (`:221,532`) and interpolate between two points. If the first point is already above the border or are adjacent to the phase identical, it may be extrapolation or division by zero. Appropriately bounded transition and non-zero phase difference must be requested.
3. **Nonrelativistic output at non-uniform cutoffes can read uninitialized data.** Performance PSTAB/RETSS/IMTSS uses local `NLQ(IQAT(1,IT))` (`:131`), writing but uses global NL (`:196,275–280`). The field doesn't zero before. The definition draws attention to explicit non-uniform SITES.NL in nonrelativistic mode; the potential content at omitted NL is not controlled here.
4. **Nonrelativistic logdrv is unavailable.** It is `IF (NL.LT.0)` (`:293`), marked TODO, and does not even calculate in this branch LDTAB. It is not appropriate to advertise logarithmic derivations as available in all modes.
5. **Relativist legends have a single digit format.** `WRITE(C2J,'(I1)')` (`:470,572`) cannot write 2*j>=10. From NL=6 the j=11/2 is created and the legend contains an asterisk. The definition draws attention to NL>5, the calculation does not reject due to a cosmetic error.

Definition of source numerical error does not fix and does not claim that valid input guarantees correct resonance split.

## Tests

15 new cases: finding task, written out and omitted defaults, file/text round trip, pickle/copy, all four relativistic modes, relative limits suppressing absolute defaults, rejection of correlations and inexplicable NE, ImE/SPLITSS/second path, ignored SOC, non-uniform cutoff and legends.

Whole set `src/ase2sprkkr/input_parameters/test`: **123 passed**, to default revision plus this separate patch (without previous patches). `git diff --check` and reversed application of patch has passed. SPR-KKR has not been started: it is an input and audit test for Fortran, not a numerical validation of phase shifts.
