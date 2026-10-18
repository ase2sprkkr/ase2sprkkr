# SOCPAR audit and definition InputParameters

Background: attached SPR-KKR 9.7, PUB9.7_260420_18h13; current ASE2SPRKKR `refactor` verified fetch, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`. SOCPAR does not have a definition here. Single patch adds `definitions/socpar.py` and `test/test_socpar.py`; does not require previous tasks.

## Entry and Real Behavior

`initvar.f:170–171` recognizes TASK SOCPAR; `gen.f:1419–1424` calls the same-name routine in KKRGEN. `socpar.f` contains no SECTION_SET/SECTION_FIND: own parameters TASK are therefore not read. There is no selection ITSEL, L, PRINT or user solver switch. All atomic types and orbitals l=1..NL-1 (`socpar.f:76, 135–136`) shall be processed. PRINT is taken from a common CONTROL.PRINT.

The definition contains CONTROL, TAU, MODE, STRCONST, ENERGY, TASK and SITES. `executable=kkrgen`, `mpi=False`: The custom loops of routines are not divided between MPI cuts and files have common names. Numerical or MPI startup was not verified here.

| Option |Behavior SPR-KKR|Definition ASE2SPRKKR|
|---|---|---|
| ENERGY.GRID |SOCPAR takes the fair part ETAB first trip|explicit GRID=1: direct equidistant real path|
| ENERGY.NE |Uses NETAB(1)|explicit NE=300; positive number|
| ENERGY.ImE |init_mod_energy converts to 0|default 0; non-zero notification|
| ENERGY.EMIN |Common native EMINDEFAULT|explicit -0.2 Ry, common option ASE2SPRKKR|
| ENERGY.EMAX |the missing limit is derived from EFERMI|None, is omitted from the input; preserved native default|
| ENERGY.EMINEV/EMAXEV |limit relative to EF in eV|common pair interface ENERGY|
| ENERGY.SPLITSS |overwrites to false|default false; when warning is switched on|
| SITES.NL |context cutoff from potential|None, is omitted from the input; explicit global maximum controlled 2..51|
| MODE.SOC |rewrites the routine SOCTL for all types/orbitals to 1e-3|When explicitly entering an alert|
| solver |set the routine BS-SOC|no invented parameter TASK|

`init_mod_energy.f:315–320` turns off SPLITSS and sets EIMAG/EILOW=0. `socpar.f:111` converts ETAB to real energy. The second path shall not be used. The definition draws attention to these choices; the more complex contour does not claim to be natively forbidden, but draws attention to the projection of its energy on the real axis.

GRID and NE are to be entered: the energy module starts with NETAB={0,0} (`mod_energy.f:43–44`) and the init_mod_energy is not independent SOCPAR default number of points. Therefore, it makes no sense to take zero native number instead of applicable explicit default. GRID=1/NE=300 is a graph option, NE can be reduced for simple checking. Listed EMIN=-0.2 Ry can be adjusted according to the desired window. Unwritten EMAX and NL delegate the actual context SPR-KKR.

## Restrictions

FULLPOT is explicitly rejected both `check_input.f:75–77` and `gen.f:1421–1422`. Input MODE does not have a FULLPOT switch and the full potential is recognized from the potential file; the validation of the definition itself cannot reveal it without the read potential. The prohibition is therefore documented and enforced SPR-KKR. The definition itself is not proof of compatibility of any potential.

DMFT and LDAU rejects GEN (`gen.f:211–222`); the definition also rejects them. SOCPAR uses relativistic wave-function channels and solve BS-SOC. To NREL/SREL/SP-SREL points out as unverified combinations: check_input for SOCPAR does not contain explicit requirement IREL=3, therefore no deemed native ban was introduced.

NL=1 produces no curve: the orbital loop starts at IL=2. LEG has 50 elements and the number of curves is NL-1, i.e. global NL>51 exceeds its capacity (`socpar.f:28, 46, 184–189`). Definition controls explicit maximum SITES.NL. U missed cutoff or details of specific potential will not perform this check itself.

## Mistakes found SPR-KKR

1. **Unattainable correction of the magnetic field sign.** Routine first turns BT and then tests `ABS(DREAL(NORM)).LT.0` (`socpar.f:81–96`). The absolute value can't be negative, so there's never gonna be another reverse sign. Probable typo in the intended test sign; precise correct correction is necessary to confirm by physical convention. It may affect the sign of the exchange split.
2. **Nesoulad of exchange output units.** SOCME contains RY_EV and has eV units (`:159`), but EXCME is `2*DREAL(XI)` without transfer (`:163–165`). The graph goes with the factor 1.0 (`:227`) and with the axis marked eV (`:216–220`). For normal internal Ry units BT it is a factor non-compliance RY_EV. Strong conclusion from the source calculation, not a comparison of numerical run with reference. The definition draws attention to it and no output data does not fix.
3. **User MODE.SOC does not control this calculation.** SOCTL=1e-3, BT=0 and SOLVER=BS-SOC are deliberately fixed directly in the routine (`:93–107`). This is not the missing parameter ASE2SPRKKR; it should be avoided that the general option SOC remains effective here.

Output files are created over XMGRHEAD for `soc` and `exc`, for individual types. The definition uses a general output reader; a specialized parser graph is not part of this patch.

## Verification

17 new tests: finding task, explicit and omitted defaults, file and text round trip, pickle/copy, relative energy limits, invalid NE/cutoffs, correlation, reminder of relativistic mode, transcription SOC/ImE and ignored second route/SPLITSS.

Full set `src/ase2sprkkr/input_parameters/test`: **125 passed** to specified default revision plus this separate patch. Number does not include previous separate patch. `git diff --check` also checked reverse patch application passed. SPR-KKR not started; tests verify input definition, not physical accuracy of calculation.
