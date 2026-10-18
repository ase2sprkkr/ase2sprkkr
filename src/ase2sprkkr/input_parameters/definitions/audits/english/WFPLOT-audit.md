# Task 1: WFPLOT

Comparison of public resources SPR-KKR 9.7 from the Annex
`PUB9.7_260420_18h13.tar.gz` and ASE2SPRKKR, branch `refactor`, commit
`8549af97ae106b08749ac1f2dc4b819a72b988e1` z 30 September 2026.
The check took place 1 October 2026. Without triggering binary SPR-KKR.

## Status and delivered files

WFPLOT does not have a default module by default checkout. ASE2SPRKKR has modules
ARPES, BSF, DOS, GILBERT, JXC, PHAGEN, SCF, TORQUE and XAS; some additional
cover aliases or variants. Deficiencies of module does not mean that cannot be manually
assemble your own InputParametersDefinition.

New Files:

- `src/ase2sprkkr/input_parameters/definitions/wfplot.py`
- `src/ase2sprkkr/input_parameters/test/test_wfplot.py`

Module appears in `InputParameters.definition_modules` automatic
search via pkgutil; the register does not need to be changed. Application:

```python
from ase2sprkkr import InputParameters

ip = InputParameters.create("WFPLOT")
ip.CONTROL.POTFIL = "Fe.pot"
ip.TASK.ITSEL = 1
ip.TASK.L = 2
ip.TASK.MJ = "+3/2"
ip.ENERGY.EMIN = 0.5
ip.save_to_file("Fe-wfplot.in", atoms=False)
```

For core functions:

```python
ip.TASK.STATE = "CORE"
ip.TASK.L.clear()  # L is not read in CORE
ip.TASK.CL = "2P"
ip.TASK.MJ = "1/2"
```

## Parameters and actual defaults

| Parameter |SPR-KKR missed|New definition| Serialization default |
|---|---|---|---|
| ITSEL |Mandatory; initial 1 not applicable fallback| 1 | Yes |
| STATE | BAND | BAND | No |
| TASK.PRINT |Current IPRINT of CONTROL.PRINT| Dynamicky CONTROL.PRINT | No |
| L, BAND | NLQ(IQAT(1,ITSEL))-1 |None means keeping the decision to SPR-KKR| No |
| CL, CORE | 2P | 2P | No |
| MJ | +1/2 |1/2, same physical value| No |
| GRID |For WFPLOT there is no separate applicable task-specific default| 3 | Yes |
| NE |NETAB begins with zero; omitted NE may end with NEMAX=0| 1 | Yes |
| EMIN | -0.2 Ry | -0.2 Ry | Yes |
| EMAX |Generally EF if not set|Same as EMIN| Yes |
| ImE | 0.01 Ry | 0.01 Ry | Yes |
| SITES.NL |3 or 4 for selected heavy elements; additional changes, if applicable|Shared ASE default 3| Yes |

Unwritten represented defaults corresponds to readers SPR-KKR.
`L=None` is not a claim that numerical SPR-KKR default is None: means
undefined value delegated to context dependent Fortran default.
It is not correct to insert omitted fixed `L=2` for each atomic type.

`GRID=3, NE=1` is more suitable than common ASE default `GRID=5, NE=32`:
WFPLOT in BAND the branch uses only `ETAB(1,1)` and the other energy will not pass.
For GRID=3 and NE=1 this point is exactly `EMIN+i*ImE`. Explicit EMAX=EMIN
does not change this point; it also allows you to automatically write both relative limits
EMINEV/EMAXEV which INIT_MOD_ENERGY requires together. Upper limits
are accessible over `ip.ENERGY["EMAX"]` and `["EMAXEV"]` as shared
ENERGY(emax="emin").

SPR-KKR default SITES.NL is better than fixed 3 for heavy elements. New Module
it maintains a shared, explicitly written policy ASE2SPRKKR; this
the general policy does not change. For f-figures it is necessary to set NL at least 4.

## Syntax, restriction and validation

- `TASK` contains a separate keyword `WFPLOT`, not TASK=WFPLOT.
- `ITSEL` is an atomic type index, not a crystallographic site index.
It is positive and always with writes. WFPLOT calls
SECTION_SET_INTEGER('ITSEL',ITSEL,9999,1), thus ends without it.
- `MJ` with writes as a semi-total fraction: `1/2`, `+3/2`, `-5/2`.
TXTMJ2NUM accepts odd numerators not more than 61; decimal `0.5` does not read.
- At abs(MJ)>l+1/2 Fortran limits the value. The definition will draw attention to this,
the value does not automatically overwrite.
- `STATE` accepts BAND or CORE. CL accepts n=1..6 and s/p/d/f;
additional physical control l<n rejects non-existent shells (e.g. 1F).
This additional check is more stringent than the reader itself SPR-KKR.
- L must be non-negative; for single SITES.NL, L<NL is checked.
For different NL different places without mapping type to places
doesn't pretend to know the right cutoff.
- TASK.PRINT inherited CONTROL.PRINT. Exactly different value with writes.
- MODE.MODE=NREL/SREL/SP-SREL leads to a warning: INITVAR for WFPLOT
then force IREL=3.
- GEN rejects WFPLOT for DMFT and LDA+U; definition rejects these options at:
final validation.
- FULLPOT results directly in WFPLOT to STOP_MESSAGE. The definition documents this;
does not itself open the specified POTFIL or verify its FULLPOT status.
- MPI must be disabled. GEN for WFPLOT in MPI the branch ends before induction
  WFPLOT. The new definition has mpi=False.
- NE>1 causes a warning because BAND only uses the first point.
- STATE=CORE set L and STATE=BAND set CL call
an alert to an ignored parameter.

## Source evidence

- `initvar.f:149–218`: detection TASK WFPLOT, switching off normal energy
loops and TAU files.
- `initvar.f:782`: enforcement IREL=3.
- `wfplot.f:67–70`: FULLPOT not implemented.
- `wfplot.f:98–114`: ITSEL, local PRINT and STATE.
- `wfplot.f:128–130`: reading CL and MJ for CORE.
- `wfplot.f:292–300`: only ETAB(1,1), default L and MJ; if there is energy
too close to zero, replacement for (0.5,0.01) Ry.
- `scan_input.f:203–280`: required sentinel 9999 and STOP_KEY=1.
- `scan_input.f:1916–1983`: TXTMJ2NUM, text format MJ.
- `scan_input.f:2138–2201`: CL default 2P and limits n and l.
- `scan_input.f:2203–2280`: default MJ=0.5 and limit to l+0.5.
- `gen.f:210–222`: refusal DMFT/LDA+U for WFPLOT.
- `gen.f:1357–1370`: MPI guard; `gen.f:1452–1453`: call WFPLOT.
- `mod_energy.f:47–54`, `init_mod_energy.f:73–76, 133–180, 310–391`:
Energy initialization, limits and EPATH.
- `epath.f:132–144`: first point GRID=3.
- `kkrmain.f:410–430`: atomic dependent defaults NL.

## Verification and open points

Full directory `input_parameters/test`: 128 of successful tests. Z toho 20
new cases WFPLOT: discovery, writing and readback via from_file,
copy/pickle, defaults, fractional syntax MJ, positive indices, CORE/BAND,
Limit L, relative energy, warning and rejection DMFT/LDA+U.

No calculation was started SPR-KKR. Thus, no numerical result is verified.
produced .agr files. The common DefaultOutputReader shall apply;
specialized parser of resulting wave functions is not part of the change.
The existence of the selected core shells for a specific element is not checked,
upper limit ITSEL according to actual potential, FULLPOT in actual
POTFIL or cutoff for individual types at non-single NL.

Other individual task by INITVAR: MEPLOT. The inventory of all tasks is not
claims that everyone in the public version is actually applicable; for example
PHONONS has a special track in the public version. Everyone needs to go through separately.
