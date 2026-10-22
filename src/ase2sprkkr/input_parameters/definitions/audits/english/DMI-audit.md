# DMI — current shared JXC/DMI definition audit

DMI is supported by the shared `jxc.py` definition, selecting the native magnetic E-INTEGRAL DMI workflow; there is no separate `dmi.py` in the integrated tree. Native physical restrictions below are based on SPR-KKR 9.7.

## Current implementation

The shared `jxc.py` schema exposes `JXC` and `DMI` selectors; a DMI-only run selects the corresponding native branch.

## Native output naming

Native `initvar.f:536–545` activates `CALC_DMI` and the internal `TASK=E-INTEGRAL`. `dmivector.f:416–418,605–606` uses `TRIM(DATSET)//'_DMIVEC_Dij.dat'` and the matching `DMIVEC_Dij_[xyz]` plot families. In `kkrmain.f:239–245`, any nonempty `ADSI` is appended to `DATSET`, so native output names include ADSI. Current `jxc.py::output_files()` uses `{dataset}_{adsi}` to match native DATSET; the separator is omitted for empty ADSI.

## Real involvement

- `initvar.f:534–553`: keyword DMI sets CALC_DMI and TASK=E-INTEGRAL;
QREPONLY limits central locations. DMI is not a separate value of internal TASK.
- `gen.f:355–393`: periodic branch prepares off-diagonal TAUIJ and requires
  IBZINT=2 (POINTS). It is not CLUSTER or WEYL calculation.
- `gen.f:975–1000`: magnetic FREL calls DMIVECTOR if not MULTISPIN.
SP-SREL instead enters the branch of the scalar J_ij and does not create DMI vectors!
- `dmivector.f:70–76`: NONMAG, KMROT other than 0, NEPATH other than 1,
GRID other than 5, IREL<=2 and NKM other than NKMMAX lead to STOP.
- `gen.f:1193–1252` and `dmivector.f:108–129`: completion of energy MPI loop,
reduction of fully initialized DMI_IJINT and output to the cut 0 are implemented.
The setting `mpi=True` has a source justification; the actual MPI run has not been executed.

Definition rejects NREL/SREL/SP-SREL, NONMAG, other GRID, multiple contours,
SPLITSS and other BZINT. KMROT and the size of the base depend on potential;
Python doesn't guarantee entry without it.

## Inputs

TASK: DMI, CLURAD, NSHLCLU, SELECT, EXCLUDE, QREPONLY.
`init_mod_tauij_star.f:123–162` reads lists of elements (maximum 100) and radius (s).
SELECT and EXCLUDE are mutually exclusive. They select places according to elements and casts,
not individual types of the resulting matrix. NSHLCLU has priority ** already with its presence**, even for zero.
`init_mod_tauij_star.f:591–635` describes different logic SELECT and EXCLUDE in mixed cast.

MULTISPIN starts a different routine XCPLDMI_FOURSPIN_3SITE; it is not part of this definition.
Embedded-cluster branch calling CLUXCPLTENSOR, is not the equivalent of periodic DMIVECTOR.
Heat fluctuations are in native routine, but their separate inputs and correctness
have not been verified by this audit. These branches are not marked off.

## Defaults and writing

| Parameter |Definitions|Native behavior without parameter|Assessment|
|---|---|---|---|
| GRID |5, writes|5 for E-INTEGRAL|the same|
| NE |30, writes|30 for E-INTEGRAL|native default taken over instead of common ASE default 32|
| ImE |0.01 Ry, writes| 0.01 Ry |Same; `mod_energy.f:43`|
| EMIN |-0.2 Ry, writes| -0.2 Ry |the same|
| EMAX |is omitted from the input, not exposed|EFERMI after conversion of the sentinel value|the same effective behaviour; `init_mod_energy.f:182`; the definition does not offer a direct change to the upper limit|
| CLURAD |2.2 ALAT, writes|none applicable default; without NSHLCLU and CLURAD is STOP|explicit option taken from JXC; radius convergence must be verified by the user|
| NSHLCLU |not entered, is omitted from the input|Internally 0 but absence switches to reading CLURAD|presence and absence are not identical|
| SELECT / EXCLUDE |not entered, omitted from the input|blank lists|the same|
| QREPONLY | false, is omitted from the input | false |Same; `mod_sites.f:73`|

The defaults that are not written do not mask another native selection.

## Significant ALE and results

1. Dz is not calculated by the same procedure as Dx and Dy.
`dmivector.f:321–395` uses weight (E-EF), for xy/yx
a different expression for these terms. `dmivector.f:483–499` additionally introduces a factor of 1/2 for Dz.
The comment explicitly states that Dz is not defined in the same way for M ©z.
This is not a documented programming error, but an important physical interpretation:
three columns are not automatically equivalent components obtained by one procedure.
2. **QREPONLY /SELECT /EXCLUDE changes the micromagnetic sum.**
Sum only runs over selected pairs; `dmivector.f:520–543` has no explicit
supplementing the number of representatives. Native output itself draws attention to conversion
the result per unit cell (`99014`). Limited sum not to consider without control
for the full sum of the whole cell. The same applies to the final CLURAD.
3. **Symmetrical equivalents require numerical control.**
BAR_M_T is created for type z IQAT(1,IT); type nuts are used on steam.
In This pair loop is not an explicit rotation within an equivalent neighboring location.
It is not possible to declare the result generally incorrect from the reading itself, because part
Transformation may be included in input TAUIJ. The audit therefore records:
non-verified area, not some error similar to previously audited JXC.

Output: DATASET_DMIVEC_Dij.dat and graphs of individual folders; pair folders in meV,
Micromagnetic totals printed as meV*A/f.u. Indexing DMI_ITJT in writing i
read is consistent; error from BSDOSIJ is not documented here. On-site members are in
leave the output at a distance less than 1e-6 ALAT.
