# 44. SPLEED separate definitions and audits

ASE already offers SPLEED as TASK option in definition ARPES but missing a separate InputParameters.create("SPLEED"). Separate definitions sharing spectroscopic sections are added; executive kkrspec. Preliminary **works but has major caveats***: this audited procedure is a fixed impact angle and energy scan, not all native polarization and geometrical branches. Numerical calculation has not been initiated.

## Different energy semantics

SPEC_RSLAB recognizes SPLEED and sets IBLOCH=1, NSPLEE=0. SPEC_INPUT then for NE>1 reads EMINEV/EMAXEV again and uses them as **conquering energy in eV**, not as relative energy to EF. Take over ARPES −8..+5 eV would therefore be wrong. Definitions writes 11..60 eV, NE=200; descriptions of both EV parameters explicitly indicate this symantics.

Absolute EMIN/EMAX are rejected in this definition at a modern input because native SPLEED the branch is not read in this place instead of EV alternatives. NE=1 is rejected: reader in this case skips reload EV borders and EFEV then takes EMIN from previous general initialization; single point calculation therefore may not have the same meaning as multipoint eV scan. For the real one-point workflow it is necessary to solve this native string, not only set NE=1.

SPEC_PH deleted: native for SPLEED lists photon geometers as obsolete. It is not meaningful to promise dependence LEED to EPHOT or photon polarization.

## Defaults and writing

| Parameter |Separate definition|Native / Meaning|
|---|---|---|
| EMINEV/EMAXEV |11 / 60 eV, written|corresponds to SPLEED start range; direct eV|
| NE |200, entered| native SPLEED start EMESH=200 |
| THETA |45..45 degrees| native SPLEED default |
| PHI |270..270 degrees| native SPLEED default |
| NT/NP |1/1, written down|fixed angles|
| SPOL |1, entered| native IBLOCH=1 default; **unpolarized** |
| EWORK_EV |4.2 eV, written|generic explicit workfunction; native may use potential|
| IMV_INI_EV/IMV_FIN_EV |0.2/0.2 eV, written| SPLEED native start VIL/VIH |
| N_LAYER/NLAT_G_VEC |1/30, written down| SPLEED native start |
| N_LAYDBL |10/10, written down|explicitly sets both variables despite native to double LANZ1 assignment|
| POL_E | PZ |existing explicit convention; native default also PZ|
| CRYS_VECS |false, not written|preserves native use of primitive crystal vectors|

SPOL=1 is a non-polarized calculation even with a draw name SPLEED. This is not the hidden difference defaults; it is written and corresponds to the reader. For spin response, the user must select the appropriate native SPOL procedure and interpret results. All spin output branches are not numerically certified here.

Output work 4.2 can be worse than the value of a specific potential; because with writes it is a real explicit overwrite. The defaults cannot be evaluated only by equality with native. The dampening is in eV and is not electron temperature.

## Other errors and restrictions

- Taken ARPES option CRYS_VEC has the wrong native name. SPEC_RSLAB reads CRYS_VECS and sets USE_CRY_PRIM_VECS=false at its presence. The separate definition contains the correct key and description of shutting down primitive vectors. This is also relevant for the original ARPES definition; this separate patch does not change it.
- The fixed-angle range does not consciously permit angular scans, a moment scans or POL0/POL0L. Native some of the following procedures have, but polarization depends on the presence of keys: cannot be used simultaneously POL_E and nonzero POL0/POL0L. Blind sharing of all ARPES defaults could end up STOP.
- For INPVER=0 it is an external old input.inp; validation of modern energy keys is not performed in this branch. Existence and content of external files shall not be verified.
- IQ_AT_SURF, MILLER_HKL and the surface structure shall correspond to potential. The definition does not pretend to control this compliance from just input options. Full-potential and other special branches were not audited throughout the chain.
- mpi=False is a conservative choice of new definition; no error of whole MPI has been demonstrated nor has a numerical comparison been made.

## Tests

11 new tests; whole set InputParameters **119 passed**, two existing NumPy DeprecationWarning messages. Discovery verified, independent task and defaults, absence of photon section, incident EV border, round-trip/pickle, serial executive, rejection of invalid borders/absolute variants/single point energy, workfunction, angular and moment scans and conflict polarization alternative. A separate round-trip test checks the correct key CRYS_VECS.

Command: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

The sourcee: spec_rslab.f:160–225; spec_inputs.f:1198–1241,1301–1328,1358–1389,1549–1648,1764–1783. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Single patch only for SPLEED; original ARPES patch audit is not a prerequisite. No numerical solver.
