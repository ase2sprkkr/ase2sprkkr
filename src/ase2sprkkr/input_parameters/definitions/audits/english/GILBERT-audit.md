# 35. GILBERT audit SPR-KKR 9.7 and ASE2SPRKKR

Already supported task; corrected definition for standard periodic magnetic FREL calculation. Preliminary classification: ** works but has major caveats**. It's not a joint final audit.

## Corrections to the definition

- SETFLUCT=MLIN has not been properly described: support is linear M(T)=1−DMAGDT*T and combination of spin fluctuations with vibrations. Completed DMAGDT, TDEBYE and INITSVEC.
- M_T loads the temperature network and magnetization from FLUCTFIL; now a file is required and rejected by ignored TMPMIN/TMPMAX/NTMP.
- Added standard MODE, CPA and STRCONST. The validation limits this audited procedure to magnetic FREL and BZINT=POINTS. This is the scope of the definition, not the claim that all other native branches explicitly end STOP.
- The number of points is positive integers. The final non-negative temperatures, their sequence, binding DMAGDT to MLIN and completeness of the temperature network for vibrations without SETFLUCT are checked.
- GRID=3, NE=1 nahrazeno GRID=0, NE=1. The single point GRID=3 results in an initialisation of the integration weight; the GIL itself sets ERYD=EFERMI and does not use the energy mesh for damping. EMIN/EMAX removed from the definition. SPLITSS refused.

## Defaults: written versus omitted

GRID=0 and NE=1 are written to the input; change is intentional for safe initialization. ImE=0 is written to the input, but own damping counts to EF. SITES.NL=4 is also reported: retained existing ASE default; not hidden native default. The suitability NL depends on the system and convergence.

Temperature and fluctuation parameters do not have forced Python defaults and are not recorded on entry. Native NFTET/NFPHI start at 1, INITSVEC to 0. NVIBRA is without requirement 1 but with vibration requirement INITSVEC=0 it native overwrites it to canonical 14 directions. For example NVIBRA=26 i.e. does not mean 26 directions; for the generated network it is necessary to choose INITSVEC=1. TDEBYE=0 means the table value by composition, not automatically division by zero. Negative TDEBYE/RMSU0 procedures are not covered by this definition.

For MLIN explicit TMPMIN is required because native default initial value 9999 is not a useful common entry. Without TMPMAX native uses one point in TMPMIN; NTMP without TMPMAX would ignore it, therefore it is rejected.

## Significant restrictions native calculation

1. GIL requires KMROT=0 and IBZINT=2 or 6 (gil.f:89–93). The definition makes the standard POINTS accessible. KMROT and magnetic moments depend on potential and this validation is does not verify.
2. gilsum.f:180–204 counts GEFF=2*(1+Morb/Mspin) and other non-zero-time normalizations. The calculation also takes place before the condition of printing Z>0. Empty or nonmagnetic components and compensated total moment can cause invalid results. CONTROL.NONMAG=False does not prevent it itself.
3. INITVAR turns on RDTAUDIR; a compatible initialization TAU file may be needed. The potential and external data have not been numerically verified here.
4. Native contains MPI branch in gil.f; existing ASE parallelisation setting has not been modified or numerically certified. MCS, RDLM and NLCPA are not part of this audit.

## Accessories

SPR-KKR 9.7: initvar.f:376–382; kkrchi.f:1425–1427; rsp.f:913–925; gil.f:89–121,266–320; thermal_init.f:79–81,105–119,170–242,322–368,548–558; mod_thermal.f:16; thermal_rmsdisp.f; gilsum.f:180–204. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee, single patch only for GILBERT.

## Tests

14 new tests: discovery/defaults/pickle/executable, round-trip MLIN with DMAGDT and vibrations, set M_T and ignored parameters, completeness of the temperature network, positive numbers, rejection of unaudited modes, NONMAG, incorrect energy mesh and missing MLIN temperature. Full set InputParameters: **122 passed**, two pre-existing NumPy DeprecationWarning messages. This is not a numerical test SPR-KKR.

Command: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.
