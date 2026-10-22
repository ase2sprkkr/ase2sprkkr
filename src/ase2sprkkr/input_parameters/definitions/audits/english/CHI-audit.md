# 37. CHI audit SPR-KKR 9.7 and ASE2SPRKKR

Missing task so far. Added definition for serial periodic ASA FREL magnetic susceptibility calculation. **Previously: works but has major caveats iteration of spin magnetization is disabled in this version of sources.** Numerical correctness has not been verified. This is an individual audit, not a joint final audit.

## Mapping and Range

TASK CHI recognizes INITVAR and switches off RDTAUDIR/WRTAUDIR; RSP for ASA calls SUSC. Executive is a krachi. The definition contains CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES and minimum optional SCF section for validation FULLPOT.

The audited procedure is limited to FREL, BZINT=POINTS, one GRID=5 contour without SPLITSS, positive NE and OP single/NONE/BROOKS. This is the definition of the normal procedure, not the explicit argument STOP for all other native variants. NONMAG is allowed: paramagnetic susceptibility can be calculated from nonmagnetic base state. The potential, rotational states and other structure-specific settings have not been numerically verified.

TASK offers CHIPRINT and NOSYMCHI. MODE.PRINT is read again inside SUSC and can overwrite CHIPRINT taken from RSP. MODE.NCHIITER is mandatory writes as 1 and other values are rejected; MODE section is mandatory due to this writing.

## Errors and major restrictions native code

- **Iteration does not occur:** SUSC reads MODE.NCHIITER (native initial value 100), MIX and TOL but to susc.f:690 necessarily sets NCHIITER=1. The following branch IF(NCHIITER.GT.1), mixing and controlling convergence are so unattainable. The result cannot be presented as a converged self-consistent response calculation. Definition MIX/TOL does not offer because they would promise ineffective iteration management in this procedure.
- **MPI is not available:** susc.f:139–147 at MPI calls STOP_MESSAGE with information on developed implementation. mpi=False is required default; explicit enforcement MPI in calculator will not pass this STOP.
- **Public FULLPOT branch not available:** rsp.f:813–828 at PUBLIC_VERSION ends STOP_REGULAR. It is not enough to enter STAT. FULLPOT is therefore rejected; if the potential activates it, the control is decisive native.
- The automatic DOS to EF is expected to be fixed Im(E)=0.001 Ry (susc.f:269), independently of ENERGY.ImE. This is a native option, not an erroneously transmitted value ASE.
- Advanced RDABXC/ITRCHI0 restarts, manual NEF for types SIGNAF and cluster/surface/weak procedures are not available. RDABXC/MIX/TOL/NCHIITER are read from MODE, not from TASK relevant for any future extension.

## Defaults and writing

**Written defaults:** GRID=5, NE=30, ImE=0.01, EMIN=−0.2, MODE.NCHIITER=1. The number of energy and k-points requires convergence. NCHIITER=1 differs from the initial native 100 but corresponds to the actual forced number of passageways and does not pretend iteration. This explicit default is more appropriate than native initial value.

**Unwritten:** EMAX is left native EF. FREL with writes absence of switch key. NOSYMCHI=False is omitted from the input and agrees native by searching for the presence of the key. FULLPOT=False is omitted from the input and does not force the conversion of full potential to ASA. MODE.OP without is omitted from the input: general ORBPOL and local CHIORBPOL are not identical variables. SUSC when missing OP sets CHIORBPOL to BROOKS to describe/choose response; it cannot be automatically derived that the underlying potential uses Brooks correction. Therefore, the existing common MODE.OP was not artificially set to BROOKS.

MODE.PRINT and TASK.CHIPRINT do not have forced defaults; native uses input IPRINT/CHIPRINT. MIX=0.5 and TOL=1e−5 are native initial values but their iterative effect is suppressed by the error mentioned above.

## Tests

**11 new tests**, full set InputParameters **119 passed**, two existing NumPy DeprecationWarning messages. Tests verify discovery, serial execution/defaults, actual writing NCHIITER=1, parse/write round-trip, pickle, both print options and NOSYMCHI, OP=BROOKS, permissible NONMAG and refusal FULLPOT/unaudited modes, incorrect networks and calculus and ineffective iterations.

Command: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

The sourcee: initvar.f:412–418; rsp.f:195–201,811–835; susc.f:139–207,240–304,334–343,690–781. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Single patch only for CHI. Tests do not start SPR-KKR.
