# CHI — current implementation audit

**Status:** Static susceptibility in the audited serial ASA/FREL path; native NCHIITER is effectively 1. This audit describes the implementation on `tasks-workflow` against SPR-KKR 9.7 (`PUB9.7_260420_18h13`); it is not numerical certification.

## Mapping and Range

TASK CHI recognizes INITVAR and switches off RDTAUDIR/WRTAUDIR; RSP for ASA calls SUSC. Executive is a krachi. the definition contains CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES and minimum optional SCF section for validation FULLPOT.

The audited procedure is limited to FREL, BZINT=POINTS, one GRID=5 contour without SPLITSS, positive NE and OP single/NONE/BROOKS. This is the definition of the normal procedure, not the explicit argument STOP for all other native variants. NONMAG is allowed: paramagnetic susceptibility can be calculated from nonmagnetic base state. The potential, rotational states and other structure-specific settings have not been numerically verified.

TASK offers CHIPRINT and NOSYMCHI. MODE.PRINT is read again inside SUSC and can overwrite CHIPRINT taken from RSP. MODE.NCHIITER is mandatory writes as 1 and other values are rejected; MODE section is mandatory due to this writing.

## Errors and major restrictions native code

- **Iteration does not occur:** SUSC reads MODE.NCHIITER (native initial value 100), MIX and TOL but to susc.f:690 necessarily sets NCHIITER=1. The following branch IF(NCHIITER.GT.1), mixing and controlling convergence are so unattainable. The result cannot be presented as a converged self-consistent response calculation. Definition MIX/TOL does not offer because they would promise ineffective iteration management in this procedure.
- **MPI is not available:** susc.f:139–147 at MPI calls STOP_MESSAGE with information on developed implementation. mpi=False is required default; explicit enforcement MPI in calculator will not pass this STOP.
- **Public FULLPOT branch not available:** rsp.f:813–828 at PUBLIC_VERSION ends STOP_REGULAR. It is not enough to enter STAT. FULLPOT is therefore rejected; if the potential activates it, the control is decisive native.
- The automatic DOS to EF is expected to be fixed Im(E)=0.001 Ry (susc.f:269), independently of ENERGY.ImE. This is a native option, not an erroneously transmitted value ASE.
- Advanced RDABXC/ITRCHI0 restarts, manual NEF for types SIGNAF and cluster/surface/weak procedures are not available. RDABXC/MIX/TOL/NCHIITER are read from MODE, not from TASK relevant for any future extension.

## Defaults and interpretation

**Written defaults:** GRID=5, NE=30, ImE=0.01, EMIN=−0.2, MODE.NCHIITER=1. The number of energy and k-points requires convergence. NCHIITER=1 differs from the initial native 100 but corresponds to the actual forced number of passageways and does not pretend iteration. This explicit default is more appropriate than native initial value.

**Unwritten:** EMAX is left native EF. FREL with writes absence of switch key. NOSYMCHI=False is omitted from the input and agrees native by searching for the presence of the key. FULLPOT=False is omitted from the input and does not force the conversion of full potential to ASA. MODE.OP without is omitted from the input: general ORBPOL and local CHIORBPOL are not identical variables. SUSC when missing OP sets CHIORBPOL to BROOKS to describe/choose response; it cannot be automatically derived that the underlying potential uses Brooks correction. Therefore, the existing common MODE.OP was not artificially set to BROOKS.

MODE.PRINT and TASK.CHIPRINT do not have forced defaults; native uses input IPRINT/CHIPRINT. MIX=0.5 and TOL=1e−5 are native initial values but their iterative effect is suppressed by the error mentioned above.

## Current implementation contract

- File roles are defined by the current `TaskDefinition` plus any conditional schema methods; required inputs, reusable cache and produced outputs have distinct semantics.
- `output_file_templates`: `( "{dataset}*_dens.agr", "{dataset}*_chifmag.agr", "{dataset}_*_GAMMA-OP.dat", "{dataset}_*_CHI-r.dat", "{dataset}_*_*_I-spin.dat", )`
- `cache_file_templates`: `("{dataset}.tau", "TKTK_tables")`
- Schema source: `src/ase2sprkkr/input_parameters/definitions/chi.py`.

## Verification

Definition and file-contract tests live in `src/ase2sprkkr/input_parameters/test/`. No new full-suite or numerical SPR-KKR execution is claimed here.
