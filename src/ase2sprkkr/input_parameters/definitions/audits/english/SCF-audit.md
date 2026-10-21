# SCF — current implementation audit

**Status:** Supported algorithms BROYDEN2, ANDERSON, TCHEBY; TCHEBY excludes FULLPOT. This audit describes the implementation on `tasks-workflow` against SPR-KKR 9.7 (`PUB9.7_260420_18h13`); it is not numerical certification.

## Interface, validation and native constraints

1. ASE offered ALG=PULAY but scf.f:427–437 accepted only BROYDEN2, ANDERSON and TCHEBY. PULAY leads to SCFALG not found. Task copy ALG now offers only these three names; PULAY is not automatically mapped to ANDERSON because such an alias native reader does not have.
2. TCHEBY s FULLPOT native explicitly rejects; validation added. FULLPOT can also be activated by potential, which input-only validation itself does not detect.
3. NITER and ISTBRY are positive integers. ITDEPT limited to 3..40: SCFBROYDEN converts values below 3 to 40 and limits values above 40 to capacity. the definition now does not pretend that these requirements will be implemented accurately.
4. MIX, MIXOP and TOL must be final positive values; EFGUESS final energy. The upper limit MIX was not arbitraryly set at 1 because overrelaxation is not merely a type error.
6. Description EFGUESS specified: This is the initial estimate EF in Ry and the omission of the initial search when creating starting potential. It does not prohibit all subsequent updates EF in SCF.

Corrections shall be made through the drawing copies SCF and ENERGY to avoid changing other definitions using the common sections. The common section ALG may still contain PULAY in other draw definitions; this is not a confirmation of its native support.

## Defaults and interpretation

Normal values NITER=200, MIX=0.2, TOL=1e−5, ISTBRY=1 and ITDEPT=40 are written to the input and correspond to initialisation in scf.f:360–371. BROYDEN2 corresponds to MOD_SCF default and is default of the draw definition. VXC options, their default and metadata have not been modified; MOD_SCF initializes VWN, but specific potential and reader can continue to influence the status input. Energy default GRID=5/NE=30 corresponds to SCF to INIT_MOD_ENERGY; EMIN=−0.2 is written ASE option EMAX is omitted from the input and leaves native EF.

**MIXOP is an important unwritten case:** scf.f first sets SCFMIX=0.2 and SCFMIXOP=SCFMIX, then reads MIX and MIXOP. At MIX=0.5 and single MIXOP it is therefore not possible to claim that MIXOP will be 0.5; for normal non-export procedure, 0.2 remains. ASE does not impose another hidden default, MIXOP is single and is omitted from the input. The description of this dependence is now explicitly highlighted. To complete explicit MIXOP is better if the user wants the same value as changed MIX; automatic change default would change native behavior.

EFGUESS is not recorded without input; native is determined by the presence of the key, not by any replacement energy. FULLPOT=False is omitted from the input; it is not a promise that the full potential will change to ASA. TAU/CPA/STRCONST/SITES/MODE are common sections, here without interference.

## Other restrictions and branches

- Native EXPORT_HARTREE_POTENTIAL sets NITER=1, MIX=1 and skips reading VXC/NITER/MIX. The table of normal defaults is not transferable to this branch.
- ANDERSON mapuje IMIX=7, BROYDEN2 IMIX=4. TCHEBY uses a separate procedure; the algorithm name is not synonymous with the same numerical mix.
- Native offers further management as SCFRESTART, MIXRHO, SIM, MAGCOUPL, SEMICORE and 2D Layer charge parameters. They are not completed here and their entire workflow is not covered by this audit.
- mpi=True preserved. SCF has parallel solutions and initial search; numerical match MPI/serial in this audit has not been tested.

## Current implementation contract

- File roles are defined by the current `TaskDefinition` plus any conditional schema methods; required inputs, reusable cache and produced outputs have distinct semantics.
- Conditional file declarations: `cache_files`, `output_files` in the schema.
- Schema source: `src/ase2sprkkr/input_parameters/definitions/scf.py`.

## Verification

Definition and file-contract tests live in `src/ase2sprkkr/input_parameters/test/`. No new full-suite or numerical SPR-KKR execution is claimed here.
