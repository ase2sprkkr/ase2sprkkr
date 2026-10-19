# BIS — SPR-KKR 9.7

Category for future joint audit: **broken**, not off placeholder. ASE2SPRKKR baseline 267512c4 does not have a separate BIS definition. V of this version native of the program cannot be created input of the definition that would lead to functional BIS of the calculation.

## The Path Shown
1. initvar.f:264–271 detects TASK BIS, sets TASK to a ten-character string `BIS       `, RDTAUDIR=true, IHYPER=0 and ICORE=0.
2. check_input.f:91–93 connects BIS to VBPES and requires IREL=3. Even valid fully relativistic input does not solve the following error.
3. gen.f:1445–1448 transmits BIS unchanged as argument TASK to VBPES.
4. vbpes.f:174 controls `TASK(3:5).NE.'PES'` before calculation and calls STOP_MESSAGE with text `TASK ???`. For BIS this substring is `S  `; the condition is always fulfilled. The calculation may fail earlier on another condition, but it will necessarily stop in achieving this check.

Searching vbpes.f, vbpesint.f and vbpes_aux.f did not find a separate branch BIS that would circumvent that stop. The mere removal of the check is not documented by the repair of the physics inverse photoemission. Rewrite TASK to VBPES would change the required calculation to photoemis and must not be presented as support BIS.

## Definition, defaults and tests
A separate bis.py or an alias to VBPES with does not create: would only provide input to a proven stopping native branch, not supported task. Cannot recommend numerical defaults as functional BIS configuration. The finding is supported by the source route; does not create with a Python test that would only copy the Fortran condition and fake verification of the solver. Numerical SPR-KKR run not started.

## What Is Needed to Support
First native correction and verification of actual BIS implementation (energy of initial/final states and transition probability), then InputParameters definition by corrected reader and numerical regression case. The current resources are not enough to claim that it is only a faulty guard.
