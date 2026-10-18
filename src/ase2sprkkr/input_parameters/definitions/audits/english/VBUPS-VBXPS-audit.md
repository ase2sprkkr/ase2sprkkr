# VBUPS / VBXPS — aliasy VBPES v SPR-KKR 9.7

This is not a new separate physical tasks. initvar.f:224–249 seeks VBXPS/VBPES/ARPES/VBUPS and according to ANGLE_INTEGRATED selects internal VBPES or ARPES. Separate VBUPS and VBXPS without other switches sets TASK=VBPES. GEN then starts the already audited VBPES route.

## Important error in expected defaults
vbpes.f:197 reads TASK.EPHOT with default1253.6 eV, for both aliases as well. Not branching according to VBUPS/VBXPS. NINT(EPHOT)=1 means1253.6 and NINT(EPHOT)=2 means1486.6 eV. VBUPS does not set itself21.2 or40.8 eV. Photon energy for UV photoemisi must be specified by the user explicitly in the definition VBPES. Adding an absent UV default which is not written would be incorrect; explicitly written21.2 eV would be possible, but would be a new Python convention, not native default alias.

## Kombinace s ARPES
VBXPS then independently controls ARPES and changes ANGLE_INTEGRATED depending on its presence. It is not a guarantee that any combination of tokens will maintain an integrated calculation. We do not represent new aliases as compensation KKRSPEC ARPES/AIPES; executive and native implementations are essential.

## Recommended use and audit output
Use existing definition VBPES and set TASK.EPHOT in eV by source. No duplicate stand-alone input definition or regression test is created, copying the same native mapping. This check complements earlier audit VBPES; does not verify numerical PES results.

Category for future joint audit: functional aliases VBPES, with significant ALE for the meaning of the name VBUPS and its photon default. Not even a new independent implementation is switched off. Numerical native calculation not started.
