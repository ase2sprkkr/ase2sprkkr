# SIGMA-XAS — SPR-KKR 9.7

Category for future joint audit: **broken from the point of view of accessibility by standard input**. Not supported by generally disabled placeholder: there is real implementation, but its selection in parser is missing.

## Imported Problem
rsp.f:16 describes an electric field induced XMCD; lines 941–943 contain a dispatch TASK='SIGMA-XAS ' to SIG_XAS. sigme.f:499 contains an appropriate selection matrix elements. SIG_XAS in sig_xas.f is not just a tube.

However, INITVAR for task conductivity reads only the whole keyword SIGMA (initvar.f:351–355) and sets TASK='SIGMA '. All delivered init files lack selection or assignment SIGMA-XAS. SECTION_FIND_KEYWORD (scan_input.f:144–201) seeks accurate substring with a gap before and after KW. SIGMA-XAS is therefore not recognised as SIGMA prefix. Adding a separate token SIGMA selects a standard SIGMA, not SIGMA-XAS. CONTROL.ADSI is not a direct assignment of internal TASK and the problem cannot be bypassed.

KKRCHI transmits internal TASK to RSP (kkrchi.f:1425–1427). The internal branch SIGMA-XAS is therefore not accessible by the usual input of this version. Adding the option to ASE2SPRKKR alone would not mean functional support.

## What Routine really expects (indicative map, not certified definition)
In TASK reads SIGSPROJ, GRID (default3, allowed3/6), NEA50, EMIN defaultEF, mandatory EMAX, IME .01, NEB1, DEB .05, EPS .02, IPOLEFIELD2, IT or IQ, TAUCORE. These options are in the section TASK, not automatically in the usual ENERGY. sig_xas.f:293 requires the first current operator1 or9; line 299 requires IBZINT2/4/6. The choice of core shell and lifetime must be independently audited before real support. This defaults is a source read, not a functional input recommendation.

The mere addition of the parser would allow you to enter the existing routine, but does not confirm its numerical correctness or sufficiency of initialization. We do not remove guards and does not createme fake alias to SIGMA/XAS.

## Output of this check
The default file and the Python regression test with does not create: input-only support is not available in the given solver. The finding is documented by parser and dispatch path, numerical calculation has not been started. K implementation is first needed native selection SIGMA-XAS, audit of its initialization and numerical reference case.
