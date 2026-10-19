# CLPES — SPR-KKR 9.7 / ASE2SPRKKR

Nová definice: TASK=CLPES, executable kkrspec, konzervativně serial. Skutečná větev IBLOCH=3: spec_rslab.f:172–173,226,239–241,862; driver spec_clpes.f.

## Parametry a defaults
SPEC_CL je skutečný název sekce (spec_inputs.f:1877–1921). Implementováno IT (explicitně 1), CL (explicitně 2P), IE_ICST (default 1). CL je vždy zapsáno: bez přítomného CL reader ponechá 1S, přestože SECTION_SET_STRING uvádí fallback 2P. Volitelné CE není v této definici exponováno; IT určuje potenciálový typ a native výběr zahrnuje také stejný atomový počet Z (spec_clpes.f:213–216), nejde o izolovaný site selector.

NE=300 zapsáno, single mesh bez SPLITSS. Driver přepisuje energetické okno rozsahem core energií rozšířeným o +/-3 eV (spec_clpes.f:301–306); EMIN/EMAX nejsou uživatelské spektrální meze. NE=1 vybírá IE_ICST a validuje se počet 4*l+2. Core energie jsou vztažené k EF (ř.224).

EPHOT=1000 eV se explicitně zapisuje místo zděděných 25 eV: vhodnější pro core PES, ale uživatel musí ověřit práh konkrétního shellu. Nejde o záruku dostupnosti každé core hladiny. Zdělěné explicitní WORK=4.2, IMV_INI_EV=.05, IMV_FIN_EV=2, NL=4, N_LAYER=50, NLAT_G_VEC=57 jsou zachovány; nejsou deklarovány jako nativní defaults. CRYS_VEC opraveno na skutečné CRYS_VECS; jeho přítomnost vypíná primitive crystal vectors. NPOL nepřítomné odpovídá native 1; 0 a 3 vyvolají STOP a validátor je odmítá.

## Velká ALE
- NCPA>0 vede k STOP (spec_clpes.f:157). Definice neověřuje obsazení potenciálu a nemůže zaručit ordered systém.
- Auditovaný rozsah je pevná konečná elektronová geometrie. Angular a momentum scans jsou odmítnuty; starý externí input INPVER=0 jejich kontrolu obchází.
- NE=1 používá ECOR/XNORM po smyčce atomů (spec_clpes.f:310–319); pro více typů/site může odpovídat poslední zpracované sadě, ne požadovanému typu. Tuto větev nelze považovat za obecně ověřenou.
- Vybraný shell musí v potenciálu existovat jako core shell, IT musí existovat a EPHOT musí překročit práh. Python validace bez potenciálu tyto podmínky nezaručí.
- MPI nebylo certifikováno, automatické použití je vypnuté; to není tvrzení, že native MPI chybí.

Kategorie pro budoucí společný audit: funguje, ale má velká ALE. Nejde o vypnutý placeholder.

## Ověření
8 nových regresních testů: zápis/čtení, selector a shell, odmítnutí nepodporované polarizace a geometrie, core state bounds. Celá input_parameters sada: 116 passed. Numerický SPR-KKR výpočet nebyl spuštěn.
