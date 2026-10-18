# 35. GILBERT — audit SPR-KKR 9.7 a ASE2SPRKKR

Již podporovaný task; opravená definice pro standardní periodický magnetický FREL výpočet. Předběžné zařazení: **funguje, ale má velké ALE**. Nejde o společný závěrečný audit.

## Opravy definice

- SETFLUCT=MLIN nebyl správně popsán: podporuje lineární M(T)=1−DMAGDT*T a kombinaci spinových fluktuací s vibracemi. Doplněny DMAGDT, TDEBYE a INITSVEC.
- M_T načítá teplotní síť a magnetizaci z FLUCTFIL; nyní se vyžaduje soubor a odmítají ignorované TMPMIN/TMPMAX/NTMP.
- Přidány standardní MODE, CPA a STRCONST. Validace omezuje tento auditovaný postup na magnetický FREL a BZINT=POINTS. To je rozsah definice, nikoli tvrzení, že všechny ostatní nativní větve explicitně končí STOP.
- Počty bodů jsou kladná celá čísla. Kontrolují se konečné nezáporné teploty, jejich pořadí, vazba DMAGDT na MLIN a úplnost teplotní sítě pro vibrace bez SETFLUCT.
- GRID=3, NE=1 nahrazeno GRID=0, NE=1. Jednobodový GRID=3 vede při inicializaci k nedefinované integrační váze; samotný GIL nastavuje ERYD=EFERMI a energetickou síť pro damping nepoužívá. EMIN/EMAX z definice odstraněny. SPLITSS odmítnuto.

## Defaults: vypsané versus nevypsané

GRID=0 a NE=1 se vypisují; změna je záměrná pro bezpečnou inicializaci. ImE=0 se vypisuje, vlastní damping však počítá na EF. SITES.NL=4 se rovněž vypisuje: zachován existující ASE default; nejde o skrytý příslib nativního defaultu. Vhodnost NL závisí na systému a konvergenci.

Teplotní a fluktuační parametry nemají vnucené Python defaults a při nezadání se nevypisují. Native NFTET/NFPHI začínají na 1, INITSVEC na 0. NVIBRA je bez požadavku 1, ale při požadavku vibrací a INITSVEC=0 jej native přepíše na kanonických 14 směrů. Například NVIBRA=26 tedy neznamená 26 směrů; pro generovanou síť je nutné zvolit INITSVEC=1. TDEBYE=0 znamená tabulkovou hodnotu podle složení, nikoli automaticky dělení nulou. Záporné TDEBYE/RMSU0 postupy nejsou touto definicí pokryty.

Pro MLIN je vyžadován explicitní TMPMIN, protože nativní výchozí počáteční hodnota 9999 není užitečným běžným zadáním. Bez TMPMAX native používá jeden bod v TMPMIN; NTMP bez TMPMAX by ignoroval, proto se odmítá.

## Významná omezení nativního výpočtu

1. GIL vyžaduje KMROT=0 a IBZINT=2 nebo 6 (gil.f:89–93). Definice zpřístupňuje standardní POINTS. KMROT a magnetické momenty závisí na potenciálu a tato validace je z něj neověřuje.
2. gilsum.f:180–204 počítá GEFF=2*(1+Morb/Mspin) a další normalizaci bez ochrany před nulovým momentem. Výpočet probíhá i před podmínkou tisku Z>0. Prázdné nebo nemagnetické komponenty a kompenzovaný celkový moment mohou způsobit neplatné výsledky. CONTROL.NONMAG=False tomu samo nezabrání.
3. INITVAR zapíná RDTAUDIR; může být potřeba kompatibilní inicializační TAU soubor. Potenciál ani externí data zde nebyly numericky ověřeny.
4. Native obsahuje MPI větev v gil.f; existující ASE nastavení paralelizace nebylo měněno ani numericky certifikováno. MCS, RDLM a NLCPA nejsou součástí tohoto auditu.

## Opora ve zdrojích

SPR-KKR 9.7: initvar.f:376–382; kkrchi.f:1425–1427; rsp.f:913–925; gil.f:89–121,266–320; thermal_init.f:79–81,105–119,170–242,322–368,548–558; mod_thermal.f:16; thermal_rmsdisp.f; gilsum.f:180–204. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee, samostatný patch pouze pro GILBERT.

## Testy

14 nových testů: discovery/defaults/pickle/executable, round-trip MLIN s DMAGDT a vibracemi, soubor M_T a ignorované parametry, úplnost teplotní sítě, kladné počty, odmítnutí neauditovaných režimů, NONMAG, nesprávné energetické sítě a chybějící MLIN teploty. Celá sada InputParameters: **122 passed**, dvě předem existující NumPy deprecation warnings. Nejde o numerický test SPR-KKR.

Příkaz: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.
