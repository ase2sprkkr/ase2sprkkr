# 36. CHILANDAU — SPR-KKR 9.7 / ASE2SPRKKR

Dosud chybějící task pro Landauovu orbitální magnetickou susceptibilitu. Přidána definice pro **sériový periodický ASA FREL postup**. Předběžně: sériová větev funguje se zásadními omezeními; **MPI větev je rozbitá**. Numerický výpočet nebyl proveden. Nejde o společný závěrečný audit.

## Definice a validace

Executable kkrchi; nativní INITVAR rozpoznává CHILANDAU, RSP volá SUSC_LANDAU. Přidány CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK a SITES; minimální volitelná SCF sekce umožňuje odmítnout FULLPOT, který RSP explicitně zakazuje.

Vyžadován BZINT=POINTS (native IBZINT=2), jeden GRID=5 kontur bez SPLITSS a kladný počet energií. SUSC_LANCALC jinak při NEPATH!=1 nebo GRID!=5 tiše vrací bez výpočtu. Stejně vrací při KMROT!=0; KMROT je vlastnost nativního stavu/potenciálu a tato definice jej nezjišťuje. FREL je omezení auditovaného postupu: angular a wavefunction rutiny používají relativistické CGC a spinorové kanály, nejde o nalezený explicitní STOP pro všechny ostatní režimy. NONMAG je přípustný — nejde o spinovou susceptibilitu vyžadující magnetický základní stav.

TASK obsahuje CLURAD, NSHLCLU, SELECT, EXCLUDE, QREPONLY a PRINT. SELECT/EXCLUDE jsou vzájemně výlučné; radius je konečný nezáporný. NSHLCLU včetně nuly má přednost před radiem. Definice nenabízí všechny možné nativní modifikace výpočtu.

## Defaults

**Vypsané:** GRID=5, NE=30, ImE=0.01, EMIN=−0.2 a CLURAD=2.2. Jde o explicitní praktické volby ASE; počet bodů a radius vyžadují konvergenci. Nativní INIT_MOD_TAUIJ_STAR nemá použitelný default radiusu: bez CLURAD i NSHLCLU končí STOP. Proto je vypsaný radius lepší než jeho vynechání, ale hodnota 2.2 není fyzikální zárukou konvergence.

**Nevypsané:** EMAX je skutečně ponecháno nativnímu EF; QREPONLY=False souhlasí s MOD_SITES.TREAT_QREP_ONLY=False. NSHLCLU není vnuceno na nulu, protože jeho pouhá přítomnost mění výběr clusteru. PRINT je nezadaný: SUSC_LANDAU po inicializaci nastavuje IPRINT=0 a z něj odvozuje CHIPRINT. FULLPOT=False není vypsán; pokud potenciál sám aktivuje FULLPOT, native RSP jej i tak odmítne. Python default sám nepřepíná full-potential potenciál na ASA.

## Konkrétní chyby a omezení nativních zdrojů

1. **MPI: neinicializovaná pole.** SUSC_LANCALC alokuje SAVE angular pole pouze při IECURR=1. SUSC_LANDAU přitom dělí energie mezi ranky. Rank počítající jiné energie nemá tato pole alokovaná (susc_lancalc.f:66–88; susc_landau.f:193–214).
2. **MPI: redukce za koncem pole.** CHILANLT má tvar (0:NLMAX,NTMAX,2), ale MPI_SUM_R začíná na CHILANLT(0,1,2) a předává délku obou bloků. Od začátku druhého bloku zbývá pouze polovina požadované délky (susc_landau.f, blok finálního MPI součtu). To je nezávislá paměťová chyba. mpi=False je proto záměrný default definice; explicitní vynucení MPI přes calculator tuto chybu neopraví.
3. **Průběžné výstupy.** Podmínka IECURR.LT.-NETAB(1) v SUSC_LANCALC je pro kladné energetické indexy vždy nepravdivá. Výstup susceptibilit se tiskne při každé energii, nikoli až na konci. V sériovém běhu je nutné použít poslední blok; dřívější bloky jsou pouze částečné integrály.
4. **Čtecí TAU větev.** RDTAUFIL otevře a přečte hlavičku textového souboru, ale v energetické smyčce chybí odpovídající čtení TAU a výpočet TAU je při RDTAUFIL přeskočen. Tento reuse postup není touto definicí nabízen ani certifikován.
5. Výpočet používá konečný součet sousedů a konkrétní kartézské složky. Konvergenci, geometrický rozsah, násobnosti při QREPONLY a interpretaci pro obecnou anisotropii je potřeba ověřit. Nižší dimensionalita a nestandardní rotační stavy nejsou tímto auditem pokryty.


## Ověření

12 nových testů: discovery/defaults/serial executable, round-trip, pickle, NSHLCLU=0 a selektory, NONMAG, odmítnutí neauditovaného režimu, FULLPOT, BZINT, GRID, SPLITSS, neplatného radiusu a NE, vzájemná výlučnost SELECT/EXCLUDE. Celá sada InputParameters: **120 passed**, dvě existující NumPy deprecation warnings.

Příkaz: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

Zdroje: initvar.f:453–454; rsp.f:134–147; susc_landau.f; susc_lancalc.f; susc_laname.f; susc_lanme.f; init_mod_tauij_star.f:150–180; mod_sites.f:73. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Samostatný patch pouze pro CHILANDAU.
