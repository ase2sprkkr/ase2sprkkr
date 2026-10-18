# NUCSSCPL — audit a experimentální definice

Porovnání dodaných zdrojů SPR-KKR 9.7 s fetchovanou větví ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. NUCSSCPL vlastní definici nemá. Samostatný patch 0022-nucsscpl.patch přidává nucsscpl.py a test_nucsscpl.py proti tomuto základu.

## Účel a vstup

NUCSSCPL počítá tensor jaderné spin-spin vazby A_ij a zvlášť Rudermanův–Kittelův s–s příspěvek R_ij. INITVAR:491–496 nastaví TASK=NUCSSCPL, RUNELOOP=true a ICORE=0. Nejde o přejmenování na E-INTEGRAL.

GEN:358–386 připravuje off-diagonal TAUIJ přes INIT_MOD_TAUIJ_STAR a vyžaduje IBZINT=2 (BZINT=POINTS). TASK proto obsahuje CLURAD, NSHLCLU, SELECT a EXCLUDE ze společné čtečky. Pokud je zadán NSHLCLU, má přednost před CLURAD; oba selektory současně vedou k STOP. Seznamy mají kapacitu100 a obsahují chemické symboly.

TASK.QREPONLY nebyl přidán: INITVAR jej čte uvnitř větve CALC_FMAG/JXC/DMI/JALF/JORB/BSDOSIJ, kam samotný NUCSSCPL nevstupuje. Zkopírovat sem stejnou nabídku jako pro BSDOSIJ by bylo chybné.

NUCSSCPL:77–83 požaduje NONMAG, IREL=3, KMROT=0, jednu energetickou cestu, bodové jádro a GRID=5. GEN:1043–1047 zastavuje pro FULLPOT. Definice validuje NONMAG, MODE, CONTROL.FINNUC, energetickou cestu, BZINT, nekladné NE a výlučnost selektorů. FULLPOT z potenciálu a KMROT se musí ověřit při běhu; Python definice je nezná.

Izotopy se nenastavují vstupem tasku. TABMUQIN (table.f) vybírá pro Z tabulkový izotop s nenulovým magnetickým momentem; tabulka pokrývá Z=0..110 a obsahuje také nulové položky. Definice proto nevymýšlí volby ISOTOPE, MNUC nebo INUC.

## Defaulty a jejich zápis

| Volba | SPR-KKR při nepřítomnosti | Definice |
| --- | --- | --- |
| CONTROL.NONMAG | běžně false, tento task potom STOP | true, explicitně vypisováno |
| CONTROL.FINNUC | false | false, vynecháno |
| MODE | běžně FREL, IREL=3 | None, vynecháno |
| TAU.BZINT | běžná obecná volba | POINTS, explicitně |
| TASK.CLURAD | bez NSHLCLU nemá použitelný default, STOP | 2,2 ALAT, explicitně |
| TASK.NSHLCLU | nepřítomnost přesměruje na CLURAD | vynecháno |
| TASK.SELECT/EXCLUDE | prázdné seznamy | vynecháno |
| ENERGY.GRID | NUCSSCPL nemá vlastní inicializační větev; rutina vyžaduje5 | explicitně5 |
| ENERGY.NE | MOD_ENERGY inicializuje0; NEMAX=0 vede k STOP | explicitně30 |
| ENERGY.EMIN | -0,2 Ry | explicitně -0,2 Ry |
| ENERGY.EMAX | obecná náhrada sentinelu za EF, INIT_MOD_ENERGY:182 | vynecháno, None |
| ENERGY.ImE | 0,01 Ry | explicitně0,01 Ry |

Explicitní rozdíly NONMAG, GRID, NE a CLURAD připravují požadovaný vstup; vynechané FINNUC, MODE a EMAX odpovídají nativnímu chování. Poloměr2,2 není prokazatelně konvergovaný a má se zkontrolovat proti větším clusterům. NE=30 je běžná integrační volba, numerická konvergence neověřena. Pro GRID=5 EPATH tvoří Gaussovu kvadraturu na komplexním oblouku a ImE v této větvi nepoužívá; zde tedy nejde o fyzikální rozdíl proti společnému explicitnímu ImE=0 v ASE2SPRKKR.

## Chyby ve Fortranu

1. **Odalkování nealokovaných polí:** na každé energii kromě poslední se provede DEALLOCATE(AITJT,RITJT,DITJT,NITJT,DVECITJT,XVEC,YVEC) na řádku341. Jejich první ALLOCATE je až v dokončovací větvi366–368. NE=30 tedy narazí na odalokování nepřipravených polí už po první energii. Změna NE na1 není věcná oprava integračního workflow.
2. **Ztráta matic předchozích typů:** smyčka DO IT=1,NT na řádku115 volá CINIT(M*M*3*NT,ME) na133. Nuluje celé ME pro všechny typy při každém IT; po smyčce zůstane naplněn pouze poslední typ. Párové příspěvky ostatních typů se tím ztrácejí.
3. **Přepsaný radiální mesh:** IM=IMT(IT) je ihned nahrazeno IM=1 (118–120). Čtení a integrace všech typů tak používají JRWS(1), nezávisle na jejich přiřazení mesh.
4. **Nebezpečné diagnostické dělení:** CHECK je pevně true. RAT dělí MEZZ/ME pro nenulové XA+XB, což nezaručuje nenulový jmenovatel ME (198–200). Výstupní filtry AIJ/A11 a RIJ/R11 (440 a454) zase dělí nulovým referenčním diagonálním členem.
5. **Chybné měřítko grafů:** po smyčce přes páry se pro grafy znovu používá FTP naposledy nastavené pro konkrétní dvojici. Více izotopových momentů vyžaduje faktor podle daného IT,JT. Navíc YMIN=YMAX*FTP a následně YMAX=YMIN*FTP (547–548) přepisují dolní mez horní a horní mez násobí dvakrát.
6. **MPI není dokončeno:** AIJINT/RIJINT se alokují pouze při IECURR=1, bez MPI redukce a bez rozlišení procesu ve výstupu. Energie rozdělené mezi procesy nemají společnou inicializaci ani úplný integrál. Definice používá mpi=False.

Hlavní soubor se jmenuje DATASET_A_ij.dat; další výstupy jsou pevné RK.dat a tab_aij.dat a grafy Aij/Rij. Není implementován Python parser výsledků. Kvůli výše uvedeným chybám podpora vstupu neopravňuje k tvrzení, že nativní numerický task funguje.

## Testy

14 nových testovacích případů: serializace, načtení a pickle; explicitní NONMAG a vynechané FINNUC/EMAX; odmítnutí magnetického systému, jiné relativity, FINNUC, jiné kontury a split cesty, jiné BZINT; seznamy selektorů a shell count; nekladné NE. Celá sada input_parameters: **122 passed**. Numerický SPR-KKR nebyl spuštěn. Patch prošel git diff --check a kontrolou zpětného aplikování.

Při kontrole společné čtečky byl opraven také předchozí audit BONDING: vynechané EMAX se tam přebírá z EF, nikoli ze sentinelu. Jeho explicitní EMAX=1 Ry zůstává zvoleným scanem.
