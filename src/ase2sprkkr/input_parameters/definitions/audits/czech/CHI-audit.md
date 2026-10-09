# 37. CHI — audit SPR-KKR 9.7 a ASE2SPRKKR

Dosud chybějící task. Přidána definice pro sériový periodický ASA FREL výpočet magnetické susceptibility. **Předběžně: funguje, ale má velké ALE — iterace spinové magnetizace je v této verzi zdrojů vypnutá.** Numerická správnost nebyla ověřena. Toto je individuální audit, nikoli společný závěrečný audit.

## Mapování a rozsah

TASK CHI rozpoznává INITVAR a vypíná RDTAUDIR/WRTAUDIR; RSP pro ASA volá SUSC. Executable je kkrchi. Definice obsahuje CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES a minimální volitelnou SCF sekci pro validaci FULLPOT.

Auditovaný postup omezen na FREL, BZINT=POINTS, jeden GRID=5 kontur bez SPLITSS, kladný NE a OP nezadané/NONE/BROOKS. Je to vymezení běžného postupu, nikoli tvrzení o explicitním STOP pro všechny jiné nativní varianty. NONMAG je povolen: paramagnetická susceptibilita může být počítána z nemagnetického základního stavu. Potenciál, rotační stavy a další nastavení specifická pro strukturu nebyla numericky ověřena.

TASK nabízí CHIPRINT a NOSYMCHI. MODE.PRINT se čte znovu uvnitř SUSC a může přepsat CHIPRINT převzatý z RSP. MODE.NCHIITER se povinně zapisuje jako 1 a jiné hodnoty jsou odmítány; MODE sekce je kvůli tomuto zápisu povinná.

## Chyby a zásadní omezení nativního kódu

- **Iterace neproběhne:** SUSC čte MODE.NCHIITER (native počáteční hodnota 100), MIX a TOL, ale na susc.f:690 bezpodmínečně nastavuje NCHIITER=1. Následující větev IF(NCHIITER.GT.1), míchání a kontrola konvergence jsou tak nedosažitelné. Výsledek nelze prezentovat jako konvergovaný samokonzistentní výpočet odezvy. Definice MIX/TOL nenabízí, protože by v tomto postupu slibovaly neúčinné řízení iterace.
- **MPI není dostupné:** susc.f:139–147 při MPI volá STOP_MESSAGE s informací o rozpracované implementaci. mpi=False je nutný default; explicitní vynucení MPI v calculatoru tento STOP neobejde.
- **Veřejná FULLPOT větev není dostupná:** rsp.f:813–828 při PUBLIC_VERSION končí STOP_REGULAR. Nestačí pouze zadat STAT. FULLPOT se proto odmítá; pokud jej aktivuje potenciál, rozhodující je nativní kontrola.
- Automatický DOS na EF se počítá s pevně daným Im(E)=0.001 Ry (susc.f:269), nezávisle na ENERGY.ImE. To je nativní volba, nikoli chybně předávaná hodnota ASE.
- Pokročilé RDABXC/ITRCHI0 restarty, ruční NEF pro typy, SIGNAF a cluster/surface/slab postupy nejsou zpřístupněny. RDABXC/MIX/TOL/NCHIITER se čtou z MODE, ne z TASK — důležité pro případné budoucí rozšíření.

## Defaults a zápis

**Vypsané defaults:** GRID=5, NE=30, ImE=0.01, EMIN=−0.2, MODE.NCHIITER=1. Počet energií a k-bodů vyžaduje konvergenci. NCHIITER=1 se liší od počátečního native 100, ale odpovídá skutečnému vynucenému počtu průchodů a nepředstírá iteraci. Tento explicitní default je vhodnější než nativní počáteční hodnota.

**Nevypsané:** EMAX je ponecháno native EF. FREL se zapisuje absencí přepínacího klíče. NOSYMCHI=False se nevypisuje a souhlasí s nativním hledáním přítomnosti klíče. FULLPOT=False se nevypisuje a nenutí převod full-potential potenciálu na ASA. MODE.OP bez hodnoty se nevypisuje: obecný ORBPOL a lokální CHIORBPOL nejsou totožné proměnné. SUSC při chybějícím OP nastavuje CHIORBPOL na BROOKS pro popis/volbu odezvy; nelze z toho automaticky odvozovat, že základní potenciál používá Brooksovu korekci. Existující společný MODE.OP proto nebyl uměle nastaven na BROOKS.

MODE.PRINT a TASK.CHIPRINT nemají vnucené defaults; native používá vstupní IPRINT/CHIPRINT. MIX=0.5 a TOL=1e−5 jsou nativní počáteční hodnoty, ale jejich iterační účinek je potlačen výše uvedenou chybou.

## Testy

**11 nových testů**, celá sada InputParameters **119 passed**, dvě existující NumPy deprecation warnings. Testy ověřují discovery, sériový executable/defaults, skutečný zápis NCHIITER=1, parse/write round-trip, pickle, obě print volby a NOSYMCHI, OP=BROOKS, přípustný NONMAG a odmítnutí FULLPOT/neauditovaných režimů, nesprávných sítí a počtů i neúčinné iterace.

Příkaz: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

Zdroje: initvar.f:412–418; rsp.f:195–201,811–835; susc.f:139–207,240–304,334–343,690–781. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Samostatný patch pouze pro CHI. Testy nespouštějí SPR-KKR.
