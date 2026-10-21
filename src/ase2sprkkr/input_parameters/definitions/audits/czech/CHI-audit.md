# CHI — audit současné implementace

**Stav:** Statická susceptibilita v auditovaném ASA/FREL sériovém režimu; nativní NCHIITER je efektivně 1. Audit vychází ze současné implementace na větvi `tasks-workflow` a nativního zdroje SPR-KKR 9.7 (`PUB9.7_260420_18h13`). Hodnocení není numerickou certifikací.

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

## Výchozí hodnoty a význam

**Vypsané defaults:** GRID=5, NE=30, ImE=0.01, EMIN=−0.2, MODE.NCHIITER=1. Počet energií a k-bodů vyžaduje konvergenci. NCHIITER=1 se liší od počátečního native 100, ale odpovídá skutečnému vynucenému počtu průchodů a nepředstírá iteraci. Tento explicitní default je vhodnější než nativní počáteční hodnota.

**Nevypsané:** EMAX je ponecháno native EF. FREL se zapisuje absencí přepínacího klíče. NOSYMCHI=False se nevypisuje a souhlasí s nativním hledáním přítomnosti klíče. FULLPOT=False se nevypisuje a nenutí převod full-potential potenciálu na ASA. MODE.OP bez hodnoty se nevypisuje: obecný ORBPOL a lokální CHIORBPOL nejsou totožné proměnné. SUSC při chybějícím OP nastavuje CHIORBPOL na BROOKS pro popis/volbu odezvy; nelze z toho automaticky odvozovat, že základní potenciál používá Brooksovu korekci. Existující společný MODE.OP proto nebyl uměle nastaven na BROOKS.

MODE.PRINT a TASK.CHIPRINT nemají vnucené defaults; native používá vstupní IPRINT/CHIPRINT. MIX=0.5 a TOL=1e−5 jsou nativní počáteční hodnoty, ale jejich iterační účinek je potlačen výše uvedenou chybou.

## Implementační kontrakt
- Autoritou pro soubory je aktuální `TaskDefinition` a případné dynamické metody schématu. Výstupy a cache se nesměšují s povinnými vstupy.
- `output_file_templates`: `( "{dataset}*_dens.agr", "{dataset}*_chifmag.agr", "{dataset}_*_GAMMA-OP.dat", "{dataset}_*_CHI-r.dat", "{dataset}_*_*_I-spin.dat", )`
- `cache_file_templates`: `("{dataset}.tau", "TKTK_tables")`
- Zdroj: `src/ase2sprkkr/input_parameters/definitions/chi.py`.

## Ověření

Testy definic a deklarací souborů jsou v `src/ase2sprkkr/input_parameters/test/`. Tento dokument netvrdí nové spuštění úplné sady ani numerický výpočet SPR-KKR.
