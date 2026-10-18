# BONDING — audit a nová experimentální definice

Zdroj: dodané SPR-KKR 9.7, ASE2SPRKKR refactor po fetchi, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Task dosud nemá definici. Patch 0021-bonding.patch přidává bonding.py a test_bonding.py samostatně proti tomuto základu.

## Vstup

INITVAR:505–512 nastavuje TASK=BONDING, vypíná čtení/zápis TAU directory a zapíná core/hyperfine. GEN:1013–1025 povoluje CLUBONDING pouze pro IBZINT=0. TAU.BZINT=CLUSTER nastaví tento režim; INIT_MOD_KMESH:717 přepne běžný výpočet na REAL-SPACE. CLUBONDING samo žádné TASK parametry nečte.

Definice vystavuje konzervativní molekulární variantu s explicitním TAU.MOL. TAU_REAL_SPACE pro MOL pracuje s IQ1=1 až IQ2=NQ; jeho čtečka i v tomto režimu žádá NSHLCLU nebo CLURAD. Definice proto vypisuje NSHLCLU=0. Scratch výstup vyžaduje jednu clusterovou kopii každého site (TAU_REAL_SPACE:663 a748); běžný periodický cluster nemusí tuto podmínku splnit. Rozšíření na něj není tímto patchem ověřeno.

Plná relativita je konzervativní omezení definice: CLUBONDING sice při čtení má větev IREL<=2, ale následně bez podmínky používá transformace REL>CLM a dvouspinové indexování. Funkčnost skalárního workflow z této větve nelze doložit. KMROT!=0 a NEPATH!=1 způsobují tichý RETURN (53–54). Jedna energetická cesta se ověřuje validátorem, KMROT závisí na potenciálu.

## Defaulty a zápis

| Volba | Nativní chování | Definice |
| --- | --- | --- |
| TAU.BZINT | běžná integrace nevyhovuje BONDING | explicitně CLUSTER |
| TAU.MOL | false při nepřítomnosti | explicitně true |
| TAU.NSHLCLU | pokud chybí, vyžaduje CLURAD | explicitně 0 |
| ENERGY.GRID | BONDING nemá vlastní větev defaultu | explicitně 3, energetický scan |
| ENERGY.NE | MOD_ENERGY inicializuje 0; při NEMAX=0 STOP | explicitně 300 |
| ENERGY.ImE | 0,01 Ry | explicitně 0,01 Ry |
| ENERGY.EMIN | -0,2 Ry | explicitně -0,2 Ry |
| ENERGY.EMAX | sentinel -99,99 Ry se v INIT_MOD_ENERGY:182 přepíše na EF | explicitně 1 Ry |
| MODE | běžně FREL | vynecháno, None |

MOL a NSHLCLU mají v obecné TAU definici expert režim, který defaulty ne vždy zapisuje. V tomto modulu je jejich zápis vynucen: vynechané MOL by neodpovídalo Python defaultu true a vynechané NSHLCLU by vedlo ke čtení chybějícího CLURAD. Regresní test kontroluje přítomnost obou voleb.

EMAX se záměrně vypisuje; při vynechání by se převzalo EF. Původní audit přehlédl obecné přepsání sentinelu v INIT_MOD_ENERGY:182; toto je opravená verze. Explicitní rozsah odpovídá běžnému DOS scanu ASE2SPRKKR; jeho fyzikální vhodnost a konvergence nejsou ověřeny. Nativní NE=0 není použitelný default. Nativní horní mez EF by naopak byla vhodnější pro analýzu omezenou na obsazené stavy; definice volí scan i nad EF. NE musí být alespoň 2, protože EPATH GRID=3 používá ESTEP i při NE=1, kdy nebyl přiřazen.

## Chyby a omezení nativního workflow

1. **Chybějící pomocné soubory:** TAU_DRIVE:70 nastavuje KTAUIJ=0. TAU_REAL_SPACE vytváří scratch jednotku18 jen pro KTAUIJ>=2 a jednotku19 jen pro KTAUIJ!=0 (657–674,741–766). CLUBONDING je přesto rewinduje a čte (78–108). Běžný tok tedy nepřipraví vstupy této rutiny. Oprava vyžaduje zásah do Fortranu; definiční soubor ji neřeší.
2. **Centrální site je pevně 1:** CLUBONDING:228 nastavuje IQ=1 a zpracovává pouze jeho typy a sousedy. Není zde čtečka TASK.IQ ani obecný výběr centra analýzy.
3. **Pevná výstupní jednotka22:** řádky355–358 zapisují bez místního OPEN a bez DATASET názvu. V běžném Fortran runtime může vzniknout fort.22. Výstup neposkytuje identifikátor typu ani hlavičku jednotlivých řádků.
4. **Shell indexy po vložení:** při vložení nového poloměru do DSH se posunou hodnoty DSH, ale ISHQ již přiřazených site se neaktualizují (164–173). Výpis shellů tak může být nesprávný. Výpočet příspěvků níže prochází přímo JQ, proto nelze bez další analýzy tvrdit, že tato diagnostická chyba mění všechny DOS hodnoty.
5. **Rozměr lokální matice:** N je přiřazeno v předchozí smyčce přes typy podle jejich IQAT. Při následném pevném IQ=1 se nepřepočítá. Různé NKMQ mohou znamenat nesprávný rozsah operací.
6. **MPI:** rutiny používají pevné scratch jednotky, lokální uložené shell tabulky inicializované pouze při IECURR=1 a výstup22 bez redukce. Modul proto volí mpi=False; paralelní workflow nebylo doloženo.

Podpora je výhradně definice vstupu. SPR-KKR 9.7 BONDING nelze na základě tohoto auditu prohlásit za funkční numerický task.

## Testy

12 nových případů: serializace a načtení, pickle, explicitní zápis nutných clusterových voleb, odmítnutí jiné integrace/MOL/relativity, více cest a NE<2. Celá sada input_parameters: **120 passed**. Numerický SPR-KKR nebyl spuštěn. Patch prošel git diff --check a kontrolou zpětného aplikování.
