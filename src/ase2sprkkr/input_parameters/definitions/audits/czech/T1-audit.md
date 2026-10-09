# T1 — audit a nová definice InputParameters

Podklad: přiložené SPR-KKR 9.7 (PUB9.7_260420_18h13) a aktuální ASE2SPRKKR `refactor`, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`, ověřený fetch. T1 zde nemá definiční modul. Samostatný patch přidává `definitions/t1.py` a `test/test_t1.py`; není závislý na předchozích patchích.

## Vstup a režimy

`initvar.f:331–334` rozpoznává T1 a vypíná RUNELOOP. GEN při T1 odmítá FULLPOT a volá NONRELT1 pro IREL<=1, jinak RELT1 (`gen.f:1485–1496`). Rutiny nečtou task-specific parametry. Nevytvářím proto neúčinné volby ITSEL, teploty, izotopu nebo magnetického pole. Jaderná data pocházejí z interních tabulek TABGAMMAN; další tabulky v nefunkční NONRELT1 nejsou uživatelské input parametry.

Definice používá CONTROL, TAU, MODE, STRCONST, ENERGY, TASK, SITES; executable=kkrgen a mpi=False. Postprocessing RELT1 není distribuovaný mezi ranky a používá lokální TAUT. Numerické ani MPI spuštění nebylo ověřeno.

- FREL/None: vstup pro RELT1; výchozí režim.
- NREL a SREL: NONRELT1 obsahuje `IF (NT.GT.0) STOP_MESSAGE('not available')` (`nonrelt1.f:83–84`). Definice tyto režimy odmítá, protože pro skutečný systém NT>0 nejsou dostupné.
- SP-SREL: GEN jej posílá do RELT1, ale rutina pracuje s relativistickými kappa kanály. Definice upozorňuje na neověřenou kombinaci, nevymýšlí explicitní nativní zákaz.
- DMFT/LDAU: zakazuje GEN (`gen.f:211–222`); definice také.
- FULLPOT: zakazuje check_input i GEN. Definice zákaz dokumentuje; samotný input bez načtení potenciálu režim FULLPOT spolehlivě neurčí.

## Energie a defaults

`init_mod_energy.f:265–275` pro T1 natvrdo nastavuje GRID=3, NETAB(1)=1, EMIN=EMAX=EF a ETAB(1,1)=EF+i*EIMAG. EIMAG začíná pro tento task na 0 a může být přepsáno ENERGY.ImE. RELT1 si navíc vytvoří CEFERMI=DCMPLX(EFERMI,EIMAG) pro SSITE (`relt1.f:159–171`). Energetické okno ani spektrální sken zde nejsou účinnými parametry.

| Volba | Definice | Vztah k SPR-KKR |
|---|---|---|
| ENERGY.GRID | explicitní 3 | stejně se přepisuje na 3 |
| ENERGY.NE | explicitní 1 | stejně se přepisuje na 1; nenulová inicializace nutná před resetem |
| ENERGY.ImE | 0 | shodný nativní default tasku; uživatelské nenulové ImE se skutečně čte |
| ENERGY.EMIN/EMAX | nejsou v definici | T1 je vždy přepíše na EF; není důvod přidávat zavádějící interval |
| SITES.NL | None, nevypisuje se | cutoff určuje potenciál; explicitní NL>5 odmítnuto kvůli kapacitě CHK |
| CONTROL.PRINT | explicitní 1 | záměrný rozdíl od obecného 0 kvůli inicializační chybě NN/MM |
| CONTROL.RUNELOOP | explicitní flag true | přepisuje nativní vypnutí smyčky, aby vznikla TAUT |

GRID/NE odlišné od {3}/{1} definice pouze varuje jako přepsané. ENERGY.NE stále musí být kladné: INIT_MOD_ENERGY alokuje ETAB a kontroluje NEMAX>0 ještě před task-specific T1 resetem.

**PRINT=1 a RUNELOOP musí být vypsané.** Nejde o omitted defaults, které by měly odpovídat nativním: jsou to explicitní volby ASE2SPRKKR obcházející dva níže uvedené inicializační problémy. V tomto případě nativní defaults nejsou lepší. Definice při PRINT=0 nebo vypnutí RUNELOOP odmítne uložení. Vstup pro využití již připravených TAU dat s jiným řízením smyček tato definice neřeší.

## Závažné chyby a omezení SPR-KKR

1. **Neinicializované NN/MM při PRINT=0.** RELT1 nastaví NN=2*NL² a MM=NKMMAX pouze uvnitř `IF (IPRINT.GE.1)` (`relt1.f:150–154`), ale později je vždy předá do CMATSTRUCT (`:187`). Obecný PRINT=0 proto může vést k chybnému přístupu k paměti při diagnostickém výpisu. Explicitní PRINT=1 v definici zajišťuje inicializaci; RELT1 jej následně samo nastaví na 0.
2. **T1 standardně nevypočte TAUT, které pak používá.** INITVAR vypne RUNELOOP pro T1 (`:331–334`). RELT1 počítá SSITE znovu, ale to řeší single-site TSST, nikoli TAUT; používá `W2(:,:)=TAUT(:,:,IT)` (`relt1.f:206`). TAUT se získává v energetické smyčce GEN, která je za `IF (RUNELOOP)` (`gen.f:534`). Pro samostatný T1 běh bez připravené TAU matice je to vadná inicializační cesta. Explicitní CONTROL.RUNELOOP aktivuje pozdější nativní override (`initvar.f:689–693`) a umožní smyčku na jednom bodě EF. Tento postup je ověřen čtením zdrojů, ne numerickým během. Při ručně přidaném RDTAU a existujících TAU souborech může hlavní program měnit flags znovu; tento režim nebyl prověřen.
3. **Záporné kappa jako horní mez smyčky.** V součtu relaxačních příspěvků se počítá podepsané KAP0/KAP1/... a používají se smyčky `DO M0=1,2*KAP0` a obdobné pro M1/M2/M3 (`relt1.f:215–245`). Záporné kappa kanály nemají žádnou iteraci; například s1/2 s kappa=-1 se vynechá. Počet m_j stavů má souviset s absolutní hodnotou kappa. Nejde jen o chybějící příspěvky: posouvané LAM indexy pak neodpovídají celé kappa bázi. Samotná definice inputu tuto numerickou chybu neopraví.
4. **Omezená kapacita popisků.** CHK má 10 položek (`relt1.f:39,54–57`), ale výpisy používají NKQ=2*NLQ-1 (`:338–347`). Pro NLQ>=6 překročí pole. Definice odmítá explicitní max(SITES.NL)>5; při omitted NL musí být cutoff zkontrolován v potenciálu. Neprohlašuji implicitní kontext za bezpečný bez jeho načtení.
5. **NONRELT1 je záměrně zablokované.** Jeho rozsáhlý kód za STOP není důkaz dostupného nerelativistického výpočtu. Definice neprezentuje tyto větve jako podporované.

Před použitím relaxačních rychlostí pro fyzikální závěry je zejména nutné opravit podepsané kappa smyčky a porovnat numerické výsledky. Patch přidává vstupní rozhraní a inicializační ochrany, nikoli validovaný řešič T1.

## Výstup a testy

RELT1 vypisuje tabulky a diagnostické matice do standardního výstupu. Specializovaný parser těchto tabulek součástí patchu není; používá se obecný output reader.

10 nových testů: discovery, explicitní PRINT/RUNELOOP, file/text round trip, pickle/copy, nezavádění neúčinných energetických mezí, nedostupné NREL/SREL, zakázané korelace, PRINT=0 a vypnutá smyčka, přepis GRID/NE, kapacita kappa popisků a účinné ImE.

Celá sada `src/ase2sprkkr/input_parameters/test`: **118 passed** na uvedené výchozí revizi plus tento samostatný patch (bez předchozích patchů). `git diff --check` i obrácená aplikovatelnost patchu prošly. Numerický běh SPR-KKR nebyl proveden; testy ověřují input definici.
