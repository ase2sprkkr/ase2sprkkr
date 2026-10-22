# PLOTPRS — audit a nová definice InputParameters

Podklad: přiložené SPR-KKR 9.7 (PUB9.7_260420_18h13) a aktuální ASE2SPRKKR `refactor`, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`, ověřený fetch. PLOTPRS zde nemá definici. Samostatný patch přidává `definitions/plotprs.py` a `test/test_plotprs.py`; nevyžaduje předchozí tasky.

## Rozhraní a rozsah

`initvar.f:187–190` při keywordu PLOTPRS zapne PLOTPRS(0:3), tedy všechny tři skupiny grafů. `gen.f:1460–1463` volá FPPLOTPRS. Samotná rutina `fpplotprs.f` nečte task-specific parametry. Definice tedy nepřidává výběr atomu, orbitalu, radial mesh, interval r ani přepínače skupin, které by zdroj nečetl.

Použití: `InputParameters.create('PLOTPRS')`. Definice používá CONTROL, TAU, MODE, STRCONST, ENERGY, TASK, SITES; executable=kkrgen, mpi=False. Rutina zapisuje společné soubory bez distribuce práce mezi ranky; MPI ani numerické spuštění zde ověřeno nebylo.

| Grafy | ASA | FULLPOT |
|---|---|---|
| Potenciál | sférické r*VT a BT pro každý atomový typ | navíc VNST a BNST podle nesférických komponent |
| Hustoty | sférická charge density, spin/orbital density při magnetickém režimu | navíc nesférické RHO2NS |
| Shape functions | žádný výstup SFN | FLMSF pro každou radiální mesh IM=1..NM |

Sférická charge křivka je RHOCHR*r²/(4*pi), nikoli prostá RHOCHR (`fpplotprs.f:254`). Radius je v a.u.; potenciálové popisky rozlišují r*V a B. Shape functions jsou koeficienty existující radiální reprezentace, nikoli 3D vizualizace tvaru buňky. NONMAG se bere ze společného CONTROL a omezuje počet panelů.

Účelem tasku je vykreslení dat načtených do radiálních polí, nikoli SCF výpočet nové hustoty. FPPLOTPRS samo shape functions nevytváří. Inicializace FULLPOT má samostatnou cestu čtení/vytváření SFN; dostupnost správné reprezentace závisí na potenciálu, mesh a souboru SFN, nikoli na novém parametru této rutiny. PLOTPRS funguje i v ASA, kde podmínka `PLOTPRS(3).AND.FULLPOT` shape-function grafy přeskočí (`:392`). Není proto důvod zavádět obecný zákaz ASA.

## Defaults a inicializace

| Volba | Definice | Zdůvodnění |
|---|---|---|
| TASK | explicitní PLOTPRS | nativně zapíná všechny tři skupiny |
| ENERGY.GRID | explicitní 1 | jednoduchá inicializační cesta; neurčuje radiální grafy |
| ENERGY.NE | explicitní 1 | INIT_MOD_ENERGY vyžaduje NEMAX>0 i při vypnutém energy loop |
| ENERGY.ImE | 0 | není použito jako smearing radiálních křivek |
| ENERGY.EMIN | explicitní -0.2 Ry | shodné s nativním EMINDEFAULT; pro grafy neúčinné |
| ENERGY.EMAX | None, nevypisuje se | nativní fallback na EF |
| SITES.NL | None, nevypisuje se | kontextový cutoff z potenciálu |
| CONTROL.SFNFIL | společný default DATASET.sfn | lze explicitně zadat vlastní filename |

`initvar.f:211–218` vypíná INITELOOP a RUNELOOP. Přesto hlavní KKRGEN volá INIT_MOD_ENERGY; ten stanoví NEMAX z NETAB a při nule ukončí program (`init_mod_energy.f:192–197`). Zcela vynechat ENERGY by proto nebylo bezpečné. NE=1 je zapsaná technická volba ASE2SPRKKR, ne tvrzení o nativním omitted defaultu. Vyšší počet bodů nepřidává body radiálního grafu; definice na něj upozorňuje. Default NE=300 používaný v jiných grafových taskech by zde nebyl lepší.

Nevypisované EMAX a NL zachovávají nativní kontext. CONTROL, TAU a STRCONST jsou sdílené definice; task při grafování nevyhodnocuje scattering-path operator. Explicitní obecné volby proto nelze interpretovat jako další parametry FPPLOTPRS.

## Omezení a selektivní varianty

GEN zakazuje DMFT a LDAU i pro PLOTPRS (`gen.f:211–222`); definice je odmítá. Ve zdroji není obecný zákaz FULLPOT ani požadavek IREL=3, proto definice připouští FREL, NREL, SREL i SP-SREL. Numerická kombinace všech režimů ověřena nebyla.

TASK.HFF může později v INITVAR nastavit IHYPER, ale pro PLOTPRS neběží obvyklý energy loop ani výpočet hyperfine výsledků. Definice upozorňuje na neúčinný požadavek; netvrdí, že se flag vůbec nenačte.

SPR-KKR rozpoznává také PLOTPOT, PLOTRHO a PLOTSFN (`initvar.f:192–196`). Tyto selektivní varianty **nejsou registrovány tímto patchem**: task je PLOTPRS se všemi skupinami. Je to výslovná hranice podpory, nikoli zamlčená ekvivalence všech variant. Samostatný PLOTRHO je navíc ve zdroji chybný, viz níže. PLOTPRS nelze kombinací falešných flags přimět k vypnutí skupin: jeho keyword je všechny zapíná.

## Nalezené chyby a rizika SPR-KKR

1. **Standalone PLOTRHO používá nealokované pole.** RHORSQ se alokuje jen uvnitř PLOTPRS(1), při kreslení potenciálu (`fpplotprs.f:62`). Hustotní větev alokuje pouze Y (`:225`), ale zapisuje RHORSQ (`:254`) a předává ho XMGRTABLE (`:280`). Pokud je zapnutá jen hustota, RHORSQ nemá alokaci. Běžné PLOTPRS tuto chybu obchází, protože první větev běží také; neznamená to, že selektivní varianta je správná.
2. **Spinová a orbitální hustota používají stejný dataset grafu.** Obě volání XMGRTABLE mají panel=1 a curve index=0 (`:282–287`). Formát tedy dostává dvě tabulky se stejným cílem místo samostatných křivek. Konkrétní efekt při importu závisí na zobrazovači; data nejsou jednoznačně označena jako dvě odlišné křivky.
3. **Limity druhého hustotního panelu ignorují orbitální hustotu.** Y1MIN/Y1MAX se určují pouze z RHOSPN (`:257–258`), ač se zapisuje také RHOORB. Orbitální křivka může přesahovat nastavený rozsah.
4. **Velké L se v některých nesférických legendách nevejde do formátu.** Potenciálová a hustotní větev používají `99001 FORMAT(I1,',',I2)` pro l,m. Pro l>=10 vzniká hvězdička místo l; SFN větev naproti tomu volí širší formát podle velikosti l,m. Obecný zákaz vysokého cutoffu kvůli kosmetické chybě nebyl zaveden.

RHORSQ je lokální allocatable; chybějící explicitní DEALLOCATE na konci není samo o sobě důkaz úniku paměti, protože Fortran lokální allocatables při návratu uvolňuje. Audit proto tento detail neoznačuje za další chybu.

Výstup používá DATSET0, odvozené od CONTROL.DATASET v hlavním programu (`kkrgen.f:239`). Soubory vznikají přes XMGRHEAD jako V, VNS, RHO, RNS a SFN_M podle skupiny. Patch neobsahuje specializovaný output parser grafů.

## Testy

12 nových případů: automatické nalezení tasku, základní textový/souborový round trip, pickle/copy, omitted EMAX/NL, nenulový počet energií pro inicializaci, všechny čtyři relativistické režimy bez vymyšleného zákazu, DMFT/LDAU, upozornění na HFF a počet energií, vlastní SFNFIL.

Celá sada `src/ase2sprkkr/input_parameters/test`: **120 passed** na výchozí revizi plus tento samostatný patch (bez předchozích patchů). `git diff --check` a obrácená aplikovatelnost patchu prošly. SPR-KKR nebyl spuštěn; testy ověřují definici inputu, nikoli numerické grafy a chování jejich importu.
