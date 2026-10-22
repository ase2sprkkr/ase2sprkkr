# PLOT2DPR — audit a nová definice InputParameters

Podklad: přiložené SPR-KKR 9.7 (PUB9.7_260420_18h13) a aktuální ASE2SPRKKR `refactor`, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`, ověřený fetch. PLOT2DPR zde nemá definici. Patch samostatně přidává `definitions/plot2dpr.py` a `test/test_plot2dpr.py`; předchozí patche nejsou nutné.

## Vstup a varianty

`initvar.f:199–205` rozpoznává PLOT2DPOT, PLOT2DRHO a PLOT2DPR. Poslední keyword zapne obě skupiny. GEN volá ve FULLPOT FPPLOT2DPR, v ASA PLOT2DPR_ASA (`gen.f:1465–1479`).

Vytvoření: `InputParameters.create('PLOT2DPR')`. Selektivní variantu lze nastavit přes `p.TASK.TASK = 'PLOT2DPOT'` nebo `'PLOT2DRHO'`. Tento patch neregistruje samostatné moduly pro `create('PLOT2DPOT')` a `create('PLOT2DRHO')`; podporuje jejich keyword v definici PLOT2DPR a načtení vstupu přes `from_file`.

| Parametr TASK | Default při vynechání | Definice |
|---|---|---|
| AVEC | {2,0,0} | tři reálné složky, default se nevypisuje |
| BVEC | {0,2,0} | tři reálné složky, default se nevypisuje |
| DVEC | {0,0,0} | tři reálné složky, default se nevypisuje |
| NA | 50 | kladný počet intervalů, default se nevypisuje |
| NB | NINT(NA*length(BVEC)/length(AVEC)) | None, nevypisuje se; výpočet delegován SPR-KKR |

Vektory jsou kartézské v jednotkách ALAT, **nikoli zlomkové souřadnice v bázi mřížkových vektorů**. AVEC a BVEC tvoří hrany obdélníkového okna; DVEC jej posouvá od společného počátku. Ve starší pravidelné FULLPOT větvi NA/NB označují intervaly a mřížka má (NA+1)*(NB+1) vrcholů. Skutečně volané PLY rutiny mají také vlastní radiální/úhlové sampling parametry pevně v kódu: NA/NB proto nelze vydávat za úplnou kontrolu hustoty výsledné PLY triangulace.

Zdroj defaults: `plot2dpr_asa.f:53–54,64–85`, `plot2dpr_fp_ply.f:76–77,126–149`, `fpplot2dpr.f:63,102–125`.

## Geometrická validace

Obě větve při nenulovém skalárním součinu upraví BVEC projekcí kolmo na AVEC a zachovají jeho původní délku. Definice na tuto změnu upozorní, původní input sama nepřepočítává. Kontroluje konečnost vektorů a odmítá nulové nebo rovnoběžné hrany: zdroj by při normalizaci dělil nulou. Téměř rovnoběžné hrany v relativní toleranci 1e-12 jsou odmítnuty jako numericky degenerované.

Pokud NB není vypsáno, zdroj může pro velmi úzké okno zaokrouhlit NB na 0 a následně dělit nulou. Definice před zápisem ověří nativní výpočet pro kladný poměr a vyžaduje jeho výsledek >=1; případně je nutné explicitně nastavit NB. NB nelze nahradit nevypisovaným konstantním defaultem 50: při jiných délkách vektorů by byl takový default nesprávný.

Fortranový array reader dovoluje kratší pole a kopíruje poslední složku do zbývajících pozic; definice vyžaduje explicitně tři složky, aby geometrie nebyla skrytě doplňována.

## Společné defaults a režimy

Definice obsahuje CONTROL, TAU, MODE, STRCONST, ENERGY, TASK a SITES. Používá `kkrgen`, `mpi=False`: rutiny zapisují sdílené pevné názvy souborů a nemají vlastní distribuci práce mezi ranky. Numerický ani MPI běh zde ověřen nebyl.

ENERGY používá explicitní GRID=1, NE=1, ImE=0, EMIN=-0.2 Ry. Nejde o parametry rozlišení prostorového grafu, ale o inicializační volby: INIT_MOD_ENERGY i pro task s vypnutým energy loop vyžaduje NEMAX>0 (`init_mod_energy.f:192–197`). Vyšší NE nezlepšuje prostorové rozlišení; definice na něj upozorňuje. EMIN odpovídá nativnímu defaultu. Nevypisované EMAX=None deleguje EF a SITES.NL=None kontext potenciálu. Vypsané defaults lze měnit; není důvod přebírat nulový počáteční NETAB místo použitelné inicializace.

DMFT a LDAU explicitně zakazuje GEN (`gen.f:211–222`); definice je odmítá. ASA ani FULLPOT není obecně zakázáno: mají samostatné rutiny s odlišnými výstupy. Režim se určuje z potenciálu, nikoli přidaným TASK přepínačem. Zdroj nezavádí obecný požadavek IREL=3.

## Skutečné chování a problémy SPR-KKR

1. **FULLPOT končí v experimentální PLY rutině.** FPPLOT2DPR bezpodmínečně volá PLOT2DPR_FP_PLY (`fpplot2dpr.f:71`). Ta při aktivním potenciálovém/hustotním tasku zapisuje PLY soubory a volá STOP_REGULAR s textem `run of test version completed` (`plot2dpr_fp_ply.f:926`). Následující pravidelná mřížka FPPLOT2DPR se za této cesty neprovede. Nelze tedy slibovat její soubory `_V`, `_B`, `_RHOCHR` jako běžný výstup této revize.
2. **Selektivní flags v experimentální FULLPOT cestě nerozlišují jednotlivé skupiny.** PLOT2DPR(1) a (2) používá PLY rutina pouze pro společné PLOT_V_RHO (`:101`). Další výpočty potenciálu i hustoty nejsou podmíněné jednotlivými flags. PLOT2DPOT/PLOT2DRHO tak nemají totéž selektivní chování jako následující starší větev. Je třeba hodnotit konkrétní PLY soubory, ne dokumentaci staršího grid writeru.
3. **ASA ignoruje výběr hustoty.** PLOT2DPR_ASA nepoužívá PLOT2DPR flags a importuje VT, nikoli hustotní pole. Zapisuje `PLY-V-ASA-*.ply` a `POT-ASA-LINES.ply`. I PLOT2DRHO tedy v ASA vyrobí potenciálové soubory. Definice při této variantě varuje; protože input sám nerozpoznává obsah potenciálu, varování popisuje obě větve.
4. **Starší FULLPOT větev obsahuje další předčasný STOP.** Po zápisu `FULLPOT-pot.ply` následuje prostý STOP (`fpplot2dpr.f:671`); pokud by se deaktivovalo ukončení v PLY rutině, potenciálový task by stále nedošel k pozdějším tabulkám a hustotám. Jde o druhou překážku ve zdroji, nikoli nyní běžně dosaženou cestu.
5. **Nekontrolovaná degenerovaná geometrie.** Orthogonalizace rovnoběžných hran normalizuje nulový zbytek; nulové hrany a NA/NB=0 vedou k dalším dělením nulou. Definice chrání explicitní input geometrickou validací.
6. **Pevné výstupní názvy.** PLY rutiny používají názvy nezávislé na CONTROL.DATASET. Opakované běhy ve stejném adresáři mohou kolidovat; nelze z DATASET odvozovat izolaci těchto výstupů. Specializovaný output parser ani změna výstupních názvů součástí patchu není.

Výstupy jsou založené na načtených radiálních datech a interpolaci geometrického okna. Počet bodů a detaily experimentálního FULLPOT modelu mají další pevné konstanty, které se nevydávají za uživatelské parametry. Numerická správnost PLY výstupu nebyla ověřena.

## Ověření

16 nových případů: discovery, omitted geometrické defaults, souborový/textový round trip, pickle/copy, vlastní okno, selektivní keywords včetně autodetekce při načtení, upozornění na ASA hustotu, nulové/rovnoběžné vektory, orthogonalizace, zaokrouhlení NB na nulu, kladné NA/NB, rozměr vektorů a zakázané korelace.

Celá sada `src/ase2sprkkr/input_parameters/test`: **124 passed**, na výchozí revizi plus tento samostatný patch (bez předchozích patchů). Jedno očekávané varování při načtení PLOT2DRHO připomíná chování ASA. `git diff --check` a obrácená aplikovatelnost patchu prošly. SPR-KKR nebyl spuštěn; testy ověřují definici vstupu.
