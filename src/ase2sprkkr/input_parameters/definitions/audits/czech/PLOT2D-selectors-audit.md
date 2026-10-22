# PLOT2DPOT / PLOT2DRHO — SPR-KKR 9.7

Selektivní vstupní tokeny, nikoli nové nezávislé solvery. Doplněny společně se sdílenou definicí _plot2d_selectors.py. Kategorie pro budoucí společný audit: fungují, ale mají velká ALE; selective FULLPOT chování a ASA density jsou rozbité/neimplementované.

## Dispatch
initvar.f:201–205 nastavuje PLOT2DPR(1) podle PLOT2DPOT a flag2 podle PLOT2DRHO; oba mapují interní TASK=PLOT2DPR. GEN (gen.f:1466–1478) vybírá FPPLOT2DPR pro FULLPOT a PLOT2DPR_ASA pro ASA.

## Zásadní skutečné chování
ASA routine plot2dpr_asa.f neobsahuje čtení density arrays ani podmíněné rozlišení PLOT2DPR flag1/2. Vytváří potential PLY a geometrické line výstupy (pevná jména PLY-V-ASA-... a POT-ASA-LINES.ply). PLOT2DRHO proto na ASA není density plot. Definice odmítá explicitní FULLPOT=false pro PLOT2DRHO; skutečný potenciál musí být FULLPOT, flag jej nekonvertuje.

FULLPOT routine fpplot2dpr.f:70 nejprve bezpodmínečně volá PLOT2DPR_FP_PLY. Ten prochází potential i density arrays a nevětví výstupy podle jednotlivých PLOT2DPR flagů; na ř.926 volá STOP_REGULAR s textem run of test version completed. Následující běžný plane-grid výpočet je v této cestě nedosažitelný. Native selectory tedy neoddělují PLY potential/density produkty tak, jak jejich názvy slibují. Nejde jen o input validation problém.

Další STOP v fpplot2dpr.f:671 po zápisu FULLPOT-pot.ply je až v následné nedosažitelné rutině; nesmí být zaměňován za první ukončení. Oprava pouze tohoto STOP by nestačila. V této následné rutině se navíc V<-1 ořezává na -1 před zápisem PLY; nepovažujeme ji za aktuálně dosažitelný výstup.

## Geometrie a defaults
TASK.AVEC/BVEC/DVEC jsou tříprvkové vektory v a.u.; oba readers defaults A=(2,0,0), B=(0,2,0), D=0. B se při nenulovém A.B ortogonalizuje a zachovává svou původní délku. Nulové/paralelní vektory způsobují dělení nulou, definice je odmítá, stejně jako nefinite geometrie.

NA=50 explicitní default odpovídá readeru. NB=50 je explicitně zapsáno; při nepřítomném NB native dopočítává nint(NA*|B|/|A|), což pro velmi malý B může dát0. Explicitní50 je pro předvídatelnost bezpečnější, ale mění aspect-ratio sampling při různě dlouhých vektorech. NA/NB>=1.

PLOT2DPR_FP_PLY má vlastní řadu pevných sampling parametrů: NA/NB nejsou záruka běžného rectangular grid výstupu. Podpora input parametrů neznamená, že se v test writeru každý použije stejně jako v pozdějším plane-grid writeru.

GRID0/NE1/ImE0 jsou explicitní startup hodnoty; plot task nevykonává energetickou smyčku. Nepřítomné FULLPOT pro PLOT2DPOT nechává native rozpoznat formát potenciálu; PLOT2DRHO zapisuje true. MPI vypnuto, protože GEN guard plot tasky nepovoluje. DMFT/LDAU plot větve native rovněž odmítá; nové definice tyto MODE parametry nenabízejí.

## Testy a dodání
7 nových případů: oba roundtripy, mesh metadata, nulové/paralelní/nefinite vektory, NA0 a odmítnutí ASA density. Celá input_parameters sada115 passed, dvě původní NumPy warnings. Numerický solver nebyl spuštěn, výstupy nejsou numericky certifikované.

Patch0056 nezávislý proti baseline267512c4, zahrnuje obě definice, společný helper a testy. Instalovat všechny3 definitions soubory společně. Jejich info výslovně popisuje native test-writer omezení.
