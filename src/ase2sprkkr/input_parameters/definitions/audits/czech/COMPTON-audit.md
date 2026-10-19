# COMPTON — audit SPR-KKR 9.7 / ASE2SPRKKR

Pořadí 30. SPR-KKR: PUB9.7_260420_18h13.tar.gz.
ASE2SPRKKR: refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee.
Patch je samostatný proti tomuto commitu, nezahrnuje předchozí tasky.

## Stav a rozsah

COMPTON dosud nemá definici InputParameters. Přidána `definitions/compton.py`
pro periodický FREL výpočet s dvojicí kontur GRID=5 a SPLITSS, přes kkrgen.
Nativní implementace je aktivní. Předběžná kategorie: **funguje, ale má významné ALE**;
konkrétní větev EMD2D s nenulovým podélným posunem obsahuje doloženou chybu.
Numerické výsledky nebyly ověřeny během tohoto auditu.

NONMAG je povoleno, skalárně relativistické a nerelativistické režimy nikoli.
Základní výstupy jsou DATASET_BND.mcp, DATASET_COR.mcp a DATASET_TOT.mcp,
s celkovou a spinovou rozdílovou složkou (`compton.f:1315–1380`). MD2D vytváří
projekci, EMD2D řez hybnostní hustotou. Nejde o zaměnitelné volby.

## Zdrojové zapojení a vstupy

- `initvar.f:480–489`: COMPTON přepne RDTAUDIR na true a běžnou energetickou
  smyčku na false. `gen.f:1572–1577` volá COMPTON.
- `compton.f:167–170`: vyžaduje IREL=3 a kompatibilní konstantu NSYMMAX.
- `kkrgen.f:264–274`: implicitní TAUDIR je DATASET.tau.dir; lze zadat TAUDIR
  a TAUPATH v sekci TAU. Oba parametry jsou doplněny do definice.
- `init_mod_energy.f:284–302`: kontroluje počet energií a kompatibilitu TAU
  adresáře; některé nesoulady zapínají přepočet, jiný počet energií skončí STOP.
  **Definice nezajišťuje vytvoření nebo obsah cache.** Výpočet potřebuje vhodné
  TAU podklady z předchozího kroku. Python testy jejich kompatibilitu neověřují.
- `compton.f:203–228`: USE_ZJOVZZ, SPLITSS, MEMODE a NPMAX.
- `init_mod_compton_pmesh.f`: PNVEC, PPVEC1, PPVEC2, PNMAX, PPMAX,
  NPN, NPP, MD2D, EMD2D, SHFT, SHFT_PP, SHFT_PN, INSIDE_SPHERE, USE_PEQUIV.
- `compton.f:263`: SQQ_RESOLVED. Doplněno jako volba rozlišení podle párů míst.

Definice záměrně pokrývá standardní split-contour postup. ONLY_CORE a USE_ZJOVZZ
mají jiné větve, zde nejsou vystaveny ani prohlášeny za vypnuté. GRID=5/5 a POINTS
jsou omezení auditovaného rozsahu, nikoli tvrzení, že všechny ostatní nativní
energetické mřížky jsou zakázané. Počty NE musí být kladné.

## Defaults: zapisované versus vynechané

| Parametr | Definice a zápis | Nativní chování | Posouzení |
|---|---|---|---|
| GRID, NE | zapisuje 5/5, 30/50 | stejné v INIT_MOD_ENERGY pro COMPTON | převzaty nativní defaults |
| ImE, EMIN | zapisuje 0.01, -0.2 Ry | stejné | nativní defaults |
| EMAX | nevystaven a nevypisuje se | nastaví se na EF pro COMPTON | stejné účinné chování |
| ENERGY.SPLITSS | true, vždy zapíše keyword | inicializace pro COMPTON nastaví true | zachováno |
| TASK.SPLITSS | true, vždy zapíše keyword | jeho absence v COMPTON přepíše SPLITSS na false | výslovná volba zachovávající dvoukonturový postup |
| MEMODE | 1, zapisuje se | 1 | stejné |
| NPMAX | 250, zapisuje se | 250; hodnoty <=0 vrací na 250, >10000 ořízne | stejné default, validace 2..10000 zabrání NP=1 |
| PNMAX, PPMAX | 0, 10 a.u., zapisují se | stejné | PNMAX=0 vede k jedinému podélnému bodu, validátor upozorní |
| NPN, NPP | 50, 50, zapisují se | stejné; pro PNMAX=0 se NPN efektivně změní na 1 | parametr a efektivní velikost se liší |
| USE_PEQUIV | 0, zapisuje se | 0; mimo 0..2 vrátí na 0 | default stejný, neplatné hodnoty odmítnuty |
| směrové vektory | nezadané, nevypisují se | normála z, příčné osy x/y | stejné |
| MD2D, EMD2D, shift flags, INSIDE_SPHERE, SQQ_RESOLVED | false, nevypisují se | false při absenci keywordů | stejné |
| TAUDIR, TAUPATH | nezadané, nevypisují se | DATASET.tau.dir, aktuální adresář | stejné |

PNMAX=0 není rozumný default pro celý profil, ale může sloužit pro jediný bod.
Nativní default proto ponechán s upozorněním; uživatel má výslovně zadat rozsah.
Nativní NPMAX=250 je vhodný start, potřebuje kontrolu konvergence; MEMODE=2
umožňuje porovnat výpočet bez interpolace. Default není zárukou konvergence.

## Omyly a omezení nativního kódu

### Dvě nezávislá čtení SPLITSS

INIT_MOD_ENERGY nastaví SPLITSS pro COMPTON, ale `compton.f:210` jej znovu
načte **z TASK**, kde absence keywordu znamená false. ENERGY.SPLITSS samo
nestačí. Definice vždy zapisuje oba keywordy a validátor požaduje oba true.
Není to automatická synchronizace dvou editovatelných hodnot: vypnutí kterékoliv
odmítne při validaci. USE_ZJOVZZ by navíc SPLITSS vynutilo false (`:212`),
proto je tato alternativa mimo rozsah definice.

### Směrové vektory se neortogonalizují

`init_mod_compton_pmesh.f:49–121` normalizuje zadané vektory a vytvoří třetí
vektor křížovým součinem. Neodstraní podélnou složku druhého zadaného vektoru.
Nulové či rovnoběžné vektory mohou dělit nulou; nekolmé vytvoří deformovaný základ.
Definice kontroluje konečné nenulové vektory délky 3 a kolmost páru.
Při PNVEC se PPVEC2 ignoruje; definice tuto zavádějící kombinaci odmítá.
Bez PNVEC musí být zadány oba příčné vektory nebo žádný.

### Sudé posunuté mřížky ztrácejí krajní body

`init_mod_compton_pmesh.f:206–245` používá ±INT((N-1)/2), ale krok původní
mřížky je MAX/(N-1). Pro N=50 se místo 50 bodů získá 49 a rozsah nedosáhne
obou zamýšlených krajů. Pro N=51 je konstrukce konzistentní.
Definice požaduje lichý počet v posouvaném směru a nejvýše jeden shift flag.
To je konzervativní ochrana rozsahu, nikoli tvrzení, že nativní sudý vstup vždy havaruje.

### EMD2D: podélné škálování proběhne dvakrát

`init_mod_compton_pmesh.f:183–185` již násobí PNVEC hodnotou DPN.
V EMD2D pak `:273–279` nastaví DPN=PNMAX a tentýž vektor znovu škáluje.
Řez se tak obecně nevytvoří na deklarovaném PNMAX. Například pro NPN=50,
PNMAX=2 bez shift flagů je první škála 2/49 a výsledná 4/49 namísto 2.
Definice dovoluje EMD2D jen pro PNMAX=0, kde tento konkrétní problém neovlivní
polohu řezu. Ostatní EMD2D vstupy odmítá. Nativní kód nebyl měněn.

### Další významná ALE

- `NPMAX=1` vede k dělení NP-1 při interpolaci (`compton.f:446`,
  `compton_mod1_me0.f:144`). Definice vyžaduje nejméně 2.
- `compton.f:185–194` při GF_CONV_RH přepne konvenci na ZJ a SOLVER_FP na BS.
  I FULLPOT tedy mění konvenci; bezprostřední kompatibilitu s cache je třeba kontrolovat.
- S INSIDE_SPHERE je PFEMAX=min(PNMAX,PPMAX); při nativním PNMAX=0 se ořízne
  na 1e-10. Takto zadaná koule neposkytuje obvyklý profil.
- Konkrétní systémy vyžadují kompatibilní velikosti NKMQ; nativní kontroly
  `compton.f:678` mohou skončit STOP. Definice bez potentialu toto nezaručí.

## MPI

`compton.f:581–592` vyžaduje NE na každé aktivní kontuře alespoň NPROCS.
Energetická smyčka používá MPI_DISTRIBUTE a výsledky redukuje MPI_SUM_R
(`:1273–1296`). `mpi_distribute.f` rozděluje energie pozpátku, takže poslední
bod zpracuje rank 0; odpovídá podmínce výpočtu core contribution (`compton.f:626`).
Proto je mpi=True zdrojově opodstatněné. Nebyl proveden skutečný MPI běh.

## Testy

22 nových testů: zápis obou SPLITSS, defaults, automatické rozpoznání souboru,
pickle, TAU cesty a směrové vektory, NONMAG, zakázané režimy, vypnutí každého
SPLITSS, chybné směry, validní příčný pár, shift/parita, EMD2D s nulovým a
nenulovým posunem, počty/interpolace, NE a auditované kontury/BZINT.

Celá sada InputParameters: **130 passed**, Python 3.12. Čtyři warnings:
dvě očekávaná upozornění PNMAX=0 při obecném testu defaults a dvě stávající
NumPy 2.5 DeprecationWarning. `git diff --check` a zpětná kontrola patche prošly.
Numerický SPR-KKR nebyl spuštěn; testy neověřují fyzikální hodnoty, TAU cache,
konvergenci, MPI výsledky ani správnost alternativ ONLY_CORE/USE_ZJOVZZ.

## Soubory

- src/ase2sprkkr/input_parameters/definitions/compton.py
- src/ase2sprkkr/input_parameters/test/test_compton.py
- 0030-compton.patch
