# COMPTON — audit současné definice

Definice `compton.py` podporuje periodický FREL výpočet Comptonových hybnostních hustot pomocí dvou kontur GRID=5 a SPLITSS. Aktivní nativní větev má dokumentovaná omezení.

## Stav a rozsah

Současná definice pokrývá periodický FREL výpočet přes `kkrgen`; zvláštní nativní větve nejsou tímto kontraktem numericky ověřeny.

## Výstupy a názvy souborů podle SPR-KKR

Nativní `compton.f:1314–1354` zapisuje `DATSET_BND.mcp`, `DATSET_COR.mcp` a `DATSET_TOT.mcp`; při `PRINT>0` také `DATSET_data.mcp`. `MD2D` vybírá příponu `_mcp2d.dat` a `EMD2D` příponu `_emd2d.dat`, vždy pro `BND`, `COR`, `TOT` (`compton.f:1439–1458`). Jde o valenční, core a celkovou složku; `MD2D` je projekce, `EMD2D` dvourozměrný řez hybnostní hustotou.

`kkrmain.f:239–245` nejprve nastaví `DATSET0=DATSET`, potom k `DATSET` připojí neprázdné `ADSI`. Výstupy COMPTON proto zahrnují ADSI. Současný `compton.py::output_files()` používá stejný prefix `{dataset}_{adsi}`; při prázdném ADSI formatter vynechá podtržítko. Manifest proto odpovídá nativnímu pojmenování, nejde o změnu Fortranu.

`init_mod_compton_pmesh.f:183–185,269–279` nejprve škáluje `PNVEC` krokem `DPN` a při `EMD2D` znovu hodnotou `PNMAX`. Zdroj tím potvrzuje omezení nenulového podélného řezu uvedené níže.

## Zdrojové zapojení a vstupy

- `initvar.f:480–489`: COMPTON přepne RDTAUDIR na true a běžnou energetickou
  smyčku na false. `gen.f:1572–1577` volá COMPTON.
- `compton.f:167–170`: vyžaduje IREL=3 a kompatibilní konstantu NSYMMAX.
- `kkrmain.f:264–274`: implicitní TAUDIR je DATSET.tau.dir; lze zadat TAUDIR
  a TAUPATH v sekci TAU. Oba parametry jsou doplněny do definice.
- `init_mod_energy.f:284–302`: kontroluje počet energií a kompatibilitu TAU
  adresáře; některé nesoulady zapínají přepočet, jiný počet energií skončí STOP.
  **Definice nezaručuje vytvoření ani obsah cache.** Chybějící TAU může podle `kkrmain.f:932–960` vyvolat přepočet; jeho úspěšnost ani fyzikální shodu existující cache Python testy neověřují.
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
| TAUDIR, TAUPATH | nezadané, nevypisují se | DATSET.tau.dir (DATASET_ADSI při neprázdném ADSI), aktuální adresář | stejné |

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
