# POSANITAU — audit SPR-KKR 9.7 / ASE2SPRKKR

Pořadí 33. SPR-KKR: PUB9.7_260420_18h13.tar.gz.
ASE2SPRKKR refactor: 267512c4b8ccbbfd37af7ede898febc066b7e3ee.
Samostatný patch 0033-posanitau.patch proti tomuto commitu.

## Stav a význam

Chybějící definice doplněna pro standardní periodický 3D ASA/FREL výpočet
pozitronové anihilační rychlosti a životnosti přes kkrgen. Kategorie pro společný
audit: **funguje, ale má významné ALE**. Nativní task je aktivní, numerická
správnost není ověřena tímto auditem.

POSANITAU není totéž jako POSANI ani POSANIPREP. Nečte DATASET_pos_dat.pan,
nevytváří jeho TAUQQ export a nepřebírá proto chybu NT místo NQ z tohoto exportu.
Pracuje s elektronovým potentialem/hustotou, sestaví pozitronový potential,
integruje pozitronovou hustotu a vyhodnotí její překryv s elektronovou hustotou.

## Zapojení, vstupy a výstupy

`gen.f:414–426` rozpoznává keywordy až po počáteční inicializaci energií;
POSANI a POSANIPREP mají prioritu před POSANITAU. Definice píše pouze POSANITAU.
`gen.f:431–449` volá POSANIPOT a zakazuje aktualizaci EF.
`gen.f:458–476` vyžaduje GRID=5 a uloží kopii elektronové hustoty před její
pozdější náhradou pozitronovou hustotou.

Definice vystavuje:
- TASK.POSVC: BM / DR / IPM, volba ASA electron-positron correlation potential.
- ENERGY: jednu GRID=5 konturu, NE=30, ImE=0.01 Ry, absolutní EMIN/EMAX.
- TAU: standardní POINTS a společné parametry k-mřížky.
- MODE: FREL; NONMAG je povoleno.
- SCF.FULLPOT: výslovné true se odmítá, protože tento audit pokrývá ASA.
  Nativní FULLPOT větev existuje a není označena za vypnutou. Režim odvozený
  přímo z potentialu není samotnou validací InputParameters ověřitelný.

POSANILAMBDA nemá vlastní keywordy pro volbu zesilovacího faktoru, angularní
mesh nebo normalizaci. Definice takové parametry nevymýšlí.

`gen.f:1508–1564` dokončí korekce hustoty, zapíše DATASET_positron.rho,
vynutí grafy PLOTPRS a pro FULLPOT také PLOT2DPR, potom volá POSANILAMBDA.
Rychlosti a životnosti se tisknou na standardní výstup: rate [1/s], lifetime [s]
a celková lifetime také [ps]. Pole RATET mají rozlišení podle typů.

## Defaults a zápis

| Parametr | Definice a zápis | Nativní účinné chování | Posouzení |
|---|---|---|---|
| GRID | zapisuje 5 | GEN vyžaduje 5, ale keyword POSANITAU nepřidává vlastní časnou inicializaci | explicitní nutná volba |
| NE | zapisuje 30 | není samostatný POSANITAU default před GEN | explicitní startovní integrace, konvergenci ověřit |
| ImE | zapisuje 0.01 Ry | obecný default 0.01 Ry | stejné |
| EMIN | zapisuje -0.2 Ry | obecný default -0.2 Ry | stejné číslo, vhodnost pro pozitronový pás nezaručena |
| EMAX | nezadáno, nevypisuje se; lze výslovně změnit | sentinel převeden na EF v INIT_MOD_ENERGY:182 | stejné účinné chování; EF z elektronového vstupu nemusí být vhodná mez pozitronového pásu |
| POSVC | BM, zapisuje se | BM | převzat nativní default |
| SPLITSS | false, nevypisuje se | false | stejné; zapnutí je mimo definovaný workflow |
| MODE | bez keywordu, FREL | FREL | stejné |
| SCF.FULLPOT | false, nevypisuje se | může vycházet z potentialu | input default nativní režim nepřepisuje; použití ASA je podmínkou definice |
| BZINT / NKTAB | POINTS / společný ASE default 250, zapisují se | společná nativní inicializace | výslovná volba standardního workflow, nikoli záruka konvergence |

Vynechaný EMAX neskrývá jiný nativní default. Uživatel jej může nastavit explicitně;
na fyzikální vhodnost integračních mezí upozorňuje validátor. Relativní EMINEV
nebo EMAXEV definice odmítá jako neauditovaný výběr vztažený k elektronové EF;
nejde o obecný nativní zákaz relativních energetických parametrů.

## Významná ALE a doložené omyly

### POSVC nemění zesilovací faktor životnosti

POSVC řídí korelační potential v POSANIPOT (`posanipot.f:240–252`).
Anihilační rychlost však vždy používá POSANI_GFACTOR
(`posanilambda.f:92`, pro FULLPOT `:250`). Tato rutina má jediný pevný polynom
v rs (`posani_gfactor.f:38–39`) a vůbec nečte POSVC.
**POSVC=IPM proto neznamená výpočet životnosti se zesilovacím faktorem 1.**
Mění potential, zatímco lifetime enhancement zůstává pevně dané.
To je významné omezení fyzikální interpretace, ne chyba mapování keywordu.

### Normalizace může maskovat nevhodnou integrační oblast

ASA větev POSANILAMBDA nejprve integruje pozitronový náboj s NAT a CONC
(`posanilambda.f:60–71`) a normalizuje hustotu na celkový náboj 1 (`:74–78`).
Nulový/malý náboj vede k STOP, ale nenulový neověřuje správnou oblast pásu.
Normalizace sama nepotvrzuje vhodnost EMIN/EMAX ani konvergenci NE.

### Cache hustoty potlačuje nový energetický výpočet

`gen.f:500–504` čte DATASET_positron.rho a při úspěchu vypne RUNELOOP.
Změna EMIN, EMAX, NE či POSVC tak nemusí znamenat novou hustotu: lze znovu
vyhodnotit starý soubor. POSANIRHOIO neověřuje fyzikální původ cache ani shodu
s aktuálními parametry; u FULLPOT také přeskakuje hlavičky bez ověření jejich obsahu.
Definice tento soubor nemaže a neřeší automaticky jeho invalidaci.

### Inverze nulové rychlosti jednotlivého typu

Na konci POSANILAMBDA se testuje pouze celková RATE>=1e-8. Potom se pro každý
aktivní typ tiskne 1/RATET(IT), bez samostatné kontroly RATET(IT).
Celková kladná rychlost nezabraňuje nulové rychlosti jednoho typu; individuální
výstup pak může obsahovat nekonečno nebo vyvolat floating-point výjimku.
Celková lifetime má kontrolu RATE, problém je v typovém výpisu.

### Rozsah konverze a potential

POSANIPOT vyžaduje 3D a obsahuje větve závislé na SCFVXC; nepodporovaný
funkcionál může skončit STOP. Tato definice nepřepisuje elektronové VXC a
audit neověřuje všechny jeho varianty. Použití kompatibilního elektronového
potentialu a hustoty zůstává podmínkou.
FULLPOT výpočet používá pevnou úhlovou kvadraturu NGAU=20×20; celý jeho
výpočet tvarových funkcí a konvergence nejsou v tomto tasku auditovány.

## MPI

Definice mpi=False: auditovaná cesta je sériová. GEN provádí závěrečnou korekci
a zápis nové hustoty pod podmínkou MPI_ID=0, ale následné ploty a POSANILAMBDA
nejsou v tomto bloku chráněny stejnou podmínkou. Keyword POSANITAU navíc může
ponechat interní TASK=NONE, který projde obecnou závěrečnou MPI kontrolou.
Bez ověření redukcí hustoty a následného lifecycle nelze slibovat bezpečný
paralelní výpočet. Nejde o explicitní nativní zákaz MPI.

## Testy

15 nových testů: defaults a vynechaný EMAX, automatické rozpoznání tasku,
pickle, BM/DR/IPM, explicitní EMAX, NREL/SREL/SP-SREL, FULLPOT, jiný GRID,
více kontur NE, relativní energetické meze, CLUSTER a povolené NONMAG.
Kontrolují upozornění na integrační oblast a význam POSVC.

Celá sada InputParameters: **123 passed**, Python 3.12. Čtyři warnings:
dvě očekávaná upozornění POSANITAU v obecném testu defaults a dvě stávající
NumPy 2.5 DeprecationWarning. Kontrola diff a zpětné aplikace patche prošla.
Numerický SPR-KKR nebyl spuštěn. Testy neověřují lifetime, hustotu, konvergenci,
cache, elektron-pozitronové modely ani numerické MPI výsledky.

## Soubory

- src/ase2sprkkr/input_parameters/definitions/posanitau.py
- src/ase2sprkkr/input_parameters/test/test_posanitau.py
- 0033-posanitau.patch
