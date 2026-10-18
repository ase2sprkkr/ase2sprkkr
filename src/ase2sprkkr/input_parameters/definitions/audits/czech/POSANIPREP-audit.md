# POSANIPREP — audit SPR-KKR 9.7 / ASE2SPRKKR

Pořadí 32. Zdroj: PUB9.7_260420_18h13.tar.gz.
ASE2SPRKKR refactor: 267512c4b8ccbbfd37af7ede898febc066b7e3ee.
Samostatný patch 0032-posaniprep.patch proti tomuto commitu.

## Stav a rozsah

POSANIPREP dosud nemá definici. Doplněna definice pro přípravu pozitronových
podkladů POSANI přes kkrgen, pro standardní 3D ASA/FREL systém, jednu energii
GRID=0, NE=1 a elektronový vstupní potential.

Kategorie: **funguje, ale má významné ALE** v omezeném workflow NQ=NT;
**rozbité** v exportu TAUQQ pro obecné NQ≠NT. Task je aktivní. Numerický
SPR-KKR výpočet nebyl spuštěn.

## Zapojení a skutečný význam

`gen.f:414–426` rozpozná POSANIPREP až po prvotní inicializaci energií;
POSANI má přednost před POSANIPREP. Definice nezapisuje oba keywordy.
`gen.f:431–455` vytvoří potřebná pole, volá POSANIPOT a zapne pravidelné
i nepravidelné vlnové funkce i maticové elementy.
`gen.f:1566–1571` pak volá POSANIPREPARE.

Příprava není samostatná minimalizace či hledání dna pozitronového pásu.
`posaniprep.f:86–88` pouze vyžaduje NETAB(1)=1 a konzistenci aktuální energie
s ETAB(1,1). Pro GRID=0, NE=1 nastaví EPATH energii přímo na EMIN
(`epath.f:81–95`). **Uživatel musí vybrat energii pozitronového stavu**, typicky
u dna pásu; -0.2 Ry není nalezený výsledek. Definice při validaci upozorňuje.

Výstup DATASET_pos_dat.pan obsahuje vybrané pozitronové vlnové funkce a
scattering data. Následující POSANI používá elektronový potential a tento
soubor; POSANIPREP naproti tomu interně převádí elektronový potential na
pozitronový (`posanipot.f:124–125` a navazující korekce).

## Parametry

- TASK.POSVC = BM / DR / IPM. `posanipot.f:240–252` čte keyword ze současné
  sekce TASK, vybírá VCEPBM, VCEPDR nebo nulovou korelaci pro IPM. Je vystaven
  pouze v ASA rozsahu. Není to parametr konverze finálního POSANI spektra.
- ENERGY.GRID=0, NE=1, EMIN zvolená absolutní energie v Ry.
- TAU.NKTAB a další společné parametry k-mřížky. POSANIPREPARE znovu čte
  NKTAB (`posaniprep.f:116`) a vytváří TAUQQ přes KMESHS/KLOOPSIJ.
- MODE musí být FREL; NONMAG je povoleno. Nativně je ještě nutné KMROT=0
  (`posaniprep.f:72–77`) a 3D (`posanipot.f:88–89`).
- SCF.FULLPOT je vystaven pouze pro odmítnutí neauditovaného FULLPOT workflow.
  False se nevypisuje a prázdná sekce SCF se nevytváří. Nativní FULLPOT větev
  existuje; není prohlášena za vypnutou. FULLPOT nastavené potentialem bez
  keywordu nelze samotnou validací vstupu vyloučit.

## Defaults a zápis

| Parametr | Definice a zápis | Nativní účinné chování | Posouzení |
|---|---|---|---|
| GRID, NE | zapisuje 0, 1 | POSANIPREP vyžaduje jedinou energii, ale nemá vlastní inicializaci mesh defaults před GEN | explicitní kompatibilní volba |
| EMIN | zapisuje -0.2 Ry, uživatel mění | obecný default -0.2 Ry | hodnota shodná, fyzikálně nemusí být dnem pozitronového pásu |
| ImE | zapisuje 0 Ry | obecný default 0.01 Ry, ale GRID=0 má reálné energie | odlišný explicitní parametr, v tomto mesh nemění energii |
| EMAX | člen není vystaven, nevypisuje se | sentinel se převede na EF | stejné nativní chování; GRID=0/NE=1 používá EMIN, nikoli EMAX |
| POSVC | BM, zapisuje se | BM | převzat nativní default |
| SPLITSS | false, nevypisuje se | false | stejné; zapnutí odmítne |
| MODE | FREL bez keywordu | FREL | stejné |
| SCF.FULLPOT | false, nevypisuje se | skutečný režim může vycházet z potentialu | absence nativní režim nepřepisuje; povinnost použít ASA potential |
| TAU.BZINT | POINTS, zapisuje se | společný výběr závisí na systému | explicitní omezení standardního periodického workflow |
| TAU.NKTAB | společný ASE default 250, zapisuje se | společná nativní inicializace k-mřížky | explicitní volba, ne záruka konvergence |

Přebíráme nativní BM a obecné -0.2 Ry; nepřidáváme smyšlený default energie
dna pásu. EMINEV je odmítnuto, aby energie nebyla vztažena k elektronové EF,
což by pro tento vstup mohlo být zavádějící. Alternativní relativní energie není
prohlášena za nativně zakázanou; jde o omezení auditované definice.

## Doložené chyby a významná ALE

1. **NT místo NQ v exportu TAUQQ.** `posaniprep.f:65–69` alokuje TAUQQ
   jako pole míst NQMAX×NQMAX, ale `:268–272` zapisuje v cyklech 1..NT.
   Následující POSANIREAD používá stejný chybný rozsah. Pro NQ>NT část
   fyzických párů chybí; POSANI je čte z neinicializovaného pole. Pro NT>NQMAX
   může export přistupovat mimo pole. NQ=NT eliminuje tento konkrétní nesoulad,
   nepředstavuje obecný důkaz správnosti. Podmínka závisí na potentialu a není
   automaticky ověřena validací InputParameters.
2. **Omezení na s_1/2 pozitronovou vlnu.** `posaniprep.f:91–106` vybírá
   ZG/ZF(:,1,1), ne celý relativistický stav. Pozdější POSANIREAD navíc potlačí
   imaginární část ZGPOS. Jde o nativní aproximace, které vyžadují fyzikální
   interpretaci; nejsou to chyby defaults definice.
3. **Přítomnost DATASET_positron.rho může potlačit přípravu vlnových funkcí.**
   `gen.f:500–504` načte existující hustotu a nastaví RUNELOOP=false. Přípravná
   rutina přitom dál vyžaduje aktuální vlnové funkce a energii. Samotná existence
   staré hustoty není záruka potřebných podkladů. Tato vazba musí být při
   skutečném běhu kontrolována; definice soubory nemaže a cache automaticky neřeší.
4. **Stejná velikost NKMQ na všech typech/místech.** `posaniprep.f:174–185`
   kontroluje shodné dimenze, jinak STOP. Také požaduje KMROT=0 a 3D. Bez
   potentialu nejsou tyto podmínky všechny rozhodnutelné.
5. **Konverze potentialu závisí na SCFVXC a hustotě.** POSANIPOT má nativní
   větve NONE/VBH/MJW/VWN/PBE, jiné varianty mohou skončit STOP. Definice zde
   nepřepisuje VXC ani nepředstírá audit celé konverze všech funkcionálů;
   konzistence s elektronovým výpočtem zůstává podmínkou použití.

## MPI

Definice má mpi=False. POSANIPREPARE nemá ochranu exportu podle MPI_ID a
otevírá/zapisuje společný DATASET_pos_dat.pan. Keyword POSANIPREP nemusí změnit
interní TASK z NONE, takže běžná závěrečná MPI kontrola GEN nemusí tento postup
ani zastavit. Bez samostatného numerického ověření lifecycle a kolektivních
výpočtů nelze paralelní přípravu považovat za bezpečnou. Nejde o tvrzení, že
nativní kód obsahuje explicitní zákaz MPI.

## Testy

13 nových testů: defaults, souborové rozpoznání tasku, pickle, absolutní energie,
BM/DR/IPM roundtrip, NREL/SREL/SP-SREL, FULLPOT, jiný GRID/NE, relativní energie,
CLUSTER a povolené NONMAG. Kontrolují také upozornění na energii a omezení
3D/KMROT/NQ=NT v textu validátoru.

Celá sada InputParameters: **121 passed**, Python 3.12. Čtyři warnings:
dvě očekávaná upozornění POSANIPREP v obecném testu defaults a dvě stávající
NumPy 2.5 DeprecationWarning. `git diff --check` a zpětná kontrola patche prošly.
Numerický SPR-KKR nebyl spuštěn. Testy nepotvrzují vhodnou energii, konvergenci,
pozitronový potential, exportovaný obsah ani správnost MPI.

## Soubory

- src/ase2sprkkr/input_parameters/definitions/posaniprep.py
- src/ase2sprkkr/input_parameters/test/test_posaniprep.py
- 0032-posaniprep.patch
