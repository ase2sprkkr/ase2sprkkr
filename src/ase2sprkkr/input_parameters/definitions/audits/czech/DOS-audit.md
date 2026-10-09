# DOS — audit SPR-KKR 9.7 / ASE2SPRKKR

Pořadí 34. Zdroj: PUB9.7_260420_18h13.tar.gz.
ASE2SPRKKR refactor: 267512c4b8ccbbfd37af7ede898febc066b7e3ee.
Samostatný patch 0034-dos.patch proti tomuto commitu.

## Stav a opravy definice

DOS již podporován. Standardní workflow má **plný smysl**; numerickou správnost
konkrétního systému tento audit nepotvrzuje. DOSREP a další doplňkové výstupy
mají vlastní omezení a nejsou tímto patchem vystaveny.

Původní definice měla CONTROL, TAU, ENERGY, TASK a SITES, ale chyběly jí MODE,
CPA a STRCONST. Doplněny společné definice těchto sekcí; lze řádně zadat NREL,
SREL, SP-SREL a společná nastavení CPA/strukturních konstant. FREL zůstává
nativním režimem bez keywordu. Původní zapisované energetické defaults zachovány.

Přidány kontroly kladného NE, alespoň dvou bodů pro GRID=2/3, konečných
energetických mezí a konečného nezáporného ImE. ImE=0 na GRID=3 vyvolá upozornění,
nikoli zákaz. Kontroly se neuplatňují při jednotlivém nastavování souvisejících
parametrů, ale při finální validaci/zápisu.

## Nativní zapojení a MPI

INITVAR nevybírá samostatnou interní hodnotu TASK='DOS' podle keywordu DOS;
běžný kkrgen postup zůstává TASK='NONE' (`mod_calcmode.f:95`, `initvar.f`).
DOS zde označuje standardní výpočet Greenovy funkce a výpis hustoty stavů,
ne samostatnou koncovou rutinu. CONTROL.ADSI=DOS není příkaz, který by sám
provedl přepnutí interního tasku.

Energetická smyčka GEN volá CALCDOS_WRITE (`gen.f:1149–1153`).
`calcdos_write.f:48–66` redukuje DOBS_LTEX přes MPI_SUM_R a nechá výstup
na ranku 0. `mpi=True` je pro standardní DOS opodstatněné. Závěrečná kontrola
GEN povoluje TASK=NONE; nelze chybně odvodit, že DOS skončí v MPI jen proto,
že v jejím seznamu není řetězec DOS. Skutečný MPI běh nebyl proveden.

## Defaults a jejich zápis

| Parametr | ASE2SPRKKR | Nativní chování bez parametru | Posouzení |
|---|---|---|---|
| GRID | 3, zapisuje se | není samostatný časný DOS override v INIT_MOD_ENERGY | explicitní vhodná volba pro energeticky rozlišený DOS |
| NE | 300, zapisuje se | není samostatný časný DOS override; NETAB má obecný počáteční stav 0 | explicitní volba rozlišení, nutno ověřit konvergenci |
| ImE | 0.01 Ry, zapisuje se | obecný default 0.01 Ry | stejné |
| EMIN | -0.2 Ry, zapisuje se | -0.2 Ry | stejné |
| EMAX | 1.0 Ry, zapisuje se | sentinel se převede na elektronovou EF (`init_mod_energy.f:182`) | odlišný zapisovaný default, záměrně umožňuje i neobsazené stavy |
| EMINEV/EMAXEV | nezadané, nevypisují se | žádný relativní interval | stejné; při zadání se předávají obě relativní meze |
| SPLITSS | false, nevypisuje se | false | stejné |
| MODE | bez keywordu, FREL | FREL | stejné |
| TAU.NKTAB | společný ASE default 250, zapisuje se | společná inicializace k-mřížky | explicitní volba, nikoli záruka konvergence |
| doplněné CPA/STRCONST | běžné defaults společných definic | společné nativní čtení | žádný DOS-specific override doplněn |

EMAX=EF by byl vhodnější při úloze omezené na obsazené stavy, ale jako výchozí
DOS scan by vynechal neobsazenou část spektra. Zapisovaný EMAX=1.0 proto ponechán.
Je to absolutní energie vůči referenci potentialu, **nikoli 1 Ry nad EF**.
Pro interval kolem EF slouží dvojice EMINEV/EMAXEV. ImE=0.01 Ry odpovídá
širokému rozmazání; jeho vhodnost je třeba ověřit konvergencí, ne zaměnit za
numericky zanedbatelnou veličinu.

Vynechané defaults nevedou k jinému nativnímu energetickému intervalu:
EMAX se zde na rozdíl od některých jiných tasků skutečně zapisuje.

## Doložené omyly a omezení

### GRID=3 a NE=1

`epath.f:136` nastaví ESTEP jen pro NE>1. Následně jej ale používá pro váhy
WETAB i při NE=1 (`:143–145`). Jde o neinicializovaný integrační krok.
Definice nyní tuto kombinaci odmítá. GRID=2 dělí NE-1 (`:115`), proto také
vyžaduje alespoň dva body. To je nutná kontrola, nikoli úplná záruka vhodné
geometrie obdélníkové kontury pro všechny počty bodů a rozsahy energií.
GRID=0, NE=1 je platný samostatný bod a zůstává povolený.

### DOSREP není běžný l-resolved DOS

`kkrgen.f:1169–1196` čte v CONTROL pouze při FREL DOSREP=RLM/CLM/REL,
DOSCORSYS=LOC/GLO, WRKAPDOS a WRPOLAR. Původní definice je nevystavovala;
tento patch je nepřidává, protože jejich správnost vyžaduje samostatný audit.

WRDOSREP ukládá energicky rozlišené pole a vypisuje je při IECURR=NETAB(1)
(`wrdosrep.f:152–185`). Nemá MPI redukci a část pole DOS se před použitím
kompletně nenuluje. Při rozdělení energií mezi ranky má poslední rank jen část
energií. Běžná redukce CALCDOS_WRITE **nezahrnuje toto lokální pole**.
Doplňkový DOSREP výstup tedy nelze v energetickém MPI považovat za správně
sestavený jen proto, že standardní DOS funguje. Jde o doložený problém této
vedlejší větve, nikoli standardního DATASET.dos.

Komentář WRDOSREP o kontrolách kubického paramagnetického systému nepředstavuje
vynucené omezení: CHECK je compile-time false (`wrdosrep.f:32–35`).

### Meze auditu

- MODE.FREL a doplněné skalární/nerelativistické vstupy jsou podporované formátem;
  jejich konkrétní výsledek závisí na potentialu a není zde numericky porovnán.
- FULLPOT může být určeno potentialem; patch nepřidává celou SCF konfiguraci.
- Společné ENERGY poskytuje další kontury a SPLITSS. Tento audit ověřuje
  standardní GRID=3 scan a uvedené nutné kontroly, ne všechny kombinace
  pokročilých kontur, jejich počty bodů a numerické váhy.
- ImE=0 je formálně dovoleno, ale výpočet přímo na reálné ose vyžaduje kontrolu
  pólů, hustoty k-mřížky a konvergence.
- Nejsou zde auditovány tepelné fluktuace, NLCPA, impurity/embedded cluster,
  self-energy ani jiné pokročilé větve. Nejsou označeny za vypnuté.

## Testy

12 nových testů: původní defaults, automatické načtení DOS souboru, pickle,
NREL/SREL/SP-SREL roundtrip, GRID=2/3 s NE=1, platný GRID=0/NE=1,
nekladné NE, upozornění ImE=0, záporné ImE a relativní interval vůči EF.

Celá sada InputParameters: **120 passed**, Python 3.12. Dvě stávající
NumPy 2.5 DeprecationWarning. `git diff --check` a zpětná kontrola patche prošly.
Numerický SPR-KKR nebyl spuštěn. Testy neověřují hodnoty DOS, konvergenci,
MPI výsledky ani doplňkové DOSREP/WRKAPDOS/WRPOLAR výstupy.

## Soubory

- src/ase2sprkkr/input_parameters/definitions/dos.py
- src/ase2sprkkr/input_parameters/test/test_dos.py
- 0034-dos.patch
