# Task 2: MEPLOT

Prověřeno 2. října 2026 proti veřejným zdrojům SPR-KKR 9.7 z přílohy
PUB9.7_260420_18h13.tar.gz a aktuální větvi ASE2SPRKKR refactor:
267512c4b8ccbbfd37af7ede898febc066b7e3ee. MEPLOT v této větvi chybí.
Nový modul se registruje automaticky přes existující pkgutil discovery.
Patch je samostatný a nevyžaduje předchozí WFPLOT patch.

## Co task dělá

ME_PLOT počítá a vykresluje dipólové maticové elementy. Výstupní křivky
obsahují absolutní hodnoty normalizovaných maticových elementů, nikoli
celé absorpční spektrum. Běží přes kkrgen, mpi=False.

- INITIALSTATE=BAND, NEPHOT=1: elektronový scan; ENERGY určuje počáteční
  energii a finální energie je počáteční energie + EPHOT/RY_EV.
- INITIALSTATE=BAND, NEPHOT>1: fotonový scan při jediné počáteční energii.
- INITIALSTATE=CORE: ENERGY určuje finální elektronovou energii, počáteční
  jádrový stav vybírá CL. Fotonové vstupní parametry se v této větvi nečtou.

Počáteční m_j=+1/2, polarizace IPOL=2 (q=0, z) a forma maticových elementů
MEFORM='ADA' jsou ve zdroji pevné. Parametry MJ, MJI, IPOL či MEFORM proto
nevymýšlíme a nevystavujeme jako podporované vstupy.

## Defaulty a zápis

| Parametr | Efektivní SPR-KKR default | Nová definice | Zapisuje se default? |
|---|---|---|---|
| ITSEL | 1 | 1 | Ne |
| IPRINT | Aktuální CONTROL.PRINT | Dynamicky CONTROL.PRINT | Ne |
| INITIALSTATE | BAND | BAND | Ne |
| L | Všechny počáteční band orbitaly; LISEL=-1 | None = bez výběru l | Ne |
| CL | 2P | 2P | Ne |
| TASK.NONMAG | False | False | Ne, flag se zapisuje jen pro True |
| SEFIN | 0 | 0 | Ne |
| EPHOT | 100 eV | 100 eV | Ne |
| NEPHOT | 1 | 1 | Ne, ani po explicitním nastavení na 1 |
| EPHOTMIN/EPHOTMAX | Nenastavené lokální proměnné | None | Ne |
| ENERGY.EMIN | -0.2 Ry | -0.2 Ry | Ano |
| ENERGY.EMAX | EF při chybějícím explicitním limitu | None = ponechat skutečné EF | Ne |
| ENERGY.ImE | 0.01 Ry | 0.01 Ry | Ano |
| ENERGY.GRID | Nemá samostatný použitelný MEPLOT default | 3 | Ano |
| ENERGY.NE | NETAB inicializované na 0; bez explicitního NE může výpočet skončit | 300 | Ano |
| SITES.NL | Kontextově 3 nebo 4 | Sdílená současná definice: None, delegace do SPR-KKR | Ne |

NE=300 a GRID=3 jsou vědomě zapisované ASE volby pro hladkou energetickou
křivku. Neprohlašujeme je za SPR-KKR default. GRID=3 je pro plot vhodnější
než obecný oblouk GRID=5: vede podél reálné osy ve vzdálenosti ImE.
Při fotonovém scanu musí uživatel nastavit ENERGY.NE=1.

Default EMAX nebyl nahrazen pevnou hodnotou: skutečné EF z potenciálu je
zde vhodnější. None v API znamená neurčený limit, nikoli numerické EF.
Totéž platí pro L=None, které znamená ponechat původní výběr všech orbitalů.
Současný refactor již vynechává SITES.NL, a zachovává tak vhodnější
atomově závislý SPR-KKR default; dřívější pevné ASE NL=3 zde neplatí.

## Podstatné rozdíly proti WFPLOT

- ITSEL není povinné klíčové slovo: STOP_KEY=0 zachová počáteční hodnotu 1.
- Klíč pro lokální tisk je IPRINT, nikoli PRINT.
- Výběr stavu je INITIALSTATE, nikoli STATE.
- Chybějící L neznamená nejvyšší orbital. Čtečka sice nejprve použije
  NLQ-1, ale při FOUND_INTEGER=False následně vrátí LISEL=-1 a vybírá vše.
- ME_PLOT skutečně prochází NETAB(1) energií, ne jen ETAB(1,1).
- FULLPOT má implementované FPSSITE větve. Není zde stop jako ve WFPLOT;
  numerická správnost těchto větví nebyla ověřena.
- IREL=3 se pro MEPLOT nevynucuje. Výstavba kanálů předpokládá relativistické
  kappa/j/m_j indexování. Pro jiné MODE vydáváme upozornění, nikoli tvrzení,
  že je Fortran explicitně odmítá. Jejich bezpečnost a fyzikální význam
  vyžadují samostatné ověření.

## Chyby a omezení původních zdrojů

1. ME_PLOT rozlišuje režim fotonového scanu podle samotné přítomnosti NEPHOT
   či některé meze, nikoli pouze podle NEPHOT>1. Explicitní NEPHOT=1 bez
   mezí proto vstoupí do range větve a může tisknout neinicializované
   EPHOTMIN/EPHOTMAX; vynechá také čtení EPHOT. Nový modul default 1 nikdy
   nezapisuje. Tím chrání i kombinaci NEPHOT=1, EPHOT=150.
2. Kontrola neúplných mezí ve Fortranu zachytí jen část neplatných kombinací.
   Definice proto vyžaduje oba limity pro NEPHOT>1 a odmítá meze v režimu
   NEPHOT=1. Pro jedinou energii se používá EPHOT.
3. Při NEPHOT>1 se YY indexuje přes IEPHOT, nikoli přes dvojici IE/IEPHOT.
   S více elektronovými energiemi se tak stejné fotonové křivky opakovaně
   přepisují; nakonec zbude poslední elektronová energie. Definice tuto
   kombinaci odmítá a vyžaduje ENERGY.NE=1. Jde o záměrně přísnější ochranu
   než ve Fortranu, který kombinaci nezastavuje.
4. SEFIN se při NEPHOT>1 automaticky vypíná. Definice na to upozorňuje;
   hodnotu automaticky nepřepisuje.
5. Pro CORE se L a fotonové parametry nečtou. Výslovné nastavení vyvolá
   upozornění. SEFIN se však může použít i v této větvi.
6. SEFIN dostává ERYD, nikoli ERYDF. U BAND je ERYD počáteční energie,
   zatímco finální stav má energii ERYD+EPHOT/RY_EV. Vzhledem k popisu
   SEFIN jako self-energy finálního stavu je to podezřelé; bez numerického
   ověření to neoznačujeme za dokázanou chybu a zdroj neměníme.
7. GEN má explicitní zákaz DMFT/LDA+U pro některé tasky včetně WFPLOT,
   ale MEPLOT v tomto seznamu není. Definice takový zákaz nevymýšlí;
   správné použití těchto metod s MEPLOT zůstává neověřené.

## Validace

ITSEL>=1, L>=0, NEPHOT>=1, SEFIN=0..3, kladné počty elektronových energií.
INITIALSTATE je BAND/CORE. CL má n=1..6, s/p/d/f; fyzikální kontrola l<n
je přísnější než samotná Fortranová čtečka. Při výslovném jednotném NL
ověřujeme l<NL. Bez potenciálu neověřujeme horní mez ITSEL, skutečnou
existenci jádrové slupky ani mapování nejednotného NL na atomové typy.
Neúplné editační stavy jsou povolené při set; konečná validace nastává
při save/parse. Použita je existující validace páru EMINEV/EMAXEV.

## Příklady

```python
from ase2sprkkr import InputParameters

ip = InputParameters.create("MEPLOT")
ip.CONTROL.POTFIL = "Fe.pot"
ip.TASK.EPHOT = 150.0  # eV
ip.ENERGY.set({"EMINEV": -10.0, "EMAXEV": 0.0})
ip.save_to_file("Fe-meplot.in", atoms=False)
```

Fotonový scan:

```python
ip.ENERGY.set({"NE": 1, "EMINEV": -2.0, "EMAXEV": -2.0})
ip.TASK.EPHOT.clear()
ip.TASK.set({"NEPHOT": 20, "EPHOTMIN": 10.0, "EPHOTMAX": 200.0})
```

Pro CORE se nastaví INITIALSTATE=CORE, CL a finální energetický interval;
předchozí fotonové parametry mají být vyčištěny, protože se ignorují.

## Zdrojové důkazy a ověření

- me_plot.f:132–146: pevné MJI/IPOL a inicializace defaultů.
- me_plot.f:151–211: skutečné názvy parametrů, L fallback, fotonová čtečka,
  ignorování SEFIN při více fotonových energiích.
- me_plot.f:297–306, 349–356, 564–571: dimenze a přepisování fotonových křivek.
- me_plot.f:465–476: energie finálního stavu a argument ERYD pro SEFIN.
- me_plot.f:379–387, 499–507: FPSSITE větve pro FULLPOT.
- sefin.f:25–28, 107–108: význam a limity SEFIN.
- initvar.f:159–160, 208–218: rozpoznání MEPLOT a vypnutí běžné ELOOP.
- gen.f:210–222: seznam zákazů DMFT/LDA+U, bez MEPLOT.
- gen.f:1357–1370, 1480–1483: MPI guard a vlastní dispatching.

Nových 20 testovacích případů ověřuje zápis, discovery/from_file, copy/pickle,
nevypisované defaulty, IPRINT, explicitní NEPHOT=1, kompletní a neúplné
fotonové scany, přepisování křivek, SEFIN, CORE, cutoffy, relativní energie
a rozdíl TASK.NONMAG versus CONTROL.NONMAG.

Celý adresář input_parameters/test: 128 úspěšných testů, včetně 20 nových.
Patch prošel kontrolou git diff --check a git apply --reverse --check.

SPR-KKR nebyl spuštěn; .agr data nemají nový specializovaný parser a použije
se existující DefaultOutputReader. Testy neprokazují numerickou správnost
maticových elementů. Další task v pořadí INITVAR: MECHECK.
