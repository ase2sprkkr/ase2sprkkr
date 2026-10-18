# Task 1: WFPLOT

Porovnání veřejných zdrojů SPR-KKR 9.7 z přílohy
`PUB9.7_260420_18h13.tar.gz` a ASE2SPRKKR, větev `refactor`, commit
`8549af97ae106b08749ac1f2dc4b819a72b988e1` z 30. září 2026.
Kontrola proběhla 1. října 2026. Bez spuštění binárního SPR-KKR.

## Stav a dodané soubory

WFPLOT ve výchozím checkoutu nemá definiční modul. ASE2SPRKKR má moduly
ARPES, BSF, DOS, GILBERT, JXC, PHAGEN, SCF, TORQUE a XAS; některé navíc
pokrývají aliasy nebo varianty. Absence modulu neznamená, že nelze ručně
sestavit vlastní InputParametersDefinition.

Nové soubory:

- `src/ase2sprkkr/input_parameters/definitions/wfplot.py`
- `src/ase2sprkkr/input_parameters/test/test_wfplot.py`

Modul se objeví v `InputParameters.definition_modules` automatickým
vyhledáváním přes pkgutil; registr není nutné měnit. Použití:

```python
from ase2sprkkr import InputParameters

ip = InputParameters.create("WFPLOT")
ip.CONTROL.POTFIL = "Fe.pot"
ip.TASK.ITSEL = 1
ip.TASK.L = 2
ip.TASK.MJ = "+3/2"
ip.ENERGY.EMIN = 0.5
ip.save_to_file("Fe-wfplot.in", atoms=False)
```

Pro jádrové funkce:

```python
ip.TASK.STATE = "CORE"
ip.TASK.L.clear()  # L se v CORE nečte
ip.TASK.CL = "2P"
ip.TASK.MJ = "1/2"
```

## Parametry a skutečné defaulty

| Parametr | SPR-KKR při vynechání | Nová definice | Zápis defaultu |
|---|---|---|---|
| ITSEL | Povinný; počáteční 1 není použitelný fallback | 1 | Ano |
| STATE | BAND | BAND | Ne |
| TASK.PRINT | Aktuální IPRINT z CONTROL.PRINT | Dynamicky CONTROL.PRINT | Ne |
| L, BAND | NLQ(IQAT(1,ITSEL))-1 | None znamená ponechat rozhodnutí na SPR-KKR | Ne |
| CL, CORE | 2P | 2P | Ne |
| MJ | +1/2 | 1/2, stejná fyzikální hodnota | Ne |
| GRID | Pro WFPLOT se nenastavuje samostatný použitelný taskový default | 3 | Ano |
| NE | NETAB začíná nulou; vynechané NE může skončit na NEMAX=0 | 1 | Ano |
| EMIN | -0.2 Ry | -0.2 Ry | Ano |
| EMAX | Zpravidla EF, pokud nebylo nastaveno | Stejné jako EMIN | Ano |
| ImE | 0.01 Ry | 0.01 Ry | Ano |
| SITES.NL | 3, nebo 4 pro vybrané těžké prvky; případně další změny | Sdílený ASE default 3 | Ano |

Nevypisované reprezentované defaulty odpovídají čtečkám SPR-KKR.
`L=None` není tvrzení, že numerický SPR-KKR default je None: znamená
neurčenou hodnotu delegovanou na kontextově závislý Fortranový default.
Není správné vložit nevypisované pevné `L=2` pro každý atomový typ.

`GRID=3, NE=1` je vhodnější než obecný ASE default `GRID=5, NE=32`:
WFPLOT v BAND větvi používá pouze `ETAB(1,1)` a ostatní energie neprojde.
Pro GRID=3 a NE=1 je tento bod přesně `EMIN+i*ImE`. Explicitní EMAX=EMIN
nemění tento bod; umožňuje též automaticky napsat oba relativní limity
EMINEV/EMAXEV, které INIT_MOD_ENERGY vyžaduje společně. Horní limity
jsou přístupné přes `ip.ENERGY["EMAX"]` a `["EMAXEV"]`, jako ve sdíleném
ENERGY(emax="emin").

SPR-KKR default SITES.NL je pro těžké prvky lepší než pevné 3. Nový modul
zde zachovává sdílenou, výslovně zapisovanou politiku ASE2SPRKKR; tuto
obecnou politiku nemění. Pro f-stavy je nutné nastavit NL alespoň 4.

## Syntaxe, omezení a validace

- `TASK` obsahuje samostatné klíčové slovo `WFPLOT`, nikoli TASK=WFPLOT.
- `ITSEL` je index atomového typu, nikoli index krystalografického místa.
  Je kladný a vždy se zapisuje. WFPLOT volá
  SECTION_SET_INTEGER('ITSEL',ITSEL,9999,1), tedy bez něj skončí.
- `MJ` se zapisuje jako polocelý zlomek: `1/2`, `+3/2`, `-5/2`.
  TXTMJ2NUM přijímá liché čitatele nejvýše 61; desetinné `0.5` nečte.
- Při abs(MJ)>l+1/2 Fortran hodnotu omezí. Definice na to upozorní,
  hodnotu automaticky nepřepisuje.
- `STATE` přijímá BAND nebo CORE. CL přijímá n=1..6 a s/p/d/f;
  dodatečná fyzikální kontrola l<n odmítá neexistující slupky (např. 1F).
  Tato dodatečná kontrola je přísnější než samotná čtečka SPR-KKR.
- L musí být nezáporné; při jednotném SITES.NL se kontroluje L<NL.
  Při různých NL jednotlivých míst se bez mapování typu na místa
  nepředstírá znalost správného cutoffu.
- TASK.PRINT dědí CONTROL.PRINT. Výslovně jiná hodnota se zapisuje.
- MODE.MODE=NREL/SREL/SP-SREL vede k upozornění: INITVAR pro WFPLOT
  následně vynutí IREL=3.
- GEN odmítá WFPLOT pro DMFT a LDA+U; definice tyto volby odmítá při
  konečné validaci.
- FULLPOT vede přímo ve WFPLOT k STOP_MESSAGE. Definice to dokumentuje;
  sama neotevírá zadaný POTFIL a neověřuje jeho FULLPOT stav.
- MPI musí být vypnuto. GEN pro WFPLOT v MPI větvi skončí před vyvoláním
  WFPLOT. Nová definice má mpi=False.
- NE>1 vyvolá upozornění, protože BAND použije jen první bod.
- STATE=CORE s nastaveným L a STATE=BAND s nastaveným CL vyvolají
  upozornění na ignorovaný parametr.

## Zdrojové důkazy

- `initvar.f:149–218`: rozpoznání TASK WFPLOT, vypnutí běžné energetické
  smyčky a TAU souborů.
- `initvar.f:782`: vynucení IREL=3.
- `wfplot.f:67–70`: FULLPOT není implementováno.
- `wfplot.f:98–114`: ITSEL, lokální PRINT a STATE.
- `wfplot.f:128–130`: čtení CL a MJ pro CORE.
- `wfplot.f:292–300`: pouze ETAB(1,1), default L a MJ; pokud je energie
  příliš blízko nule, náhrada za (0.5,0.01) Ry.
- `scan_input.f:203–280`: required sentinel 9999 a STOP_KEY=1.
- `scan_input.f:1916–1983`: TXTMJ2NUM, textový formát MJ.
- `scan_input.f:2138–2201`: CL default 2P a limity n a l.
- `scan_input.f:2203–2280`: default MJ=0.5 a omezení na l+0.5.
- `gen.f:210–222`: odmítnutí DMFT/LDA+U pro WFPLOT.
- `gen.f:1357–1370`: MPI guard; `gen.f:1452–1453`: volání WFPLOT.
- `mod_energy.f:47–54`, `init_mod_energy.f:73–76, 133–180, 310–391`:
  inicializace energií, limity a EPATH.
- `epath.f:132–144`: první bod GRID=3.
- `kkrmain.f:410–430`: atomově závislé defaulty NL.

## Ověření a otevřené body

Celý adresář `input_parameters/test`: 128 úspěšných testů. Z toho 20
nových případů WFPLOT: discovery, zápis a zpětné načtení přes from_file,
copy/pickle, defaulty, frakční syntaxe MJ, kladné indexy, CORE/BAND,
omezení L, relativní energie, varování a odmítnutí DMFT/LDA+U.

Nebyl spuštěn výpočet SPR-KKR. Není tak ověřen číselný výsledek ani
produkované .agr soubory. Použije se obecný DefaultOutputReader;
specializovaný parser výsledných vlnových funkcí není součástí změny.
Neověřuje se existence zvolené jádrové slupky u konkrétního prvku,
horní mez ITSEL podle skutečného potenciálu, FULLPOT ve skutečném
POTFIL ani cutoff pro jednotlivé typy při nejednotném NL.

Další jednotlivý task podle INITVAR: MEPLOT. Inventář všech tasků není
tvrzením, že je každý ve veřejné verzi skutečně použitelný; například
PHONONS má ve veřejné verzi zvláštní stop. Je nutné každý projít zvlášť.
