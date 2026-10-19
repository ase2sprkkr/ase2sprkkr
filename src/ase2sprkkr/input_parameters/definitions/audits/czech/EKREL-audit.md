# EKREL — audit a nová definice InputParameters

Podklad: přiložené SPR-KKR 9.7 (PUB9.7_260420_18h13) a aktuální ASE2SPRKKR refactor, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`, ověřený fetch. EKREL nemá definici. Patch samostatně přidává `definitions/ekrel.py` a `test/test_ekrel.py`.

## Parametry TASK a defaults

EKREL počítá E(k) pro uspořádané systémy, případně hledá Fermiho povrch. INITVAR rozpoznává keyword EKREL a vypíná běžný RUNELOOP (`initvar.f:459–467`); GEN pak volá vlastní EKREL. Rutina má vlastní energetickou a k smyčku.

| Parametr | Nativní default / interpretace | Definice |
|---|---|---|
| NK | 1000, požadavek počtu k bodů | >=2, default se nevypisuje |
| EMIN | 0.01 absolutních Ry | default se nevypisuje |
| EMAX | 1.0 absolutních Ry | default se nevypisuje |
| NE | 501 hledacích energií | >=2, default se nevypisuje |
| F-SURF | false | flag |
| ESLICE | EF, pouze při F-SURF | None a nevypsání deleguje EF; explicitní hodnota v Ry |
| KPATH | 1 | sdílené popisy cest BSF/KDIRTAB; default se nevypisuje |
| NKDIR | 0, použít KDIRTAB | 0..9; nenulové znamená manuální cestu |
| KA1..KA9 / KE1..KE9 | při manuální cestě oba konce povinné | tři složky v kartézských jednotkách 2*pi/a, vyplňují se explicitně |

Zdroj: `ekrel.f:92–154`. EMIN, EMAX a NE se čtou z **TASK**, nikoli ENERGY. Meze jsou absolutní Ry. EMINEV/EMAXEV v obecné ENERGY sekci nejsou náhradou těchto TASK mezí.

Ruční cesta se zadává například TASK.NKDIR=1, TASK.KA1=[0.,0.,0.], TASK.KE1=[1.,0.,0.]. Konkrétní očíslované options jsou vytvořené krátkou deklarativní smyčkou pro nativní rozsah 1..9. Validace kontroluje úplnost, konečnost, nenulové segmenty a nepoužité indexy. Zde nedochází ke skrytému zkrácení NKDIR>9: definice jej odmítá, zatímco nativní rutina jej ořízne na 9.

KPATH používá stejnou tabulku kdirs jako BSF, proto definice přebírá existující slovník popisů. Dostupnost konkrétního čísla závisí na Bravaisově mřížce; samotná validace inputu bez struktury ji neověřuje. Rozsah podporovaných mřížek není rozšířen jen přidáním této definice.

## Fermiho povrch

F-SURF vytváří celobrilouinovský k grid pro BXSF výstup. TASK.NK je výchozí požadavek hustoty, nikoli záruka přesného výsledného počtu bodů: zdroj jej přepočítává, včetně zdvojených hraničních bodů (`ekrel.f:270–305`). Ruční k-path zde neřídí výsledné vzorkování. KDIRTAB se ale volá před větví F-SURF, a proto i pro Fermiho povrch mohou platit jeho omezení mřížky.

ESLICE při vynechání zůstává EF, což může ležet mimo výchozí hledací interval 0.01..1 Ry. Bez potenciálu nelze tento vztah ověřit. Explicitní ESLICE validátor kontroluje proti TASK mezím. Zdroj ve skutečnosti připouští rovnost s mezí, přestože jeho chybová hláška tvrdí striktní nerovnost (`ekrel.f:132–134`); definice odpovídá podmínce `EMIN <= ESLICE <= EMAX`.

## Společné sekce a režimy

Definice má CONTROL, TAU, MODE, STRCONST, ENERGY, TASK a SITES; executable=kkrgen. `mpi=True` odpovídá skutečné distribuci vlastní k smyčky pomocí MPI_SET_WORKERS_ARRAY a IPROCK (`ekrel.f:435–448`). MPI ani numerické spuštění zde ověřeno nebylo.

ENERGY má explicitní GRID=1, NE=1, ImE=0, EMIN=-0.2 Ry; EMAX=None se nevypisuje a deleguje EF. ENERGY slouží společné inicializaci a nepředstavuje hledací interval EKREL. INIT_MOD_ENERGY potřebuje nenulový počáteční počet bodů před pozdějšími task-specific úpravami. Na vícenásobné ENERGY.NE definice upozorňuje; počet vlastních hledacích energií řídí TASK.NE. SITES.NL=None se nevypisuje a zachovává kontext potenciálu.

Nevypisované TASK defaults odpovídají nativní rutině. Pozitivní EMIN=0.01 Ry je nativní volba, nikoli volba vynucená ASE2SPRKKR; kvůli pólům a chování u E=0 ji ponechávám. Pro konkrétní energetické okno je nutné obě TASK meze nastavit. Nevypisované ESLICE musí zůstat kontextové, nelze jej nahradit konstantním číselným EF.

**CPA slitiny jsou nativně odmítnuté** (`ekrel.f:81–82`); zdroj doporučuje Bloch spectral function. Vstupní definice bez potenciálu nemůže stanovit NCPA, proto omezení dokumentuje a spoléhá na nativní kontrolu. Nevytváří zbytečnou CPA sekci. FULLPOT nemá v check_input obecný zákaz pro EKREL; stejně tak GEN nezakazuje DMFT/LDAU u tohoto tasku. EKREL má vlastní čtení DMFT self-energy (`:384–399`). Definice proto nevymýšlí jejich zákaz, ale numerická kompatibilita zde nebyla ověřena.

## Omyly a rizika ve zdroji

1. **NE>=2 se kontroluje jen pro F-SURF, přestože dělení NE-1 je obecné.** `ETOL=(EMAX-EMIN)/(NEK-1)` (`ekrel.f:179`) se provádí také pro obyčejné E(k). Definice proto vyžaduje TASK.NE>=2 vždy.
2. **Zdroj obecně nekontroluje platné energetické okno.** Při EMIN=EMAX vzniká ETOL=0 a následná kontrola průchodu E=0 může dělit nulou; obrácené meze navíc porušují předpoklad pro hledání. Definice vyžaduje konečné EMIN<EMAX.
3. **Chybová hláška ESLICE neodpovídá testu.** Text uvádí EMIN<ESLICE<EMAX, ale podmínka odmítá jen hodnoty striktně mimo interval. Rovnost není nativně zakázaná a ani definice ji nezakazuje.
4. **NK není přesná velikost výsledné cesty.** Segmenty se dělí podle délek a každý dostává minimálně dva body (`:204–227`); F-SURF má vlastní grid. Definice popisuje NK jako požadavek, nikoli garantovaný počet.
5. **Ruční nulové segmenty nebo celková nulová délka nejsou řádně ošetřené.** Dělení podle SUMD ve zdroji může dělit nulou. Definice kontroluje konkrétní manuální segmenty; dostupnost a geometrii preset cest ověřuje SPR-KKR podle mřížky.

Přechod přes E=0 zdroj zkouší řešit snížením NEK o 1 (`ekrel.f:181–187`). Skutečný počet energií tedy nemusí odpovídat vypsanému NE. Definice toto chování neopravuje.

## Testy

9 nových případů: discovery a MPI metadata, omitted TASK defaults, textový/souborový round trip, pickle, manuální očíslovaná cesta, chybějící nebo nulový segment, počty NE/NK/NKDIR, neplatné energetické okno, F-SURF a ESLICE.

Celá sada `src/ase2sprkkr/input_parameters/test`: **117 passed** na výchozí revizi plus tento samostatný patch (bez předchozích patchů). `git diff --check` a obrácená aplikovatelnost patchu prošly. SPR-KKR nebyl spuštěn; testy ověřují definici vstupu, nikoli numerické pásy nebo Fermiho povrch. Specializovaný output parser není součástí patchu.
