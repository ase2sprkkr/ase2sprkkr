# SOCPAR — audit a definice InputParameters

Podklad: přiložené SPR-KKR 9.7, PUB9.7_260420_18h13; aktuální ASE2SPRKKR `refactor` ověřený fetch, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`. SOCPAR zde nemá definici. Samostatný patch přidává `definitions/socpar.py` a `test/test_socpar.py`; nevyžaduje předchozí tasky.

## Vstup a skutečné chování

`initvar.f:170–171` rozpoznává TASK SOCPAR; `gen.f:1419–1424` volá stejnojmennou rutinu v KKRGEN. `socpar.f` neobsahuje žádné SECTION_SET/SECTION_FIND: vlastní parametry TASK se tedy nečtou. Není zde výběr ITSEL, L, PRINT ani uživatelský přepínač solveru. Zpracují se všechny atomové typy a orbitaly l=1..NL-1 (`socpar.f:76, 135–136`). PRINT se bere ze společného CONTROL.PRINT.

Definice obsahuje CONTROL, TAU, MODE, STRCONST, ENERGY, TASK a SITES. `executable=kkrgen`, `mpi=False`: vlastní smyčky rutiny nejsou rozdělené mezi MPI ranky a soubory mají společná jména. Numerické ani MPI spuštění zde ověřeno nebylo.

| Volba | Chování SPR-KKR | Definice ASE2SPRKKR |
|---|---|---|
| ENERGY.GRID | SOCPAR bere reálnou část ETAB první cesty | explicitní GRID=1: přímá ekvidistantní reálná cesta |
| ENERGY.NE | používá NETAB(1) | explicitní NE=300; kladný počet |
| ENERGY.ImE | init_mod_energy přepisuje na 0 | default 0; upozornění při nenulové volbě |
| ENERGY.EMIN | obecný nativní EMINDEFAULT | explicitní -0.2 Ry, společná volba ASE2SPRKKR |
| ENERGY.EMAX | chybějící mez se odvozuje od EFERMI | None, nevypisuje se; zachován nativní default |
| ENERGY.EMINEV/EMAXEV | meze relativně k EF v eV | zachováno společné párové rozhraní ENERGY |
| ENERGY.SPLITSS | přepisuje na false | default false; při zapnutí upozornění |
| SITES.NL | kontextový cutoff z potenciálu | None, nevypisuje se; explicitní globální maximum kontrolováno 2..51 |
| MODE.SOC | rutina přepisuje SOCTL pro všechny typy/orbitaly na 1e-3 | při explicitním zadání upozornění |
| solver | rutina nastaví BS-SOC | žádný vymyšlený parametr TASK |

`init_mod_energy.f:315–320` vypíná SPLITSS a nastavuje EIMAG/EILOW=0. `socpar.f:111` navíc převádí ETAB na reálnou energii. Druhá cesta se nepoužije. Definice na tyto volby upozorňuje; složitější contour neprohlašuje za nativně zakázaný, ale upozorní na projekci jeho energií na reálnou osu.

GRID a NE mají být zadány: modul energie začíná s NETAB={0,0} (`mod_energy.f:43–44`) a v init_mod_energy není samostatný SOCPAR default počtu bodů. Nemá proto smysl přebírat nulový nativní počet místo použitelného explicitního defaultu. GRID=1/NE=300 je volba pro graf, NE lze snížit pro jednoduchou kontrolu. Vypsané EMIN=-0.2 Ry může být upraveno podle požadovaného okna. Nevypisované EMAX a NL delegují skutečný kontext SPR-KKR.

## Omezení

FULLPOT je explicitně odmítnut jak `check_input.f:75–77`, tak `gen.f:1421–1422`. Vstupní MODE nemá FULLPOT přepínač a plný potenciál se rozpoznává z potenciálového souboru; samotná validace definice jej bez načteného potenciálu nemůže odhalit. Zákaz je proto dokumentován a vynucuje jej SPR-KKR. Definice sama není důkazem kompatibility libovolného potenciálu.

DMFT a LDAU odmítá GEN (`gen.f:211–222`); definice je odmítá také. SOCPAR používá relativistické wave-function kanály a solver BS-SOC. Na NREL/SREL/SP-SREL upozorňuje jako na neověřené kombinace: check_input u SOCPAR neobsahuje explicitní požadavek IREL=3, proto nebyl zaveden domnělý nativní zákaz.

NL=1 neprodukuje žádnou křivku: orbitalová smyčka začíná na IL=2. LEG má 50 prvků a počet křivek je NL-1, tedy globální NL>51 překračuje jeho kapacitu (`socpar.f:28, 46, 184–189`). Definice kontroluje explicitní maximum SITES.NL. U vynechaného cutoffu nebo detailů konkrétního potenciálu tuto kontrolu vstup sám neprovede.

## Nalezené omyly ve SPR-KKR

1. **Nedosažitelná korekce znaménka magnetického pole.** Rutina nejprve obrátí BT a následně testuje `ABS(DREAL(NORM)).LT.0` (`socpar.f:81–96`). Absolutní hodnota nemůže být záporná, takže další obrácení znaménka nikdy nenastane. Pravděpodobný překlep v zamýšleném testu znaménka; přesnou správnou opravu je nutné potvrdit fyzikální konvencí. Může ovlivnit znaménko výměnného štěpení.
2. **Nesoulad jednotek výměnného výstupu.** SOCME obsahuje RY_EV a má jednotky eV (`:159`), ale EXCME je `2*DREAL(XI)` bez převodu (`:163–165`). Do grafu jde s faktorem 1.0 (`:227`) a s osou označenou eV (`:216–220`). Při obvyklých interních Ry jednotkách BT jde o nesoulad faktoru RY_EV. Silný závěr ze zdrojového výpočtu, nikoli porovnání numerického běhu s referencí. Definice na něj upozorňuje a žádná výstupní data neopravuje.
3. **Uživatelské MODE.SOC neřídí tento výpočet.** SOCTL=1e-3, BT=0 a SOLVER=BS-SOC jsou záměrně pevně nastavené přímo v rutině (`:93–107`). Nejde o chybějící parametr ASE2SPRKKR; je třeba zabránit dojmu, že obecná volba SOC zde zůstává účinná.

Výstupní soubory se vytvářejí přes XMGRHEAD pro `soc` a `exc`, pro jednotlivé typy. Definice používá obecný output reader; specializovaný parser grafů součástí tohoto patchu není.

## Ověření

17 nových testů: nalezení tasku, explicitní a omitted defaults, souborový a textový round trip, pickle/copy, relativní energetické meze, neplatné NE/cutoffy, korelace, upozornění na relativistický režim, přepis SOC/ImE a ignorovanou druhou cestu/SPLITSS.

Celá sada `src/ase2sprkkr/input_parameters/test`: **125 passed** na uvedené výchozí revizi plus tento samostatný patch. Počet nezahrnuje předchozí samostatné patche. `git diff --check` i kontrola obrácené aplikace patchu prošly. SPR-KKR nebyl spuštěn; testy ověřují input definici, nikoli fyzikální správnost výpočtu.
