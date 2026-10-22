# PSHIFT — audit a nová definice InputParameters

Podklad: přiložené SPR-KKR 9.7 (PUB9.7_260420_18h13) a aktuální ASE2SPRKKR `refactor`, ověřený fetch, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`. PSHIFT v této revizi nemá definiční modul. Samostatný patch přidává `definitions/pshift.py` a `test/test_pshift.py`; nepotřebuje předchozí patche.

## Vstup

`initvar.f:173–174` rozpoznává TASK PSHIFT, `gen.f:1426–1431` volá rutinu PSHIFT v KKRGEN. `pshift.f` nečte žádné task-specific parametry; není zde ITSEL, L, PRINT ani přepínač výstupů. Rutina počítá všechny atomové typy, používá CONTROL.PRINT a cutoff z potenciálu/SITES.NL. Definice přidává společné CONTROL, TAU, MODE, STRCONST, ENERGY, TASK a SITES.

Používá `kkrgen`, `mpi=False`: smyčky PSHIFT neobsahují distribuci práce mezi ranky a výstupní názvy jsou sdílené. MPI ani numerický běh ověřen nebyl.

| Volba | Chování / default SPR-KKR | Definice |
|---|---|---|
| TASK | pouze keyword PSHIFT | žádné vymyšlené vlastní parametry |
| MODE.MODE | vynechání znamená plnou relativitu | FREL/None, NREL, SREL a SP-SREL; všechny mají zdrojové větve |
| ENERGY.GRID | bere reálnou část první ETAB cesty | explicitní 1, přímá ekvidistantní reálná cesta |
| ENERGY.NE | NETAB(1), globální počáteční hodnota 0 bez task-specific defaultu | explicitní 300, kladný počet |
| ENERGY.ImE | reset na 0 | default 0; upozornění při nenulové volbě |
| ENERGY.EMIN | -0.2 Ry (`init_mod_energy.f:49,74`) | explicitní -0.2 Ry, shodný default |
| ENERGY.EMAX | při vynechání skutečná EF | explicitní 1.0 Ry, viz níže |
| ENERGY.EMINEV/EMAXEV | meze relativně k EF v eV | společné párové rozhraní ENERGY zachováno |
| ENERGY.SPLITSS | reset na false | false, při zapnutí upozornění |
| SITES.NL | cutoff z potenciálu | None, nevypisuje se; kontext zůstává SPR-KKR |

ENERGY.ImE a SPLITSS přepisuje `init_mod_energy.f:315–320`. PSHIFT pracuje jen s první cestou a reálnou částí energie (`pshift.f:66,79,123,368`). Na druhou cestu definice upozorňuje. Jiná GRID není nativně zakázaná, ale může vytvářet opakované nebo obráceně uspořádané reálné energie; to je nevhodné pro postupné rozbalování fází. Definice proto doporučuje GRID=1 a jiné cesty pouze varuje.

### Rozdíly defaults

Vypsaná horní mez **EMAX=1.0 Ry** se záměrně liší od nativního vynechání, které používá EF. PSHIFT hledá rezonance pod absolutní energií 1 Ry; takové okno je pro diagnostiku užitečnější než zastavení na EF, které může vynechat nezaplněné rezonance. Proto zde nativní default nepovažuji automaticky za lepší. Jde o explicitní volbu ASE2SPRKKR, nikoli o nesprávně modelovaný omitted default. Pro jiné okno uživatel nastaví absolutní nebo relativní meze.

NE=300 a GRID=1 jsou rovněž explicitně zapsané praktické defaults. Nativní NETAB začíná na nule (`mod_energy.f:43–44`), init_mod_energy nemá zvláštní PSHIFT default počtu bodů; tento stav není vhodné přebírat. Vynechané SITES.NL zůstává kontextové, a zachovává proto nativní chování.

## Režimy a výstupy

- FREL: přímý DIRBS výpočet phase shifts a logaritmických derivací; výstupní panely zahrnují nemagnetické kanály a vybrané magnetické kanály.
- NREL/SREL: jeden spinový kanál, phase shifts a reálná/imaginární část t-matice.
- SP-SREL: dva spinové kanály a navíc odhad výměnného rezonančního štěpení.
- Výstupní energetická osa je **absolutní eV**, `Re(ETAB)*RY_EV`, nikoli `(E-EF)*RY_EV` (`pshift.f:79`). Ani zadání EMINEV/EMAXEV tento význam výstupní osy nemění.
- Relativistická větev volá DIRBS přímo bez SOCTL (`:396–400`); obecné MODE.SOC její SOC neškáluje. Definice při explicitní volbě v FREL upozorňuje. Nezaměňuje tuto větev za SOCPAR, který SOCTL aktivně přepisuje.
- Výstupní grafy zůstávají bez specializovaného parseru; patch přidává InputParameters a používá obecný output reader.

FULLPOT explicitně odmítá `check_input.f:80–82` i `gen.f:1428–1429`. Rozpoznává se z potenciálu, a samotný input bez načtení potenciálu jej nemůže spolehlivě validovat; omezení je dokumentované a vynucuje je SPR-KKR. DMFT a LDAU odmítá GEN (`gen.f:211–222`); definice také. Není zaveden neexistující zákaz NREL/SREL/SP-SREL.

## Omyly a rizika ve zdroji

1. **Rezonanční štěpení může použít chybějící rezonanci.** V nerelativistické větvi je RESONANCE společný příznak pro oba spiny (`:213–232`). Stačí najít jednu rezonanci; druhá ERES může zůstat -9999.90 Ry a přesto se odečíst. Stejný problém je v relativistické větvi (`:519–545`): společný příznak pro různé m_j a DESO navíc používá ERES předchozího kappa kanálu bez kontroly jeho nalezení. Číselné štěpení pak nemusí být fyzikální výsledek.
2. **Detekce rezonance nekontroluje skutečné překročení pi/2.** Testuje jen aktuální fázi >pi/2 (`:221,532`) a interpoluje mezi dvěma body. Pokud už první bod leží nad hranicí nebo jsou sousední fáze shodné, může jít o extrapolaci nebo dělení nulou. Je nutné požadovat vhodně ohraničený přechod a nenulový rozdíl fází.
3. **Nerelativistický výstup při neuniformních cutoffech může číst neinicializovaná data.** Plnění PSTAB/RETSS/IMTSS používá lokální `NLQ(IQAT(1,IT))` (`:131`), zápis ale používá globální NL (`:196,275–280`). Pole se před tím nenuluje. Definice upozorňuje na explicitně neuniformní SITES.NL v nerelativistickém režimu; obsah potenciálu při omitted NL zde kontrolován není.
4. **Nerelativistický logdrv je nedostupný.** Je za `IF (NL.LT.0)` (`:293`), označen TODO, a v této větvi se LDTAB ani nevypočítává. Není vhodné inzerovat logaritmické derivace jako dostupné ve všech režimech.
5. **Relativistické legendy mají jednociferný formát.** `WRITE(C2J,'(I1)')` (`:470,572`) nedokáže zapsat 2*j>=10. Od NL=6 vzniká j=11/2 a legenda obsahuje hvězdičku. Definice upozorňuje na NL>5, výpočet neodmítá kvůli kosmetické chybě.

Definice zdrojové numerické chyby neopravuje a netvrdí, že validní input zaručuje správné rezonanční štěpení.

## Testy

15 nových případů: nalezení tasku, vypsané a omitted defaults, souborový/textový round trip, pickle/copy, všechny čtyři relativistické režimy, relativní meze potlačující absolutní defaults, odmítnutí korelací a nekladného NE, upozornění na ImE/SPLITSS/druhou cestu, ignorované SOC, neuniformní cutoff a legendy.

Celá sada `src/ase2sprkkr/input_parameters/test`: **123 passed**, na výchozí revizi plus tento samostatný patch (bez předchozích patchů). `git diff --check` a obrácená aplikovatelnost patchu prošly. SPR-KKR nebyl spuštěn: jde o testy definice vstupu a audit Fortranu, nikoli numerickou validaci phase shifts.
