# XES — nová definice a audit

Dodaný SPR-KKR 9.7 porovnán s fetchovanou větví ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. XES nemá vlastní definici. Patch 0024-xes.patch přidává xes.py a test_xes.py samostatně proti tomuto základu.

## Rozhraní

INITVAR rozpoznává TASK.XES a nastavuje rentgenový workflow IXRAY=1, ICORE=0, RDTAUDIR=true a INITELOOP=true. GEN:1378–1385 volá XRAYSPEC, která interně přejmenuje XES na RXES. CHECK_INPUT:43–47 vyžaduje IREL=3. Definice proto odmítá ostatní MODE.

XRAYSPEC:376–380 vyžaduje ENERGY.NE a jeho rovnost počtu bodů z prvního TAU contouru. Samotná definice nečte obsah TAU directory: její mesh, počet bodů, typy a potenciál musí odpovídat skutečnému výpočtu. Při dostupné TAU se čte, jinak hlavní workflow může připravit smyčku; numerické propojení zde nebylo spuštěno.

TASK vystavuje IT, CL, ME, OUTPUT, POLOBS, DIPOL, GFCONVRH, WRREGIRR, NONDIPE2, NONDIPM1, TAUCORE a Eulerovy úhly FRAMETET/FRAMEPHI/FRAMEGAM. Nejde o rozhraní ARPES ani VBPES. Frame úhly otáčejí společný souřadný rámec, nikoli prostý jediný vektor elektrického pole; init_mrot_frame.f:48 tvoří matici GETMROT(-gamma,-theta,-phi).

## Defaulty a zápis

| Volba | Nativní efekt | Definice |
| --- | --- | --- |
| IT | 1 | 1, vynecháno |
| CL | 2P ze společné čtečky | 2P, vynecháno |
| ME | ADA | ADA, vynecháno |
| OUTPUT | BARN; běžný výstup v Mb | BARN, vynecháno |
| POLOBS | S | S, vynecháno |
| DIPOL | 0; pro NAB se přepíše na1 | 0, vynecháno |
| GFCONVRH | false při nepřítomnosti keywordu | false, vynecháno |
| WRREGIRR | DATA uvádí true, ale SECTION_FIND_KEYWORD je přepíše podle přítomnosti | false, vynecháno |
| TAUCORE | dvě nuly, hodnoty<1e-6 přebírají tabulku WCOREHOLE | None, vynecháno |
| Eulerovy úhly | 0/0/0 | None, vynecháno; nepřítomnost nezapíná CHANGE_FRAME |
| GRID/NE | XES nemá vlastní inicializační default; NE=0 vede k neplatné kapacitě | explicitně3/300 |
| ImE | 0,01 Ry | explicitně0,01 Ry |
| EMIN | -0,2 Ry | explicitně -0,2 Ry |
| EMAX | XES výslovně přepíše na EF | None, vynecháno |

GRID3/300 je explicitní volba energetického scanu obsazených stavů, nikoli nevypsaný údaj předstírající shodu s nativním defaultem. Konvergence EMIN/NE/ImE není ověřena. Nativní NE=0 není vhodnější; ostatní základní hodnoty jsou zachované. Vynechané EMAX, efektivní WRREGIRR=false a Eulerovy úhly odpovídají skutečnému nativnímu chování.

Definice omezuje n<=5, protože SHELL obsahuje jen K,L,M,N,O, a l<n. IT musí být kladné; jeho horní mez a existence core slupky vyžadují potenciál. NE>=2 chrání GRID3 před neinicializovaným ESTEP při jednom bodu. Split a dvě cesty se odmítají kvůli použití prvního TAU contouru. ME=NAB s DIPOL=0 vyvolá upozornění na přepsání; NONDIPE2 má přednost před NONDIPM1.

## Velká ALE a nativní omyly

1. **MPI spektrum vůbec nespustí.** GEN:1355–1368 ve výčtu pokračujících tasků nemá XES; po energetické smyčce STOP_REGULAR ukončí běh před XRAYSPEC. Definice má mpi=False. Obdobné riziko se týká existující XAS definice s mpi=True; úplný XAS audit bude samostatný krok.
2. **RXES sdílí vlastní rate calculation s RXAS.** Jediný výskyt RXES mimo komentář je při přejmenování tasku. Výpočtová a výstupní větev rozlišuje XMO/XRS proti všem ostatním, takže XES vstoupí do společné absorption-rate větve. Výsledek proto nelze bez fyzikálního porovnání s referencí prohlásit za ověřené emisní spektrum. Nemá prostý zástupný STOP jako FORCETEN, takže task nebyl přeskočen, ale bude mít výrazné ALE ve společném auditu.
3. **Výstup je stále označen XAS.** Název .rat i text SPEC používají stejnou core shell konvenci jako XAS; text uvádí XAS spectrum. XAS a XES pro stejný DATASET, typ a slupku mohou zapsat tentýž soubor.
4. **MESPLIT není přepínač:** po načtení keywordu se K_ME_SPLIT vždy nastaví true (297–298). Definice jej nevystavuje jako účinnou volbu. DYNAMICAL se čte pouze pro RXAS/XMO, proto rovněž není nabídnuto pro XES.
5. **TAUCORE se čte, ale jeho použití pro RXES rate není doloženo.** Hodnoty se převádějí z eV do Ry, popřípadě přebírají WCOREHOLE. Následné použití je v rozšířených XMO/XRS větvích; standardní matrix-element call nastaví CECOR=0. Parametr není označen jako účinné řízení emisního broadeningu. Změna TAUCORE sama nemusí RXES spektrum změnit.
6. **ME není důkaz tří plných gauge výpočtů.** MEFORM se načítá; NAB vynutí DIPOL=1, GRV připraví radiální derivace. Vlastní běžné maticové elementy jsou ale počítány ME_ALF_ALF, bez předání MEFORM. Pouhé přijetí ADA/NAB/GRV čtečkou nelze zaměnit za ověřenou nezávislou implementaci všech tří forem.
7. **Non-dipole varianty závisí na formalismu.** ME_ALF_ALF_NONDIP se volá jen pro .NOT.FULLPOT a DIPOL=0 (805–812). Pro FULLPOT tedy stejné flags automaticky nezaručují tytéž příspěvky. Potenciálový formalismus není čten čistým InputParameters validátorem.
8. **Jen jeden atomový typ:** NTXRAY=1 a IT se načítá jako jediné číslo; nejde o seznam typů. Výstup má pro každý core stav a energii tři polarizační záznamy, nikoli jediný sloupec intenzity. POLOBS=S dává (+,-,z), C dává (x,y,z).
9. **Stávající XAS má neúčinný default EMAX v Pythonu:** ENERGY(defaults={"EMAX":4.0,...}) tento klíč nepoužívá. XAS sice nativně nastaví EMAX=4 Ry, ale deklarace nevystavuje danou mez jako skutečný člen. Nový XES modul nepřebírá tuto zavádějící konstrukci. Další rozšíření XAS je samostatné.

## Testy

16 nových případů: serializace, opětovné načtení, pickle; nevypisované defaulty; volby typu/core shell/polarizace/šířek/rámce; omezení relativity, mesh a split contour; neplatné IT/DIPOL; fyzikální slupky a limit n=5; NAB override. Celá sada input_parameters: **124 passed**. Numerický SPR-KKR nebyl spuštěn, nebyly ověřeny emise ani .rat data. Patch prošel git diff --check a kontrolou zpětného aplikování.
