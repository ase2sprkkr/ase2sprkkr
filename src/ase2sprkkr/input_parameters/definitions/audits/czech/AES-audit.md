# AES — audit současné definice

AES má na `tasks-workflow` vlastní definici `aes.py`. Nativní úloha je plně relativistická Augerova spektroskopie: přednostně čte existující TAU, ale `kkrmain.f:932–960` při jeho nedostupnosti může spustit přepočet. ASE2SPRKKR jej proto deklaruje jako volitelnou cache, nikoli povinný vstup. Závěry vycházejí z dodaných zdrojů SPR-KKR 9.7 (`PUB9.7_260420_18h13`), nikoli z numerické validace.


## Workflow a vstup

Plně relativistická Augerova elektronová spektroskopie v kkrgen. initvar.f:299–308 nastaví INITELOOP=true, RUNELOOP=false, RDTAUDIR=true, WRTAUDIR=false a USELPLS1=true. Při dostupné cache rutina AES čte počáteční TAU adresář; nativní fallback může data nejprve dopočítat. Není to postprocessing DOS souboru jako NRAES. Výchozí název TAU v ASE2SPRKKR je `{dataset}_{adsi}.tau.dir` a explicitní `TAUDIR`/`TAUPATH` jej mohou změnit.

TASK čte IT, CL, NEME, PRINTME a PRINTAME (aes.f:128–144). Volby CALCLS, OUTGNU či GAMMA z nerelativistických větví sem nepatří. Zdroj pracuje pouze s NQ=1 (aes.f:230). Definice toto omezení dokumentuje; obecný validátor InputParameters nemá metadata potenciálu, a proto nemůže počet míst ověřit. NQ=1 neznamená nutně NT=1 u slitiny. Explicitní IT má být index typu existujícího v potenciálu.

check_input.f zakazuje FULLPOT a vyžaduje IREL=3. GEN zakazuje DMFT a LDA+U. Definice odmítá nekompatibilní MODE a obě korelace; FULLPOT a NQ musí ověřit SPR-KKR podle potenciálu. MPI je vypnuto, protože AES není v seznamu tasků pokračujících po MPI energetické smyčce GEN.

## Defaulty

| Volba | Nativní chování | Definice | Zápis |
| --- | --- | --- | --- |
| TASK.IT | 1 | 1 | default vynechán |
| TASK.CL | 2P | 2P | default vynechán |
| TASK.NEME | 10 | 10 | default vynechán |
| PRINTME / PRINTAME | false při nepřítomnosti | false | vynechány |
| ENERGY.GRID / NE / ImE | bez specifického AES meshe; počáteční mesh musí odpovídat TAU | 3 / 300 / 0,01 Ry | explicitní |
| ENERGY.EMIN | -0,2 Ry | -0,2 Ry | explicitní |
| ENERGY.EMAX | EF | None, delegace na zdroj | vynechán |

Vynechané defaulty TASK odpovídají zdroji. GRID=3 a NE=300 jsou explicitní volbou definice pro počáteční mesh a kapacitu, nikoliv tvrzením o nativních AES defaultech. Při použití existujícího TAU musí vstupní mesh odpovídat souboru; fixní 300 není detekcí jeho obsahu. Lepší řešení by bylo odvodit mesh a kapacitu z metadat konkrétního TAU. Vynechané EMAX správně ponechává EF z potenciálu.

NEME je hrubá síť matrix elements, skutečně použitý počet je min(TASK.NEME, NE z TAU). Definice povoluje NEME>=3, protože interpolace používá sousední body na obě strany. ENERGY.NE>=3 ověřuje jen alokační kapacitu, nezaručuje skutečné NE>=3 v souboru.

## Nálezy zdroje

1. **NEME=1 či 2 není bezpečné.** DEME=(EFERMI-real(EBOT))/(NEME-1) (aes.f:292–293) pro 1 dělí nulou. Interpolace clampuje index na NEME-1 a používá index-1 i index+1 (řádky 513–540). Pro 2 se může dostat na index 0. Proto definice odmítá TASK.NEME<3. Stejný problém nastane, pokud má existující TAU méně než 3 energie, i při nativním NEME=10.
2. **F core slupky mohou překročit tabulku cache názvů.** STRME má jen ME1..ME12 (aes.f:79–80), zatímco f slupka dává NCST=14. Smyčka přes core stavy použije STRME(ICST) při nepřítomnosti společného souboru ME (řádek 363). Definice f slupky odmítá, i když výjimečná větev s existujícím společným ME tabulku obejde. To je konzervativní omezení nové definice do opravy Fortranu.
3. **Názvy hlavních slupek končí u n=5.** SHELL obsahuje pouze K,L,M,N,O (aes.f:81). Definice povoluje n<=5 a kontroluje l<n.
4. **Automatická cache ME je bez ověření původu.** Rutina sama hledá soubor ME a následně MEi, čte matrix elements a neověřuje shodu potenciálu ani energetické sítě. Existující cache také vypne PRINTME. Stejný pracovní adresář s dřívějšími výpočty může tedy změnit výsledek i bez vstupního přepínače. Čtečka používá indexy z cache přímo, bez kontroly rozsahu. Definice tuto nativní cache nespravuje.
5. **Kapacita TAU není lokálně ověřena.** NE se přebírá z TAU metadat po alokaci podle NEMAX. ETAB, IAES a ETABFIN jsou potom používány podle tohoto NE (řádky 226–236, 484–486). Je třeba NE<=NEMAX a soulad vstupu s TAU; test serializace to neověří.
6. **Mesh musí být rovnoměrný a rostoucí.** DE se bere pouze z prvních dvou ETAB energií a celá konvoluce používá tento krok. Matrix-element síť končí u EF, zatímco EBOT je první TAU energie. Pokud EBOT=EF, DEME=0 i pro jinak platné NEME. Zdroj neověřuje tyto předpoklady přímo v AES.
7. **Výstup není jednoduchá tabulka NE řádků.** Hlavička je vytvořena s NE počátečních energií, ale pro každý core stav se zapisuje 2*NE-1 konečných energií. Každý bod má dva řádky a intenzity jsou komplexní hodnoty. Specializovaný parser musí tuto strukturu respektovat.
8. **Výstupní energie zůstávají v Ry.** aes.f:657–659 vypisuje real(ETABFIN) a relativní rozdíl bez RY_EV; nelze sem přenést interpretaci eV ze spektrálního výstupu NRAES. Konečné energie se sestavují z 2*EF-ECOR s broadeningem převzatým z první počáteční energie (řádky 483–487).

Numerická správnost AES, chování cache a konkrétní TAU soubor nebyly ověřeny běžícím SPR-KKR. Definice nemění nativní implementaci Fortranu.
