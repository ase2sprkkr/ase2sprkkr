# TORQUE — audit existující definice a oprava validace

Dodaný SPR-KKR 9.7 porovnán s ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. TORQUE již podporován je. Patch 0023-torque.patch upravuje existující torque.py a přidává test_torque.py; je samostatný proti uvedenému základu. FORCETEN byl přeskočen, protože forcetensor.f obsahuje pouze výpis nedostupnosti a STOP. Bude uveden až na úplném konci závěrečného společného auditu.

## Co znamenají parametry

THETAQ/PHIQ v TASK nejsou rotace magnetizace. Určují jednotkový vektor u, osu vyhodnocovaného torquu. TORQUE:74–75 nejprve vyžaduje MOMENTS_ROTATED; KKRGEN:846 nastaví tento příznak pro KMROT=1 nebo2. MROTQ pochází z natočení momentů. TORQUE:168–175 kontroluje m dot u = 0 s tolerancí1e-8. Úhly TASK samy tuto podmínku ani rotaci momentů nevytvoří.

Pořadí hodnot je podle site IQ, nikoli atomového typu IT. SECTION_SET_REAL_ARRAY (scan_input.f:701–707) opakuje poslední dodaný prvek do zbývajících site. Jednoprvkový seznam tedy skutečně platí pro všechny site. Delší pole nesmí překročit NQ; definice bez potenciálu nezná NQ a tuto mez nekontroluje. Délky obou seznamů nemusejí být stejné, protože každý se rozšiřuje samostatně.

## Defaulty

| Parametr | SPR-KKR při vynechání | Existující ASE2SPRKKR a patch |
| --- | --- | --- |
| THETAQ/PHIQ | když chybí oba:45°/0°; pokud je zadán jen jeden, druhé pole zůstane0 | explicitně90°/90°, ponecháno |
| ENERGY.GRID | TORQUE nemá vlastní inicializační default | explicitně8, ponecháno |
| ENERGY.NE | bez přepsání MOD_ENERGY začíná0 | explicitně36, ponecháno |
| ENERGY.ImE | 0,01 Ry | explicitně0 Ry, ponecháno |
| ENERGY.EMIN | -0,2 Ry | explicitně -0,2 Ry, ponecháno |
| ENERGY.EMAX | obecná náhrada sentinelu za EF | vůbec není členem této definice, nativně EF |
| MODE | běžně FREL | None, vynecháno |

90°/90° dává u podél y. Nativní45°/0° dává směr v rovině xz. Ani jedna volba není univerzálně lepší: u musí být kolmé ke skutečné magnetizaci. Default ASE2SPRKKR je explicitně zapisovaný, takže jeho odlišnost sama není omyl. Případnou změnu je vhodné volit společně s konvencí natočení momentů; tento patch default nemění.

GRID=8 používá logaritmickou integrační konturu. EPATH:249–264 nastavuje pevné EPS=0,050 a předává ho do EGRID8_IME, místo aby tento parametr odvodil z ENERGY.ImE. Samotný rozdíl explicitního ImE=0 vůči nativnímu0,01 tak neprokazuje rozdílnou konturu GRID8. Vhodnost GRID8/36 proti jiné integrační kontuře vyžaduje konvergenční srovnání.

Původní ENERGY(defaults={"EMIN": -0.2,...}) předával EMIN do slovníku, který jej nepoužívá. Faktická hodnota byla správně -0,2 ze standardního argumentu emin. Patch odstraňuje tento zavádějící neúčinný klíč; hodnotu nemění.

## Opravy definice

- THETAQ/PHIQ jsou neprázdné seznamy reálných hodnot; při konečné validaci se odmítá NaN/Inf. Původní definice tyto fyzikálně nepoužitelné směry nehlídala.
- ENERGY.NE musí být kladné.
- Odmítá se SPLITSS a více cest GRID/NE: TORQUE čte vždy WETAB(IECURR,1) a ETAB(IECURR,1), nikoli obecnou cestu. Bez ochrany mohl vstup vzbudit dojem výpočtu obou příspěvků.
- Dokumentace parametrů rozlišuje osu torquu od natočení magnetizace a popisuje rozšiřování seznamů.

NREL/SREL a NONMAG již mají varování v původní definici; patch je ponechává. Jejich text tvrdí, že výpočet selže, ale vlastní TORQUE nemá kontrolu NONMAG/IREL. Používá plně relativistické operátory J=L+S a závisí na okolním workflow. Skutečnou dostupnost jiných režimů nelze rozhodnout pouhým tímto varováním; do definitivního nativního zákazu se nepřevádí. SP-SREL původní varování nepokrývá. Toto zůstává nedořešená otázka podpory režimů, nikoli doložený funkční scalar-relativistic TORQUE.

## Nálezy a omezení nativního výpočtu

1. Požadavek MOMENTS_ROTATED a kolmosti u vůči m není ze samotných InputParameters ověřitelný. Existující default vstup se standardním nenatočeným potenciálem může skončit STOP.
2. Rutina nastavuje IPRINT=1 bez ohledu na CONTROL.PRINT a produkuje rozsáhlé matice i v tichém vstupu.
3. J a UJ se plní do globálního NKM, ale lokální násobení používá NKMQ(IQ). WKM4 pak kopíruje rozdíl do 1:NKM, tedy i mimo lokálně právě naplněné submatice při rozdílných cutoffs. Konkrétní dopad na trace je nutné ověřit; není to dokázaná chyba všech výsledků.
4. Při WGT_OCC se bere pouze ITOQ(1,IQ). Site s více komponentami tak vyžaduje kontrolu významu vážení podle prvního typu.
5. Hlavní výsledky jsou textový výpis site torquů a součtu; jednotlivé T(IQ) jsou v interních Ry, součet se také vypíše přepočtený na meV. TASK nemá vlastní volbu výstupního souboru ani nový Python parser.

MPI=true je zde podložené: TORQUE:186–199 redukuje T přes MPI_SUM_R a na procesu0 dokončí výpis; GEN:1264 má dodatečné volání pro energetickou paralelizaci. To odlišuje TORQUE od nedokončených MPI větví některých předchozích tasků. Numerické MPI výsledky nebyly ověřeny.

## Testy

11 nových případů: serializace, načtení a pickle, seznamy různých site, odmítnutí prázdných a nekonečných úhlů, split contour a nekladného NE. Celá sada input_parameters: **119 passed**. SPR-KKR nebyl spuštěn; testy neověřují natočení momentů, kolmost ani numerický torque. Patch prošel git diff --check a kontrolou zpětného aplikování.
