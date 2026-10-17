# 41. BSF — audit společné definice BSFEK/BSFKK

Již podporovaný task; opravena společná definice. Předběžně: **má plný smysl**, s omezeními defaults a nativní geometrie. Audit se zaměřuje na volbu E-K/K-K, počty a geometrii; nejde o numerické ověření ani úplnou kontrolu všech tabulek KPATH pro všechny Bravaisovy mřížky.

## Opravy

- NK1/NK2 jsou nyní Integer(min=2). Native BLOCHSF při parametrizaci K-K mapy dělí NK1−1 a NK2−1 (blochsf.f:939–940); předchozí prostý int připouštěl dělení nulou. NK pro E-K cestu má minimum 2, NE minimum 1.
- GRID default je dynamický: pro E-K zůstává 3, pro K-K je 0. Dřívější GRID=3/NE=1 při nativní inicializaci vede k nedefinované jednobodové integrační váze. Změna se skutečně vypisuje. Po explicitním nastavení GRID se uživatel samozřejmě může od defaultu odchýlit.
- K-K validace odmítá nekonečné/NaN směry a nulovou či kolineární rovinu. Kontrola probíhá vůči KA, nikoli pouze podle samotných K1/K2.
- Opraven popis K1/K2: jde o koncové body; native od nich odečítá KA. Pokud nejsou relativní směry kolmé, druhý z nich ortogonalizuje. Takové zadání tedy nemusí vykreslit původně zamýšlený rovnoběžník (blochsf.f:438–455).

## Defaults: vypsané versus nevypsané

E-K: zachovány vypsané NE=200, GRID=3, NK=300, ImE=0.001 Ry, EMIN=−0.2 Ry a EMAX=−1 Ry. **EMAX=−1 je skutečná energie, nikoli sentinel pro EF.** Native init_mod_energy používá náhradu EF až pro EMAX<−99. Výchozí rozsah −0.2→−1 je tedy klesající a pro běžné spektrum překvapivý. Toto je upozornění na explicitní ASE default, nikoli oprava parametru ignorovaného native readerem. Je vhodné nastavit vlastní rozsah, případně v budoucnu default změnit na fyzikálně srozumitelný rozsah vůči EF; bez potenciálu nelze určit vhodnou absolutní horní energii. Default nyní zachován podle zadání auditu defaults.

Native NK při čtení E-K je 800, ASE zapisuje 300. Jde o platnou explicitní odchylku; 800 může být lepší pro jemné cesty, ale konečný počet se rozděluje mezi segmenty a může se upravit podle kapacity k-tabulek. Ani 300 ani 800 není zárukou konvergence.

K-K: vypsaná energie EMIN=0.7 Ry je explicitní ASE volba, EMAX kopíruje EMIN; nejde automaticky o EF. Vypsané NE=1 a nově GRID=0 odpovídají jednomu bodu. Pro hledání Fermiho plochy je třeba nastavit správnou energii či relativní EMINEV. Nevypsané KA native doplní nulovým počátkem; Python validace používá tentýž počátek. NK1/NK2/K1/K2 nemají vymyšlené defaults a musí je dodat uživatel.

Existující mpi=True ponecháno. BLOCHSF má paralelní energetický postup a GEN jej připouští; tento audit však neověřuje numerickou shodu MPI výsledku. MODE/CPA/STRCONST/SITES zůstávají společné sekce beze změn.

## Další omezení

E-K se vybírá skutečným počtem energií NE>1, K-K při NE=1; názvy BSFEK/BSFKK jsou ASE aliasy a hinty, nikoli dva rozdílné native dispatch klíče. TASK se zapisuje jako BSF, CONTROL ADSI jako BLOCHSF. Existující validace varuje při rozporu požadovaného aliasu a NE.

V E-K režimu native rozděluje body podle délky segmentů, s minimem dvou bodů na segment. Degenerované explicitní segmenty a nulová celková délka mohou vést k problematickému rozdělení. Jejich úplná kontrola, spolu s geometrií každé KPATH/Bravais kombinace, není v tomto patchi řešena.

Pro 2D native v K-K větvi nuluje z-složky K1/K2 před odečtením KA. Input-only kontrola nezná dimensionalitu potenciálu, takže rovina platná podle zadaných 3D vektorů může po nativním přepisu degenerovat. Pro 2D je nutné explicitně konzistentní rovinné zadání.

Při native čtení TAU cache může BLOCHSF převzít počet energií z cache; konzistence externích souborů zde není ověřena. Tento audit nedoplňuje advanced BSFSPN/BSFORB/BSFETD/BSFALF ani EPHOT.

## Ověření

11 nových testů, celá sada InputParameters **119 passed**, dvě existující NumPy deprecation warnings. Testy kontrolují K-K GRID=0, rovnost EMIN/EMAX, round-trip/pickle, zachování E-K defaults, nenulové KA s endpointy, minimální počty, degeneraci, kolinearitu a NaN. Existující testy nadále ověřují sdílení aliasů, změnu režimu a základní gramatiku.

Příkaz: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

Zdroje: blochsf.f:220–268,310–357,419–455,939–940; init_mod_energy.f:182; epath.f (GRID=3 jednobodová větev). ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Samostatný patch pouze pro BSF. Žádný numerický SPR-KKR běh.

## Aktuální integrační stav (2026-10-17)

Source audit výše zůstává historickým záznamem původního baseline definice. Na větvi `tasks` nyní společná definice `bsf` modeluje TAU directory jako volitelnou reusable cache a explicitně manifestuje nativní persistentní výsledky `{dataset}_BLOCHSF.bsf` a `{dataset}_BLOCHSF_spol.bsf` pro `BSF`, `BSFEK` i `BSFKK`. Aliasy nemají vlastní odlišná jména souborů, protože všechny tři volí nativní `ADSI=BLOCHSF`.

Tato integrace souborového kontraktu je přítomná v commitu `032d932b2e05aa3a540ce232c98522b1fb7065c6`. Historický počet testů výše nepřepisuji; pozdější regresní test manifestu pokrývá všechna tři veřejná jména tasku.
