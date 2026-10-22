# MECHECK — audit SPR-KKR 9.7 a nová definice ASE2SPRKKR

Podklad: přiložené zdroje PUB9.7_260420_18h13; ASE2SPRKKR `refactor`, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`. V této revizi chybí definice MECHECK. Patch přidává samostatný modul `definitions/mecheck.py` a `test/test_mecheck.py`; nezávisí na předchozích patchích WFPLOT ani MEPLOT.

## Rozhraní

MECHECK je diagnostika maticových elementů, kterou `gen.f:1455–1459` předává do `ME_TEST_11`. `initvar.f:162–167` nastavuje ITEST=11 a vyžaduje KKRGEN. Definice používá `kkrgen`, MPI=False: rutina nemá rozdělení své práce mezi ranky a zapisuje sdílené názvy souborů. Numerické spuštění ani chování MPI nebylo ověřeno.

| Parametr TASK | Default při vynechání / chování | Definice |
|---|---|---|
| ITSEL | 1; index atomového typu, nikoli site | kladné celé číslo, 1 se nevypisuje |
| PRINT | zachová CONTROL.PRINT | dynamický default, shodná hodnota se nevypisuje; zde skutečně PRINT, nikoli MEPLOT.IPRINT |
| INITIALSTATE | BAND | BAND nebo CORE; default se nevypisuje |
| CL | 2P pro CORE | shell n=1..6, S/P/D/F; validace l<n; default se nevypisuje |
| ERYDI | BAND: efektivně {0.5,0.01} Ry | přesně dvě reálné složky, default se nevypisuje |
| ERYDF | BAND: efektivně {0.8,0.01} Ry | přesně dvě reálné složky, default se nevypisuje |
| NONMAG | false | TASK flag potlačuje BT vybraného typu; odlišný od CONTROL.NONMAG |
| NMEFORM | reader začíná na 3, ale efektivní hodnota je vždy 2 | reader default 3 se nevypisuje; explicitní nastavení vyvolá upozornění, že se přepíše |

Zdroj: `me_test_11.f:131–177, 297–346`. BAND vyhodnocuje jen jednu dvojici komplexních energií; CORE počítá počáteční stav z CL a finální stavy z první mřížky ENERGY. ERYDF se pro CORE sice načte, ale přepíše ETAB; ERYDI se v CORE ani nenačítá. Na explicitní ignorované volby definice upozorňuje.

Pokud je absolutní hodnota BAND energie menší než 1e-6, zdroj ji nahradí výše uvedeným defaultem. Požadovat přesně nulu tedy tímto vstupem nelze. Definice upozorňuje a vstup nepřepisuje. Fortranový reader přijímá i jednu reálnou složku a kopíruje ji do druhé (`scan_input.f:698–707`): skalár zde neznamená komplexní energii s nulovou imaginární částí. Definice vyžaduje explicitní dvojici, aby tato interpretace nebyla skrytá.

## Omezení

- `check_input.f:99–100` vyžaduje IREL=3. NREL/SREL/SP-SREL definice odmítá.
- `gen.f:211–222` explicitně odmítá DMFT a LDAU i pro MECHECK. Stejné omezení je v definici.
- FULLPOT není pro MECHECK ve zdroji obecně zakázán; rutina volá FPSSITE a má odpovídající větve. Definice proto nevymýšlí zákaz. Jeho numerická správnost zde ověřena není; režim se odvozuje z potenciálu.
- SITES.NL zůstává None a nevypisuje se: cutoff určí SPR-KKR/potenciál. U explicitního uniformního cutoffu se kontroluje dostupnost CORE orbitalu. Neuniformní cutoff a dostupnost shellu konkrétního atomového typu bez potenciálu ověřit nelze.

## Problémy ve zdrojích SPR-KKR

1. **NMEFORM nemá účinek.** Po načtení a případném FULLPOT omezení následuje `NMEFORM=3` a při pevně nastaveném `CMP_NAB_VS_NAB=.FALSE.` další přepsání na 2 (`me_test_11.f:96, 172–177`). Běží NAB a ADA; GRV nelze vstupem zapnout.
2. **Druhá kontrola popsaná v komentáři neběží.** Kompilační konstanta `K_CALC_ME=2` nastaví NTEST=1 (`:49–50, 370–374`). Nahrazení J pomocí Z a následná kontrola nepravidelných členů jsou nedostupné; nejde o chybějící input parametr.
3. **Souhrnné statistiky se nevypíšou.** Smyčka `DO IMEFORM = 2, - NMEFORM` (`:636`) s implicitním krokem +1 nemá při NMEFORM=2 žádnou iteraci. Pravděpodobný překlep v horní mezi; dílčí porovnání před ní běží.
4. **Vedlejší zápis na unit 601.** NAB větev volá CMATSTRUCT i s jednotkou 601 (`:478–480`), pro niž tato rutina neobsahuje OPEN. Může vznikat implicitní soubor `fort.601` podle Fortran runtime. Je to podezřelý zbytek diagnostiky, nikoli parametr vstupu; runtime dopad nebyl ověřen.

## Defaults ENERGY

ASE2SPRKKR zde explicitně zapisuje GRID=3, NE=1, ImE=0.01 a EMIN=-0.2 Ry. Jsou to záměrné volby jednoduché diagnostiky, nikoli tvrzení o nativních defaults. Pro BAND je rutina při výběru obou stavových energií ignoruje. Pro CORE skutečně definují finální energetickou cestu; pro spektrum je třeba zvýšit NE a nastavit vhodné meze. NE=1 je vhodné pro základní kontrolu, pro CORE sken by vyšší default byl praktičtější; změna je návrh uživatelského chování, nikoli oprava parseru. EMAX=None se nevypisuje a deleguje skutečnou Fermiho energii SPR-KKR. Pro konkrétní okno lze použít párové EMINEV/EMAXEV.

## Ověření

18 nových testů: automatické nalezení tasku, omitted defaults, zápis/načtení komplexních dvojic a souboru, pickle/copy, dědičnost PRINT, CORE sken, neplatné shelly a cutoff, near-zero energie, ignorované NMEFORM/ENERGY, relativita a nepodporované korelace, oddělený TASK.NONMAG.

Celá sada `src/ase2sprkkr/input_parameters/test`: **126 passed**. Varování při zpětném načtení BAND vstupu jsou očekávaná: explicitně zapsaná ENERGY se při výběru BAND energií nepoužívá. `git diff --check` a obrácená aplikovatelnost patchu prošly. Nebyl spuštěn SPR-KKR; ověřena je definice vstupu a její shoda s přečtenými zdroji, nikoli numerické výsledky. Výstupní numerické soubory zatím nemají specializovaný parser.
