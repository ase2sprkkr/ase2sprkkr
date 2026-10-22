# POSANI — audit SPR-KKR 9.7 / ASE2SPRKKR

**Zdroje:** nativní SPR-KKR 9.7 (`PUB9.7_260420_18h13`) a aktuální definice `ase2sprkkr/input_parameters/definitions/posani.py`. Nativní výpočetní tok a Python souborový kontrakt jsou posuzovány odděleně; fyzikální výsledky nebyly numericky certifikovány.

## Stav

ASE2SPRKKR deklaruje standardní periodický FREL výpočet pozitronové anihilační hybnostní hustoty přes `kkrgen`, s jednou energetickou konturou `GRID=5` a elektronovým potenciálem. Nativní cesta je použitelná s významnými omezeními pro `NQ=NT`; pro obecné `NQ≠NT` obsahuje doložený problém indexování pozitronové `TAUQQ`.

## Skutečný workflow a vstupy

`gen.f:399–426` rozlišuje POSANI, POSANIPREP, POSANITAU a MODE.POSITRON.
Priorita je POSANI > POSANIPREP > POSANITAU > POSITRON. Jsou to odlišné kroky:

1. Elektronový potential a elektronové TAU podklady musí odpovídat systému,
   energiím a použité konvenci vlnových funkcí.
2. POSANIPREP připraví DATSET_pos_dat.pan s pozitronovým stavem. Tento task
   a POSANITAU budou posouzeny samostatně; tato definice je nespouští.
3. POSANI čte **elektronový potential**, TAU adresář a tento připravený soubor.
   `gen.f:431–449` pro POSANI nevolá POSANIPOT. POSVC proto není vstupem tohoto
   kroku; patří do přípravy pozitronového potentialu.

`compton.f:424–426` volá POSANIREAD. Název pozitronového souboru je pevně
DATSET_pos_dat.pan (`posaniread.f:38`), nikoli cesta odvozená z TAUDIR.
Změna DATASET tedy mění očekávané jméno pozitronového souboru; TAUPATH jej nepřesměruje.
Definice záměrně nevymýšlí neexistující keyword pro jeho libovolnou cestu.

Společné vstupy hybnostní mřížky v TASK:
MEMODE=1/2, NPMAX, PNMAX, PPMAX, NPN, NPP, USE_PEQUIV,
PNVEC, PPVEC1, PPVEC2, MD2D, EMD2D, SHFT, SHFT_PP, SHFT_PN,
INSIDE_SPHERE a SQQ_RESOLVED. Čte je COMPTON a INIT_MOD_COMPTON_PMESH.
TAU.TAUDIR a TAU.TAUPATH nastavují elektronové TAU podklady.

Definice povoluje NONMAG, vyžaduje FREL a POINTS. GRID=5 a jedna kontura jsou
omezení auditovaného workflow, nikoli tvrzení, že nativní rutina zakazuje všechny
jiné kontury. ONLY_CORE a USE_ZJOVZZ nejsou vystaveny; vyžadují další audit.

## Důležitý rozdíl proti COMPTON

INITVAR rozpozná COMPTON před INIT_MOD_ENERGY, ale POSANI se hledá až
v GEN (`gen.f:417`); TASK se na POSANI přepne v `gen.f:522`.
INIT_MOD_ENERGY je již voláno v `kkrgen.f:1039`.
Proto samotné TASK POSANI **neznamená automatické převzetí** COMPTON
GRID=5/5, NE=30/50 a SPLITSS=true.

Definice nastaví CONTROL.RDTAU=true a vždy zapíše keyword; INITVAR pak vypne
běžnou energetickou smyčku (`initvar.f:687`). Krok využívá kompatibilní elektronové TAU podklady. Nativní `kkrmain.f:932–960` při jejich nedostupnosti umožňuje přepočet, ale úspěch konkrétní větve není zaručen. ASE2SPRKKR deklaruje TAU jako volitelnou cache a vyžaduje `RDTAU=true` v podporovaném workflow.

COMPTON následně znovu čte SPLITSS z TASK (`compton.f:210`). Definice pro
POSANI vystavuje ENERGY.SPLITSS i TASK.SPLITSS=false a jejich zapnutí odmítne,
aby příprava energií a pozdější výpočet byly konzistentní. Na rozdíl od předchozího
COMPTON tasku se zde žádný SPLITSS keyword nezapisuje.

## Defaults a zápis

| Parametr | Definice | Nativní účinné chování bez parametru | Posouzení |
|---|---|---|---|
| RDTAU | true, vždy zapisuje | INITVAR default false; POSANI jej sám nenastaví | výslovná volba workflow, nikoli shoda nativního defaultu |
| GRID, NE | 5, 30, zapisují se | není POSANI-specific default před inicializací energií; záleží na obecném stavu a TAU podkladech | explicitní volba, musí odpovídat elektronové cache |
| ImE | 0.01 Ry, zapisuje se | 0.01 Ry | stejné |
| EMIN | -0.2 Ry, zapisuje se | -0.2 Ry | stejné |
| EMAX | nevystaven, nevypisuje se | sentinel se převede na EF v INIT_MOD_ENERGY:182 | stejné účinné chování |
| SPLITSS v ENERGY i TASK | false, nevypisuje se | false v tomto postupu | stejné; zapnutí odmítnuto jako jiný workflow |
| MEMODE | 1, zapisuje se | 1 | stejné |
| NPMAX | 250, zapisuje se | 250; interní limit 10000 | stejné default; definice vyžaduje 2..10000 |
| PNMAX, PPMAX | 0, 10 a.u., zapisují se | stejné | PNMAX=0 dá jeden podélný bod a vyvolá upozornění |
| NPN, NPP | 50, 50, zapisují se | stejné; při PNMAX=0 se efektivní NPN změní na 1 | stejný vstupní default |
| USE_PEQUIV | 0, zapisuje se | 0 | stejné, validace omezuje 0..2 |
| směrové vektory | nezadané, nevypisují se | normála z, příčné x/y | stejné |
| ostatní flagy hybnostní mřížky | false, nevypisují se | false | stejné |
| TAUDIR / TAUPATH | nezadané, nevypisují se | DATASET.tau.dir / aktuální adresář | stejné |

Defaults, které se nevypisují, nemaskují odlišný nativní výběr.
NE=30 je výslovný startovní počet, nikoli automatický nativní POSANI default.
Zásadnější než jeho izolovaná hodnota je shoda s elektronovým TAU adresářem.
Nativní PNMAX=0 není celý profil; ponechán s upozorněním, rozsah musí uživatel
zvolit. Konvergence NE, hybnostní mřížky a NPMAX není testy potvrzena.

## Doložená chyba: počty míst versus typů

`posaniprep.f:65–69` alokuje TAUQQ s dimenzemi NQMAX,NQMAX. Připravený soubor
ale zapisuje TAUQQ v cyklech **IQA/IQB=1..NT** (`:268–272`).
`posaniread.f:88–94` obdobně čte TAUQQPOS v cyklech 1..NT, přestože dummy
pole má dimenze NQMAX,NQMAX.

Následná anihilační rutina běží přes fyzická místa IQA/IQB=1..NQ a používá
TAUQQPOS(IQA,IQB) (`compton.f:1043–1079`). Toto pole před čtením není
kompletně inicializováno.

- NQ>NT: část použitých prvků se vůbec neuloží ani nenačte a zůstane neinicializovaná.
- NT>NQ: příprava/čtení zacházejí s dalšími indexy jako s místy; při NT>NQMAX
  jde o přístup mimo pole. Zvětšení NQMAX samo neopraví význam těchto prvků.
- NQ=NT odstraňuje tento konkrétní nesoulad počtů; nezaručuje všechny další
  předpoklady ani numerickou správnost.

Běžné smíšené obsazení nebo více ekvivalentních míst jednoho typu může mít NQ≠NT.
Python vstup bez potentialu tento vztah neurčí. Omezení je proto explicitně
v popisu definice a auditu; nelze slibovat, že validace InputParameters systém
s chybnými počty automaticky zastaví. Nativní Fortran nebyl opraven.

## Další významná ALE

- POSANIPREP používá jen s_1/2 pozitronovou parciální vlnu
  (`posaniprep.f:91–93`); POSANIREAD potlačí imaginární část ZGPOS
  (`posaniread.f:67–70`). Jde o nativní aproximace, nikoli defaults vstupu.
- Reader čte řádky v pořadí typů/míst a neověřuje shodu struktury, radialního
  meshe či pořadí typů se současným potentialem. Záměna připravených podkladů
  není zachycena testy definice.
- Pro směrové vektory nativní kód normalizuje, ale neortogonalizuje zadaný pár.
  Definice vyžaduje konečné nenulové tříkomponentní vektory a kolmost.
  Ignorované PPVEC2 při PNVEC odmítne; bez PNVEC vyžaduje oba příčné vektory.
- Shift flagy se vyhodnocují s prioritou; definice dovoluje nejvýše jeden.
  Posunuté sudé mřížky používají ±INT((N-1)/2) a ztrácejí kraje; v posouvaném
  směru požadována lichá velikost.
- EMD2D násobí podélný vektor dvakrát (`init_mod_compton_pmesh.f:183–185`
  a `:273–279`). Povoleno pouze PNMAX=0, ostatní řezy odmítnuty.
- SQQ_RESOLVED není implementováno pro MD2D/EMD2D
  (`compton.f:354–360`); definice kombinaci odmítá místo tichého vypnutí.
- Při INSIDE_SPHERE a PNMAX=0 je PFEMAX prakticky nulové (1e-10).
- MEMODE=1 vyžaduje nejméně dva interpolační body kvůli dělení NP-1.
- Rutina přepíná GF_CONV_RH na ZJ a SOLVER_FP na BS i při FULLPOT
  (`compton.f:185–194`). Podklady musí odpovídat používané konvenci.

## MPI a výstup

POSANI je povoleno v závěrečné MPI kontrole GEN. Společná COMPTON rutina má
rozdělení energií a redukce výsledků (`compton.f:581–619`, `:1273–1296`).
Poslední energii dostane rank 0 v MPI_DISTRIBUTE, což odpovídá výpočtu core
příspěvku na ranku 0. Na každé aktivní kontuře musí být NE>=NPROCS.
mpi=True má zdrojové opodstatnění, skutečný MPI běh nebyl spuštěn.

Výstupy jsou DATSET_BND.mpa, DATSET_COR.mpa, DATSET_TOT.mpa;
2D výstupy používají příponu _pan2D. Soubory rozlišují sumu a spinový rozdíl.
Interpretaci aproximací a numerickou normalizaci musí ověřit následný výpočet.

## Souborový kontrakt ASE2SPRKKR

Nativní `DATSET` vzniká jako `DATASET_ADSI` (při prázdném ADSI bez podtržítka), `DATSET0` je původní `DATASET`. Šablony používají samostatné `{dataset}` a `{adsi}`; při výchozím `ADSI=POSANI` je povinným vstupem `{dataset}_{adsi}_pos_dat.pan`, připravený taskem POSANIPREP. Elektronový TAU adresář `{dataset}_{adsi}.tau.dir` je volitelná cache s explicitním přesměrováním `TAUDIR` a `TAUPATH`. Zbylé výsledné soubory se určují podle `MD2D`/`EMD2D` a používají nativní `DATSET`.

## Rozsah ověření

Python testy pokrývají parametry, serializaci, podmíněné výstupy a základní souborové názvy. Neprokazují kompatibilitu elektronové cache, správnost připraveného pozitronového souboru, fyzikální hodnoty ani MPI numeriku. Nativní SPR-KKR nebyl v tomto auditu numericky spuštěn.
