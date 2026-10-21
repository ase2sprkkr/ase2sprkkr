# XAS — audit aktuální implementace

**Stav:** Rentgenová absorpční spektroskopie v nativně podporovaných režimech; pouze sériový výpočet. Auditovaná implementace je na `tasks-workflow`, nativní referencí je SPR-KKR 9.7 (`PUB9.7_260420_18h13`). Nejde o numerickou certifikaci.

## Doložené chyby původní definice

4. **CL jako volný str nehlídal fyzikální slupku ani hranici tabulky.** XRAYSPEC má SHELL(5) s K,L,M,N,O a přímo indexuje NCXRAY na491. n=6 by překročilo pole. Nově jsou povoleny n=1..5, s/p/d/f a l<n. Existence slupky v konkrétním potenciálu se ověřuje až mimo tento validátor.
6. **Frame popisky zaměňovaly rotaci rámce za jeden vektor světla.** INIT_MROT_FRAME čte tři Eulerovy úhly a tvoří GETMROT(-gamma,-theta,-phi). Nejde o prostý parametr jednoho elektrického pole; XRAYSPEC počítá tři polarizační kanály. Popisky byly opraveny a chybějící FRAMEGAM doplněn.

## Parametry a validace

Doplněny ME, POLOBS, DIPOL, GFCONVRH, WRREGIRR, NONDIPE2 a NONDIPM1 podle XRAYSPEC:295–329,442–465. POLOBS=S vypisuje (+,-,z), C vypisuje (x,y,z). ME=NAB přepíše DIPOL na1. Pokud jsou přítomné oba non-dipole flags, NONDIPE2 má přednost. Definice na obě kombinace upozorňuje.

Validace odmítá jiné MODE než FREL, neplatné core slupky, nečíselné rámcové úhly, nekladné/jednobodové NE a více energetických cest. CHECK_INPUT:43–47 přímo požaduje IREL=3. XRAYSPEC:376–380 požaduje shodu ENERGY.NE s počtem bodů prvního TAU contouru; obsah TAU a potenciálu modul nečte, takže tuto shodu nevydává za ověřenou.

Nepřidány MESPLIT (po čtení je K_ME_SPLIT vždy true), TAUCORE jako údajně účinné řízení XAS broadeningu ani DYNAMICAL bez auditu jeho samostatné rutiny. TAUCORE se společnou čtečkou načte, ale jeho účinné použití je ve větvích XMO/XRS; standardní call ME_ALF_ALF nastavuje CECOR=0. DYNAMICAL pro RXAS/XMO skutečně přepíná na XRAYSPECDYN, což není tento auditovaný statický workflow.

## Defaulty: zapisované a vynechané

| Volba | Nativní XAS | Opravená definice |
| --- | --- | --- |
| ENERGY.GRID | 6 | 6 explicitně, beze změny |
| ENERGY.NE | 180 | 180 explicitně, beze změny |
| ENERGY.ImE | 0,01 Ry | 0,01 explicitně, beze změny |
| ENERGY.EMIN | po obecné čtečce XAS default EF | -0,2 Ry explicitně, ponecháno |
| ENERGY.EMAX | 4 Ry | 4 Ry explicitně, nově skutečný parametr |
| TASK.IT / CL | 1 / 2P | 1 / 2P explicitně, ponecháno |
| OUTPUT | BARN, výstup v Mb | None vynecháno; přijímá BARN, MBARN, SIGMA |
| POLOBS | S | S vynecháno |
| ME / DIPOL | ADA / 0 | ADA / 0 vynecháno |
| GFCONVRH | false při nepřítomnosti keywordu | false vynecháno |
| WRREGIRR | efektivně false: SECTION_FIND_KEYWORD přepíše DATA true | false vynecháno |
| Frame úhly | 0/0/0, bez keywordu CHANGE_FRAME=false | stejné nuly vynecháno |
| MODE | FREL | None vynecháno |

EMIN=-0,2 Ry není chybně vynechaný default: zapisuje se. Nativní EF je pro standardní absorpční interval zpravidla přirozenější a vyhne se části pod EF. Neměním ho automaticky, protože existující ASE2SPRKKR chování může sloužit záměrnému scanu širšího intervalu. Konkrétní vhodnost závisí na EF a požadovaném spektru; konvergence nebyla testována.

Výstup SIGMA je v10^15/s, ostatní OUTPUT používá Mb. Původní DefKeyword neměl explicitní základní default; vynechání OUTPUT správně nechávalo nativní BARN. Nový Keyword uchovává toto vynechání a doplňuje přímo nativní název BARN vedle kompatibilního MBARN.

## Další nativní omezení

- MEFORM se načítá, NAB ovlivňuje DIPOL a GRV připravuje derivace. Standardní výpočet přesto volá ME_ALF_ALF bez předání MEFORM. Nelze tvrdit, že tři přijaté názvy znamenají tři nezávislé, numericky ověřené gauge implementace.
- Non-dipole call je jen pro .NOT.FULLPOT a DIPOL=0. FULLPOT proto nemusí vyhodnotit stejné příspěvky při stejných flags.
- Výpočet vybere jediný typ IT (NTXRAY=1). Nevystavuje seznam typů.
- Počet bodů, mesh a druh potenciálu se musí shodovat s TAU directory. Pouhá úspěšná serializace inputu tuto podmínku netestuje.
- XAS a XES se stejným DATASET, typem a slupkou mohou používat stejný .rat název. Je nutné rozlišovat výstupy při práci s oběma tasky.

## Aktuální souborový kontrakt

Vstupy, volitelná cache a výstupy jsou určovány aktuální definicí tasku, nikoli starým odděleným patchem. Zdroj: `src/ase2sprkkr/input_parameters/definitions/xas.py` a související task registry.

## Ověření

Testy definic jsou v `src/ase2sprkkr/input_parameters/test/`. Nový kompletní pytest ani numerický SPR-KKR běh se zde netvrdí.
