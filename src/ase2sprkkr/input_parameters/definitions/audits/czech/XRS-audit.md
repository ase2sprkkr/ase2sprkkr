# XRS — experimentální definice a audit

SPR-KKR 9.7 z dodaných zdrojů versus ASE2SPRKKR refactor po fetchi, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. XRS zatím nemá definici. Samostatný patch 0026-xrs.patch přidává xrs.py a test_xrs.py; nevyžaduje XAS/XES patche.

## Skutečné rozhraní

XRS počítá komplexní tensor rezonanční rentgenové rozptylové amplitudy ve společné rutině XRAYSPEC. Na rozdíl od RXES má vlastní činitele, rezonanční jmenovatele a výstup SIG. Vstupní TASK je XRS. CHECK_INPUT vyžaduje FREL. INITVAR aktivuje IXRAY, čtení TAU a energetickou smyčku.

Pro XRS se používají tři energetické integrace: první s TAU a odčtením single-site části, druhá se single-site částí na reálných energiích první sítě a třetí přes EPATH GRID7 od ENERGY.EMAX do TASK.EMAX3 s TASK.NE3 body (xrayspec.f:684–702,1213–1220). TSELECT=0 zahrnuje oba členy; 1 vynuluje nerezonanční WTERM2, 2 rezonanční WTERM1 (413–419).

TASK obsahuje IT, CL, ME, POLOBS, DIPOL, GFCONVRH, WRREGIRR, non-dipole flags, tři Eulerovy frame úhly a TAUCORE. Šířky TAUCORE jsou v eV, nejvýše dva prvky pro core kappa větve; společná čtečka opakuje poslední hodnotu. Hodnoty pod1e-6 přebírají WCOREHOLE. Na rozdíl od XAS/XES zde mají prokazatelné použití v rezonančních jmenovatelích (1259,1264).

OUTPUT není nabídnuto: XRS jej přepisuje na SIGMA bez ohledu na vstup. MESPLIT se přepisuje na true. DYNAMICAL se pro XRS nečte. Tyto volby by nepředstavovaly funkční přepínače tasku.

## Defaulty a zápis

| Volba | Nativní efekt při vynechání | Definice |
| --- | --- | --- |
| IT / CL | 1 / 2P | stejné hodnoty, vynecháno |
| ME / POLOBS / DIPOL | ADA / S / 0 | stejné hodnoty, vynecháno |
| NE3 / EMAX3 / TSELECT | 40 / 8 Ry / 0 | stejné hodnoty, vynecháno |
| TAUCORE | [0,0], delegace na WCOREHOLE | None, vynecháno |
| Frame úhly | 0/0/0, bez přítomnosti CHANGE_FRAME=false | None, vynecháno |
| GFCONVRH / WRREGIRR | false podle nepřítomnosti keywordu | false, vynecháno |
| ENERGY.GRID / NE | XRS nemá vlastní default větev; NE začíná0 | explicitně6/180 |
| ENERGY.EMIN | XRS následně nastaví na EF | vynecháno, sekce nemá dolní mez |
| ENERGY.EMAX | obecné přepsání sentinelu na EF | explicitně4 Ry |
| ENERGY.ImE | 0,01 Ry | explicitně0,01 Ry |

GRID6/180 a horní mez4 Ry představují explicitně zvolený rozsah neobsazených stavů, obdobný nativnímu XAS scanu. XRS sám nemá nativní EMAX4; jeho obecné EF by vedlo k prvnímu intervalu EF..EF, proto není pro tento vystavený scan vhodnější. Vynechaná EMIN přebírá skutečné EF a neprezentuje neúčinný absolutní default. Konvergence prvního i třetího intervalu nebyla numericky ověřena. ImE odpovídá nativní hodnotě.

Validator kontroluje FREL, fyzikální core shell n<=5/l<n, kladné IT/NE3, NE>=2, jednu energetickou cestu, finite úhly a šířky, nezáporné TAUCORE a EMAX3>ENERGY.EMAX při známé absolutní mezi. Horní mez IT, existence slupky a shoda ENERGY.NE s prvním TAU contourem vyžadují skutečný potenciál a TAU directory.

## Závažné chyby a omezení

1. **SIG je inicializováno jen částečně:** alokace SIG(NEMAX+50,NCSTMAX,3,3) na636, ale CINIT má délku NEMAX*NCSTMAX*3*3 na639. Následně se sčítá přes NESIG=50+NE. Kvůli sloupcovému uložení jde o různé oblasti core/polarizačních kanálů, nikoli pouze posledních50 energií. Výsledek může obsahovat neinicializovaná data; oprava vyžaduje Fortran.
2. **Dodatečná síť není monotónní:** před prvním energetickým bodem se generuje50 bodů výrazem (NESIG-IE-1)^3. IE=49 už odpovídá první energii a IE=50 leží lehce nad ní, pak se znovu připojí první energie TAU. Metadata počtu i třídění osy musí čtečka výsledků respektovat.
3. **MPI nedojde ke spektru:** GEN:1355–1368 pro XRS skončí před voláním XRAYSPEC; definice má mpi=False.
4. **Název a hlavička .rat zůstávají rentgenové absorption konvence:** pro stejný DATASET, typ a slupku může XRS kolidovat s XAS/XES. Výstup XRS zapisuje pro každou energii a core stav9 komplexních komponent SIG. Nelze jej interpretovat jako3 absorpční polarizační intenzity.
5. **Výstupní počet bodů není prosté ENERGY.NE:** NESIG začíná50 a připojí NE; současně pozdější třetí integrace změní NE na NE3. Čtečka nesmí odvozovat počet řádků z aktuálního NE nebo jeho původní hodnoty bez znalosti struktury.
6. **Non-dipole flags nezaručují non-dipole XRS tensor:** jejich pole se plní, ale vlastní XRS SIG součet používá MZAZB/MZBZA a jejich conjugate varianty. Vystavené flags odpovídají společné čtečce; účinek na plný rozptylový tensor vyžaduje další fyzikální audit. ME také nelze vydávat za ověřenou nezávislou implementaci tří gauges.
7. **Rezerva pro třetí integraci:** pokud NE3>NEMAX, XRAYSPEC přealokuje ETAB/WETAB/ETAB0. Další dříve alokovaná pole tím automaticky nevyrostou. Přístup k těmto polím ve třetí větvi je nutné ověřit před používáním NE3 nad původní kapacitou; tento patch chyby alokací neopravuje.

Na rozdíl od FORCETEN nejde o vypnutou placeholder rutinu, proto byl vstup doplněn. Kvůli částečné inicializaci SIG patří v plánovaném společném auditu mezi rozbité tasky do opravy Fortranu. Žádný numerický výsledek není touto definicí certifikován.

## Testy

16 nových případů: serializace/načtení/pickle, nevypisované nativní defaulty a explicitní EMAX, rezonanční parametry a TAUCORE, relativity, core slupky, split contour, IT/NE3/TSELECT, finite a pořadí třetí meze. Celá sada input_parameters: **124 passed**. SPR-KKR nebyl spuštěn. Patch prošel git diff --check a kontrolou zpětného aplikování.
