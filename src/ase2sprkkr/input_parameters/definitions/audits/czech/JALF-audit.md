# JALF — audit současného tasku

`jalf.py` definuje sériový task `kkrgen` pro relativistickou valenční proudovou hustotu z Diracových matic alfa. Následující omezení vycházejí z nativních zdrojů; numerická správnost nebyla ověřena.

## Výpočet a vstup

JALF používá kkrgen a rutinu CURDNS: valenční proudová hustota z Diracových matic alpha, rozvinutá do reálných sférických harmonik s kartézskými komponentami. Není to paramagnetický gradientový příspěvek JORB. CURDNS současně vyhodnocuje magnetický moment z proudu.

TASK JALF aktivuje CALC_JALF a přejmenuje interní task na E-INTEGRAL (initvar.f:538–550). GEN zavolá AMEBI a alokuje/nuluje JALF_LMCT (gen.f:331–336). Při každé energii volá CURDNS s integrační vahou WETAB, po smyčce finální CURDNS s EFCORRECT=true (1091–1092 a1580–1582). CURDNS nemá vlastní čtečku TASK; definice nepřidává neexistující parametry síťování, vektorového potenciálu nebo rotační symetrie.

## Defaulty

| Volba | Nativní chování | Definice | Zápis |
| --- | --- | --- | --- |
| ENERGY.GRID | 5 pro E-INTEGRAL | 5 | explicitní |
| ENERGY.NE | 30 pro E-INTEGRAL | 30 | explicitní |
| ENERGY.ImE | obecný default 0,01 Ry | 0 Ry | explicitní |
| ENERGY.EMIN | -0,2 Ry | -0,2 Ry | explicitní |
| ENERGY.EMAX | EF při nepřítomnosti | None, delegace | vynechán |
| MODE | plná relativita při běžné nepřítomnosti | None, FREL | vynechán |
| TASK | JALF příznak | JALF | explicitní |

GRID/NE odpovídají init_mod_energy.f:89–92. ImE=0 Ry je explicitní volba ASE2SPRKKR, stejná jako ve společné integrační definici ENERGY. Liší se od nativních 0,01 Ry; nejedná se tedy o chybný vynechaný default. Bez numerického konvergenčního srovnání nelze říci, že nativní hodnota je lepší. Vynechané EMAX správně přebírá EF.

Definice požaduje kladné NE; varuje pro GRID odlišný od výchozího oblouku 5. Nativní zdroj obecně dovoluje jiné kontury, takže není zaveden blanket zákaz. SPLITSS a dvě položky GRID/NE se odmítají: GEN pro E-INTEGRAL nastaví NEPATH=1. Lokální větve SPLITSS v CURDNS samy neprokazují funkční split-contour workflow.

## Nálezy a omezení

1. **MPI vynechá dokončení.** GEN:1355–1368 pro samotný CALC_JALF zastaví před finální větví1580–1582. CURDNS sama nemá MPI redukci JALF_LMCT. Definice proto používá mpi=False. Pouhá paralelní energetická smyčka JALF nedokončí hlavní soubor.
2. **Diagnostika překročí pole při malém cutoffu.** CURDNS alokuje JDNST(NRMAX,NLABIMAX,-1:+1,NTMAX), ale pevná smyčka LR=1,3 zapisuje JDNST(IR,LR,LAM,1) (curdns.f:368–375). KKRMAIN nastavuje NLABIMAX=NLMAX (471). Pro NLMAX<3 tedy diagnostika přistupuje mimo pole. Definice odmítá explicitní SITES.NL s max(NL)<3; cutoff převzatý z potenciálu musí dodatečně zkontrolovat běhový workflow.
3. **Diagnostické curdns soubory směšují typy/meshe.** Soubor používá R(IR,1) a JDNST(...,1), ale IRTOP zůstane po smyčce přes typy hodnotou posledního typu. Mesh 1 navíc nemusí být mesh typu1, protože IMT(1) může být jiný. Tato chyba je v pomocných curdns1..3; hlavní `{dataset}_JALF.dat` vybírá IM=IMT(IT) a IRTOP pro každý typ.
4. **Diagnostická jména nejsou podle DATASET.** curdns1, curdns2 a curdns3 jsou pevné názvy, takže různá data v jednom adresáři kolidují. Do typu1 diagnostiky se navíc může číst nenaplněný typ, pokud ITBOT>1.
5. **Rotational symmetry je compile-time volba.** FORCE_ROTATIONAL_SYMMETRY je PARAMETER=false (curdns.f:45–46). Není to vstupní default, který může uživatel změnit zadáním TASK. Aktivní výstup uchovává obecné harmonické; redukované JDNST slouží magnetickému momentu a diagnostice.
6. **Vektorový potenciál se standardně nepočítá.** CALCVECPOT je při každém volání nastaveno false (curdns.f:93). Po výpisu proudu rutina vrátí na řádku379; následující rozsáhlá větev pro A není dosažitelná bez změny zdroje. Definice nepředstírá podporu TASK.CALCVECPOT.
7. **Relativistický workflow.** CURDNS používá WAVFUN_READ_REL a AME_RLMSIGC bez vlastního scalar-relativistic rozlišení. Definice konzervativně odmítá NREL/SREL/SP-SREL; nejde o doložený samostatný native check JALF. Numerická správnost neplně relativistických vstupů nebyla prokázána.
8. **Korekce EF a renormalizace jsou stav programu.** CURDNS pro EFCORRECT použije SHFTEF jako váhu a při RENORMALIZE násobí uložený proud SCLNOS. Tyto hodnoty předává GEN, nejsou načítány jako TASK volby JALF. Definice je nevystavuje jako fiktivní parametry.

FULLPOT ani ASA nejsou pro JALF ve zdroji explicitně zakázány; zavádět zákaz podle analogie s FMAG by bylo nesprávné. Numerická validace obou režimů zde neproběhla. Totéž platí pro kombinace s DMFT/LDA+U; rutina rozlišuje levé a pravé řešení, ale jejich konzistence nebyla numericky ověřena a definice nezavádí nepodložený zákaz.

## Výstup

Hlavní soubor `{dataset}_JALF.dat` má hlavičku J_alpha a následné TYPE IT bloky s IRTOP, NLAMAX-1 a NLMAMAX. Pro každou harmoniku a radiální bod zapisuje tři kartézské komponenty. Hlavičkové NETAB(1) je počet integračních energií, nikoliv počet řádků hustoty. Proud při akumulaci obsahuje faktor C_AU (curdns.f:255–257). Magnetický moment je v komentáři deklarován v Bohr magnetonech a používá prefaktor PRE a UNITS=0,5 (342–364); numerické jednotky a znaménko momentu nebyly nezávisle validovány.

## Nativní souborový kontrakt

Nativní `kkrmain.f:239–245` nastaví `DATSET0=DATSET` před připojením ADSI k `DATSET`; `{dataset}_JALF.dat` proto ADSI neobsahuje. `curdns.f:279–283,369–377` zapisuje tento soubor a `curdns1`–`curdns3`. `gen.f:1357–1367` ukončí nepodporovaný MPI task před finálním voláním `CURDNS`.

## Rozsah ověření

Regresní testy definice existují; během této revize nebyl spuštěn nový pytest ani numerický SPR-KKR.
