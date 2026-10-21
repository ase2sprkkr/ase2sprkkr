# JORB — audit současného tasku

`jorb.py` definuje sériový task `kkrgen` pro relativistickou paramagnetickou orbitální proudovou hustotu. Audit pokrývá nativní zdrojovou větev; numerická správnost nebyla ověřena.

## Význam a vstup

JORB aktivuje CALC_JORB. initvar.f:539–550 přepne interní TASK na E-INTEGRAL a nastaví NOWRDOS=true, RDTAUDIR=false a WRTAUDIR=false. Nejde o čtení existujícího DOS ani TAU directory. GEN pro E-INTEGRAL nastaví zápis regulárních i irregulárních vlnových funkcí, CALCINT=true, GETIRRSOL=true, ICORE=1 a NEPATH=1.

V každé energii volá GEN ORBCURDNS s WETAB(IE,IEPATH) (gen.f:1098–1099); po integraci následuje finální ORBCURDNS s EFCORRECT=true (1586–1588). Rutina popisuje relativistickou paramagnetickou proudovou hustotu valenčních elektronů v atomic units, ve full-potential formalismu. Výsledkem nejsou Stonerův parametr ani magnetický form factor.

ORBCURDNS nemá vlastní TASK čtečku. Definice proto nepřidává neexistující volby pro síť nebo komponenty. Obecné QREPONLY čte společná větev E-INTEGRAL, ale není doloženo jako samostatná parametrizace JORB výstupu; nevystavuji jej zde jako JORB volbu.

## Defaulty a jejich zápis

| Volba | Nativní chování | Definice | Zápis |
| --- | --- | --- | --- |
| ENERGY.GRID | 5 pro E-INTEGRAL | 5 | explicitní |
| ENERGY.NE | 30 pro E-INTEGRAL | 30 | explicitní |
| ENERGY.ImE | obecný default 0,01 Ry, bez JORB override | 0 Ry | explicitní |
| ENERGY.EMIN | -0,2 Ry | -0,2 Ry | explicitní |
| ENERGY.EMAX | EF při nepřítomnosti | None | vynechán |
| MODE | běžně plná relativita při nepřítomnosti | None, FREL | vynechán |
| TASK | JORB příznak | JORB | explicitní |

GRID=5 a NE=30 přímo odpovídají init_mod_energy.f:89–92. ImE=0 je explicitní volba odpovídající společné integrační definici ENERGY, liší se od obecného nativního EIMAG=0,01. Proto nejde o chybný vynechaný default. Nativní 0,01 není bez konvergenčního srovnání jednoznačně lepší; pro komplexní oblouk záleží na samotné kontuře. Vynechané EMAX správně přebírá EF z potenciálu.

Jiný GRID vyvolá upozornění na potřebný obsazený integrační interval, nikoliv tvrdý zákaz všech alternativních kontur. NE má být kladné. SPLITSS nebo dvě položky GRID/NE jsou odmítnuty, protože E-INTEGRAL nastaví jednu cestu. ORBCURDNS sice obsahuje větve SPLITSS, ale samy nezajišťují dokončený dvoukonturový workflow v této větvi GEN.

## Omezení a nálezy zdroje

1. **MPI nedokončí finální JORB výstup.** Po energetické smyčce GEN:1355–1368 pokračuje v MPI jen pro vyjmenované tasky či CALC_JXC/CALC_DMI. Interní E-INTEGRAL s pouhým CALC_JORB podmínku nesplní, a zastaví se před finálním voláním na řádku1588. ORBCURDNS navíc sama neredukuje JORB_LMSCT mezi procesy. Definice má mpi=False. Přítomnost společné MPI energetické smyčky proto neznamená funkční paralelní JORB.
2. **Relativistické předpoklady.** ORBCURDNS používá WAVFUN_READ_REL a relativistické indexy/couplings; nemá vlastní IREL větev pro scalar-relativistic řešení. Definice konzervativně odmítá NREL/SREL/SP-SREL. Nejde o doložený samostatný nativní check JORB, ale o omezení vystaveného workflow podle použité implementace.
3. **FULLPOT je zamýšlený formalismus, nikoliv doložený zákaz ASA.** Hlavička ORBCURDNS uvádí full-potential formalism. CHECK_INPUT nemá JORB zákaz ASA či FULLPOT. Definice proto nevymýšlí zákaz podle pouhého komentáře. Numerická správnost pro ASA ani FULLPOT nebyla ověřena; potenciál musí poskytovat potřebnou angular reprezentaci.
4. **Angular cutoff má vnitřní podmínku.** AME_NAB_RLMC zastaví NKM>UBOUND(CGC,1) a NLAMAX<2*NLMAX (orbcurdns.f:462–466). To jsou odvozené interní rozměry; samotné SITES.NL nestačí k jejich ověření. Nevytvářím neexistující TASK.NLAMAX.
5. **Automatický diagnostický výstup fort.*.** Pro liché l a abs(m)=1 rutina zapisuje na unit 1000+LMR bez explicitního OPEN s názvem (řádky 374–386). Fortran runtime typicky vytvoří fort.<unit>. Tyto soubory nejsou rozlišeny podle DATASET ani typu; při více typech se do stejného unit přidávají data. Popisují diagnostické harmonické, nejsou jednoznačným hlavním výstupem.
6. **Finální korekce není vstupní konstanta JORB.** WE při EFCORRECT=true je SHFTEF; při běžné smyčce WEINP. RENORMALIZE případně násobí uložené JORB_LMSCT hodnotou SCLNOS. GEN inicializuje SHFTEF=0 a RENORMALIZE=false, případně je upraví dalšími větvemi (např. Lloyd). Definice nevystavuje SHFTEF/SCLNOS jako nedoložené TASK volby.
7. **Výstupní NE není počet řádků hustoty.** WRHEAD dostane NETAB(1), ale data jsou zapsána po typech, LM a radiálních bodech, se třemi kartézskými komponentami. Parser nemůže považovat hlavičkové NE za počet datových řádků.

Hlavní soubory jsou `{dataset}_JORB.dat`, `{dataset}_JORB_UP.dat` a `{dataset}_JORB_DN.dat` (nativní `DATSET0`). Celkový výstup je součet dvou spinových složek. Pro každý typ je uveden IRTOP, NLAMAX-1 a NLMAMAX, následované LM bloky radiálních hodnot. Nejde o jednotlivé atomové proudy bez harmonického rozvoje.

DMFT/LDA+U nejsou v definici bez přímého zdrojového zákazu odmítnuty; ORBCURDNS rozlišuje odlišné left/right solutions. Konzistence těchto kombinací nebyla numericky ověřena. Definice vstupu nemění nativní Fortran ani nedokládá specializovaný parser proudové hustoty.

## Nativní souborový kontrakt

`orbcurdns.f:346–352` vytváří `{dataset}_JORB.dat`, `{dataset}_JORB_UP.dat` a `{dataset}_JORB_DN.dat` přes nativní `DATSET0`, tedy bez ADSI připojeného v `kkrmain.f:239–245`. Diagnostické zápisy na unit `1000+LMR` nemají pojmenovaný `OPEN` (`orbcurdns.f:374–386`) a nejsou kanonickými výsledky. `gen.f:1357–1367` ukončí nepodporovaný MPI běh před finálním `ORBCURDNS`.

## Rozsah ověření

Regresní testy definice existují; během této revize nebyl spuštěn nový pytest ani numerický SPR-KKR.
