# APS — audit současné definice

APS má na `tasks-workflow` vlastní definici `aps.py`; závěry o fyzice a parametrech jsou vztažené k dodanému zdroji SPR-KKR 9.7.

## Vstup a defaulty

APS používá kkrgen a NRAPS. initvar.f:281–287 nastaví INITELOOP=false, RUNELOOP=false, USELPLS1=true a RDDOS=true. kkrmain.f:264 sestaví DOSFIL=DATASET.dos a při RDDOS otevře tento soubor. NRAPS očekává starý spinově rozlišený formát DOS s konkrétní strukturou hlavičky. Jméno libovolného DOS souboru tato větev nečte z TASK.

| Parametr | Nativní chování | Definice | Zápis |
| --- | --- | --- | --- |
| TASK.IT | 1 | 1 | default vynechán |
| TASK.CL | 2P | 2P | default vynechán |
| TASK.CALCLS | false | false | vynechán |
| TASK.OUTGNU | false | false | vynechán |
| TASK.GAMMA | lokálně -0,1; nepřítomnost vypne broadening | None | vynechán |
| MODE | IREL=min(1,IREL), při běžném původním IREL=3 tedy 1 | SREL | explicitně zapsán |
| ENERGY.GRID/NE/ImE | task nevolí vlastní mesh, energie se později čtou z DOS | 1 / 300 / 0 Ry | explicitně zapsány |
| ENERGY.EMIN/EMAX | EMIN se pro APS nastaví na EF; tabulka se potom přepíše z DOS | nejsou v této definici | vynechány |

GAMMA=None vyjadřuje nepřítomnost a vypnutý broadening, nikoliv fyzikální šířku -0,1. Vynechané defaulty IT, CL a přepínačů odpovídají zdroji. SREL explicitně vyjadřuje výsledek nativního přepnutí, NREL zůstává povoleno. Při FREL nebo SP-SREL definice varuje, že zdroj režim změní.

GRID=1, NE=300 a ImE=0 jsou explicitní volby definice, nikoliv tvrzení o nativních APS defaultech. NE určuje alokovanou kapacitu NEMAX; skutečný počet energií a ETAB se přebírají z DOS. Lepší než pevný default 300 by bylo převzít kapacitu z hlavičky konkrétního souboru, případně včetně rezervy broadeningu. Tento patch automatické čtení DOS metadat nezavádí.

GAMMA je zdrojem vypisováno jako eV a převáděno dělením RY_EV. Jeho skutečné použití není konstantní Lorentzova šířka: v NRAPS násobí absolutní reálnou energii, viz níže. CALCLS spouští další (l,s)-rozlišené spektrální výstupy, OUTGNU další gnuplot data.

## Nálezy ve zdroji

1. **Nechráněná kapacita DOS.** NRAPS nejprve alokuje podle NEMAX z ENERGY.NE. Později přečte NE z DOS hlavičky a zapisuje ETAB(IE,1) bez kontroly NE<=NEMAX (nraps.f:835–865). I bez broadeningu musí ENERGY.NE stačit celému souboru. IT, NT a NL rovněž musejí odpovídat potenciálu a struktuře DOS; pouze vstupní definice jejich shodu neověří.
2. **Nejméně dvě energie a rovnoměrný rostoucí mesh.** Rutina bez podmínky čte ETAB(2)-ETAB(1) a dále pracuje s jediným krokem DE. DOS s jednou energií není použitelný; neuniformní nebo sestupný mesh není adekvátně ověřován. Definice odmítá ENERGY.NE<2, ale skutečný počet a mesh ze souboru musí ověřit čtečka DOS.
3. **Broadening může překročit více polí.** Pro GAMMA>0 se přidává 5 eV pod EF. NADD=NINT((5/RY_EV)/DELE)+1; zdroj kontroluje pouze NADD+1<=101. Později nastaví NE=NE+NADD a rozšíří ETAB, které má jen NEMAX řádků (nraps.f:947–954). Navíc NEBRD=2*NADD+NE+1 musí vejít do NEMAX+101. Kontrola samotného NADD není dostatečná. Nastavení NE musí zahrnout rezervu; zvolení větší kapacity nenahrazuje kontrolu zdroje.
4. **Šířka je vztažena k absolutní energii.** Pro E>=EFIN se používá WLORTAB=GAMMAV*DREAL(EBRD), nikoliv GAMMAV*(E-EFIN) (nraps.f:929–931). Energie pod EFIN dostanou nulu. Závislost výsledku na absolutním energetickém počátku a jednotky koeficientu zaslouží fyzikální opravu či vysvětlení; parametr nelze poctivě popsat jako běžnou konstantní šířku v eV.
5. **Režim se přepíná natvrdo.** initvar.f:779 provede IREL=min(1,IREL); NRAPS navíc zastaví IREL>1 (nraps.f:270). FULLPOT zakazuje check_input.f:65–68. DMFT a LDA+U zakazuje GEN. Definice odmítá obě korelace; FULLPOT z potenciálu musí zachytit SPR-KKR.
6. **MPI větev APS nedokončí.** APS není v seznamu tasků pokračujících po MPI energetické smyčce v gen.f:1355–1368. Definice má mpi=False.

Definice validuje fyzikální shell l<n a kladný jedničkový IT. Horní mez IT závisí na potenciálu. GAMMA musí být konečné; hodnota <=0 vyvolá upozornění, že broadening vypne. Bez GAMMA nejsou zaváděna omezení odvozená od broadeningu.
