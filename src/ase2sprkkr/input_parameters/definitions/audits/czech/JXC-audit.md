# JXC — audit existující definice

SPR-KKR 9.7 z dodaných zdrojů versus ASE2SPRKKR refactor po fetchi, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Patch 0028-jxc.patch upravuje existující jxc.py a přidává test_jxc.py samostatně proti tomuto základu.

## Rozsah a opravy

JXC je již podporován. INITVAR čte klíčové slovo JXC, nastaví CALC_JXC a interní task E-INTEGRAL. GEN:930–990 rozlišuje SP-SREL (IREL=2) s XCPLJ0/XCPLJIJ a FREL (IREL=3) s XCPLTENSOR. NREL/SREL nejsou stejné jako SP-SREL; původní warning_condition ostatně zakazovala pouze první dva režimy. Nové testy výslovně ověřují, že SP-SREL zůstává přijat.

Doplněny skutečně čtené TASK.NSHLCLU, SELECT, EXCLUDE a QREPONLY. INIT_MOD_TAUIJ_STAR:125–161 čte seznamy chemických symbolů, zakazuje jejich současné použití a při přítomném NSHLCLU dává shell count přednost před CLURAD. QREPONLY čte INITVAR uvnitř E-INTEGRAL větve; ovlivňuje seznam centrálních site v INIT_MOD_TAUIJ_STAR. CLURAD je v jednotkách ALAT, nikoli angstromů.

Nová konečná validace odmítá NREL/SREL/NONMAG, nekladné NE, více energetických cest, neauditované integrace a současné SELECT/EXCLUDE. Původní warning_condition pro režim a NONMAG zůstává; platný vstup nemění. Pro standardní periodický výpočet je BZINT=POINTS explicitní. GEN:358–386 pro off-diagonal TAUIJ vyžaduje IBZINT=2.

Rozsah tohoto patchu je běžný periodický JXC na GRID5. Ve zdroji existují další větve embedded-cluster, occupation-weighted GRID3 (WGT_OCC), MULTISPIN, SLC_JXC a SLC_TRQ. Nemají zde kompletní audit; validátor proto není tvrzením, že SPR-KKR obecně nikdy neumí GRID3 nebo cluster. Obecná předchozí definice je ani nedokládala správnými task-specific parametry a ověřeným tokem. Pro zpřístupnění rozšířených větví je nutný samostatný krok.

## Defaulty a jejich zápis

| Volba | SPR-KKR | Definice |
| --- | --- | --- |
| TASK.CLURAD | bez NSHLCLU se vyžaduje hodnota; obecný default není | explicitně2,2 ALAT, ponecháno |
| NSHLCLU | při nepřítomnosti volí CLURAD | None, vynecháno |
| SELECT/EXCLUDE | prázdné seznamy | None, vynecháno |
| QREPONLY | false při nepřítomnosti | false, vynecháno |
| ENERGY.GRID | 5 pro E-INTEGRAL | explicitně5, ponecháno |
| ENERGY.NE | 30 pro E-INTEGRAL | explicitně32, ponecháno |
| ENERGY.ImE | obecně0,01 Ry | explicitně0 Ry, ponecháno |
| ENERGY.EMIN | -0,2 Ry | explicitně -0,2 Ry, ponecháno |
| ENERGY.EMAX | obecná náhrada sentinelu za EF | vynecháno, přebírá se EF |
| MODE | běžně FREL | None, vynecháno; SP-SREL se zapisuje při zadání |
| TAU.BZINT | obecná volba závisí na workflow | explicitněPOINTS pro doložený off-diagonal TAU tok |

NE32 versus30 je zapisovaný rozdíl, není to chybný omitted default. Nativní30 by o něco zkrátilo integraci; bez konvergenčního srovnání není důvod měnit dosavadních32. GRID5 tvoří komplexní oblouk vlastní geometrie, takže ImE0 versus0,01 není prostý rozdíl jeho výšky. Rádius2,2 je pohodlný start, nikoli záruka konvergovaného součtu J_ij. Konvergence NE, spodní meze, k-mesh a clusteru zůstává na numerickém ověření.

## Závažné ALE ve SPR-KKR

1. **Vypnutá rotace podle symetrie ve FREL:** XCPLTENSOR:312 má JQ.EQ.IQREPQ(JQ) .OR. 1.EQ.1. Druhá větev s rotací ME_ROT_MREG tedy nikdy neproběhne. Pro nereprezentativní sousední site může být použitá matice ve špatném rámci. QREPONLY omezuje centrální site, neřeší automaticky všechny sousedy. Toto neoznačuji za rozbití každého jednoatomového kolineárního běhu; jde ale o důležitou chybu obecného tensorového workflow.
2. **FREL má potenciálové podmínky:** NONMAG=false, KMROT=0, NKM=NKMMAX a reprezentativní IQAT(1,IT) jsou přímo kontrolovány XCPLTENSOR:74–80,155. Python zná explicitní NONMAG, ostatní podmínky vyžadují skutečný potenciál. Nelze zadáním QREPONLY certifikovat všechny jeho vlastnosti.
3. **Více cest není zpracováno obecně:** rutiny používají ETAB/WETAB první cesty a mají NEPATH=1 guard. Pozdější vynucení E-INTEGRAL tuto podmínku nezmění na podporu split contouru.
4. **Clusterový SP-SREL není běžný periodický tok:** CLUXCPLJIJ čte scratch unit19, ale TAU_DRIVE běžně volá TAU_REAL_SPACE s KTAUIJ=0, který její export vypne. Bez auditu alternativního připravení unit19 není správné slíbit tento workflow jen přepnutím BZINT=CLUSTER. FREL rozlišuje embedded-cluster zvlášť.
5. **Vyšší spinové a spin-lattice větve mají vlastní interface:** GEN načítá MULTISPIN/SLC_JXC/SLC_TRQ a volí jiné rutiny. Tento patch neprezentuje tyto keywordy jako plně podporované pouhým přidáním názvu.

Na rozdíl od BSDOSIJ používá XCPLTENSOR pro JC správný offset BAR_M_JT_JQ(1,1,JC) a výstupní J_ITJT má konzistentní párové indexování. Chyby z BSDOSIJ nelze mechanicky připsat JXC jen kvůli podobnému zdroji.

## MPI a výstup

mpi=True je zachováno a podložené: J_IJINT je při prvním volání celé nulováno, MPI_SUM_C provede redukci a GEN:1193–1241 má dodatečné dokončovací volání. SP-SREL má vlastní XCPLJIJ redukci. Numerická shoda serial/MPI nebyla měřena.

FREL zapisuje DATASET_XCPLTEN_Jij.dat, tensor v meV (DIMAG(J_IJINT)/pi*RY_EV*1000), a další výstupy; SP-SREL má jinou strukturu a scalar Jij. Tento audit nemění stávající output parser ani konvenci Hamiltoniánu/faktory převzaté dalším softwarem.

## Testy

14 nových případů: serializace/načtení/pickle, původní defaulty, zachování SP-SREL, nové shell/element selectors/QREPONLY, jejich výlučnost, energetická cesta, neauditované integrace, nekladné NE a zakázané nemagnetické režimy. Celá sada input_parameters: **122 passed**. Numerický SPR-KKR nebyl spuštěn. Patch prošel git diff --check a kontrolou zpětného aplikování.
