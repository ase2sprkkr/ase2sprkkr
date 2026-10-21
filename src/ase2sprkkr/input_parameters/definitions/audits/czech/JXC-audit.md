# JXC — audit současného tasku

`jxc.py` je sdílená definice JXC/DMI. Nativní `initvar.f:535–545` volí magnetický `E-INTEGRAL`; `gen.f:930–996,1193–1255` rozlišuje SP-SREL skalární výměnu a FREL výměnný tenzor. Následující omezení platí pro podporovaný periodický postup, nikoli jako numerická certifikace.

## Rozsah a opravy

JXC je již podporován. INITVAR čte klíčové slovo JXC, nastaví CALC_JXC a interní task E-INTEGRAL. GEN:930–990 rozlišuje SP-SREL (IREL=2) s XCPLJ0/XCPLJIJ a FREL (IREL=3) s XCPLTENSOR. NREL/SREL nejsou stejné jako SP-SREL; původní warning_condition ostatně zakazovala pouze první dva režimy. Nové testy výslovně ověřují, že SP-SREL zůstává přijat.

Integrovaná definice zpřístupňuje nativně čtené TASK.NSHLCLU, SELECT, EXCLUDE a QREPONLY. INIT_MOD_TAUIJ_STAR:125–161 čte seznamy chemických symbolů, zakazuje jejich současné použití a při přítomném NSHLCLU dává shell count přednost před CLURAD. QREPONLY čte INITVAR uvnitř E-INTEGRAL větve; ovlivňuje seznam centrálních site v INIT_MOD_TAUIJ_STAR. CLURAD je v jednotkách ALAT, nikoli angstromů.

Současná validace odmítá NREL/SREL/NONMAG, nekladné NE, více energetických cest, neauditované integrace a současné SELECT/EXCLUDE. Validaci doplňují upozornění pro MODE a NONMAG. Pro standardní periodický výpočet je BZINT=POINTS explicitní. GEN:358–386 pro off-diagonal TAUIJ vyžaduje IBZINT=2.

Rozsahem auditu je běžný periodický JXC s GRID=5. Ve zdroji existují další větve embedded-cluster, occupation-weighted GRID3 (WGT_OCC), MULTISPIN, SLC_JXC a SLC_TRQ. Nemají zde kompletní audit; validátor proto není tvrzením, že SPR-KKR obecně nikdy neumí GRID3 nebo cluster. Rozšířené nativní větve nejsou součástí současného vstupního kontraktu.

## Defaulty a jejich zápis

| Volba | SPR-KKR | Definice |
| --- | --- | --- |
| TASK.CLURAD | bez NSHLCLU se vyžaduje hodnota; obecný default není | explicitně2,2 ALAT, ponecháno |
| NSHLCLU | při nepřítomnosti volí CLURAD | None, vynecháno |
| SELECT/EXCLUDE | prázdné seznamy | None, vynecháno |
| QREPONLY | false při nepřítomnosti | false, vynecháno |
| ENERGY.GRID | 5 pro E-INTEGRAL | explicitně5, ponecháno |
| ENERGY.NE | 30 pro E-INTEGRAL | explicitně 30 |
| ENERGY.ImE | obecně 0,01 Ry | explicitně 0,01 Ry |
| ENERGY.EMIN | -0,2 Ry | explicitně -0,2 Ry, ponecháno |
| ENERGY.EMAX | obecná náhrada sentinelu za EF | vynecháno, přebírá se EF |
| MODE | běžně FREL | None, vynecháno; SP-SREL se zapisuje při zadání |
| TAU.BZINT | obecná volba závisí na workflow | explicitněPOINTS pro doložený off-diagonal TAU tok |

Integrovaná definice používá NE=30 a ImE=0,01 Ry, což odpovídá těmto nativním výchozím hodnotám. Pro GRID=5 je nadále nutné ověřit fyzikální konvergenci; numerické srovnání nebylo provedeno. Rádius2,2 je pohodlný start, nikoli záruka konvergovaného součtu J_ij. Konvergence NE, spodní meze, k-mesh a clusteru zůstává na numerickém ověření.

## Závažné ALE ve SPR-KKR

1. **Vypnutá rotace podle symetrie ve FREL:** XCPLTENSOR:312 má JQ.EQ.IQREPQ(JQ) .OR. 1.EQ.1. Druhá větev s rotací ME_ROT_MREG tedy nikdy neproběhne. Pro nereprezentativní sousední site může být použitá matice ve špatném rámci. QREPONLY omezuje centrální site, neřeší automaticky všechny sousedy. Toto neoznačuji za rozbití každého jednoatomového kolineárního běhu; jde ale o důležitou chybu obecného tensorového workflow.
2. **FREL má potenciálové podmínky:** NONMAG=false, KMROT=0, NKM=NKMMAX a reprezentativní IQAT(1,IT) jsou přímo kontrolovány XCPLTENSOR:74–80,155. Python zná explicitní NONMAG, ostatní podmínky vyžadují skutečný potenciál. Nelze zadáním QREPONLY certifikovat všechny jeho vlastnosti.
3. **Více cest není zpracováno obecně:** rutiny používají ETAB/WETAB první cesty a mají NEPATH=1 guard. Pozdější vynucení E-INTEGRAL tuto podmínku nezmění na podporu split contouru.
4. **Clusterový SP-SREL není běžný periodický tok:** CLUXCPLJIJ čte scratch unit19, ale TAU_DRIVE běžně volá TAU_REAL_SPACE s KTAUIJ=0, který její export vypne. Bez auditu alternativního připravení unit19 není správné slíbit tento workflow jen přepnutím BZINT=CLUSTER. FREL rozlišuje embedded-cluster zvlášť.
5. **Vyšší spinové a spin-lattice větve mají vlastní interface:** GEN načítá MULTISPIN/SLC_JXC/SLC_TRQ a volí jiné rutiny. Současné schéma tyto větve nevystavuje jako plně podporované tasky.

Na rozdíl od BSDOSIJ používá XCPLTENSOR pro JC správný offset BAR_M_JT_JQ(1,1,JC) a výstupní J_ITJT má konzistentní párové indexování. Chyby z BSDOSIJ nelze mechanicky připsat JXC jen kvůli podobnému zdroji.

## MPI a výstup

mpi=True je zachováno a podložené: J_IJINT je při prvním volání celé nulováno, MPI_SUM_C provede redukci a GEN:1193–1241 má dodatečné dokončovací volání. SP-SREL má vlastní XCPLJIJ redukci. Numerická shoda serial/MPI nebyla měřena.

Tenzorový výstup FREL je v meV (DIMAG(J_IJINT)/pi*RY_EV*1000); SP-SREL používá skalární Jij a samostatnou rodinu souborů. Numerická správnost parseru ani konvence Hamiltoniánu přebíraná jiným softwarem nebyla ověřena.

## Současný souborový kontrakt

Sdílený manifest `jxc.py` deklaruje pro SP-SREL JXC `{dataset}_J_ij.dat` a `{dataset}_Jij_*.agr`, pro FREL JXC `{dataset}_XCPLTEN_Jij.dat`, `{dataset}_XCPLTEN_Dij.dat` a odpovídající grafy. Při zapnutém `TASK.DMI` přidává soubory `DMIVEC_Dij`. Nativní `initvar.f:535–545` vybírá společný `E-INTEGRAL`, `gen.f:951–996` rozhoduje o výpočetní větvi a `dmivector.f:416–418` potvrzuje prefix `TRIM(DATSET)` pro DMI. Protože `kkrmain.f:239–245` k nativnímu `DATSET` připojuje neprázdné `ADSI`, může současný manifest, který používá pouze `CONTROL.DATASET()`, deklarovat nesprávné názvy. Jde o nesoulad manifestu, ne důkaz chybné numeriky JXC.

## Rozsah ověření

Serializace a validace vstupů mají cílené regresní testy; při této úpravě nebyly nově spuštěny testy ani numerický výpočet SPR-KKR.
