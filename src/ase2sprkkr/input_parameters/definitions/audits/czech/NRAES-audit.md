# NRAES — audit SPR-KKR 9.7 a ASE2SPRKKR

NRAES má v ASE2SPRKKR aktuální definici `definitions/nraes.py`. Nativní požadavky na starý formát DOS a omezení konvoluce vycházejí ze SPR-KKR 9.7; existence Python definice jejich správnost nezaručuje.

## Vstup

Nerelativistická Augerova elektronová spektroskopie v kkrgen. initvar.f:290–296 nastaví INITELOOP=false, RUNELOOP=false, RDDOS=true a USELPLS1=true. DOS se čte z DATSET.dos (DATASET_ADSI při neprázdném ADSI), otevřeného na unit 10. Čtečka vyžaduje starý formát spinově rozlišeného DOS; nedochází k automatickému přepočtu DOS z potenciálu.

TASK čte pouze IT, CL, CALCLS a OUTGNU. GAMMA z podobné rutiny NRAPS sem nepatří. Definice umožňuje n=1..5 a s/p/d/f při fyzikálním l<n. Horní mez IT, počet typů a cutoff musejí odpovídat potenciálu i DOS; samotná definice nemá jejich metadata.

## Defaulty a zápis

| Volba | Nativní chování | Definice | Zápis |
| --- | --- | --- | --- |
| TASK.IT | 1 | 1 | default vynechán |
| TASK.CL | 2P | 2P | default vynechán |
| CALCLS / OUTGNU | false při nepřítomnosti | false | vynechány |
| MODE | IREL=min(1,IREL) | SREL | explicitně zapsán |
| ENERGY.GRID / NE / ImE | bez specifického NRAES meshe; energie načteny z DOS | 1 / 300 / 0 Ry | explicitně zapsány |
| EMIN / EMAX | EMIN obecně -0,2 Ry; EMAX inicializováno na EF, poté se ETAB čte z DOS | bounds nejsou vystaveny | vynechány |

Vynechané TASK defaulty odpovídají SPR-KKR. Explicitní SREL vyjadřuje nativní přepnutí běžného relativistického režimu na IREL=1; NREL je také povoleno. Definice varuje při FREL/SP-SREL, protože zdroj je přepíše. FULLPOT zakazuje check_input.f:60–62. DMFT/LDA+U zakazuje GEN a odmítá je i definice.

GRID=1, NE=300 a ImE=0 jsou explicitní volby inicializace, nejsou nativními taskovými defaulty. NE určuje NEMAX a tím kapacitu polí; skutečné energie a počet NE se čtou ze souboru. Lepší než pevné NE=300 by bylo odvození z hlavičky konkrétního DOS. Definice metadata DOS automaticky nečte. ENERGY bounds nemají být prezentovány jako ovládání spektrálního rozsahu NRAES, když rutina tabulku přepíše.

## Chyby a důležité předpoklady zdroje

1. **IAES je jen částečně inicializováno.** nraes.f:99 alokuje IAES(2*NEMAX-1,0:NLMAX,0:2). Před konvolucí se však volá RINIT pouze s počtem NEMAX*(1+NLMAX)*(1+2) (řádek 814), místo (2*NEMAX-1)*(1+NLMAX)*3. Následná smyčka běží přes 2*NE-1 a sčítá do IAES. Část pole může obsahovat neinicializované hodnoty. Kvůli sloupcovému uložení nelze problém popsat jen jako chybějící nulování posledních energií: částečné nulování zasahuje různé kanály. Toto vyžaduje opravu Fortranu; vstupní definice ji nenahrazuje.
2. **Chybí kontrola kapacity čteného DOS.** Před čtením je DOS i ETAB alokováno podle NEMAX. Rutina potom načte NE, NT a NQ z DOS hlavičky a zapisuje řádky ETAB/DOS bez ověření proti kapacitě a potenciálu (řádky 211–246). ENERGY.NE musí stačit souboru; definice ověří jen minimální kapacitu 2.
3. **Mesh musí být rostoucí a rovnoměrný.** DE=real(ETAB(2)-ETAB(1)) je použit bez kontroly počtu energií, následně jako jediný krok celého spektra. Jedna energie nebo neuniformní mesh není adekvátně ošetřen. Interpolace matrix elements používá DEME=(EMAX-EBOT)/(NEME-1); EMAX=EBOT způsobí nulový krok. Tato vlastnost konkrétního DOS není ověřitelná pouze InputParameters.
4. **Omezená tabulka názvů slupek.** SHELL obsahuje pouze K,L,M,N,O (nraes.f:86) a je indexována NCXRAY. Proto n=6 přistupuje mimo tabulku. Definice jej odmítá. Stejná tabulka je v NRAPS: původní původní definice APS chybně připouštěl 6S..6F. V tomto kroku byl APS opraven na n<=5, přidán test odmítnutí 6S a obnoven patch 0014-aps.patch. Původní APS audit před touto opravou popisoval 8 případů; aktualizovaný test_aps.py má 9, všechny prošly.
5. **Energie výstupu je relativní k onsetu.** EONSET=2*EFERMI-E1 a EFIN=(2*EBOT-E1+(IE-1)*DE)-EONSET (řádky 995–998). E1 se proto z osy odečtením vyruší. Není to absolutní kinetická energie Augerova elektronu; posun na absolutní osu vyžaduje odpovídající onset.
6. **Metadata výstupu neodpovídají počtu bodů spektra.** WRHEAD dostane NE načtené z DOS (řádek 198), ale hlavní spektrum se zapisuje s 2*NE-1 body (řádky 996–999). Parser nemůže bez další kontroly předpokládat, že hlavičkové NE je počet spektrálních řádků.
7. **MPI větev nedokončí NRAES.** NRAES není v seznamu povolených tasků po MPI energetické smyčce GEN. Definice má mpi=False.

CALCLS zapíná další (l,s)-rozlišený výstup. OUTGNU zapíná další gnuplot data. Specializovaný parser NRAES nebyl vytvořen a numerická správnost nebyla ověřena běžícím SPR-KKR.

## Rozsah ověření

Aktuální Python testy kontrolují definici tasku, serializaci a vybrané validační podmínky. Historické počty testů z oddělených patchů nejsou výsledkem současné větve. Numerický referenční běh SPR-KKR nebyl součástí tohoto zdrojového auditu.
