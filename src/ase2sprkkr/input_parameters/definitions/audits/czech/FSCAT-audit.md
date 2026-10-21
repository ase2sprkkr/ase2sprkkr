# FSCAT — audit současné implementace

Integrovaná definice `fscat.py` zpřístupňuje nativní FSCAT v `kkrgen`. Níže uvedené chyby indexování ve Fortranu tím nejsou odstraněny.

## Skutečné rozhraní

`initvar.f:176–177` rozpoznává TASK FSCAT; `gen.f` volá FSCAT v KKRGEN. Rutina `fscat.f` nečte SECTION_SET/SECTION_FIND, tedy nemá task-specific parametry. Nevytvářím ITSEL, NGEO, DIRINI, DIRFIN ani přepínače výstupu, které vstup skutečně neřídí.

Počítají se všechny atomové typy. NGEO=3 je pevně nastavené; počáteční směry jsou [001], [100], [111]/sqrt(3), finální směr je vždy opačný. Jde o tři geometrie zpětného rozptylu, nikoli uživatelsky nastavitelné obecné úhly. Výstup pro každý typ a geometrii obsahuje čtyři spinové složky amplitudy a jejich fáze.

Definice používá CONTROL, TAU, MODE, STRCONST, ENERGY, TASK a SITES; `executable=kkrgen`, `mpi=False`. Vlastní smyčky rutiny nemají rozdělení mezi ranky a vytvářejí stejné názvy výstupů. MPI ani numerický běh nebyl ověřen.

| Volba | Nativní chování | Definice |
|---|---|---|
| ENERGY.GRID | první cesta ETAB | explicitní 1, přímá reálná energetická cesta |
| ENERGY.NE | NETAB(1), bez samostatného task-specific defaultu | explicitní 300, kladný počet |
| ENERGY.ImE | init_mod_energy přepisuje na 0 | default 0; při nenulové volbě upozornění |
| ENERGY.EMIN | -0.2 Ry | explicitní -0.2 Ry, shodné s nativním defaultem |
| ENERGY.EMAX | při vynechání skutečná EF | None, nevypisuje se a deleguje EF SPR-KKR |
| ENERGY.EMINEV/EMAXEV | relativní meze vůči EF v eV | společné párové rozhraní ENERGY |
| ENERGY.SPLITSS | reset na false | false; při zapnutí upozornění |
| SITES.NL | kontextový cutoff | None, nevypisuje se; explicitní hodnoty se kontrolují |
| MODE.MODE | default plná relativita | na jiné režimy upozornění, nikoli domnělý nativní zákaz |

GRID=1 a NE=300 se vypisují: jsou to praktické volby ASE2SPRKKR, nikoli tvrzení o nativních defaults. NETAB začíná na {0,0} (`mod_energy.f:43–44`) a init_mod_energy pro FSCAT nenastavuje zvláštní počet bodů; není vhodné jej ponechat na tomto počátečním stavu. Pro rychlou kontrolu lze NE snížit. EMAX a NL jsou omitted defaults a zachovávají kontext SPR-KKR. Nativní EMAX=EF je zde ponechán; pro nezaplněné stavy musí uživatel zvýšit horní mez, například relativním EMINEV/EMAXEV.

`init_mod_energy.f:315–320` přepisuje ImE a vypíná SPLITSS. FSCAT používá pouze první ETAB cestu. Na přepis a druhou cestu definice upozorňuje. Na rozdíl od PSHIFT nepřevádí samotnou energii před SSITE na její reálnou část: předává celé ETAB (`fscat.f:162–166`). Reset ImE proto nelze vydávat za záruku, že libovolný jiný contour bude mít všechny energie reálné. Default GRID=1 dává přímou reálnou cestu.

Výstupní energetická osa je **E-EF v eV** (`fscat.f:163`), na rozdíl od absolutní osy PSHIFT. Grafy se vytvářejí pomocí XMGRHEAD s GEO1/GEO2/GEO3; tento audit neověřuje specializovaný parser výstupů.

## Omezení a validace

FULLPOT explicitně odmítají `check_input.f:70–72` a příslušná FSCAT větev GEN. Rozpoznává se z potenciálu, který samotná definice vstupu nenačítá; zákaz je zdokumentovaný a vynucuje jej SPR-KKR. DMFT a LDAU zakazuje `gen.f:211–222`, definice je také odmítá.

FSCAT sestavuje relativistické kappa a spinové kanály, ale check_input neobsahuje přímý požadavek IREL=3. Na NREL/SREL/SP-SREL definice upozorňuje jako na neověřené kombinace. Pro IREL<2 se NK/NKM v init_mod_angmom nastavují jinak než pro relativistickou bázi, takže jejich použití v relativistické konstrukci je podezřelé; numerický dopad zde ověřen není.

FSCAT přepisuje NLT i NLQ na globální NL a NKMQ na NKM (`fscat.f:79–85`). Explicitně neuniformní cutoff tedy nezůstává účinný; definice upozorňuje. Pokud globální NL=1, CALC_CHPLM alokuje T(0:0,0:0), ale zapisuje T(1,1) a T(1,0); LEGENDRE navíc zapisuje P(2) do pole délky 1. Explicitní s-only cutoff proto definice odmítá. Při omitted NL bez potenciálu tuto kontrolu nemůže provést.

## Závažné chyby ve zdroji SPR-KKR

1. **Nesprávný index sférických harmonik a čtení mimo pole.** FSCAT používá `LM=L*(L+1)+L+ML+1` (`fscat.f:131`). CALC_CHPLM ukládá m=0 na `L*(L+1)+1` a m složky s posunem ±m (`calc_chplm.f:155–162`). Ve FSCAT je navíc +L. Například při NL=3, l=2, m=2 vzniká LM=11, ale harmoniky mají pouze 9 prvků. Toto není pouze kosmetická chyba: může znehodnotit amplitudy nebo vyvolat bounds error. Navíc je třeba vynechat nulové Clebsch-Gordan koeficienty před přístupem k harmonikám, protože nepřípustné m=l+1 se může objevit i v jejich nulových složkách.
2. **Výpočet fáze před ošetřením nuly.** Po `AA=ABS(BB)` se okamžitě počítá `LOG(BB/AA)` (`fscat.f:207–209`); kontrola AA>1e-8 následuje až potom. Při BB=0 jde o dělení 0/0. Pozdější přiřazení nul do tabulek nezabrání předchozí invalid operaci a PHI se ještě používá pro limity grafu. Správné pořadí vyžaduje vypočítat PHI až uvnitř nenulové větve.
3. **Omezení NL=1 není kontrolováno ve zdroji.** CALC_CHPLM a LEGENDRE zapisují prvky mimo pole pro l_max=0, jak je uvedeno výše. Definice omezuje explicitní input, ale neodstraňuje chybu knihovní rutiny ani chybu harmonického indexu pro NL>=2.
4. **Výstupní text je zavádějící.** STOP_REGULAR uvádí `phase shifts calculated` (`fscat.f:284`), přestože vypočítané výstupy jsou matice amplitud rozptylu. Lokální pomocné FT se počítá, ale dále nepoužívá; není důvod jej inzerovat jako dostupný výstup.

Definice nevydává syntakticky platný input za numericky spolehlivý výpočet. Zvláště první dva problémy vyžadují opravu SPR-KKR; nejsou řešitelné doplněním task-specific parametrů ASE2SPRKKR.

## Výstupní kontrakt ověřený ve SPR-KKR

`fscat.f:87–99,237–245` používá tři pevné geometrie `GEO1`–`GEO3` a volá `XMGRHEAD(DATSET,...)` pro grafové soubory. Python manifest deklaruje `{dataset}_GEO?_*.agr`; glob pokrývá tři nativní rodiny. Protože `kkrmain.f:239–245` k `DATSET` připojuje neprázdné `ADSI`, skutečný prefix výstupů obsahuje ADSI. Manifest vychází pouze z `CONTROL.DATASET()`, takže při neprázdném ADSI nemusí souhlasit. Numerická správnost tím není prokázána.

## Rozsah ověření

Definice vstupu a cílené testy existují v `src/ase2sprkkr/input_parameters/`. Tento audit netvrdí nové spuštění testovací sady ani numerickou validaci SPR-KKR.
