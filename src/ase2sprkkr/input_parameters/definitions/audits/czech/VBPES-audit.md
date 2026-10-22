# VBPES — audit aktuální implementace

**Stav:** Valenční fotoemise přes kkrgen, nikoli KKRSPEC; final-state TAU závisí na USETAUNN/USETFIN. Auditovaná implementace je na `tasks-workflow`, nativní referencí je SPR-KKR 9.7 (`PUB9.7_260420_18h13`). Nejde o numerickou certifikaci.

## Vstup

initvar.f:220–255 rozpoznává VBXPS, VBPES a VBUPS jako úhlově integrovaný VBPES; ARPES přepíná odlišnou větev. Nová definice registruje pouze kanonický VBPES. Aliasové názvy nejsou samostatně registrovány. check_input.f:91–94 vyžaduje IREL=3; definice odmítá NREL/SREL/SP-SREL. FULLPOT zde není zakázán, ale ve vbpes.f:224 násilně nastaví TASK.MODE=2.

TASK podporuje EPHOT, DEFEP0, ME, MODE, EWORK, GEOMETRY, QVEC, QVECPH, TETPH, PHIPH, NSPIN, ImEFIN, USETAUNN a USETFIN. Sekce ENERGY navíc podporuje ERANGE, který nastavuje EMIN=EF-ERANGE/RY_EV. EPHOT, DEFEP0, EWORK a ERANGE jsou v eV; ImEFIN a ENERGY.ImE v Ry. TETPH a PHIPH se zadávají ve stupních.

VBPES při existujícím počátečním TAU čte počet energií z jeho metadat (vbpes.f:282–291); při chybějícím TAU obsahuje `kkrmain.f:932–960` nativní přepočítávací větev, jejíž úspěch závisí na nastavení úlohy. Dostupnost a shodu souborů nelze ověřit pouze serializací vstupu. USETAUNN používá nebo dopočítá TAU konečných stavů; USETFIN jej přebije a vypne. Výpočet konečných stavů je posunut o EPHOT.

## Defaulty: vynechané proti explicitním

| Parametr | Nativní chování | Nová definice | Zápis defaultu |
| --- | --- | --- | --- |
| TASK.EPHOT | 1253,6 eV | 1253,6 | vynechán |
| TASK.DEFEP0 | 0 eV | 0 | vynechán |
| TASK.ME | ADA | ADA | vynechán |
| TASK.MODE | 2; FULLPOT vždy 2 | 2 | vynechán |
| TASK.GEOMETRY | 0 | 0 | vynechán |
| TASK.ImEFIN | Im(ETAB(NE,1)) | None, delegace na zdroj | vynechán |
| TASK.NSPIN | -QHAT po výběru geometrie | None | vynechán |
| TASK.EWORK | proměnná EWORKEV není před voláním explicitně inicializována; setter dostane sentinel 9999 | None | vynechán |
| USETAUNN/USETFIN | false při nepřítomnosti | false | vynechány |
| ENERGY.EMAX | EF při nativní inicializaci | None | vynechán |
| ENERGY.EMIN | -0,2 Ry před čtením ENERGY, případně override ERANGE | -0,2 Ry | explicitní |
| ENERGY.GRID/NE/ImE | bez speciálního taskového nastavení GRID/NE; počty se následně čtou z TAU | 3 / 300 / 0,01 Ry | explicitní |

GRID=3 a NE=300 jsou volbou definice pro běžný fotoemisní mesh, nejsou tvrzením o nativních defaultech VBPES. Jejich zápis je explicitní. Při opětovném použití existujícího TAU musí mesh odpovídat jeho metadatům; pevný počet 300 není automatická detekce obsahu TAU. Lepší budoucí řešení je získat metadata z konkrétního souboru. ImEFIN se naopak nesmí nahradit pevným vynechaným defaultem: je odvozen z poslední počáteční energie.

EWORK je v této rutině pouze použito k výpočtu vypsaného EVAC a nepředává se VBPESINT. Sentinel 9999 není dokumentovaný fyzikální default. Definice jej nevystavuje jako známou hodnotu a nepřebírá sem 4,2 eV z odlišné definice ARPES. Pokud je požadován smysluplný diagnostický EVAC, EWORK musí uživatel explicitně dodat. Neověřená lokální inicializace EWORKEV je nedostatek nativního zdroje, který by měl být opraven tam.

## Nálezy

1. **BIS je rozpoznán, ale nedostane se přes kontrolu tasku.** GEN volá VBPES pro VBPES i BIS (gen.f:1445), ale vbpes.f:174 vyžaduje TASK(3:5)=='PES'. BIS zde končí STOP. Tato definice BIS nevystavuje; jde o chybu nativní větve.
2. **GRV není implementováno.** ME_DRIVE má pro GRV zakomentované volání a STOP_MESSAGE (me_drive.f:107–116). Definice povoluje ADA a NAB, které mají skutečné větve výpočtu. Numerická správnost ani jednoho nebyla testována.
3. **EPHOT má číselné aliasy s překvapivým rozsahem.** NINT(EPHOTEV)==1/2 přepíše hodnotu na 1253,6/1486,6 eV (vbpes.f:199–200). Není to pouze přesná hodnota 1 nebo 2: přepsány budou i sousední energie zaokrouhlené na tato čísla. Definice chování dokumentuje a vstup nepřepisuje.
4. **QVEC a vestavěné geometrie se špatně skládají.** QVEC se normalizuje, ale geometrie 0/1/2 následně nastaví jedinou složku QHAT na 1 bez vynulování ostatních nebo nové normalizace (vbpes.f:208–250). Bez QVEC se nejprve nulují všechny složky a vestavěná geometrie je jednoznačná. Definice při této kombinaci varuje. Nulové a nefinite vektory QVEC, QVECPH, NSPIN odmítá kvůli dělení normou.
5. **Custom geometrie je závislá na úplnosti vstupu.** Pro GEOMETRY mimo 0/1/2 se při přítomném TETPH vyžaduje PHIPH; jinak QVECPH. Zdroj nastaví KGEO=3/4 podle nalezených hodnot. Definice nabízí 3/4, ověřuje potřebné hodnoty a varuje při přebití QVECPH úhly.
6. **DEFEP0 neprovádí deklarovatelný posun spektra.** Ve VBPES se pouze převede a vypíše; nepoužívá se při sestavení ETABFIN ani se nepředává VBPESINT. Totéž platí pro EWORK/EVAC. Tyto parametry nelze zaměnit za funkční energetické posuny v KKRSPEC.
7. **POLCONV míchá jednotky vrácených úhlů.** Ve vektorové větvi vypočítá TETPHOT=ACOS(-QHAT(3)) bez převodu na stupně, ale PHIPHOT převádí na stupně (vbpes_aux.f:282–289). Samotné QHAT se sestaví správně z vektoru; dopad na následné úhlové diagnostiky vyžaduje pozornost.
8. **MPI postprocessing se nespustí.** VBPES není v seznamu tasků povolených po MPI energetické smyčce v gen.f:1355–1368. Definice má mpi=False.
9. **Názvy výstupů mohou kolidovat.** Rutina otevírá pevné arpes.dat, arpes_spol.dat a spin.matrix.xps. Jména TAU konečných stavů obsahují celočíselnou část EPHOT, takže různé necelé energie mohou sdílet stejný název. Shodu existujícího souboru je nutné kontrolovat podle obsahu.

Korelace nejsou v definici bez podkladu zakázány: VBPES obsahuje explicitní DMFT větev, která pro konečný stav nastaví ORBPOLFIN=NONE. FULLPOT, korelace a numerická správnost spekter nebyly ověřeny během skutečného běhu SPR-KKR.

## Aktuální souborový kontrakt

Nativní `DATSET` zahrnuje případné neprázdné `ADSI` (zdroj `kkrmain.f:239–266`), `DATSET0` zůstává původním DATASET. ASE2SPRKKR používá samostatné `{dataset}` a `{adsi}`, které při prázdném ADSI netvoří podtržítko. Počáteční `{dataset}_{adsi}.tau.dir` je volitelná cache; vlastní jméno lze nastavit pomocí `TAUDIR`/`TAUPATH`. Výstupy `arpes.dat`, `arpes_spol.dat` a `spin.matrix.xps` jsou nativně pevná jména, zatímco volitelné TAU a DOS soubory konečných stavů používají nativní DATSET. Zdroj: `definitions/vbpes.py`.

## Ověření

Testy definic jsou v `src/ase2sprkkr/input_parameters/test/`. Nový kompletní pytest ani numerický SPR-KKR běh se zde netvrdí.
