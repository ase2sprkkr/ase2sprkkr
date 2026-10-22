# OPTABS — audit SPR-KKR 9.7 a ASE2SPRKKR

OPTABS má aktuální definici `definitions/optabs.py` pro absorpční část optické vodivosti. Nativní SPR-KKR 9.7 obsahuje zdokumentovaná omezení; vystavení tasku v Pythonu není numerickou certifikací.

## Definice

CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES; minimální SCF pro odmítnutí FULLPOT. RSP tuto větev v FULLPOT zakazuje. Auditovaný rozsah je FREL a BZINT=POINTS; native OPT_ABS připouští také IBZINT=4 a 6, které zde nejsou certifikovány. NONMAG je přípustný. Many-body optika není tímto auditem pokryta.

TASK nabízí NEB, NOMEGA, OMEGA1/2, TEMPEL, NKT a Eulerovy FRAMETET/FRAMEPHI/FRAMEGAM. ENERGY má GRID=0, NE=1 jen pro předchozí inicializaci a ImE pro skutečný výpočet. OPT_ABS sám vytváří dvě energetické sítě, nastavuje EMIN=EF−OMEGA2 a energetická rozpětí z ENERGY nepoužívá. NEA a TASK.IGRID se sice čtou, ale následný vlastní mesh jejich význam přepíše; definice je proto nenabízí. OMEGA2 je **v eV**, teploty v K, ImE v Ry, úhly ve stupních.

## Konkrétní chyby a omezení

1. **Neinitializované KBT:** při přítomném TASK a TEMPEL≤1 není lokální KBT přiřazeno. Následně řídí tepelnou extenzi a vstupuje do výpočtu NKT_ACTUAL. Nulová teplota tak v tomto zdroji není bezpečným zadáním. Validace vyžaduje konečné TEMPEL>1 K; skutečný T=0 postup by potřeboval opravu Fortranu.
2. **NOMEGA není skutečný počet frekvencí:** OMEGATAB se nejprve alokuje podle zadaného NOMEGA, poté se NOMEGA přepíše na NINT(OMEGA2/DELTA_E), tedy původní NEB. Zadání NOMEGA<NEB způsobuje zápis za alokované pole. Definice vyžaduje NOMEGA≥NEB a popisuje parametr jako kapacitu. NEB≥3 předchází tichému native clampu.
3. **Párování energií:** finální a počáteční sítě mají krok DELTA_E=OMEGA2/NEB. Přechody se započtou pouze při shodě ERYDA−ERYDB s OMEGATAB v toleranci 1e−5 Ry. Obecný posun OMEGA1 může vyřadit velkou část či všechny přechody. Definice omezuje běžný workflow na OMEGA1=0, kterou native nahradí DELTA_E. Výstupní síť tedy nezačíná nulou, ale OMEGA2/NEB. Speciální komensurabilní posuny nejsou nabízeny.
4. Native nahrazuje téměř nulové ImE hodnotou 0.001 Ry; definice vyžaduje kladné konečné ImE a explicitně zapisuje 0.01 Ry. Šířku je nutné zvolit podle požadovaného spektra.
5. Native obsahuje MPI distribuci energetických bodů a redukce spektrálních polí. Zde nebyly numericky ani v celém řetězci certifikovány, proto mpi=False; to není tvrzení, že native MPI větev nutně nefunguje.
6. Tepelná extenze zvětšuje obě energetické sítě, ale počet frekvencí zůstává původní NEB. Vyšší teplota, příliš jemná síť nebo velký NKT mohou významně zvětšit paměť a čas výpočtu. Konvergují se NEB, k-body i ImE. Tensorový rámec a fyzikální interpretace výsledků vyžadují kontrolu pro konkrétní strukturu.

## Defaults: skutečné versus vypsané

| Parametr | Efektivní native | ASE definice a zápis |
|---|---|---|
| NEB | 30, minimum 3 | 30, vypsáno |
| NOMEGA | počáteční kapacita 100, skutečný počet pak NEB | 100, vypsáno; validace kapacity |
| OMEGA1 | 0 eV → DELTA_E | 0, vypsáno |
| OMEGA2 | 8 eV | 8, vypsáno |
| TEMPEL | reader default 1e6 K při TASK | 300 K, vypsáno |
| NKT | lokálně 6, reader při TASK přepíše na 99 | 6, vypsáno |
| ENERGY GRID/NE | vlastní optická síť vytvořena znovu | 0/1, vypsáno jen pro inicializaci |
| ImE | téměř nula nahrazena 0.001 Ry | 0.01 Ry, vypsáno |
| EMIN/EMAX | EMIN přepočítá OPT_ABS, jeho síť určí TASK | nejsou v definici |
| Eulerovy úhly | 0, rámec změněn jen při přítomném klíči | nezadané, nevypsané |

TEMPEL=300 a NKT=6 jsou explicitní odchylky, nikoli skryté defaults. Jsou vhodnější pro běžný start než native 1e6/99; experimentální teplotu si uživatel musí nastavit. Význam readeru ověřen v scan_input.f: default odlišný od sentinel 9999 přepisuje hodnotu při chybějícím klíči. Nativní počáteční DATA hodnoty Eulerových úhlů souhlasí s jejich nevypsaným stavem.

## Rozsah ověření

Aktuální Python testy kontrolují definici tasku, serializaci a vybrané validační podmínky. Historické počty testů z oddělených patchů nejsou výsledkem současné větve. Numerický referenční běh SPR-KKR nebyl součástí tohoto zdrojového auditu.
