# SCF — audit současné implementace

**Stav:** Podporované algoritmy BROYDEN2, ANDERSON, TCHEBY; TCHEBY není pro FULLPOT. Audit vychází ze současné implementace na větvi `tasks-workflow` a nativního zdroje SPR-KKR 9.7 (`PUB9.7_260420_18h13`). Hodnocení není numerickou certifikací.

## Rozhraní, validace a nativní omezení

1. ASE nabízelo ALG=PULAY, ale scf.f:427–437 přijímá pouze BROYDEN2, ANDERSON a TCHEBY. PULAY vede do SCFALG not found. Tasková kopie ALG nyní nabízí jen tyto tři názvy; PULAY není automaticky mapováno na ANDERSON, protože takový alias native reader nemá.
2. TCHEBY s FULLPOT native explicitně odmítá; přidána validace. FULLPOT může být aktivován i potenciálem, což input-only validace sama nezjistí.
3. NITER a ISTBRY jsou kladná celá čísla. ITDEPT omezen na 3..40: SCFBROYDEN přepisuje hodnoty pod 3 na 40 a hodnoty nad 40 omezuje na kapacitu. Definice nyní nepředstírá, že tyto požadavky budou provedeny přesně.
4. MIX, MIXOP a TOL musí být konečné kladné hodnoty; EFGUESS konečná energie. Horní hranice MIX nebyla svévolně stanovena na 1, protože overrelaxace není pouhým typovým omylem.
5. NE je kladné; jednobodový GRID=2/3 odmítnut kvůli nativním integračním vadám. Ostatní pokročilé grid kombinace nejsou tímto patchem kompletně kontrolovány.
6. Popis EFGUESS upřesněn: jde o počáteční odhad EF v Ry a vynechání počátečního hledání při tvorbě startovního potenciálu. Nezakazuje všechny pozdější aktualizace EF v SCF.

Opravy se provádějí přes taskové kopie SCF a ENERGY, aby se nezměnily další definice používající společné sekce. Společná sekce ALG může v jiných taskových definicích stále obsahovat PULAY; to není potvrzení jeho nativní podpory.

## Výchozí hodnoty a význam

Běžné hodnoty NITER=200, MIX=0.2, TOL=1e−5, ISTBRY=1 a ITDEPT=40 se vypisují a odpovídají inicializaci ve scf.f:360–371. BROYDEN2 odpovídá MOD_SCF defaultu a je default taskové definice. VXC volby, jejich default a metadata nebyly měněny; MOD_SCF inicializuje VWN, ale konkrétní potenciál a reader mohou vstupní stav dále ovlivnit. Energetický default GRID=5/NE=30 odpovídá SCF větvi INIT_MOD_ENERGY; EMIN=−0.2 je vypsaná ASE volba, EMAX se nevypisuje a ponechává native EF.

**MIXOP je důležitý nevypsaný případ:** scf.f nejprve nastaví SCFMIX=0.2 a SCFMIXOP=SCFMIX, až potom čte MIX a MIXOP. Při MIX=0.5 a nezadaném MIXOP tedy nelze tvrdit, že MIXOP bude 0.5; pro běžný neexportní postup zůstává 0.2. ASE nevnucuje jiný skrytý default, MIXOP je nezadaný a nevypisuje se. Popis tuto závislost nyní výslovně upozorňuje. Doplnit explicitní MIXOP je lepší, pokud chce uživatel stejnou hodnotu jako změněný MIX; automatická změna defaultu by měnila native chování.

EFGUESS se bez zadání nevypisuje; native rozhoduje podle přítomnosti klíče, nikoli podle libovolné náhradní energie. FULLPOT=False se nevypisuje; není příslibem, že se full-potential potenciál změní na ASA. TAU/CPA/STRCONST/SITES/MODE jsou společné sekce, zde bez zásahů.

## Další omezení a větve

- Native EXPORT_HARTREE_POTENTIAL nastavuje NITER=1, MIX=1 a přeskočí čtení VXC/NITER/MIX. Tabulka běžných defaults není přenosná na tuto větev.
- ANDERSON mapuje IMIX=7, BROYDEN2 IMIX=4. TCHEBY využívá samostatný postup; název algoritmu není synonymem pro stejný numerický mix.
- Native nabízí další řízení jako SCFRESTART, MIXRHO, SIM, MAGCOUPL, SEMICORE a 2D layer charge parametry. Nejsou zde doplněny a jejich celé workflow není tímto auditem pokryto.
- mpi=True zachováno. SCF má paralelní řešení a počáteční search; numerická shoda MPI/serial v tomto auditu nebyla testována.

## Implementační kontrakt
- Autoritou pro soubory je aktuální `TaskDefinition` a případné dynamické metody schématu. Výstupy a cache se nesměšují s povinnými vstupy.
- Podmíněné deklarace souborů: `cache_files`, `output_files` v definici.
- Zdroj: `src/ase2sprkkr/input_parameters/definitions/scf.py`.

## Ověření

Testy definic a deklarací souborů jsou v `src/ase2sprkkr/input_parameters/test/`. Tento dokument netvrdí nové spuštění úplné sady ani numerický výpočet SPR-KKR.
