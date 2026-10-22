# TDOBS — SPR-KKR 9.7 / ASE2SPRKKR

Nová definice TDOBS pro kkrgen, serial. Kategorie: funguje, ale má velké ALE — není samostatnou novou analýzou oproti generic DOS/observables.

## Doložené chování
initvar.f:523–530 rozpozná TDOBS, nastaví RDTAUDIR=true, WRTAUDIR=false a zapne initialization/energy loop. kkrgen.f:932–961 ověřuje tau cache: s dostupnou cache původně vypne RUNELOOP/INITELOOP, bez ní nastaví recomputation a zápis tau. Následně kkrgen.f:971–974 pro TDOBS vždy znovu zapne obě smyčky. GEN v energy loop volá CALCDOS (gen.f:869–876) a CALCDOS_WRITE (ř.1151). Žádná další dispatch větev TDOBS ve zdrojích nalezena nebyla.

CALC_VECOBS je odlišný interní přepínač, inicializovaný false v mod_calcmode.f:116; parser TDOBS jej nenastavuje. Název TDOBS proto neprezentujeme jako záruku vektorové nebo časově závislé analýzy. Obecné DOBS arrays samy takovou analýzu nedokládají.

GEN MPI guard (gen.f:1360–1371) TDOBS nepovoluje a končí STOP_REGULAR; mpi=False je nutné. Podpora obecné DOS energetické smyčky není numerické ověření všech případů CPA/FULLPOT či opětovného použití cache.

## Vstupy a defaults
TASK=TDOBS a CONTROL.ADSI=TDOBS se zapisují. TAU.TAUDIR/TAUPATH jsou optional bez Python defaultu: native použije DATASET.tau.dir a aktuální adresář (kkrgen.f:266–273). Cache musí odpovídat potenciálu, energii a konfiguraci; definice její obsah nekontroluje.

GRID3, NE300, ImE.01, EMIN-.2, EMAX1 Ry jsou explicitní spektrální defaults převzaté z auditované DOS definice. Nemusejí odpovídat generic native initial defaults, ale nejsou tiše vynechány. EMAX1 není EF sentinel a zahrnuje neobsazené stavy; užší okno může být vhodnější podle cíle. Neexponujeme RDTAUDIR jako nepravdivě přepsatelnou volbu: jeho true je native chování tasku. Chybějící tau nezpůsobuje automaticky konec; native zvolí recomputation fallback.

Validace NE>=1, GRID2/3 nejméně2 body, konečné bounds, ImE>=0 a warning při GRID3/ImE0. Přidány MODE, CPA a STRCONST stejně jako pro DOS. Definice je nezávislá na dřívějším DOS patchi; nedědí jeho případnou neaplikovanou opravu.

## Ověření
5 nových případů: roundtrip TASK/ADSI a tau paths, MPI/executable metadata a odmítnutí neplatných mesh/bounds. Celá input_parameters sada 113 passed, dvě původní NumPy warnings. Numerický SPR-KKR nebyl spuštěn, cache reuse nebylo numericky testováno.
