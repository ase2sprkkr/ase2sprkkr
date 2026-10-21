# DMI — audit současné definice

DMI je selektor sdílené definice `jxc.py` pro magnetický nativní E-INTEGRAL. V integrované implementaci neexistuje samostatný `dmi.py`; omezení níže vycházejí ze SPR-KKR 9.7.

## Implementace

Sdílené schéma `jxc.py` poskytuje selektory `JXC` a `DMI`. Volba DMI aktivuje odpovídající nativní větev.

## Výstupní názvy podle SPR-KKR

Nativní `initvar.f:536–545` aktivuje `CALC_DMI` a interní `TASK=E-INTEGRAL`. `dmivector.f:416–418,605–606` používá `TRIM(DATSET)//'_DMIVEC_Dij.dat'` a odpovídající grafy `DMIVEC_Dij_[xyz]`. V `kkrmain.f:239–245` se k `DATSET` připojuje případné `ADSI`, takže nativní jména souborů ADSI obsahují. Současná implementace `jxc.py::output_files()` vytváří jména pouze z `CONTROL.DATASET()`; s neprázdným ADSI je proto manifest nesprávný. To nemění správnost výběru nativní větve DMI.

## Skutečné zapojení

- `initvar.f:534–553`: keyword DMI nastaví CALC_DMI a TASK=E-INTEGRAL;
  QREPONLY omezuje centrální místa. DMI není samostatná hodnota interního TASK.
- `gen.f:355–393`: periodická větev připraví off-diagonal TAUIJ a vyžaduje
  IBZINT=2 (POINTS). Nejde o CLUSTER ani WEYL výpočet.
- `gen.f:975–1000`: magnetický FREL volá DMIVECTOR, pokud není MULTISPIN.
  SP-SREL místo toho vstupuje do větve skalárních J_ij a nevytváří DMI vektory!
- `dmivector.f:70–76`: NONMAG, KMROT jiné než 0, NEPATH jiné než 1,
  GRID jiné než 5, IREL<=2 a NKM jiné než NKMMAX vedou k STOP.
- `gen.f:1193–1252` a `dmivector.f:108–129`: dokončení energetické MPI smyčky,
  redukce kompletně inicializovaného DMI_IJINT a výstup na ranku 0 jsou implementovány.
  Nastavení `mpi=True` má zdrojové opodstatnění; skutečný MPI běh nebyl proveden.

Definice odmítá NREL/SREL/SP-SREL, NONMAG, jiné GRID, více kontur,
SPLITSS a jiné BZINT. KMROT a velikosti základu závisejí na potentialu;
Python vstup bez něj jejich splnění nezaručí.

## Vstupy

TASK: DMI, CLURAD, NSHLCLU, SELECT, EXCLUDE, QREPONLY.
`init_mod_tauij_star.f:123–162` čte seznamy elementů (nejvýše 100) a poloměr/skořápky.
SELECT a EXCLUDE jsou vzájemně výlučné. Vybírají místa podle elementů a obsazení,
ne jednotlivé typy výsledné matice. NSHLCLU má přednost **už svou přítomností**, i pro nulu.
`init_mod_tauij_star.f:591–635` popisuje rozdílnou logiku SELECT a EXCLUDE při smíšeném obsazení.

MULTISPIN spouští jinou rutinu XCPLDMI_FOURSPIN_3SITE; není součástí této definice.
Embedded-cluster větev volá CLUXCPLTENSOR, není ekvivalentem periodického DMIVECTOR.
Tepelné fluktuace jsou v nativní rutině, ale jejich samostatné vstupy a správnost
nebyly tímto auditem ověřeny. Tyto větve nejsou označeny za vypnuté.

## Defaults a zápis

| Parametr | Definice | Nativní chování bez parametru | Posouzení |
|---|---|---|---|
| GRID | 5, zapisuje se | 5 pro E-INTEGRAL | stejné |
| NE | 30, zapisuje se | 30 pro E-INTEGRAL | nativní default převzat místo společného ASE defaultu 32 |
| ImE | 0.01 Ry, zapisuje se | 0.01 Ry | stejné; `mod_energy.f:43` |
| EMIN | -0.2 Ry, zapisuje se | -0.2 Ry | stejné |
| EMAX | nevypisuje se, člen není vystaven | EFERMI po převodu sentinel hodnoty | stejné účinné chování; `init_mod_energy.f:182`; definice nenabízí přímou změnu horní meze |
| CLURAD | 2.2 ALAT, zapisuje se | žádný použitelný default; bez NSHLCLU i CLURAD je STOP | výslovná volba převzatá z JXC; konvergenci poloměru musí ověřit uživatel |
| NSHLCLU | nezadán, nevypisuje se | interně 0, ale absence přepíná na čtení CLURAD | přítomnost a nepřítomnost nejsou totožné |
| SELECT / EXCLUDE | nezadány, nevypisují se | prázdné seznamy | stejné |
| QREPONLY | false, nevypisuje se | false | stejné; `mod_sites.f:73` |

Defaults, které se nevypisují, nemaskují jiný nativní výběr.

## Významná ALE a výsledky

1. **Dz není spočteno stejným postupem jako Dx a Dy.**
   `dmivector.f:321–395` používá pro příčné složky váhu (E-EF), pro xy/yx
   členy jiný výraz. `dmivector.f:483–499` navíc pro Dz přidává faktor 1/2.
   Komentář výslovně říká, že Dz není definováno stejným způsobem pro M||z.
   To není doložená programátorská chyba, ale důležitá fyzikální interpretace:
   tři sloupce nejsou automaticky rovnocenné složky získané jednou procedurou.
2. **QREPONLY / SELECT / EXCLUDE mění mikromagnetický součet.**
   Součet běží pouze přes vybrané páry; `dmivector.f:520–543` nemá explicitní
   doplnění násobností reprezentantů. Nativní výstup sám upozorňuje na přepočet
   výsledku na jednotkovou buňku (`99014`). Omezený součet nepovažovat bez kontroly
   za plný součet celé buňky. Totéž platí pro konečný CLURAD.
3. **Symetrické ekvivalenty vyžadují numerickou kontrolu.**
   BAR_M_T se vytváří pro typ z IQAT(1,IT); na páru se používají typové matice.
   V této párové smyčce není explicitní rotace na rámec ekvivalentního sousedního místa.
   Nelze ze samotného čtení prohlásit výsledek obecně chybným, protože část
   transformací může být zahrnuta ve vstupních TAUIJ. Audit tedy zaznamenává
   neověřenou oblast, nikoli jistou chybu analogickou dříve auditovanému JXC.

Výstup: DATASET_DMIVEC_Dij.dat a grafy jednotlivých složek; párové složky v meV,
mikromagnetické součty tisknuté jako meV*A/f.u. Indexování DMI_ITJT v zápisu i
čtení je konzistentní; není zde doložena chyba z BSDOSIJ. On-site členy se ve
výstupu vynechávají při vzdálenosti menší než 1e-6 ALAT.
