# BSF — audit společné definice BSF/BSFEK/BSFKK

Stav: **smysluplná úloha s dokumentovanými omezeními defaultů a geometrie**. Tento audit popisuje aktuální stav větve `tasks-workflow` a dodaný nativní zdroj SPR-KKR 9.7 (`PUB9.7_260420_18h13`). Pokrývá režimy E-K a K-K realizované společným BLOCHSF kódem; nejde o numerickou certifikaci všech kombinací KPATH/Bravais.

## Aktuální kontrakt ASE2SPRKKR

`BSF`, `BSFEK` a `BSFKK` jsou veřejná jména jedné nativní úlohy: `TASK=BSF`, `CONTROL.ADSI=BLOCHSF`. Aliasy proto sdílejí jeden vstupní a výstupní kontrakt a nemají vlastní odlišná jména souborů.

- `NK1`/`NK2` vyžadují alespoň 2 body. Nativní BLOCHSF v K-K režimu dělí `NK1-1` a `NK2-1`.
- `NK` v E-K režimu vyžaduje alespoň 2 body; `NE` alespoň 1.
- `GRID` závisí na režimu: E-K má výchozí hodnotu 3, K-K 0. Tím se zabrání neplatné jednobodové integraci s `GRID=3`, `NE=1`.
- K-K validace odmítá ne-konečné vektory a nulové nebo kolineární roviny vzhledem ke `KA`.
- `K1` a `K2` jsou koncové body, ne směrové vektory. Nativní kód od nich odečte `KA` a podle potřeby ortogonalizuje druhý relativní směr.

## Souborový kontrakt

TAU adresář používaný BLOCHSF je **volitelná znovupoužitelná cache**, nikoli povinný vstup: pokud cache chybí, nativní kód může TAU data znovu spočítat.

Persistentní výsledné artefakty jsou:

- `{dataset}_{adsi}.bsf`
- `{dataset}_{adsi}_spol.bsf`

`{dataset}_{adsi}` je nativní `DATSET=DATASET_ADSI`. Všechny aliasy mají výchozí `ADSI=BLOCHSF`, a tedy standardně vytvářejí `{dataset}_BLOCHSF.bsf` a `{dataset}_BLOCHSF_spol.bsf`. Při explicitně změněném ADSI manifest následuje nové nativní jméno.

Do result manifestu nepatří dočasné scratch soubory.

## Energie a defaulty mřížky

### E-K režim

Aktuální explicitní defaulty jsou `NE=200`, `GRID=3`, `NK=300`, `ImE=0.001 Ry`, `EMIN=-0.2 Ry` a `EMAX=-1 Ry`.

`EMAX=-1` je skutečná absolutní energie, **nikoli** sentinel pro Fermiho hladinu. Nativní `init_mod_energy` nahrazuje hodnotu `EF` až pro přibližně `EMAX<-99`. Výchozí interval tedy klesá z `-0.2` na `-1 Ry`; uživatel by měl obvykle zamýšlené energetické okno nastavit explicitně.

Nativní reader má pro `NK` default 800, ASE2SPRKKR explicitně zapisuje 300. Je to vědomý ASE-side default, nikoli tvrzení o konvergenci.

### K-K režim

Aktuální explicitní energie je `EMIN=0.7 Ry`, `EMAX=EMIN`, `NE=1` a `GRID=0`. Nejde automaticky o Fermiho energii. Pro výpočet Fermiho plochy je třeba požadovanou absolutní nebo relativní energii nastavit explicitně.

Nevyplněné `KA` nativní kód doplní nulovým počátkem; ASE2SPRKKR validace používá stejnou interpretaci. `NK1`, `NK2`, `K1` a `K2` nemají vymyšlené defaulty a tam, kde jsou potřeba, je musí dodat uživatel.

## Geometrie a nativní omezení

E-K versus K-K se vybírá podle skutečného počtu energií (`NE>1` versus `NE=1`), nikoli samostatnými nativními dispatch klíči. ASE2SPRKKR validuje rozpor mezi požadovaným veřejným aliasem a výsledným režimem.

V E-K režimu nativní kód rozděluje body podle délek segmentů, nejméně dva body na segment. Degenerované explicitní segmenty nebo nulová celková délka cesty zůstávají nativní/geometrické okrajové případy vyžadující smysluplný vstup.

U 2D systémů nativní K-K větev před odečtením `KA` nuluje z-složky `K1`/`K2`. Rovina platná pro zadané 3D vektory proto může po nativní 2D projekci degenerovat. Samotná input definice to z vektorů spolehlivě nepozná; 2D výpočty vyžadují geometricky konzistentní rovinné zadání.

BLOCHSF může znovu použít uložená TAU data. ASE2SPRKKR tuto možnost modeluje, ale netvrdí fyzikální konzistenci libovolné externě dodané cache.

Pokročilé nativní přepínače jako `BSFSPN`, `BSFORB`, `BSFETD`, `BSFALF` a `EPHOT` jsou mimo rozsah tohoto auditu.

## Rozsah ověření

Aktuální regresní pokrytí kontroluje sdílení aliasů, chování režimů E-K/K-K, K-K `GRID=0`, shodu `EMIN`/`EMAX` v jednobodovém energetickém režimu, round-trip/pickle, minimální počty bodů, konečnou a nedegenerovanou geometrii a společný BLOCHSF result manifest pro všechna tři veřejná jména tasku.

Tento audit netvrdí provedení numerického referenčního běhu SPR-KKR.

Nativní zdroje: `blochsf.f` (včetně K-K geometrie a dělení `NK1-1`/`NK2-1`), `init_mod_energy.f` a `epath.f`.
