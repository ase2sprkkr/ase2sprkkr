# FMAG — audit SPR-KKR 9.7 a nová definice

Základ: dodané zdroje SPR-KKR 9.7 a aktuálně fetchovaný ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. FMAG zde nemá vlastní definici. Samostatný patch 0017-fmag.patch přidává fmag.py a test_fmag.py, bez dřívějších patchů.

## Co skutečně spouští vstup FMAG

TASK FMAG aktivuje CALC_FMAG a initvar.f:534–550 přejmenuje interní TASK na E-INTEGRAL. Nastaví NOWRDOS=true, RDTAUDIR=false, WRTAUDIR=false, ICORE=0 a IHYPER=0. GEN pro E-INTEGRAL zapne zápis obou vlnových funkcí, CALCINT a GETIRRSOL, nastaví ICORE=1 a NEPATH=1 (gen.f:296–303). Číselné nastavení ICORE v initvar tedy není finální.

GEN volá DENSITY při každé energii a po dokončení integrace (gen.f:1078–1084,1180–1186). DENSITY akumuluje charge, spin a orbital density, vytvoří dens a fmag výstupy a ukončí program pomocí STOP_REGULAR. FMAG nemá vlastní čtečku parametrů hybnostní sítě ani specifické TASK volby pro form factor. Nová definice nevymýšlí QMAX nebo NQ. Obecný přepínač QREPONLY je sice čten ve společné větvi E-INTEGRAL, ale DENSITY jej sama nepoužívá; v této definici není vystaven jako FMAG parametr.

## Defaulty a zápis

| Parametr | Nativní chování | Definice | Zápis |
| --- | --- | --- | --- |
| ENERGY.GRID | 5 pro interní E-INTEGRAL | 5 | explicitní |
| ENERGY.NE | 30 pro interní E-INTEGRAL | 30 | explicitní |
| ENERGY.ImE | obecný default 0,01 Ry, bez FMAG override | 0 Ry | explicitní |
| ENERGY.EMIN | -0,2 Ry | -0,2 Ry | explicitní |
| ENERGY.EMAX | EF při nepřítomnosti | None, delegace na zdroj | vynechán |
| MODE | běžně FREL při nepřítomnosti | None, FREL | vynechán |
| TASK | FMAG příznak, žádné další form-factor volby | FMAG | explicitní |

GRID=5 a NE=30 odpovídají init_mod_energy.f:89–92. ImE=0 Ry je explicitní volba ASE2SPRKKR pro integraci po komplexním oblouku, stejná jako společná definice ENERGY používaná JXC. Liší se od obecného nativního EIMAG=0,01 Ry; protože je zapsaná, nejde o chybný vynechaný default. Nativní 0,01 není automaticky lepší: důležitá je konstrukce a konvergence samotné integrační kontury. Numerické porovnání obou hodnot nebylo provedeno. Vynechané EMAX naopak musí ponechat skutečné EF, nikoliv pevnou hodnotu.

Definice má kkrgen a mpi=True. DENSITY obsahuje MPI barrier a MPI_SUM_R pro hustoty a OBS (density.f:254–264), následné soubory zapisuje MPI_ID=0. Tato větev dokončí výstup ještě před obecnou postprocessing MPI kontrolou v GEN. Paralelní běh nebyl numericky ověřen.

## Nálezy zdroje

1. **FULLPOT je skutečně zakázán, i když původní check může minout task.** check_input.f:119–123 kontroluje TASK='FMAG' a GRID=5. Interní přejmenování na E-INTEGRAL může tuto větev přeskočit. GEN však přímo v CALC_FMAG větvích zastaví FULLPOT (1078–1081 a1180–1183). Definice dokumentuje ASA; FULLPOT je vlastnost potenciálu, kterou obecný validátor InputParameters nezná.
2. **Integrační kontura má být obsazený energetický interval.** GRID=5 je nativní požadavek původního FMAG checku. Definice jej zachovává a pro jiný GRID varuje, místo tvrzení, že zdroj vždy takový vstup odmítne. Konvergence závisí na spodní mezi i počtu bodů.
3. **Pouze první energetická cesta.** GEN pro E-INTEGRAL nastaví NEPATH=1 a DENSITY bere WETAB(IECURR,1). Definice odmítá SPLITSS i dvě položky GRID/NE. Jde o ochranu před ignorovanou druhou cestou; žádný specializovaný split-density workflow zde není zaveden.
4. **Normalizace dělením nulou.** density.f:397–400 dělí FAVR, FSPN a FORB jejich hodnotami při q=0 bez kontroly. Nulový celkový, spinový nebo orbitální moment může způsobit NaN či Inf. Vstupní definice nemá spočtené momenty, takže obecně nemůže problém detekovat. Pro explicitní NREL/SREL varuje před normalizací nemagnetického výsledku. Ani FREL samozřejmě nezaručuje nenulové momenty.
5. **Pevná síť hybností.** NQNEU=min(200,100)=100, QNEUMAX=4*pi*A0_ANG a výstupní osa x=q/(4*pi*A0_ANG) běží od 0 do 1 (density.f:380–403,430–431). NQNEUMAX=200 není počet skutečně vypisovaných bodů ani vstupní default. Rozsah není čten z TASK.
6. **Rozsah typu a MPI redukce nejsou stejně široké.** Hustoty se při inicializaci nulují pro ITBOT..ITTOP (density.f:83–91), ale MPI_SUM_R má délku NRMAX*NTMAX. Pro výpočet omezený jen na část typů tak může redukovat i neinicializované nepoužité části polí. Konkrétní dopad na výsledky vyžaduje běhové ověření.
7. **Relativistický reader je bezpodmínečný.** DENSITY volá WAVFUN_READ_REL, ale obsahuje také vlastní větve IREL<=1 i varování o rozměrech (density.f:107,121–156). Nezavádím proto nezdůvodněný úplný zákaz neplně relativistických režimů. Jejich numerická správnost není ověřena.
8. **Výstup jsou normalizované křivky.** FAVR je součet spinového a orbitálního příspěvku před samostatnou normalizací. Po normalizaci obecně neplatí FAVR=FSPN+FORB. Parser či další zpracování nesmí tyto normalizované křivky bez původních momentů zpětně sčítat.

Kombinace FMAG s DMFT/LDA+U nejsou v tomto patchi bez přímého nativního zákazu odmítnuty; jejich fyzikální konzistence se vstupním potenciálem a kernel výpočtu nebyla numericky ověřena. Patch nevytváří specializovaný parser XMGR souborů a neopravuje Fortran.

## Ověření

8 nových případů: registrace/defaulty/MPI deklarace, vynechané EMAX, souborový roundtrip a pickle; tři způsoby zadání nepodporované druhé kontury; upozornění pro jiný GRID; oba nemagnetické režimy; odmítnutí nulového počtu energií.

Celá sada input_parameters: **116 passed**, Python 3.12, 5,09 s.

```bash
PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 /workspace/scratch/e79ca15d9801/test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q
```

Prošly git diff --check a kontrola reverzní aplikace patche. Numerický SPR-KKR běh nebyl spuštěn.
