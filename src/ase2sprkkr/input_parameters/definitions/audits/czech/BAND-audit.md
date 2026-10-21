# BAND v KKRSPEC — samostatná definice a audit

BAND má na `tasks-workflow` samostatnou definici `band.py` dostupnou jako InputParameters.create("BAND"). Executable kkrspec. Předběžně **funguje, ale má velké ALE**: komplexní Blochovy kořeny a jejich filtr nelze zaměnit za běžný BSF nebo EKREL výstup. Numerický solver nebyl spuštěn.

## Význam native výpočtu

SPEC_RSLAB rozpoznává BAND a nastaví IBLOCH=0. SPEC_BLOCH v spec_band.f řeší zobecněný vlastní problém a vytváří komplexní kz z −i*log(eigenvalue)/BRUSEP, s další korekcí AKCP. Přijímá kořeny podle abs(Im(kz))<EIG_IME_TOL. EIG_IME proto není energie v Ry ani linewidth v eV; jde o cutoff imaginární části wavevectoru v nativní délkové konvenci.

SPEC_INPUT při NT=1 tiskne označení PROJECTED BAND STRUCTURE a při NT>1 NORMAL (along kz). Samotný SPEC_BLOCH však při NT=1 ukládá do K_EIG reálné kz, při NT>1 ukládá KPARA. Pro interpretaci k-sloupce jsou rozhodující tyto skutečně zapisované veličiny, nikoli převzetí obecného popisku z readeru. SPEC_VBRUN na konci odečítá EFERM od E_EIG.

## Samostatná definice

Sdílí strukturu a přípravu s ARPES, ale nastavuje TASK BAND, CONTROL ADSI BAND a přidává SPEC_BAND:

- NEIGMAX: kladná kapacita uložených kořenů, default 100000. Native při překročení končí STOP.
- EIG_IME: konečný kladný cutoff, explicitně 0.01. Nula by při striktní native nerovnosti nepřijala ani přesně reálný kořen.
- CALC_TRANS_PROB: volitelný příznak výpočtu TRLEED transmission probabilities; default false není vypsán a odpovídá MOD_SPEC_BAND defaultu.

Default K1={0,0}, NK1=1 je vypsán a volí explicitní momentum bod k_parallel=0. Native má platnou jednobodovou k větev. Pro cestu lze nastavit K1 a NK1>1, případně K2/NK2. Jde o posunové vektory přičítané ke KA jako v SPEC_INPUT, nikoli endpointy BLOCHSF.

Validace vyžaduje tento explicitní momentum postup, kladné počty, konečné vektory a jejich počty, jednu energetickou síť bez SPLITSS a konečné nezáporné workfunction/tlumení. Odmítá K3/K4/NK3/NK4, které reader nečte. INPVER=0 se řídí externím starým input.inp, jehož obsah zde neověřujeme.

## Defaults: explicitní odchylky a nevypsané chování

V BAND definici jsou **vypsány** IMV_INI_EV=0 a IMV_FIN_EV=0 eV. Pro neabsorpční band postup jsou vhodnější než převzaté photoemission 0.05/2 eV, ale pro komplexní/disordered band workflow si uživatel může zvolit nenulové tlumení. Nejde o elektronovou teplotu.

EIG_IME=0.01 se **vypisuje**. Native začne na 0.01, ale při nezadaném EIG_IME a VIL>1e−4 nastaví 9999, tedy prakticky vypne filtr. Python proto neslibuje nevypsaný default 0.01 a neponechává závislost filtru skrytou při změně tlumení.

NEIGMAX=100000 se vypisuje a odpovídá native inicializaci. Nulové CALC_TRANS_PROB se nevypisuje a zachová native false. CRYS_VECS=False se nevypisuje a zachová použití primitivních krystalových vektorů; přítomnost správně pojmenovaného klíče native vypne primitive vectors. Původní ARPES CRYS_VEC se nečte a v nové definici je opraven.

Vypsané NE=300, EMINEV=−8 a EMAXEV=+5 eV vůči EF, EWORK_EV=4.2 eV, NL=4, N_LAYER=50, NLAT_G_VEC=57 a N_LAYDBL={10,10} zůstaly praktickými volbami. Pro pásy zahrnující obsazené i neobsazené stavy má +5 eV na rozdíl od valenční AIPES smysl. Počet bodů a cutoffy je nutné konvergovat; nativní menší inicializační počty nejsou automaticky lepší.

SPEC_PH je ponecháno jako součást sdílené nativní spektroskopické přípravy. Nelze z toho vyvozovat, že BAND je fotoemisní intenzita nebo že EPHOT je posun všech band energií: SPEC_VBRUN pro IBLOCH=0 nastavuje EHIGH=ELOW před FINALSTATE. Tento audit necertifikuje každou závislost pomocné přípravy na photon/angle volbách.

## Konkrétní nativní chyby a ALE

1. spec_band.f:70–80 při nevhodném poměru EIGA/EIGB vypíše SPEC_BLOCH: not converged, ale neprovede CYCLE. EIGENVAL zůstává nulové a následuje CDLOG(EIG). Takto hlášené selhání nelze ignorovat; může vytvořit neplatný komplexní kz. Python definice to neopravuje.
2. Filtr imaginárního kz a neabsorpční default nejsou důkaz fyzikální správnosti nebo existence všech pásů. Numerická tolerance, energy mesh, povrchový směr a reciprocal basis ovlivňují přijaté kořeny.
3. CALC_TRANS_PROB je native dostupné, ale transmission amplitudes a jejich normalizace nebyly numericky ověřeny. Default false zachován.
4. Native BAND v SPEC_VBRUN vynucuje ATA=1; disorder/CPA výsledek proto vyžaduje vlastní interpretaci a není touto definicí obecně certifikován.
5. mpi=False je konzervativní nastavení. Native má vlastní paralelní sběr band dat přes rankové soubory; tento řetězec nebyl numericky testován. Není zde prokázána chyba všech MPI výpočtů.
6. Povrchová struktura, workfunction, FULLPOT, many-body a thermal větve nejsou tímto input-only auditem kompletně certifikovány.
