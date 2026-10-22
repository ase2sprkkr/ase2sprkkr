# TUTORIAL — úplný audit SPR-KKR 9.7 / ASE2SPRKKR

Tento dokument nahrazuje dílčí TUTORIAL audity. Celá definice má 14 voleb TUT-TASK. Implementace připravuje a validuje input; není opravou Fortran solveru ani potvrzením numerické správnosti všech větví.

## Kontrola předchozího souboru
Před posledním rozšířením tutorial.py skutečně obsahoval POT-XC a CORE-ALL, včetně info a validace FULLPOT. Pro jednoznačné stažení nyní přiložen tutorial-complete.py; jeho obsah je shodný s aktuálním definitions/tutorial.py a instaluje se pod názvem tutorial.py. Předchozí dělení implementace po větvích je ukončeno.

## Všechny volby
| TUT-TASK | Výpočet | Stav / hlavní podmínky |
|---|---|---|
| YLM | Reálná sférická harmonická jako prostorová plocha; kladná/záporná část zvlášť | Plný smysl pro výuku; pevná výstupní jména |
| JLX | Sférické j_l(x), n_l(x) | Plný smysl; LMAX<=4, NX>=2 |
| JLX-ASYMPT | Totéž proti malému/velkému x asymptotě | Plný smysl; SMALL vynucuje XMIN=.001 |
| CHR-POT | Uložené V/B, charge/spin density, FULLPOT shape functions | Plný smysl pro diagnostiku dostupných dat |
| POT-DECOMP | Jaderný, Hartreeho a MJW XC/X potenciál | Velká ALE: ASA, původní XC VBH/MJW/VWN, NONMAG bug |
| POT-XC | VBH/MJW/VWN/PBE na stejné hustotě | Velká ALE: ASA, NONMAG bug, pevný rozsah grafu |
| RWF-L | Regular radial functions napříč l při EMIN | Velká ALE; dodatečný matching graf je rozbitý |
| RWF-E | Regular radial functions napříč energií při zadaném l | Velká ALE: scalar relativity, ASA, reálná část a normalizace |
| CORE-ALL | Nemagnetické core stavy typu 1 | Velká ALE; bez core slupek rozbitá větev |
| DOS-FEG | V=.0001 Ry, B=0, single-site DOS | Smysluplná demonstrace; STANDARD KKR, není plná bulk free-electron integrace |
| DOS-SS | Single-site DOS původního potenciálu | Smysluplné v STANDARD KKR; CPA není certifikováno |
| DOS-CLU | DOS konfigurace clusteru | Smysluplné s explicitním TAU clusterem; název jej nezapne |
| DOS-BULK | DOS periodického systému | Smysluplné s periodickou BZ integrací |
| GEN | Obecná GEN větev a DOS/observables | Smysluplné, výsledek určuje konfigurace TAU/potenciálu |

Žádná z těchto 14 voleb není obecně vypnutý placeholder. Rozdělení se týká tohoto tasku; společný audit všech tasků ještě není vytvářen.

## Native dispatch a vstupy
initvar.f:700–718 vyžaduje TUT-TASK; pro DOS* a GEN zapíná energetickou smyčku. check_input.f:167–172 vyžaduje KKRGEN. tutorial.f:53–184 obsahuje celou dispatch tabulku, po dokončení STOP_REGULAR (normální ukončení). Přesné canonical hodnoty jsou deklarovány v Python; native prefixové varianty nepřidáváme jako další tasky.

YLM čte L/M (2/0), tiše ořízne M. JLX čte LMAX/NX/XMIN/XMAX (2/300/.001/5), ořízne LMAX na4, opravuje meze. JLX-ASYMPT navíc X-RANGE SMALL/LARGE. Validace odmítá neplatné meze/indexy místo tichých oprav.

RWF čte TASK.ITSEL (1), NORMALISATION (1), ENERGY.EMIN/EMAX/NE (1e-7 Ry/1 Ry/10). NORMALISATION=1 integrálně normalizuje **zobrazovanou reálnou část** radiální funkce; jakýkoliv jiný native řetězec zachová MST amplitudu, definice používá canonical MST. RWF-E navíc TASK.L (2, native clamp na NL-1). Existence typu ITSEL a skutečný NL vyžadují potenciál; L kontrolujeme vůči SITES.NL, pokud je explicitní.

CHR-POT nastaví celé PLOTPRS=true a volá FPPLOTPRS. Routine fpplotprs.f větví výstupy podle FULLPOT a dostupných fields; nemá další vlastní čtené input parametry. Tento výběr nepotřebuje ASA omezení ani omezení na původní XC.

DOS-FEG přepíše VT=.0001 a BT=0, DOS-SS zachová VT/BT; obě vynutí IBZINT=1, NPTMIN=NPTMAX=0. V tau_drive.f:88–102 se tau=t použije pouze pod STANDARD-KKR; definice proto tento režim explicitně zapisuje a odmítá cluster pro tyto dvě volby. Nelze jej zaměnit za běžnou bulk FEG integraci. DOS-CLU/BULK/GEN pouze volají GEN s DOSCORSYS=LOC; samy nepřepínají geometrii. DOS-CLU vyžaduje TAU.CLUSTER/MOL/BZINT=CLUSTER a CLURAD či NSHLCLU; DOS-BULK odmítá cluster. GEN zachovává nastavenou geometrii. GEN pro TASK=TUTORIAL vypíná core calculation (gen.f:394) a zapisuje CALCDOS_WRITE.

## Defaults: zapsané versus nepřítomné
Default volba YLM je explicitně zapsaná Python volba; native bez TUT-TASK zastaví. Matematické defaults souhlasí s native.

ENERGY defaults se mění podle aktuálního TUT-TASK, pokud je uživatel nepřepsal:
- matematika/plots/core: GRID0, NE1, ImE0; meze -.2..1 Ry jsou jen startup hodnoty, tyto větve je nevyužívají;
- RWF: GRID3, NE10, EMIN1e-7, EMAX1, ImE0; samotná routine znovu sestaví EPATH(3), s nulovým imaginary part a ignoruje vstupní GRID/ImE/eV meze. Definice odmítá relativní eV meze, více kontur/SPLITSS a NE1;
- DOS/GEN: GRID3, NE300, EMIN-.2, EMAX1, ImE.01 Ry. Vše se explicitně zapisuje; jde o praktickou spektrální volbu, nikoliv native default contour. EMAX1 umožňuje neobsazenou část; alternativou je okno do EF podle cíle výpočtu.

MODE explicitně SREL pro RWF, FREL pro ostatní (native FREL je absence NREL/SREL/SP-SREL; token FREL nemění tento výsledek). TAU.KKRMODE=STANDARD explicitně ve všech výchozích variantách; uživatel jej může změnit kromě DOS-FEG/SS. Cluster v native init_mod_kmesh nastavuje REAL-SPACE podle IBZINT0. BZINT defaults POINTS, s CLUSTER/MOL CLUSTER; tím se současně vyhýbáme původnímu wrapper problému callable None převáděného na řetězec NONE.

SCF sekce obsahuje jen FULLPOT flag, při defaultu je nepřítomná. To není záruka ASA: native přebírá formát potenciálu. Není zaveden falešný SCF.VXC override: tuto hodnotu pro POT-DECOMP čte XC-POT z potenciálu (fallback VWN v potio.f:963), POT-XC ji interně přepíná.

## RWF native chyby a limity
TUTORIAL_RWF (tutorial.f:1422–1760) používá WAVFUN_READ_SRA, proto definice dovolí jen NREL/SREL/SP-SREL, default SREL; FREL odmítá. IRTOP se nastavuje JRWS(IM), zatímco FULLPOT reader očekává JRCRI(IM): FULLPOT je z auditovaného rozsahu vyloučen. Numerická energie blízko nuly se interně mění na1e-5 Ry.

RWF-L hlavní graf vytváří l=0..NL-1 při první energii. Další matching graf má zásadní chyby: YMATCH(I,1)=YY(I,0,0), přičemž YY(:,0,:) zůstalo nulové; L po smyčce odpovídá poslednímu kanálu. Pro vnější radius se plní YMATCH(:,2), ale YMATCH(:,1) se nevyplní. Matching graf proto obsahuje nulovou vnitřní wavefunction a neinicializovaná vnější data. Samotná definice tento native problém neopraví.

RWF normalizace má 1/sqrt(norm) bez guardu pro nulovou normu. Výstup kreslí DREAL(ZG), nikoliv komplexní wavefunction nebo celou relativistickou spinorovou amplitudu. INDEX0 je v main grafu nulová referenční křivka. Rutina přeskládá ETAB/WETAB, používá NEMAX z inicializace a vyžaduje konzistentní jeden mesh; definice drží stejné NE při startupu i rereadu. Platnost potential data/ITSEL nelze garantovat bez jejich načtení.

## Paralelismus a rozsah ověření
mpi=False pro celý task. GEN MPI guard nepovoluje TUTORIAL, matematické a plot větve používají společné soubory. Žádná numerická SPR-KKR kalkulace nebyla provedena.

40 nových tutorial regresních případů; celá input_parameters sada: 148 passed, dvě původní NumPy warnings. Pokryto všech14 canonical voleb roundtripem, pickle RWF, změna závislých defaults, geometrie DOS, mesh, explicitní FULLPOT, relativita, energy reference a matematické indexy. Testy ověřují Python definici, nikoliv správnost native výstupů.

## Detailní nálezy potenciálových/core větví
## Co skutečně počítá
Pro každý typ ITBOT..ITTOP a radius R do JRWS vytvoří pět křivek: V_tot převzaté z VT, V_nuc=-2*Z/r, reziduální V_H=V_tot-V_nuc-V_xc(previous), MJW V_xc a MJW exchange-only V_x. Původní V_xc se znovu vyhodnotí z uložené RHOCHR/RHOSPN funkcionálem označeným v potenciálu. Pro magnetické typy vytváří analogické B křivky ze spinového rozdílu. Výstup je Grace/XMGR soubor vytvořený XMGRHEAD s DATSET0, klíčem Vdecomp a TXT_T, radius v a.u. a potenciál v Ry. Neprovádí SCF ani DOS.

Pozor: vykreslený MJW V_xc není vždy V_xc(previous) odečtený při konstrukci Hartreeho rezidua. Pro vstup VBH/VWN proto V_nuc+V_H+vykreslené V_xc obecně nereprodukuje V_tot. V_x je samostatný výměnný příspěvek, nikoli další člen, který se má přičíst k V_xc.

## Skutečná omezení
- FULLPOT: tutorial.f:633 provede RETURN bez výstupu. Nová expert optional sekce SCF.FULLPOT umožňuje explicitní požadavek odmítnout. FULLPOT obsažený pouze v POTFIL bez načtení potenciálu Python validátor nezjistí.
- Původní SCFVXC musí být VBH, MJW nebo VWN (ř.671–693). PBE vede ke STOP, ostatní nepodporované funkcionály také. Hodnota se čte z XC-POT v potio.f:963; tamní fallback je VWN. SCF.VXC je čteno v scf.f:384, tato tutorial větev SCF nespouští. Proto není vystaven falešný override ani vymyšlený default funkcionálu.
- Hartreeho položka je reziduum vůči uloženému celkovému potenciálu, jadernému členu a rekonstruovanému XC, nikoli nezávisle řešená Poissonova rovnice. Vyžaduje konzistentní density/potential data.
- NONMAG je sdílené stavové pole: při abs(B) pod 1e-4 se nastaví true (ř.722–723), ale v následujícím typu se nevrací na false. Smíšený systém, kde nemagnetický typ předchází magnetickému, může přijít o B grafy. Už počáteční NONMAG=true B grafy potlačí. Je to native chyba, neopravitelná samotnou input definicí.
- Serial nastavení TUTORIAL zachováno; společné výstupní soubory nejsou certifikovány pro MPI.

Kategorie do budoucího společného auditu: funguje, ale má velká ALE. FULLPOT cesta neposkytne výsledek; nejde o obecně vypnutý placeholder.



## Co porovnává
Pevně zakódovaný seznam VXCTAB (ř.847) obsahuje VBH, MJW, VWN (LDA) a PBE (GGA). Pro každý typ ITBOT..ITTOP vyhodnotí všechny čtyři funkcionály na stejné uložené RHOCHR/RHOSPN. U PBE připraví spinové hustoty, radiální derivace, Laplacián a grad rho . grad |grad rho| (ř.934–987). Výstup Vxc=(Vup+Vdown)/2 a Bxc=(Vup-Vdown)/2 zapisuje do Grace/XMGR souborů s DATSET0, klíčem Vxc a TXT_T; radius v a.u., potenciály v Ry.

Nejde o čtyři samostatná self-consistent řešení. Výchozí XC-POT potenciálu není omezeno na VBH/MJW/VWN jako u POT-DECOMP: zde se SCFVXC vnitřně přepíná mezi pevnými čtyřmi funkcionály. Nezavádíme vstupní SCF.VXC, protože tato větev jím seznam neřídí. Kvalita výsledku závisí na uložené hustotě a jejích radiálních derivacích.

## Omezení / omyly native implementace
- FULLPOT vede k tichému RETURN (ř.867). Validátor odmítá explicitní SCF.FULLPOT; full-potential formát samotného POTFIL bez načtení potenciálu neověří.
- NONMAG se při malém Bxc nastaví true (ř.1028–1029), ale pro další typ se nevrací false. Nemagnetický typ před magnetickým proto může potlačit Bxc grafy pozdějšího typu; problém je obdobný POT-DECOMP.
- Rozsah grafu Vxc se po výpočtu přepíše na -15..0 Ry (ř.1031–1032). Data mimo rozsah se zapisují, ale v běžném zobrazení mohou být skryta. Nejde o numerický cutoff.
- SCFVXC zůstane po výpočtu PBE, původní hodnota se neobnovuje. TUTORIAL potom normálně končí STOP_REGULAR, proto ve standardním standalone běhu nenásleduje SCF s nečekaným funkcionálem. Při budoucím vložení rutiny do delšího workflow by to byl významný side effect.
- PBE příprava upozorňuje na kladné radiální derivace hustoty, ale pokračuje; výrazy GDGAG jsou sestavené jako -n'*n''. Nejde o potvrzení správnosti pro libovolný profil hustoty.
- MPI zůstává vypnuté, paralelní výstupy nejsou certifikovány.

Kategorie pro budoucí společný audit: funguje, ale má velká ALE. Není vypnutý placeholder.



## Skutečný význam
Navzdory názvu ALL používá IT=1 (ř.1199), nikoli všechny atomové typy. BT se vynuluje pro všechny typy a NONMAG se nastaví true (ř.1177–1178). Přepočítává nemagnetické core stavy ve sférickém potenciálu VT prvního typu, nikoli původní spinové rozštěpení. Z jednotlivých slupek vybírá m_j=1/2; z výpočtu CORE ukládá první radiální komponentu pro každý vybraný stav.

Výběr slupek: pevná tabulka n,l a odhad podle Z; pokud NCORT(1) !=0, počet slupek se odvodí z kumulativních zaplněných slupek odpovídajících NCORT. Pokud shoda neexistuje, native STOP. Nejde o volně nastavitelný CL selector.

Výstup Grace/XMGR obsahuje VT(r), energie core hladin a radiální g/f komponenty škálované a posunuté na energii hladiny: E+0.2*g(r), E+0.2*f(r). Proto svislá poloha křivky není původní amplituda wavefunction. Radius je v Bohr, energie grafu v Ry; konzole uvádí i eV, normy a další kvantová čísla. Energie zde není odečtena od EF. Soubor je tvořen XMGRHEAD z DATSET0, klíče core-all a TXT_T(1).

## Native omyly / omezení
- FULLPOT vede k STOP_MESSAGE (ř.1153–1156), nikoli k podporované výpočetní větvi. Validátor odmítá explicitní SCF.FULLPOT; formát samotného POTFIL bez načtení potenciálu neověří.
- Bez core slupek, např. Z<=2 a NCORT=0, není NNKCOR inicializováno: nastavuje se až v cyklu nalezených core stavů, ale následně se používá pro počet křivek a další smyčky. Po skončení prázdného shell cyklu je navíc používáno MJ/LCXRAY. Tato větev je rozbitá. Definice bez čtení potenciálu nezaručuje neprázdné core stavy.
- Při zápisu VT se používá délka JRWS(IT)-IVBOT+1 místo JRWS(IM)-IVBOT+1 (ř.1366–1367), zatímco data a radius jsou indexované IM=IMT(1). Pro první typ s mesh indexem jiným než 1 může být délka chybná; při odlišných mezích hrozí čtení mimo platnou radiální oblast.
- BT, NONMAG, NCXRAY a LCXRAY jsou mutovány. V obvyklém standalone tutorial běhu rutina končí STOP_REGULAR; tyto mutace se neukládají jako nový potenciál.
- MPI zůstává vypnuté; serial výstupy nejsou certifikovány pro paralelní použití.

Kategorie pro budoucí společný audit: funguje, ale má velká ALE pro ASA potenciál typu 1 s neprázdnými konzistentními core slupkami; bez core slupek rozbitá větev. Nejde o obecně vypnutý placeholder.


## Dodání
0053-tutorial-complete.patch je úplný kumulativní patch vůči baseline 267512c4, nahrazuje 0049–0052. Nepoužívat staré tutorial patche současně. Po již aplikovaném starém patchi nahraďte tutorial.py a test_tutorial.py celými soubory.
