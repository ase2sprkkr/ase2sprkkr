# PHAGEN — audit existující definice

Dodaný SPR-KKR 9.7 versus fetchovaný ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. PHAGEN již podporován je. Samostatný patch 0027-phagen.patch upravuje phagen.py a přidává test_phagen.py proti tomuto základu.

## Co task dělá

TASK PHAGEN v INITVAR:179–184 přepne interní TASK na SCF, nastaví RUN_PHAGEN a vyžaduje KKRSCF. Hlavní program volá PHAGEN místo běžného SCF_DRIVE (kkrscf.f:1392–1400). PHAGEN pak spouští POTWRPHAGEN a pravidelně končí.

Jde o export vstupu pro externí PHAGEN/MsSpec. Neprovádí vlastní phase-shift výpočet. PHAGEN_CALCTYPE je ve phagen.f:34–36 pevně EXP, čtečka parametru je zakomentovaná. Nenabízí se tedy jako uživatelský přepínač.

POTWRPHAGEN uloží původní VT a VMTZ, připojí k POTFIL příponu .only_coulomb, nastaví EXPORT_HARTREE_POTENTIAL=true, VXC=NONE, NSCFITER=1 a MIX=1 a provede SCF_DRIVE. Až pak exportuje radiální mesh, r*V_Coulomb, r^2*rho a r*V_total. Není to čisté přepsání formátu bez výpočtu. Celkový potenciál je z původního vstupu, Coulombův potenciál a hustota z následného kroku; interpretace musí toto rozlišení respektovat.

## Opravy definice

- Nahrazen placeholder popis skutečným účelem a omezeními.
- Doplněna společná MODE sekce: umožní zadat relativitu a validovat explicitní DMFT/LDA+U.
- Doplněny chyby při SCF.FULLPOT, SCF.VXC odlišném od VWN, MODE.OP=DMFT/LDA+U, nekladném NE a split contouru.
- Zachován správný executable=kkrscf a mpi=False.
- Nejsou přidány neexistující TASK parametry pro druh exportu, phase-shift energie ani název výstupu.

Kontroly NCPA, USENLCPA, NT a THERMAL_VIBRA_FLUCT závisejí na potenciálu a pracovním stavu hlavního programu. Native POTWRPHAGEN:49–81 požaduje NT=1, žádné CPA/NLCPA, žádné thermal fluctuations, ASA, VWN a KKRSCF. Samotné InputParameters tyto vlastnosti potenciálu neznají. FULLPOT z potenciálu může být aktivní i při SCF.FULLPOT=false; konečná kontrola zůstává nativní.

## Defaulty a účinnost

| Volba | SPR-KKR | Definice |
| --- | --- | --- |
| PHAGEN_CALCTYPE | pevně EXP | není nabízen |
| ENERGY.GRID | 5, protože interní task SCF | explicitně5, ponecháno |
| ENERGY.NE | 30 pro SCF | explicitně32, ponecháno |
| ENERGY.ImE | obecně0,01 Ry | explicitně0 Ry, ponecháno |
| ENERGY.EMIN | -0,2 Ry | explicitně -0,2 Ry |
| ENERGY.EMAX | EF pro SCF contour | vynecháno, přebírá se EF |
| SCF.VXC | export přijímá pouze VWN | VWN, ponecháno |
| SCF.NITER / MIX při exportním kroku | vynuceno1 / 1 | obecné200 / 0,2 jsou stále zapisovány, ale export je přepíše |
| MODE | běžně FREL | None, vynecháno; explicitní NREL/SREL nejsou v exportní rutině zakázány |

NE32 versus30 je explicitní rozdíl. Nativní30 by o něco snížilo práci, ale bez konvergenčního srovnání není důvod měnit existující32; nejde o nesoulad vynechaného defaultu. ImE0 versus0,01 je rovněž explicitní. GRID5 v EPATH používá vlastní komplexní oblouk, nikoli tento parametr jako jednoduchou konstantní výšku; samotný rozdíl proto nedokazuje rozdílný výsledek.

SCF.NITER a MIX nepředstavují počet iterací a mixing PHAGEN exportního kroku. Tyto hodnoty nejsou znovu načítány při EXPORT_HARTREE_POTENTIAL (scf.f:361–389). Nepřepisují se v Pythonu jen pro kosmetickou shodu: audit přesně rozlišuje vstupní společné defaulty a jejich nativní přepsání. Nelze tvrdit, že jejich změnou uživatel získá konvergovaný Coulombův SCF export.

## Velká ALE ve zdroji

1. **Pouze jeden typ a bez disorder:** NT musí být1. Nejde o obecný konvertor víceprvkového potenciálu.
2. **Spinově polarizovaný export není implementován:** hlavička i writer deklarují non-spin-polarized ASA a zapisují VT/RHOCHR, nikoli výměnné pole BT/spinovou hustotu. Nativní safeguard NONMAG však chybí. Úspěšný běh s magnetickým potenciálem tedy sám nedokazuje zachování spinové informace.
3. **Chybné indexování radiální sítě:** writer zjistí IM=IMT(IT), ale v samotném WRITE používá R(IR,IT), nikoli R(IR,IM) (potwrphagen.f:158–160). Při odlišném mesh indexu se může exportovat nesprávná radiální osa a její násobky. Pro NT=1 a IMT(1)=1 je tento konkrétní problém skrytý.
4. **Hustota a potenciály nejsou všechny z téhož stavu:** původní VTOTAL se kopíruje před Hartree SCF krokem, RHOCHR se vypisuje až po něm. Nejde o jednu univerzální self-consistent sadu sloupců.
5. **MPI větev je nebezpečná:** SCF_DRIVE a některé bariéry jsou uvnitř MPI_ID=0, takže jiné procesy se tohoto volání neúčastní. DEALLOCATE(VTOTAL) je naopak mimo tuto větev, ač alokace proběhla jen na procesu0. Existující mpi=False je správně; pouhé bariéry v rutině nejsou důkaz podpory MPI.
6. **Omezení relativity je zakomentované:** safeguard IREL!=0 byl vypnut. Patch nevymýšlí zákaz všech relativistických režimů; jejich správnost pro následný externí PHAGEN není numericky potvrzena.
7. **Výstupy jsou odvozené od POTFIL:** hlavní dump POTFIL_ex_sprkkr a Coulombový běh s POTFIL.only_coulomb. CONTROL.DATASET samo není izolací těchto názvů. Writer mění globální POTFIL a VXC; následně PHAGEN celý běh ukončuje, takže se původní SCF stav neobnovuje pro další výpočet.

Task není vypnutý placeholder. V podporovaném jednoduchém ASA/VWN/NT1/mesh1 případě má věcný exportní účel, ale tyto podmínky a omezení musejí být uvedeny v závěrečném rozdělení tasků. Výstup neprošel běhovým porovnáním.

## Testy

13 nových případů: serializace/načtení/pickle, původní energetické defaulty, správný executable/MPI, nepřítomnost nefunkčního CALCTYPE, odmítnutí FULLPOT/jiného VXC/DMFT/LDA+U/split contouru/nekladného NE a absence vymyšleného plošného zákazu NREL/SREL. Celá sada input_parameters: **121 passed**. Numerický SPR-KKR ani externí PHAGEN/MsSpec nebyly spuštěny. Patch prošel git diff --check a kontrolou zpětného aplikování.
