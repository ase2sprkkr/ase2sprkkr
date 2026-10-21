# AIPES — samostatná definice a audit

AIPES má na `tasks-workflow` vlastní definici `aipes.py` a používá `kkrspec`. Posouzení nativní fyziky: **funguje, ale má velké ALE** — nativní AIPES má konkrétní lokální/nekoherentní aproximační význam, který nelze zaměnit za numerické zintegrování všech ARPES úhlů. Numerický solver nebyl spuštěn.

## Skutečná native větev

SPEC_RSLAB rozpoznává AIPES a nastavuje IBLOCH=2. Ve SPEC_TAU_PREPARE se pro IBLOCH=2 vynucuje CPAATOM=1 pro všechny layer atoms, i když materiál není chemická slitina. To je interní větev výpočtu projektorů a nekoherentních členů, nikoli změna fyzikální occupancy.

SPEC_MECALC volá ZMATINCO pro všechny tyto atomy. SPEC_SUMUP_INT pro AIPES vynuluje ROTRLA/ROTELA a skládá ROALL z ROAT+ROSUR plus ROINC. Pro ARPES naopak používá ROAT+ROTRA+ROTER+ROSUR a jinou korekci nekoherentních členů. Native tedy výslovně označuje úhlově integrovanou fotoemisi s povrchovým příspěvkem a implementuje jinou aproximační větev než prostý součet úhlové sítě.

NT/NP nejsou počet uzlů nějaké obecné hemisférické kvadratury. Zde je auditován postup s pevnou konečnou elektronovou geometrií a NT=NP=1. Definice odmítá angular/momentum scans jako mimo svůj auditovaný rozsah; nejde o tvrzení, že native všechny jejich varianty zakazuje.

## Definice a validace

Sdílí strukturální a spektroskopické sekce ARPES, ale má vlastní TASK AIPES, CONTROL ADSI, validátor a horní hranici energie. Photon sekce zůstává: na rozdíl od SPLEED mají EPHOT a polarizace fyzikální význam.

Validují se pevné konečné elektronové úhly, NT/NP, absence neauditovaných momentum scans, nezáporná konečná workfunction a tlumení, kladné konečné EPHOT a jedna kladná energetická síť bez SPLITSS. INPVER=0 používá externí starý input.inp a moderní parametry se v tomto případě nepovažují za autoritativní; obsah a existenci externích souborů tento input-only validátor neověřuje.

## Defaults a rozdíl proti převzetí ARPES

EMINEV=−8 eV a nově EMAXEV=0 eV **vůči EF** se vypisují. Pro valenční fotoemisi je rozsah do EF vhodnější než převzaté +5 eV, které zahrnuje neobsazené stavy. Změna je explicitní ASE volba, nikoli nevypsaný default, a uživatel ji může změnit. SPEC_INPUT fotoemisní větve vytváří vlastní energetický scan; nelze mechanicky porovnávat tento rozsah s defaults obecné GEN kontury.

Vypsané NE=300, EWORK_EV=4.2 eV, IMV_INI_EV=0.05 eV, IMV_FIN_EV=2 eV, EPHOT=25 eV, NL=4, N_LAYER=50, NLAT_G_VEC=57 a N_LAYDBL={10,10} zůstávají praktické volby zděděné z ARPES. Native photoemission start obsahuje EMESH=200, VIL/VIH=0.1/1 eV, N_LAYER=30 a NLAT_G_VEC=30. Jsou to vypsané odchylky: nelze z nich vyvozovat chybu mapování. Větší počet vrstev může být lepší pro konvergenci, zatímco generická workfunction může být horší než hodnota z konkrétního potenciálu.

Electron THETA/PHI=0..0 a NT/NP=1 odpovídají photoemission defaultu. SPOL=2 zůstává explicitním photoemission defaultem, oproti SPLEED native SPOL=1. Expert polarization konvence nebyly kompletně auditovány; nepředstíráme ekvivalenci vynechání klíče a vypsání stejné numerické hodnoty.

CRYS_VECS=False se **nevypisuje** a ponechává použití primitivních krystalových vektorů. Native přítomnost CRYS_VECS nastaví USE_CRY_PRIM_VECS=false. Původní název CRYS_VEC se nečte, proto samostatná definice nabízí opravený název a popis; pro kontrakt AIPES se používá samostatná definice `aipes.py`.

## Další ALE

- Lokální/nekoherentní aproximaci a povrchový člen je nutné vzít v úvahu při interpretaci a srovnávání s experimentální angle-integrated PES. Tato definice nedokládá, že každá geometrie nebo materiál dává stejný výsledek jako integrál ARPES.
- Je třeba konvergovat vrstvy, reciprocal cutoffy, angular momentum, energetickou síť a tlumení. Výstupní práci i photon energy nastavit podle materiálu a experimentu.
- STRVER/INPVER a generovaná povrchová geometrie musí odpovídat potenciálu. Full-potential, mnohočásticové a všechny spin-density-matrix větve nebyly v celém řetězci certifikovány.
- mpi=False nové definice je konzervativní volba, nikoli důkaz rozbité native MPI implementace.
- Native překlepy dvojího LANZ1 assignmentu a absolutního IMV_INI do VIH platí i zde. Vypsaný N_LAYDBL nastaví obě směrové hodnoty a používané IMV_INI_EV zapisuje správně VIL.
