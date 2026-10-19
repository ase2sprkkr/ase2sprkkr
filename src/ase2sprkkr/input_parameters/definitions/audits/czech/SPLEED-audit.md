# 44. SPLEED — samostatná definice a audit

ASE již nabízí SPLEED jako TASK volbu v definici ARPES, ale chybí samostatné InputParameters.create("SPLEED"). Doplněna samostatná definice sdílející spektroskopické sekce; executable kkrspec. Předběžně **funguje, ale má velké ALE**: tento auditovaný postup je pevný úhel dopadu a energetický scan, ne všechny native polarizační a geometrické větve. Numerický výpočet nebyl spuštěn.

## Rozdílná energetická sémantika

SPEC_RSLAB rozpoznává SPLEED a nastaví IBLOCH=1, NSPLEE=0. SPEC_INPUT pak pro NE>1 čte EMINEV/EMAXEV znovu a používá je jako **dopadající energie v eV**, nikoli jako relativní energie vůči EF. Převzít ARPES −8..+5 eV by proto bylo chybné. Definice zapisuje 11..60 eV, NE=200; popisy obou EV parametrů tuto taskovou sémantiku výslovně uvádějí.

Absolutní EMIN/EMAX jsou v této definici při moderním vstupu odmítány, protože native SPLEED větev je v tomto místě nepřečte místo EV alternativ. NE=1 je odmítnuto: reader v tomto případě přeskočí opětovné načtení EV hranic a EFEV pak převezme EMIN z předchozí obecné inicializace; jednobodový výpočet tedy nemusí mít stejný význam jako vícebodový eV scan. Pro skutečný one-point workflow je potřeba vyřešit tento nativní řetězec, ne pouze nastavit NE=1.

SPEC_PH odstraněno: native pro SPLEED uvádí photon geometry jako obsolete. Není smysluplné slibovat závislost LEED na EPHOT nebo photon polarization.

## Defaults a zápis

| Parametr | Samostatná definice | Native / význam |
|---|---|---|
| EMINEV/EMAXEV | 11 / 60 eV, vypsáno | odpovídá SPLEED start rozsahu; přímé eV |
| NE | 200, vypsáno | native SPLEED start EMESH=200 |
| THETA | 45..45 stupňů | native SPLEED default |
| PHI | 270..270 stupňů | native SPLEED default |
| NT/NP | 1/1, vypsáno | pevné úhly |
| SPOL | 1, vypsáno | native IBLOCH=1 default; **unpolarized** |
| EWORK_EV | 4.2 eV, vypsáno | generická explicitní workfunction; native může použít potenciál |
| IMV_INI_EV/IMV_FIN_EV | 0.2/0.2 eV, vypsáno | SPLEED native start VIL/VIH |
| N_LAYER/NLAT_G_VEC | 1/30, vypsáno | SPLEED native start |
| N_LAYDBL | 10/10, vypsáno | explicitně nastaví obě proměnné navzdory native dvojímu LANZ1 assignmentu |
| POL_E | PZ | existující explicitní konvence; native default také PZ |
| CRYS_VECS | false, nevypsáno | zachová native použití primitivních krystalových vektorů |

SPOL=1 je nepolarizovaný výpočet i při taskovém názvu SPLEED. Nejde o skrytý rozdíl defaults; je vypsán a odpovídá readeru. Pro spinovou odezvu musí uživatel zvolit příslušný native SPOL postup a interpretovat výsledky. Všechny spinové výstupní větve zde nejsou numericky certifikovány.

Výstupní práce 4.2 může být horší než hodnota z konkrétního potenciálu; protože se zapisuje, jde o skutečné explicitní přepsání. Defaults nelze hodnotit pouze podle rovnosti s native. Tlumení je v eV a není elektronovou teplotou.

## Další chyby a omezení

- Převzatá ARPES volba CRYS_VEC má nesprávný nativní název. SPEC_RSLAB čte CRYS_VECS a při jeho přítomnosti nastaví USE_CRY_PRIM_VECS=false. Samostatná definice obsahuje správný klíč a popis vypnutí primitive vectors. Toto je další nález relevantní i pro původní ARPES definici; tento samostatný patch ji nemění.
- Fixed-angle rozsah vědomě nepovoluje angular scans, momentum scans ani POL0/POL0L. Native některé z těchto postupů má, ale polarizace závisí na přítomnosti klíčů: nelze současně používat POL_E a nonzero POL0/POL0L. Slepé sdílení všech ARPES defaults by mohlo skončit STOP.
- Pro INPVER=0 jde o externí starý input.inp; validace moderních energetických klíčů se v této větvi neprovádí. Existence a obsah externích souborů se neověřují.
- IQ_AT_SURF, MILLER_HKL a povrchová struktura musí odpovídat potenciálu. Definice nepředstírá kontrolu tohoto souladu z pouhých input options. Full-potential a další speciální větve nebyly v celém řetězci auditovány.
- mpi=False je konzervativní volba nové definice; nebyla prokázána chyba celé MPI větve ani provedeno numerické srovnání.

## Testy

11 nových testů; celá sada InputParameters **119 passed**, dvě existující NumPy deprecation warnings. Ověřeny discovery, samostatný task a defaults, absence photon sekce, incident EV hranice, round-trip/pickle, serial executable, odmítnutí neplatných hranic/absolutní varianty/jednobodové energie, workfunction, angular a momentum skenů a konfliktní polarization alternative. Samostatný round-trip test ověřuje správný klíč CRYS_VECS.

Příkaz: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

Zdroje: spec_rslab.f:160–225; spec_inputs.f:1198–1241,1301–1328,1358–1389,1549–1648,1764–1783. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Samostatný patch pouze pro SPLEED; původní ARPES audit patch není předpokladem. Bez numerického solveru.
