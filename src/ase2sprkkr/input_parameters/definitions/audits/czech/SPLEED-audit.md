# SPLEED — audit aktuální implementace

**Stav:** Samostatný veřejný task KKRSPEC s experimentálními nativními omezeními. Auditovaná implementace je na `tasks-workflow`, nativní referencí je SPR-KKR 9.7 (`PUB9.7_260420_18h13`). Nejde o numerickou certifikaci.

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

- Fixed-angle rozsah vědomě nepovoluje angular scans, momentum scans ani POL0/POL0L. Native některé z těchto postupů má, ale polarizace závisí na přítomnosti klíčů: nelze současně používat POL_E a nonzero POL0/POL0L. Slepé sdílení všech ARPES defaults by mohlo skončit STOP.
- Pro INPVER=0 jde o externí starý input.inp; validace moderních energetických klíčů se v této větvi neprovádí. Existence a obsah externích souborů se neověřují.
- IQ_AT_SURF, MILLER_HKL a povrchová struktura musí odpovídat potenciálu. Definice nepředstírá kontrolu tohoto souladu z pouhých input options. Full-potential a další speciální větve nebyly v celém řetězci auditovány.
- mpi=False je konzervativní volba nové definice; nebyla prokázána chyba celé MPI větve ani provedeno numerické srovnání.

## Aktuální souborový kontrakt

Vstupy, volitelná cache a výstupy jsou určovány aktuální definicí tasku, nikoli starým odděleným patchem. Zdroj: `src/ase2sprkkr/input_parameters/definitions/spleed.py` a související task registry.

## Ověření

Testy definic jsou v `src/ase2sprkkr/input_parameters/test/`. Nový kompletní pytest ani numerický SPR-KKR běh se zde netvrdí.
