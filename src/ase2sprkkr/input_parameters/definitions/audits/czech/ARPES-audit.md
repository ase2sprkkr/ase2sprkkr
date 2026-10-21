# ARPES — audit současné implementace

**Stav:** Použitelný pro auditované geometrie; expert režimy a native konvence vyžadují kontrolu. Audit vychází ze současné implementace na větvi `tasks-workflow` a nativního zdroje SPR-KKR 9.7 (`PUB9.7_260420_18h13`). Hodnocení není numerickou certifikací.

## Rozhraní, validace a nativní omezení

- ENERGY.EWORK_EV správně popsáno jako **výstupní práce v eV**, nikoli vnitřní potenciál. SPEC_INPUT jej čte a převádí /RY_EV; v příslušné větvi při zadání nastavuje EWORK.
- NT, NP a NK1..NK4 jsou kladná celá čísla; NK1=1 zůstává přípustné. Reader má pro jedno k explicitní větev a dělení NT−1/NP−1 chrání podmínkou >1. Vyžadovat minimum 2 by neoprávněně vyřadilo platný native postup.
- Při K1 se vyžaduje NK1, při současném K2 také NK2. Záporné a nulové počty dříve mohly vést k prázdným či neplatně dimenzovaným sítím. Existující kontrola souladu angular ranges a NT/NP zůstala zachována.
- Pro ARPES odmítnuty K3/K4/NK3/NK4: v dodaném SPEC_INPUT readeru se tyto klíče nikde nečtou. Ve veřejné Python definici zůstávají pro zpětné parsování, ale validace nepředstírá jejich účinek.
- Konečné nezáporné EWORK_EV a obě IMV_*_EV, kladné konečné EPHOT; konečné elektronové úhly a k-vektory. Nulové tlumení není svévolně zakázáno.
- N_LAYDBL musí mít přesně dvě nezáporná celá čísla. N_LAYER a NLAT_G_VEC jsou kladná celá čísla. INPVER/STRVER omezeny na 0/1. Při INPVER=0 je nové řízení moderního readeru vynecháno, protože data pocházejí ze starého externího input.inp.

## Výchozí hodnoty a význam

| Parametr | ASE default | Srovnání a zápis |
|---|---|---|
| EWORK_EV | 4.2 eV | vypsáno; přepíše workfunction ze SCF/potenciálu |
| IMV_INI_EV / IMV_FIN_EV | 0.05 / 2 eV | vypsáno; native photoemission start VIL/VIH je 0.1/1 eV |
| NE | 300 | vypsáno; native photoemission start EMESH=200 |
| EMINEV / EMAXEV | −8 / +5 eV vůči EF | vypsáno; explicitní spektrální rozsah, ne obecná fyzikální mez |
| NL | 4 | vypsáno; explicitní angular cutoff |
| NT / NP | 1 / 1 | vypsáno, konzistentní s nulovým rozsahem úhlů |
| Electron THETA/PHI | 0..0 pro ARPES | task-specific default; SPLEED má 45/270 stupňů |
| EPHOT | 25 eV | vypsáno; explicitní experimentální volba |
| N_LAYER | 50 | vypsáno; native photoemission start 30 |
| NLAT_G_VEC | 57 | vypsáno; native start 30 |
| N_LAYDBL | 10 / 10 | vypsáno; native má překlep LANZ1=10 dvakrát místo inicializace LANZ2 |
| INPVER / STRVER | 1 / 1 | vypsáno; moderní vstup a generování struktury |

Tyto vypsané odchylky byly ponechány. N_LAYER=50 či více reciprocal vectors mohou zlepšit konvergenci za cenu času; nelze bez systému doporučit nativní menší čísla jako lepší. Výstupní práce a tlumení jsou citlivé materiálové/experimentální parametry: hodnota z potenciálu může být vhodnější než generická 4.2 eV, ale její použití vyžaduje vynechat tento explicitní klíč. Default není pouhé formátovací nastavení.

Expert photon polarization a spin-amplitude defaults se při nezadání obecně nevypisují. Tento audit nemění jejich logiku a nepředstavuje kompletní certifikaci všech kombinací POL_P/NPOL/ICIRC/IDREH ani POL_E/POL0/POL0L. Reader podle přítomnosti alternativních klíčů některé kombinace odmítá; samotná shoda numerické default hodnoty nestačí k zaměnitelnosti zápisu a vynechání.

## Konkrétní nativní chyby a omezení

1. **IMV_INI v Ry má překlep:** spec_inputs.f:1257 volá SECTION_SET_REAL('IMV_INI',VIH,...), tedy zapisuje do final-state proměnné místo VIL. ASE nabízí IMV_INI_EV, jehož větev zapisuje správně do VIL; nepřidán chybný Ry alias. VIL/VIH reader následně převádí do interního spektroskopického zadání v eV.
2. **Dvojí LANZ1 inicializace:** reader nastavuje LANZ1=10 dvakrát. Explicitní N_LAYDBL={10,10} v ASE nastaví obě proměnné, a je proto lepší než slepé vynechání pole s očekáváním obou native defaults.
3. K1/K2 v tomto readeru jsou **posunové vektory přičtené ke KA**, nikoli endpointy jako v BLOCHSF. Počet bodů určuje NT/NP a zadané vektory se použijí pro k-mesh. Názvy nelze interpretovat jednotně mezi tasky.
4. Native photon convention uvádí normal incidence THETA=180 stupňů. Nelze automaticky přenášet geometrický význam THETA z jiných tasků. Specifika krystalového rámce a sample tilt je nutné ověřit pro konkrétní povrch.
5. Validace neověřuje existenci a obsah externího input.inp/struc.inp při INPVER/STRVER=0, soulad povrchové struktury s potenciálem ani konvergenci layer doubling a reciprocal cutoffů.
6. Executable kkrspec a existující paralelní nastavení ponecháno; numerický MPI/serial výpočet nebyl proveden. Zde provedený ARPES audit není potvrzením všech ostatních TASK voleb sdílejících tutéž definici.

## Implementační kontrakt
- Autoritou pro soubory je aktuální `TaskDefinition` a případné dynamické metody schématu. Výstupy a cache se nesměšují s povinnými vstupy.
- Zdroj: `src/ase2sprkkr/input_parameters/definitions/arpes.py`.

## Ověření

Testy definic a deklarací souborů jsou v `src/ase2sprkkr/input_parameters/test/`. Tento dokument netvrdí nové spuštění úplné sady ani numerický výpočet SPR-KKR.
