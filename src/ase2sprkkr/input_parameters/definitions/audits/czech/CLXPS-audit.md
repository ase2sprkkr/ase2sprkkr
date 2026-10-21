# CLXPS — audit současné definice

CLXPS má na `tasks-workflow` vlastní definici `clxps.py`; nativní omezení vycházejí z dodaných zdrojů SPR-KKR 9.7.

## Rozsah

Core-level XPS v kkrgen, ASA a plně relativistický režim IREL=3. FULLPOT zakazuje check_input.f:85–88, DMFT/LDA+U gen.f:211–222. Definice odmítá explicitně nekompatibilní MODE a korelace; FULLPOT z potenciálu se kontroluje ve SPR-KKR. IT je jedničkový index typu, jeho horní mez závisí na potenciálu.

TASK obsahuje ME, IT, CL, EPHOT, ImE, MECHECK, LINPOL, USETSS, ANGRES a geometrie ANGQ0..10, ANGE0..10, ANGK0..10. Úhly jsou dvojice polárního a azimutálního úhlu ve stupních. EPHOT je v eV, ImE v Ry. USETSS nahrazuje scattering-path operator single-site t maticí. LINPOL přepíná z kruhové na lineární polarizaci.

## Defaulty a jejich zápis

| Parametr | SPR-KKR | Nová definice | Zápis defaultu |
| --- | --- | --- | --- |
| TASK.ME | ADA | ADA | vynechán |
| TASK.IT | 1 | 1 | vynechán |
| TASK.CL | 2P | 2P | vynechán |
| TASK.EPHOT | 1253,6 eV | 1253,6 eV | vynechán |
| TASK.ImE | nejprve 0, pak při abs(ImE)<1e-8 změna na 0,01 Ry | 0,01 Ry | vynechán |
| čtyři přepínače TASK | false při nepřítomnosti | false | vynechány |
| ENERGY.GRID / NE / ImE | 3 / 18 / 0,01 Ry pro CLXPS | stejné | explicitně zapsány |
| geometrie | interní nulové úhly, počet podle nalezených trojic | neuvedené hodnoty None | vynechány |

TASK.ImE odpovídá efektivnímu nativnímu defaultu, nikoliv mezihodnotě 0. Explicitní téměř nulová hodnota vyvolá upozornění na nativní náhradu. Zde není důvod měnit nativní defaulty. Obecné inicializační ENERGY.EMIN není fyzikálním počátkem XPS spektra; konečné energie se sestavují jako ECOR + EPHOT + i*TASK.ImE (clxps.f:203–210), nikoliv skenováním ENERGY.EMIN/EMAX.

## Nálezy ve zdrojích

1. **ME=GRV není dokončeno.** clxps.f:186 připraví derivace pro GRV, ale následné větve výpočtu podporují ADA a ostatní skončí STOP. Definice proto přijímá pouze ADA. Nejde o numerické ověření správnosti ADA.
2. **Jedna geometrie se skutečně čte s indexem 0.** Smyčka sestaví ANGQ0, ANGE0, ANGK0. Následné LL=4 nezmění již sestavené názvy předané čtečce. Nepřidávám tedy zdánlivě zamýšlené nesufixované ANGQ/ANGE/ANGK. Kompletní geometrie 0 má přednost před 1..10.
3. **Číslované geometrie musí být souvislé.** Čtení 1..10 skončí na první nekompletní trojici. Zdroj tak může další geometrii tiše ignorovat. Definice odmítá mezery a nekompletní trojice; upozorní, pokud geometrie 0 překryje další.
4. **ANGRES bez geometrie může ponechat NGEO=0.** Definice vyžaduje kompletní geometrii. Bez ANGRES zdroj nastaví jedinou geometrii; při více zadaných definice upozorní.
5. **Velikost energetické tabulky musí stačit počtu core stavů.** Rutina přepíše NETAB(1)=4*l+2 a zapíše tolik řádků ETAB. Definice zachovává nativní NE=18 a kontroluje alokovaný rozsah při jeho snížení. Pro f shell je třeba 14 stavů. Kontroluje také fyzikální l<n.
6. **Konečná energie pod EF zastaví výpočet.** Tuto podmínku nelze ověřit jen ze vstupu bez core energií potenciálu; zůstává kontrolou SPR-KKR.

MPI je v definici konzervativně vypnuto pro tento postprocessing a jeho výstupy. To není tvrzení, že zdroj explicitně zakazuje paralelní energetickou smyčku; paralelní CLXPS nebyl ověřen.
