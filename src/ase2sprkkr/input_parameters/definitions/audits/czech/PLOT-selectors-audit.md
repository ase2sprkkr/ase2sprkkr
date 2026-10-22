# PLOTPOT / PLOTRHO / PLOTSFN — SPR-KKR 9.7

Tři selektivní vstupní volby stejné interní výpočetní větve PLOTPRS, doplněné společně. Nejde o tři nezávislé implementace grafiky. Definice mají společný helper _plot_selectors.py; ZIP obsahuje všechny čtyři soubory pro definitions/.

## Native mapování
initvar.f:187–200: PLOTPRS zapne všechna výstupní flags. Pokud není PLOTPRS přítomno, PLOTPOT zapne flag1, PLOTRHO flag2, PLOTSFN flag3; libovolný aktivní flag nastaví interní TASK=PLOTPRS. Tyto selektory lze v native inputu kombinovat, nové jednotlivé definice poskytují samostatné canonical volby.

GEN dispatch gen.f:1461–1463 volá FPPLOTPRS. Routine fpplotprs.f:8–11 a příslušné větve poskytují:
- PLOTPOT: uložené VT/BT, u FULLPOT navíc VNST/BNST;
- PLOTRHO: uložené RHOCHR/RHOSPN, u FULLPOT navíc RHO2NS;
- PLOTSFN: uložené radiální koeficienty FLMSF a LMISF shape functions pro jednotlivé meshe, pouze FULLPOT (ř.392). Nejde o generátor SFN ani přímo o 3D vizualizaci buňky.

Neprovádí SCF ani DOS. Plotuje se existující potenciál/density/SFN data, kvalita a dostupnost závisí na POTFIL. Native výstupy jsou Grace/XMGR, ne obrázky vytvořené ASE2SPRKKR.

## Defaults a omezení
TASK i ADSI mají explicitní příslušný selector. GRID0, NE1, ImE0 jsou explicitní startup hodnoty, nikoliv tvrzení o native plot defaults; energetická smyčka pro PLOTPRS je vypnutá. EMIN/EMAX nejsou exponovány, tato větev je nepotřebuje.

PLOTSFN explicitně zapisuje SCF.FULLPOT=true a odmítá jeho false. Flag nepromění ASA potenciál na správný FULLPOT soubor; skutečný format a SFN data musí být přítomná v potenciálu. Pro ASA native PLOTSFN tiše nevytvoří SFN grafy. PLOTPOT a PLOTRHO nechávají nepřítomný SCF.FULLPOT na rozpoznání formátu souboru, fungují pro ASA i FULLPOT.

V GEN jsou PLOTPRS větve zakázány pro DMFT/LDAU (gen.f:215–223). Nové definice tyto MODE volby nenabízejí; přítomnost/platnost uložených korelačních dat potenciálu nebyla numericky certifikována. MPI je vypnuté, native GEN MPI guard PLOTPRS nepovoluje. Pokud metadata některých typů/meshů nemají odpovídající data, samotná input definice to neodhalí.

Kategorie do budoucího společného auditu: mají plný smysl jako selektivní diagnostický výstup; PLOTSFN podmíněno skutečným FULLPOT potenciálem. Nejde o chybějící fyzikální solver.

## Ověření
4 nové regresní případy: roundtrip a pickle každé definice, metadata serial/executable, explicitně zapsané FULLPOT pro PLOTSFN a odmítnutí false. Celá input_parameters sada 112 passed, dvě původní NumPy warnings. Numerický native běh nebyl spuštěn.

Patch0055 je nezávislý vůči baseline267512c4; nevyžaduje dřívější PLOTPRS patch. Všechny čtyři definitions soubory musí být instalovány společně.
