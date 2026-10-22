# VBUPS / VBXPS — aliasy VBPES v SPR-KKR 9.7

Nejde o nové samostatné fyzikální tasky. initvar.f:224–249 hledá VBXPS/VBPES/ARPES/VBUPS a podle ANGLE_INTEGRATED volí interní VBPES nebo ARPES. Samostatné VBUPS a VBXPS bez dalších přepínačů nastaví TASK=VBPES. GEN potom spustí již auditovanou VBPES routine.

## Důležitý omyl v předpokládaných defaults
vbpes.f:197 čte TASK.EPHOT s defaultem1253.6 eV, pro oba aliasy stejně. Není větvení podle VBUPS/VBXPS. NINT(EPHOT)=1 znamená1253.6 a NINT(EPHOT)=2 znamená1486.6 eV. VBUPS tedy samo nenastaví21.2 ani40.8 eV. Fotonovou energii pro UV fotoemisi musí uživatel zadat explicitně v definici VBPES. Přidávat nepřítomný UV default, který se nezapíše, by bylo chybné; explicitně zapsaný21.2 eV by byl možný, ale byl by novou Python konvencí, nikoliv native defaultem aliasu.

## Kombinace s ARPES
VBXPS následně samostatně kontroluje ARPES a podle jeho přítomnosti mění ANGLE_INTEGRATED. Nejde o záruku, že libovolná kombinace tokenů zachová integrovaný výpočet. Nové aliasy nepředstavujeme jako náhradu KKRSPEC ARPES/AIPES; executable a native implementace jsou podstatné.

## Doporučené použití a výstup auditu
Používat existující definici VBPES a nastavit TASK.EPHOT v eV podle zdroje. Nevytváří se duplicitní samostatná input definition ani regresní test kopírující stejné native mapování. Tato kontrola doplňuje dřívější audit VBPES; neověřuje numerické PES výsledky.

Kategorie pro budoucí společný audit: aliasy funkčního VBPES, s výrazným ALE pro význam názvu VBUPS a jeho photon default. Nejsou vypnuté ani novou nezávislou implementací. Numerický native výpočet nebyl spuštěn.
