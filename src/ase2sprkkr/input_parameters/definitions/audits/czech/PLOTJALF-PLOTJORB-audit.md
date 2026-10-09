# PLOTJALF / PLOTJORB — SPR-KKR 9.7

Kategorie pro společný audit: **rozbité**, ne obecně vypnuté placeholdery. Jde o postprocessing voleb JALF/JORB, nikoliv přepočet proudů.

## Zamýšlený vstup
CONTROL.PLOTJALF a CONTROL.PLOTJORB jsou skutečné flags (initvar.f:122–123), nikoliv samostatné TASK keywords. Každý nastaví interní TASK=PLOT2DPR (ř.205–206). Geometrie se čte z TASK.AVEC/BVEC/DVEC/NA/NB. V očekávané proudové větvi fpplot2dpr.f:194–231 se načítá DATSET0_JALF.dat nebo DATSET0_JORB.dat; při současném nastavení obou má JALF přednost. Požaduje odpovídající typy, radiální meshe a cutoff, jinak STOP data file not consistent. Vstupní data mají pocházet z dřívějšího JALF/JORB výpočtu.

## Proč definice nemůže zajistit funkčnost
GEN dispatch (gen.f:1466–1478) pro interní PLOT2DPR vybere FPPLOT2DPR na FULLPOT, PLOT2DPR_ASA na ASA.

FPPLOT2DPR nejprve bezpodmínečně volá PLOT2DPR_FP_PLY (fpplot2dpr.f:70). Tato rutina nečte PLOT_JALF/PLOT_JORB ani proudový soubor. Pracuje s potential/density a na konci končí STOP_REGULAR s textem run of test version completed (plot2dpr_fp_ply.f:926). Návrat do následné proudové větve se tedy v běžném FULLPOT spuštění nekoná, ani když PLOT2DPR(1/2) nejsou zapnuté. Podmíněný RETURN na ř.103 zde nepomůže: vyžaduje NOT FULLPOT.

ASA routine PLOT2DPR_ASA proudové flags ani data nečte, vykresluje potenciál. Proto ani změna FULLPOT flagu není obcházení problému.

Výsledkem může být jiný plot output nebo dřívější chyba, nikdy však touto doloženou cestou požadované vykreslení uloženého JALF/JORB pole. Přítomnost proudového souboru sama dispatch neopraví.

## Dodání a další práce
Nevytváří se input definice ani Python test kopírující nedosažitelnou native větev. Nejprve je nutné oddělit testovací PLY writer od proudové větve ve FPPLOT2DPR a ověřit numerický postprocessing referenčního JALF/JORB souboru. Pouhé vystavení CONTROL flagů v ASE2SPRKKR by nebylo funkční podporou. Native Fortran nebyl měněn; numerický běh neproběhl.
