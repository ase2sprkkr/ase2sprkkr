# Opravy obou review KKRunneru

Práce je ve větvi `kkrunner` všech tří repozitářů. Tento dokument nahrazuje
otevřený stav historických reportů `kkrunner-review.md` a `kkrunner-review-2.md`.
Schéma RunRecord zůstává **1**, bez migrací vývojových záznamů.

## Provedené opravy

| Nález | Oprava |
| --- | --- |
| R01: souběh submitu a postprocessingu | Claim ukládá backendovou identitu; dokončení kontroluje vlastní ID pokusu a stav místo rovnosti celého snapshotu. Běžná aktualizace ho nezneplatní. |
| R02: duplicitní ID | Identita je unikátní napříč datumovými shardy. Druhý submit je odmítnut ještě před backendem. |
| R03, R11: souborové kolize | Lokální cesty se kontrolují nezávisle na backendu. Vzdálené cesty se vyřeší před rezervací bez mkdir/upload/otevření výstupů. Přibyl explicitní filesystem_id pro aliasy a sdílené mounty. |
| R04: nechráněný save | save slouží k vytvoření nebo identickému opakování; další změny používají atomické update se zámkem podle globálního ID. |
| R05: cizí běhy v GUI | GUI vybírá vlastní namespace. Cizí platná metadata neblokují historii; poškozená vlastní metadata se stále hlásí. |
| R06: stdin mimo manifest | stdin je samostatná čtecí deklarace v manifestu, rezervaci i zachycených verzích; nezaměňuje se za primární vstup. |
| R07: vnořená validace | Execution, endpointy a manifesty odmítají neplatné vnořené typy. Dekodér ověřuje i shodu identity s umístěním záznamu a vrací CorruptRunRecord. |
| R08: streamy při chybě submitu | Nově otevřené pathname streamy mají vlastní ExecutionResources a chybová cesta je uzavře. Staging je neotevírá. |
| R09: dostupnost projektů | ASE i GUI zachovají původní definice projektů. Nové záznamy běhů jsou oddělené od starých vývojových formátů; explicitní repozitáře mají vlastní store. |
| R10: retry postprocessingu v GUI | Přibylo explicitní tlačítko a worker volající aplikační adaptér, následné obnovení historie a převzetí výsledku na GUI straně. |
| R12: povinné výstupy | Chybějící povinný výstup vede k chybě artefaktů, volitelný se přeskočí. Kontrola je i v raw lokálním a PSI/J běhu. |
| R13: backendové receipts | Dokončovací cesty sdílejí ACK po trvalém zápisu. Refresh umí opakovat neúspěšné ACK. Recoverable download failure receipt ponechá pro reattach a retry. |
| Dřívější deadlock + chybný JSON kořen | Opravy jsou zahrnuté v regresní sadě. |
| A01: callback volající wait | Monitorovací vlákno může přečíst vlastní výsledek bez čekání na sebe; ostatní čekající stále čekají i na callback a jeho případnou výjimku. |
| A02: platformy | Linuxový požadavek spravovaného lokálního backendu je zdokumentovaný a odmítnutí proběhne před spuštěním na nepodporované platformě. |
| B01, B03: obsah sdist | Explicitní manifest zahrnuje test package, conftest, Docker infrastrukturu a dokumentaci. |
| B02: Docker zakrýval wheel | Obraz sestaví a instaluje wheel, odstraní zdrojový modul a používá skutečný console entry point. |
| B04: metadata licence | SPDX metadata a minimum setuptools 77 odstraňují build deprecation. Text licence je zachován. |

Model a store jsou v `records.py`, souborové rezervace v `conflicts.py`.
Veřejné role mají jména `RunManager` a `ProjectExecutionManager`; původní názvy
ExecutionManager zůstávají přímými aliasy, bez druhé implementace.
Přibyly explicitní protokoly pro raw a spravované runnery a generický návratový
typ. Mechanicky zhuštěné testy jsou znovu čitelně naformátované.

Při nedostupném požadovaném postprocessingu `wait()` vyvolá `ProcessingPending`.
Nevrací raw ExecutionOutcome jako náhražku aplikačního výsledku. Trvalý stav je
nadále obnovitelný; aplikační historický výsledek načte příslušný adaptér.

## Ověření

- Úplná sada rozbaleného sdist proti instalovanému wheelu v čistém prostředí:
  **707 passed, 3 skipped**; samostatné Docker SSH/Slurm případy: **3 passed**.
- Čisté prostředí bez ASE/ASE2SPRKKR/PSI-J: import všech modulů, raw běh,
  persistentní obnova s dodatečným procesorem, oba console entry pointy a skutečný
  remote executable prošly bez zdrojového PYTHONPATH. `pip check` prošel.
- ASE2SPRKKR execution + sprkkr: **366 passed, 4 skipped**.
- Skutečný SPR-KKR 9.7: Al SCF se dvěma iteracemi a následný DOS, včetně
  opětovného načtení výsledků, prošly.
- Guy4ASE: **793 passed, 24 failed, 19 skipped**. Všech 24 zbývajících failure ID
  odpovídá předchozí baseline; nové testy historie a retry procházejí.
- Izolovaný PEP 517 build vytvořil sdist i wheel; všechny runtime moduly ve wheelu
  odpovídají zdroji. Sdist obsahuje potřebnou testovací infrastrukturu a docs.
- Ruff pro změněný runtime a adaptéry a `git diff --check` prošly.

Lokální runtime je Python 3.14.4; Docker používá Python 3.12. CI workflow navíc
definuje Python 3.9, ale jeho hostovaný běh nebyl v této práci spuštěn. Původní
commity se nepřepisují; opravy a úprava čitelnosti testů mají samostatné commity.
Globální Python prostředí ani uživatelská projektová data se neměnila. Lokální doklady jsou v `/tmp/kkrunner-repair/`.
