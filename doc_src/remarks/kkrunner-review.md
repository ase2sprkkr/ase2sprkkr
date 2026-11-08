> Historical AI development note. This snapshot may be outdated and is not
> the official API or user guide. See [the current execution guide](../execution.rst)
> and [KKRunner documentation](https://ase2sprkkr.github.io/kkrunner/).

# Review vyčlenění KKRunneru — 7. 11. 2026

> Historický review: opravy a aktuální výsledky jsou v [kkrunner-repair.md](kkrunner-repair.md).

Review vychází z KKRunneru `c9ce78b`, ASE2SPRKKR `914969d01` a Guy4ASE
`7b1c634`. Posuzuje hranice knihoven, lifecycle, perzistenci, vlastnictví
prostředků, veřejné API, integraci GUI a kvalitu testů. Následující dvě lokální
opravy jsou již v pracovním stromu KKRunneru; ostatní nálezy zůstávají otevřené.

**Verdikt: rozdělení odpovědností mezi knihovnami má smysl a je převážně
provedeno správně. Implementaci zatím nepovažuji za připravenou ke sloučení.**
Největší dluh je uvnitř runneru: řízení stavů nemá dostatečně soustředěné
pravomoci a identita exekutoru supluje identitu souborového prostoru.

## Závažné otevřené nálezy

### R01 / P1 — Úspěšný rychlý běh může selhat kvůli vlastnímu submitu

Místa: `kkrunner/manager.py`, `_process_outcome`, `create_run`,
`ProjectRunner.submit`, `transition`.

`PostprocessingRunner` spustí zpracovací vlákno ještě před návratem z `submit()`.
Toto vlákno může převzít postprocessing a uložit jeho claim dříve, než vnější
`ProjectRunner.submit()` uloží backend ID. Procesor je pak úspěšný, ale následný
zápis backend ID změní záznam. Dokončení porovnává celý záznam s původním claimem
(`expected_record`) a vyvolá `Run changed during postprocessing retry`.
Vnější dokončovací cesta z toho může udělat trvalé selhání postprocessingu.

Nezávislý deterministický test toto pořadí vynutí pomocí událostí a okamžitého
falešného backendu. `wait()` skutečně vyvolá uvedenou výjimku, přestože procesor
vrátil 42. Jde o **novou regresi** v zavedeném claim mechanismu.

Architektonická příčina: dokončení zapisuje `_process_outcome` i
`ProjectRunner._sync_completion`; ochrana celým snapshotem zaměňuje ztrátu
vlastnictví za běžnou aktualizaci záznamu. `processing_attempt` se ukládá, ale
nepoužívá pro kontrolu práva dokončit pokus.

Doporučení: nejprve zveřejnit backend identitu, potom povolit processing;
ukončení pokusu kontrolovat jeho identifikátorem, stavem a integritou souborů.
Jediná cesta má rozhodovat o trvalém dokončení a následném ACK. Pouhé odstranění
porovnání snapshotů by odstranilo i potřebnou ochranu před cizím zápisem.

### R02 / P1 — Opakované run_id může spustit druhý backendový výpočet

Místo: `ExecutionManager.reserve`, zejména vyloučení `active.id == record.id`.

Dvě nová volání `create_run(..., run_id='same-id')` vytvoří dva živé objekty.
Rezervace druhého ignoruje první záznam se stejným ID. Ve stejném denním shardu
`transition` převezme již existující záznam, ale `ProjectRunner.submit` přesto
spustí druhý raw runner. Nesoulad backend ID zjistí až následný zápis, tedy po
spuštění druhého procesu. Reprodukce zaznamenala dvě backendová odeslání.

Mechanismus je **převzatý z baseline**, ale odporuje slibované jediné identitě
běhu. Napříč denními shardy navíc může být stejné ID uloženo vícekrát; existující
test `test_run_store_lookup_prefers_newest_record` toto chování výslovně očekává.

Doporučení: oddělit atomické vytvoření identity od aktualizace existujícího
záznamu. Opakované vytvoření musí skončit před backendovým submit, nebo explicitně
vrátit existující běh; připojení má dál používat `attach`. Záruka musí platit pro
celý store, nikoliv jen pro jednu datumovou složku.

### R03 / P1 — Kolize souborů se ztrácejí mezi backendy

Místo: `submission_conflicts` a `_manifest_accesses`.

První podmínka vrátí `False`, když se liší `ExecutorEndpoint`. Lokální proces
a lokální PSI/J přitom mohou zapisovat stejný soubor. Dva různé vzdálené servery
mohou publikovat stdout do stejného lokálního souboru. Obě varianty prošly
kontrolou jako nekolidující. U vzdálených běhů se porovnává pouze cílový prostor,
nikoliv současně lokální prostor stahovaných a streamovaných dat.

Problém je **převzatý z baseline**. Není důvod vracet pravidlo „SCF rezervuje
celý adresář“; úplné manifesty jsou správná cesta. Chybí správné určení prostoru,
ve kterém se jejich cesty porovnávají.

Doporučení: oddělit identitu připojení/backendu od identity souborového prostoru.
Kontrolovat lokální zápisy i vzdálené přístupy. Sdílení vzdálených mountů nelze
spolehlivě odvodit jen ze jména hostu; pro takový případ dává smysl explicitní
identifikátor souborového prostoru nebo zdokumentované konzervativní pravidlo.
To je návrhová volba, kterou je vhodné uzavřít před přepisem kolizní logiky.

### R04 / P1 — RunStore.save obchází ochranu identity i souběžných změn

Místo: `RunStore.save` oproti `RunStore.update`.

`update` používá zámek souboru, načítá aktuální záznam a kontroluje neměnnost
identity, metadat a backend ID. Veřejné `save` pouze přepíše JSON pod zámkem
jedné instance. `save(replace(record, metadata=...))` prokazatelně přepíše
„neměnná“ metadata. Stejně může zastaralý snapshot přepsat novější stav od jiné
instance nebo procesu. Atomický rename chrání před roztrženým JSON, nikoliv
před ztracenou aktualizací.

Nechráněná zapisovací cesta je **převzatá**; tvrzení o neměnných aplikačních
metadatech je nová záruka, kterou tato cesta nesplňuje.

Doporučení: jedna kontrolovaná aktualizační cesta a samostatná operace vytvoření
nového záznamu. Alternativou je `save` s povinným očekávaným snapshotem/verzí.
Volné nahrazení existujícího záznamu nesmí být běžné veřejné API.

### R05 / P1 — Cizí platný běh rozbije načtení historie Guy4ASE

Místa: `guy4ase/gui/application/persistent_jobs.py`, `reload` a
`PersistentJob.from_run_record`.

GUI čte všechny záznamy sdíleného obecného store a na každý bez výběru zavolá
`application_metadata(record)`. Jediný platný generický běh bez namespace
`ase2sprkkr` pak vyvolá `UnrecoverableProjectResultError` už při konstrukci
`PersistentJobSource`. Nezávislý test toto potvrdil.

Jde o **novou chybu integrace**: dříve byla task metadata povinná pro všechny
záznamy. Doporučení: aplikační adaptér explicitně vybere záznamy svého vlastníka.
Cizí platný záznam se přeskočí; poškozený vlastní záznam se má nadále srozumitelně
ohlásit. Obecný runner nemá získat znalost namespace ASE2SPRKKR.

## Další otevřené chyby a mezery

### R06 / P2 — stdin není zahrnutý do čtecího manifestu

`FileManifest.from_spec` zahrne stdin do kontroly kolizí uvnitř jednoho specu,
ale neuloží jej do `inputs`. Druhý běh proto může rezervovat zápis do jeho stdin
souboru. Nezávislý test to potvrdil. Kód je **převzatý z baseline**.

Při nápravě je třeba zachovat rozdíl mezi primárním `input_file` a ostatními
čtenými soubory: dnes `FileManifest.input_file` znamená první prvek `inputs`.
Prosté přidání stdin na začátek by tento význam změnilo. Manifest potřebuje
úplnou množinu čtených cest a jednoznačnou identitu primárního vstupu.

### R07 / P2 — Validace vnořených dat připouští neplatné záznamy

`Execution.__post_init__` konvertuje slovníky, ale neověří výsledné typy
`endpoint` a `directory`. Záznam s `execution.endpoint=[]` se úspěšně načte;
později spadne běžný přístup k endpointu. Reprodukce toto potvrdila.
Jde o **převzatou vadu**, kterou přísnější validace RunRecord nepokryla.
Stejný princip je třeba systematicky zkontrolovat u dalších vnořených typů.

Validace má proběhnout na hranici deserializace a skončit `CorruptRunRecord`
s cestou k vadnému souboru. Oprava pouze obalením všech pozdějších výjimek by
chybu posouvala dál od zdroje.

### R08 / P2 — Nově otevřené vzdálené streamy zůstávají po chybě submit otevřené

`_RemoteRun._submit_async` materializuje pathname stdout/stderr na otevřené
handly, potom zavolá `_submit`. Při výjimce nevznikne monitor, který je běžně
uzavírá. Reprodukce se simulovanou ztracenou odpovědí ověřila otevřený handle.
Jde o **novou vadu** podpory pathname streamů.

Handly získané runnerem potřebují lokálního vlastníka a úklid při neúspěšném
předání monitoru. Nejasnost výsledku vzdáleného submitu vyžaduje zachovat
obnovitelný záznam a soubory; sama nevyžaduje držet otevřený lokální deskriptor.

### R09 / P2 — Změna defaultní cesty skryje dosavadní definice projektů

`default_projects_directory` nově používá `user_state_dir('kkrunner')`, baseline
používala `user_state_dir('ase2sprkkr')`. Formát definic projektů zůstává verze 2,
ale ASE i GUI nově otevřou prázdnou jinou složku. Dosavadní projektové UUID pak
nelze dohledat. Změna není popsána jako záměrná nekompatibilita.

Jde o **novou změnu chování**, potvrzenou porovnáním zdroje. Zákaz migrací
vývojových RunRecord není rozhodnutím zahodit dostupnost definic projektů.
Možnosti: aplikace explicitně předá původní umístění, nebo se provede jasně
popsaný přesun projektových definic. Obecný runner nemusí znát starý namespace.

### R10 / P2 — GUI nemá explicitní opakování selhaného postprocessingu

ASE poskytuje `retry_project_postprocessing`; GUI zpřístupňuje pouze retry
downloadu. Failed/postprocessing je terminální a neattachovatelný, automatická
obnova se týká jen aktivního stavu postprocessing. Kód GUI procesorový retry
nikde nevolá. Je to **převzatá funkční mezera**, ale nové zadání jeho ověření
v Guy4ASE výslovně vyžaduje. Není korektní tento bod považovat za dokončený.

## Již opraveno při review

1. **Nový deadlock při cancel + neúspěšném submitu.** `submit` držel
   `_submission_lock`, jeho chybová cesta přes `_release_after_failure` četla
   `self.state`, která po přijetí cancel znovu získávala stejný nereentrantní
   zámek. Oprava čte přímo stav obaleného runneru. Deterministický test ověřuje
   ukončení obou vláken, původní výjimku a uzavření prostředků.
2. **Převzaté chybné ošetření JSON kořene.** Pro `[]`, `null` a číslo původně
   unikala `AttributeError` při `data.get`. Dekodér nyní explicitně vyžaduje
   JSON objekt a vrací `CorruptRunRecord` s cestou souboru. Tři parametrické
   regresní případy jsou součástí opravy.

Opravy jsou v `kkrunner/src/kkrunner/base.py`, `manager.py` a novém
`kkrunner/tests/test_review_fixes.py`; jsou zatím necommitnuté.

## Architektura a kvalita kódu

**Co zachovat:** jednosměrnou závislost ASE2SPRKKR → KKRunner, TaskFile jako
obecný deklarativní typ, opaque namespaced metadata, callback dodaný klientem,
aplikační rekonstrukci výsledků, jedinou implementaci remote služby a práci
s úplnými manifesty. V produkčních zdrojích runneru nejsou aplikační importy.
Příprava výpočtu nadále zůstává na aplikační straně. Další distribuci ani
globální registr procesorů výsledků není potřeba zavádět.

**Co zpřesnit před stabilizací veřejného API:**

- `kkrunner.ExecutionManager` a `kkrunner.manager.ExecutionManager` mají stejné
  jméno a odlišný kontrakt: první vyžaduje projekt, druhý nikoliv. Doporučuji
  buď jeden manager s explicitní projektovou politikou, nebo dva výstižné názvy
  (`RunManager` a `ProjectExecutionManager`). Nejde o důvod k přepisu backendů.
- `Runner` Protocol nepopisuje velkou část schopností, které orchestrátor
  skutečně používá: `files`, `execution`, `return_code`, `start`,
  `prepare_submission`, `stop_observation`. Časté `getattr` schovává chybějící
  kontrakt. Vymezil bych minimum raw runneru a explicitní managed/attach
  schopnosti; aplikační výsledek lze typovat přes `Runner[T]`.
- `wait()` může po `ProcessingPending` vrátit raw `ExecutionOutcome`, přestože
  byl dodán procesor vracející aplikační objekt. Totéž nastane při zjištění, že
  jiný pozorovatel již dokončil záznam. Stav trvalého běhu, konec pozorování
  a dostupnost aplikačního výsledku potřebují samostatně popsaný kontrakt.
  Například explicitní výjimka při nedostupném zpracovaném výsledku lépe chrání
  stávající návratový typ kalkulátoru; alternativou je jednotný výsledek
  pozorování. Tuto změnu API jsem v review nevybral a neprovedl.
- `manager.py` soustřeďuje model, JSON store, konflikty, snapshoty, stahování,
  retry i běhový obal. Samotných přibližně 1300 řádků není chyba; problém jsou
  překrývající se pravomoci popsané v R01–R04. Rozdělení má následovat tyto
  hranice a centralizaci zápisů, samotné rozřezání souboru nic nevyřeší.
- Otevírání pathname streamů vrací pouze změněný spec a vlastnictví spoléhá
  na pozdější vedlejší efekt monitoru. Vhodnější je explicitně předávaná sada
  získaných prostředků s úklidem v chybové cestě (R08).

**Kvalita testů a údržba:**

- Přesun části testů zbytečně zhoršil čitelnost: dlouhé jednořádkové konstruktory,
  mechanické výrazy `str('DOS').upper()` a ztráta původního členění. Například
  `test_run_store_merge.py` byl v baseline výrazně přehlednější. Obnovil bych
  původní formát a komentáře a ponechal pouze potřebné změny importů a metadat.
- Některé názvy nesouhlasí s ověřovaným kontraktem: například
  `test_unknown_observation_does_not_fire_completion_callback` ve skutečnosti
  očekává zavolání callbacku. Tato nejasnost je již v baseline a ukazuje potřebu
  výslovně pojmenovat konec pozorování versus konec výpočtu.
- Široká původní sada dobře prověřuje jednotlivé mechanismy, ale nepokryla
  pořadí „procesor před návratem submitu“, opakované ID, různé backendy nad
  společným souborem ani cizí namespace v GUI. Tyto kontraktní testy jsou pro
  obecnou knihovnu důležitější než další testy pouhého přesunu importů.
- Nový projekt zatím nemá CI workflow s runtime maticí deklarovaných Pythonů
  a izolovanou instalací bez ASE. Manuální validační doklady jsou užitečné,
  ale další změny potřebují opakovatelnou automatickou kontrolu.
- První commit runneru spojuje přenos přibližně 19 tisíc řádků s behaviorálními
  změnami. Před předložením PR by byla přehlednější historie, která oddělí
  funkčně ekvivalentní přesun od změn recovery a lifecycle. Historii jsem při
  review nepřepisoval.

## Ověření tohoto review

Nezávislé sondy jsou v `/tmp/kkrunner-self-review/test_review.py` a
`test_gui_review.py`; úplný poslední výstup je
`/tmp/kkrunner-self-review/results-after-fixes.txt`. Tyto soubory jsou dočasné
review artefakty, nejsou součástí běžné testovací sady.

Před dvěma opravami selhalo všech 13 parametrických případů. Po opravách:
**4 passed, 9 failed**. Devět zbývajících selhání pokrývá osm otevřených
nálezů R01–R08 (R03 má dvě varianty). R09 a R10 byly ověřeny čtením zdroje
a porovnáním s baseline. Pořadí vláken se vynucuje událostmi; testy nepoužívají
pravý scheduler. Remote chybová cesta používá simulovanou odpověď.

```sh
PYTHONPATH=/home/logik/skola/kkr/kkrunner/src:/home/logik/skola/kkr/guy4ase:/home/logik/skola/kkr/ase2sprkkr/src \
  python -m pytest -q /tmp/kkrunner-self-review/test_review.py \
  /tmp/kkrunner-self-review/test_gui_review.py
```

Čtyři nové trvalé regresní případy a související existující testy na skutečném
pracovním stromu: **34 passed, 3 deselected**. Tři vynechané testy spouštějí
supervizované procesy; nebyly součástí tohoto cíleného běhu. Ruff pro oba
změněné produkční moduly a nový test prošel, stejně jako `git diff --check`.
Velké sady, SSH a Slurm z předchozí implementace jsem v tomto review znovu
nespouštěl; jejich tehdejší průchod nevyvrací zde reprodukované mezery.

Další postup: nejdříve uzavřít kontrakt identity, vlastnictví dokončení
a souborových prostorů (R01–R04), potom opravit aplikační integraci a zbývající
lokální chyby. Následovat má návrat k nezávislým sondám a plná regresní kontrola.
