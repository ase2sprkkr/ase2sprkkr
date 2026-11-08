> Historical AI development note. This snapshot may be outdated and is not
> the official API or user guide. See [the current execution guide](../execution.rst)
> and [KKRunner documentation](https://ase2sprkkr.github.io/kkrunner/).

# Druhé review KKRunneru včetně sestavení — 7. 11. 2026

> Historický review: opravy a aktuální výsledky jsou v [kkrunner-repair.md](kkrunner-repair.md).

Posuzován je KKRunner `c9ce78b` se dvěma necommitnutými opravami z prvního
review. ASE2SPRKKR je na `914969d01`, Guy4ASE na `7b1c634`. V tomto kole jsem
produkční kód neměnil. Předchozí report je v `kkrunner-review.md`.

**Verdikt: základní build a instalovatelný runtime fungují. Zdrojová distribuce
má neúplný obsah pro testování a dál zůstávají otevřené chyby lifecycle.
Samotný průchod buildem zatím nestačí k doporučení změnu sloučit.**

## Nově potvrzené nálezy

### B01 / P2 — Ze zveřejněného sdist nelze spustit dodanou testovací sadu

Místo nápravy: `kkrunner/pyproject.toml:22` a explicitní pravidla obsahu sdist
(nyní chybí `MANIFEST.in` nebo jejich ekvivalent).

Standardní `python -m build` vytvoří platný sdist a následně z něj wheel.
Do sdist setuptools automaticky přibalí `tests/test_*.py`, ale nepřibalí
`tests/__init__.py`, `tests/conftest.py` ani `tests/docker/`.

Po rozbalení sdist sběr testů skončí dvěma chybami:

```text
tests/test_remote_stream_protocol.py: attempted relative import with no known parent package
tests/test_remote_wildcards.py: attempted relative import with no known parent package
645 tests collected, 2 errors
```

Ze stejného zdrojového snapshotu před balením se sesbírá **685 testů bez chyb**.
Chybějící conftest navíc odstraňuje automatickou izolaci executor poolu
a fixture Docker obrazu. Samotné přidání `__init__.py` tedy nestačí.

Doporučení: explicitně přibalit celou potřebnou testovací infrastrukturu,
například pravidlem pro `tests/**/*.py`, Dockerfile, shellové skripty a `.conf`
soubory. Do automatické kontroly zařadit testy rozbaleného sdist. Nejde o chybu
runtime wheelu; tento nález se týká ověřitelnosti zdrojové distribuce.

### R11 / P1 — Výchozí SSH mapování obejde kontrolu kolizí i na stejném endpointu

Místa: `RemoteEndpointExecutor.submission_directory` (`remote_endpoint.py:458`),
`ExecutionManager.reserve` a následná příprava v `ProjectRunner.submit`.

První běh po stagingu uloží například cílovou cestu
`/home/tester/case/stdout`. Další běh se stejnou výchozí konfigurací vstupuje
do rezervace ještě s cestou `~/case/stdout`. Kontrola porovnává řetězce,
neexpanduje `~` a konflikt nenajde. Až po rezervaci se jméno přeloží na stejnou
absolutní cestu; další kontrola konfliktů už nenásleduje.

Nezávislá sonda používá skutečné sestavení specu pro `SshTransport`, žádné SSH
spojení. První záznam reprezentuje uložený stav po stagingu, druhý běžnou novou
rezervaci. Dvě rezervace pro tentýž výstup projdou. Jde o další projev R03
z prvního review, tentokrát u **stejného endpointu a výchozího remote_root**.
Podstatné mechanismy jsou převzaté z baseline.

Doporučení: oddělit read-only vyřešení mapování od stagingu. Nejdříve získat
kanonické cílové cesty, potom je atomicky rezervovat a teprve potom přenášet
nebo otevírat výstupy. Přesun celé stávající přípravy před rezervaci by dovolil
vedlejší efekty před kontrolou konfliktů a nebyl by správnou opravou.

### R12 / P1 — Chybějící povinný lokální výstup se vydává za úspěšný běh

Místo: `file_availability.py:14`, `discover_local_files`, a následné zachycení
manifestu v `ExecutionManager.complete`.

Lokální discovery bez ohledu na `FileEntry.optional` označí chybějící výstup
jako `skipped`. Zachycení verzí uloží pouze soubory, které existují. Výpočet
s návratovým kódem 0 tak může skončit `FINISHED`, přestože nevytvořil povinný
výstup; historický manifest bez zachycených souborů pak hlásí `CURRENT`.

Sonda nad managerem porovnala jinak totožnou deklaraci s `optional=False`
a `optional=True`. Obě skončily jako úspěšné. Na vzdálené straně se přitom
chybějící povinný výstup za chybu považuje (`RemoteRunner._sync_output_files`).
Je to **převzatá chyba** a nekonzistence mezi backendy.

Doporučení: jednotná validace povinných výstupů před úspěšným dokončením;
volitelný chybějící soubor může být přeskočen. Rozlišit chybu vzniklého artefaktu
od chyby přenosu, aby chybějící lokální výstup nebyl označen jako retry downloadu.
Existující test `test_missing_optional_files_are_skipped_without_creating_them`
navíc deklaruje soubor bez `optional=True` a tím nekonzistenci zakrývá.

### R13 / P2 — Dokončení přes refresh neuklidí backendový outcome receipt

Místa: `ExecutionManager.update` (`manager.py:675`), `retry_postprocessing`
a jediná cesta volající `executor.acknowledge` v `ProjectRunner._sync_completion`.

Refresh nyní správně umí uložit `FINISHED` pro běh bez procesoru. Provede však
pouze `acknowledge_submission`, nikoliv potvrzení terminálního výsledku.
Záznam se následně vyřadí z aktivních běhů. Backendový outcome receipt zůstává
bez běžné cesty k úklidu. Podobně skončí explicitní retry přerušeného
postprocessingu, které nepoužívá ProjectRunner. Ani selhání ACK v běžném
dokončení nemá jasnou následnou cestu opakování.

Sonda s backendem evidujícím ACK ověřila uložení `FINISHED` a nulový počet
potvrzení výsledku. Absence společného úklidu je převzatá; nová možnost
dokončení přes refresh ji zpřístupňuje dalšímu běžnému scénáři.

Doporučení: společný krok potvrzení výsledku až po trvalém dokončení, volaný
všemi dokončovacími cestami. Selhání ACK nesmí přepsat úspěšný RunRecord, ale
potřebuje idempotentní pozdější reconciliation. Další autoritativní lifecycle
kvůli tomu není potřeba zavádět.

## Další návrhové a build připomínky

### B02 / P2 — SSH testovací obraz může zakrýt chyby instalovaného balíčku

Dockerfile sice projekt instaluje, ale `tests/docker/entrypoint.sh:30` zároveň
vytváří wrapper `kkrunner`, který přidává `/workspace/kkrunner/src` na
`PYTHONPATH`. Při jeho použití se může spouštět zdroj místo instalovaného wheelu.
Dosavadní důkaz „bez ASE na cílovém hostu“ tím zůstává platný; důkaz správnosti
obsahu wheelu z něj ale neplyne.

Doporučení: pro distribuční smoke test vytvořit wheel v builder fázi a do
čistého runtime prostředí dát pouze wheel a deklarované závislosti. Spouštět
skutečné console entry pointy bez vloženého source PYTHONPATH. V tomto review
jsem takový test mimo Docker samostatně provedl a prošel.

### B03 / P3 — Sdist neobsahuje dokumentaci, na kterou odkazuje README

Chybí `docs/architecture.md` a `docs/validation.md`. README je přitom obsahuje
jako odkazy. Doporučuji přibalit docs do sdist a později doplnit publikované
adresy dokumentace a repozitáře do `project.urls`. Dokumentace nemusí být
součástí importovatelného Python balíčku.

### B04 / P3 — Backend upozorňuje na zastaralý formát licence

`pyproject.toml:11` používá `license = {file = 'LICENCE'}`. Setuptools 84.2.1
při skutečném buildu opakovaně vypisuje deprecation warning s požadavkem
na přechod na SPDX pole. Dnes build projde. Úpravu metadat je vhodné provést
spolu s odpovídajícím minimem setuptools; nový formát nelze prostě zavést
a ponechat deklarovaný backend `>=68` bez ověření podpory. Text licence jsem
neměnil ani právně neposuzoval.

### A01 — Kontrakt callbacku dovoluje nečekaný deadlock při reentranci

`finish_monitor` zavolá callback před nastavením `_done`. Pokud callback zavolá
`run.wait()`, čeká na dokončení sebe sama. Sonda nad skutečným LocalRunnerem
deadlock potvrdila. Jde o převzaté chování; existující testy explicitně vyžadují,
aby běžní čekající dostali výsledek až po callbacku, takže nelze jen prohodit
pořadí a ignorovat propagaci výjimek callbacku.

Tuto sondu hodnotím jako návrhovou past, nikoliv jako důkaz porušení výslovně
slíbené reentrance. Doporučení: buď bezpečná reentrance pro monitorovací vlákno,
nebo okamžitá srozumitelná výjimka a popsané omezení. Nekonečné čekání je pro
veřejné API špatně diagnostikovatelné.

### A02 — Podpora operačních systémů potřebuje přesnější vymezení

Core se zabalí jako čistý Python wheel `py3-none-any`, ale spravovaný lokální
supervizor používá `/proc/<pid>/stat`, Unix sockety, `pass_fds` a procesní skupiny.
Na systému bez Linuxového `/proc` nemůže stávající identifikace procesu fungovat.
Windowsové větve ve file locking a durability samy nezaručují přenositelnost
spravovaného lokálního backendu.

Universal wheel tag je pro čistý Python v pořádku; není důvod uměle vyžadovat
kompilaci nebo omezit vzdálené klienty na Linux. Dokumentace ale má přesně říct,
které backendy vyžadují Linux a které běží i jinde. V tomto review jsem netestoval
jiný operační systém ani plnou runtime matici Pythonů od deklarované 3.9.

## Co bylo skutečně ověřeno při sestavení

Pracoval jsem s kopií aktuálních souborů v `/tmp/kkrunner-review2/source`.
Sestavení nezapisovalo egg-info ani build adresáře do pracovního repozitáře.

| Kontrola | Výsledek |
| --- | --- |
| `python -m build --no-isolation` | sdist i wheel vytvořeny |
| Standardní izolovaný `python -m build` | sdist i wheel z tohoto sdist vytvořeny |
| Porovnání Python souborů v obou wheelech | všech 36 modulů shodných se zdrojem |
| Obsah wheelu | žádné přibalené ASE, ASE2SPRKKR ani tests |
| Metadata | core závislost platformdirs, volitelné psij a test extras |
| Instalace wheelu do nového venv | prošla, `pip check` bez porušených závislostí |
| Import všech modulů bez ASE a PSI/J | prošel |
| Raw lokální výpočet z instalovaného wheelu | prošel |
| Persistentní běh a dodatečný procesor po novém manageru | prošel |
| `kkrunner-remote` jako skutečný podproces přes LocalTransport | prošel |
| Oba console entry pointy `--help` | prošly |
| Sběr testů ze zdrojového checkoutu | 685 testů, bez chyb |
| Sběr testů z rozbaleného sdist | 645 testů, 2 chyby |

Čistý venv má pouze pip, platformdirs a KKRunner. Smoke test běžel mimo checkout,
s odstraněným `PYTHONPATH`, ověřil nepřítomnost `ase` a `ase2sprkkr` a načtení
modulu ze site-packages venv. Instalace použila lokálně dostupný wheel
platformdirs; standardní izolovaný build stáhl svůj deklarovaný backend.
Python byl 3.14.4 a setuptools 84.2.1. První síťový pokus a první pokus se socketem
omezil sandbox; po povoleném opakování uspěly. Tyto první chyby nejsou vady buildu.

Nové nezávislé sondy: **4 failed, 1 passed**. Tři selhání potvrzují R11–R13;
čtvrté ukazuje popsanou reentrantní past A01. Pro R12 prošla kontrolní varianta
s volitelným souborem. Sondy neprováděly skutečné SSH spojení ani schedulerový běh.

Znovu jsem spustil všech 13 sond prvního review: výsledek zůstává
**9 failed, 4 passed**. Předchozí opravy tedy zůstávají funkční, otevřené nálezy
nebyly mezitím odstraněny. Celou runtime sadu, skutečné SSH ani Slurm jsem v tomto
kole znovu nespouštěl.

Doklady jsou v `/tmp/kkrunner-review2/`: `isolated-build.log`, `build.log`,
`artifact-check.json`, `install.log`, `smoke.py`, `smoke.log`,
`sdist-collection.log`, `checkout-collection.log`, `test_second_review.py`,
`probes.log` a `first-review-recheck.log`. Jde o dočasné lokální artefakty.

## Doporučený výsledný build proces

1. Z čistého zdrojového stromu vytvořit sdist a wheel standardním frontendem.
2. Ověřit obsah obou artefaktů; sdist musí obsahovat to, co jeho testy a README
   potřebují. Wheel musí obsahovat všechny runtime moduly a správná metadata.
3. Instalovat wheel do nového prostředí bez source cest a bez ASE; spustit import,
   oba CLI entry pointy, lokální supervisor a obnovu běhu.
4. Z rozbaleného sdist spustit testy proti instalované distribuci; minimum je
   korektní collection, následovat má běžná sada a volitelné integrační testy.
5. Automatizovat tento řetězec v CI a přidat ověřenou runtime matici podporovaných
   Pythonů. Docker test distribuce má používat pouze instalovaný artefakt.

Změna build backendu není potřeba. Současný setuptools se src layoutem je
pro tento čistě Pythonový projekt přiměřený; schází především výslovný obsah
zdrojové distribuce a automatické testování výsledných artefaktů.
