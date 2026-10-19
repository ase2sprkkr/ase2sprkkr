# SIGMA-XAS — SPR-KKR 9.7

Kategorie pro budoucí společný audit: **rozbité z hlediska dosažitelnosti běžným inputem**. Není doložen obecně vypnutý placeholder: existuje skutečná implementace, ale chybí její výběr v parseru.

## Doložený problém
rsp.f:16 popisuje elektrickým polem indukovaný XMCD; ř.941–943 obsahují dispatch TASK='SIGMA-XAS ' do SIG_XAS. sigme.f:499 obsahuje odpovídající výběr maticových elementů. SIG_XAS v sig_xas.f není pouhý stub.

INITVAR však pro task vodivosti čte pouze celé klíčové slovo SIGMA (initvar.f:351–355) a při něm nastaví TASK='SIGMA     '. Ve všech dodaných init souborech chybí výběr nebo přiřazení SIGMA-XAS. SECTION_FIND_KEYWORD (scan_input.f:144–201) hledá přesný substring s mezerou před i za KW. SIGMA-XAS tedy není rozpoznáno jako SIGMA prefix. Přidání samostatného tokenu SIGMA vybere standardní SIGMA, nikoliv SIGMA-XAS. CONTROL.ADSI není přímé přiřazení interního TASK a problém tím nelze obejít.

KKRCHI předává interní TASK do RSP (kkrchi.f:1425–1427). Interní větev SIGMA-XAS je proto nedosažitelná obvyklým inputem této verze. Samotné přidání volby do ASE2SPRKKR by neznamenalo funkční podporu.

## Co routine skutečně očekává (orientační mapa, ne certifikovaná definice)
V TASK čte SIGSPROJ, GRID (default3, povoleno3/6), NEA50, EMIN defaultEF, povinné EMAX, IME .01, NEB1, DEB .05, EPS .02, IPOLEFIELD2, IT nebo IQ, TAUCORE. Tyto volby jsou v sekci TASK, nikoliv automaticky v obvyklé ENERGY. sig_xas.f:293 požaduje první proudový operátor1 nebo9; ř.299 vyžaduje IBZINT2/4/6. Výběr core shellu a lifetime musí být samostatně auditován před skutečnou podporou. Uvedené defaults jsou zdrojově přečtené, nikoliv doporučení funkčního inputu.

Pouhé doplnění parseru by umožnilo vstoupit do existující rutiny, ale nepotvrzuje její numerickou správnost ani dostatečnost initialization. Neodstraňujeme guardy a nevytváříme falešný alias na SIGMA/XAS.

## Výstup této kontroly
Definiční soubor ani Python regresní test se nevytváří: input-only podpora není v daném solveru dosažitelná. Nález je doložen parserem a dispatch cestou, numerický výpočet nebyl spuštěn. K implementaci je nejprve potřeba native výběr SIGMA-XAS, audit jeho initialization a numerický referenční případ.
