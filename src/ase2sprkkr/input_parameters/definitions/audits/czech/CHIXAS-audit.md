# CHIXAS — audit současné definice

Definice `chixas.py` popisuje dichroismus XAS indukovaný polem v sériovém ASA/FREL režimu. Externí fyzikální data susceptibility jsou povinná; numerická správnost nebyla ověřena.

## Mapování a definice

INITVAR rozpoznává CHIXAS; RSP volá CHIXAS (susc_xas.f). To není alias běžného XAS: SUSC_XASCALC má výslovně účel FIELD-INDUCED MCXD. Absorbující typ se zadává TASK.ITXRAY, nikoli XAS TASK.IT. Core shell se zadává CL.

Definice obsahuje CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES a minimální SCF sekci k odmítnutí FULLPOT. RSP CHIXAS v FULLPOT explicitně zakazuje. FREL, POINTS a jeden GRID=3 jsou zvolený auditovaný spektroskopický rozsah, nikoli tvrzení o explicitním zákazu každé jiné native sítě. Native CHIXAS znovu vytváří energetickou cestu a nastavuje NEPATH=1; SPLITSS proto nenabízíme. NE>=2 zabraňuje jednobodovému GRID=3 s nedefinovanou vahou.

TASK dále obsahuje NEF (číslovaný dictionary → NEF1, NEF2, ...), CHISPN, CHIVV, ABSUNIT a PRINT. Nejsou vytvořeny smyšlené fyzikální defaults. Bez NEF a obou susceptibilit validace selže; parser/skeleton samotný lze vytvořit a potom doplnit.

## Povinná externí data a validace

- NEF pro každý potenciálový typ je l-resolved DOS na EF. Native při chybějícím NEF pro kterýkoli typ končí STOP. Počet složek nesmí překročit nativní NLMAX a musí odpovídat l-kanálům i normalizaci, kterou používá SUSC_GAMMA. To input-only validace z potenciálu nezjistí.
- CHISPN a celkové CHIVV jsou per-type susceptibility v 10^-6 cm^3/mol. SUSC_XASINIT je převádí do atomových jednotek. Total CHIVV dělí dvěma a rozděluje stejně mezi spinové kanály.
- Definice požaduje shodnou délku CHISPN/CHIVV a souvislé NEF1..NEFn, kontroluje konečné hodnoty, nezáporný DOS a ITXRAY v tomto rozsahu. Skutečný NT v potenciálu se zde nečte; shoda vstupních polí navzájem není důkaz shody s potenciálem. Čísla v testech jsou pouze grammar fixtures, nikoli materiálová data.
- Spin-resolved CHIVVUP/CHIVVDN nejsou zpřístupněny, aby nevznikaly konfliktní kombinace s celkovým CHIVV. Native parser při načítání susceptibilit dočasně přepisuje globální NT počtem nalezených hodnot. Špatné délky mohou dát nulové fallbacky a jiný stav; je nutné dodat přesně všechny potenciálové typy.

## Defaults: co se vypisuje

| Parametr | Native | ASE definice |
|---|---|---|
| ITXRAY | explicitně vyžadován | 1, vypsáno; nutno zvolit správný typ |
| CL | 2P | 2P, vypsáno |
| GRID / NE | pro prefix CHI: 5 / 30 | 3 / 180, vypsáno pro energeticky rozlišené spektrum |
| EMIN | EMINDEFAULT | −0.2 Ry, vypsáno |
| EMAX | EF | 4 Ry, vypsáno; explicitní spektroskopická volba |
| ImE | 0.01 Ry | 0.01 Ry, vypsáno |
| ABSUNIT | false | false, nevypsáno |
| NEF, CHISPN, CHIVV | NEF povinné, susceptibility bez dat mohou skončit nulové | bez defaults, vyžadovány validací |
| TASK.PRINT | vstupní IPRINT | nezadáno, nevypsáno |

Vypsané GRID/NE/EMAX jsou záměrné odchylky od prefixového CHI defaultu. Jsou vhodnější jako start spektrálního mesh než integrace pouze do EF, ale rozsah, rozlišení a šířku je nutné nastavit pro konkrétní edge. EMAX=4 není důkaz, že pokrývá požadované spektrum, a −0.2..4 Ry může zbytečně zahrnovat velký interval. Native default nelze slepě doporučit jako lepší jen proto, že jde o CHI prefix.

TASK.OP se nevypisuje: CHIXAS pak nastavuje lokální CHIORBPOL=BROOKS. Tato proměnná není automaticky totéž jako obecný MODE.OP určující korekci potenciálu. Nevypsaný ABSUNIT přesně zachovává nativní chování. RDTAUFIL začíná ve zdrojích false; restart a čtení starých textových TAU/TAUTAU cache nejsou touto definicí nabízeny.

## Další ALE

Native CHIXAS není samokonzistentní odvození susceptibilit z materiálu: nejdřív potřebuje NEF, z nich počítá GAMMA a pak používá zadané CHI při indukované spektroskopii. Propojení se skutečnými externími daty a core wavefunctions nebylo numericky ověřeno.

RDABXC/LASTITRXC, TTTFMT, surface/slab, cache a spin-resolved orbital postupy nejsou auditovány. Při POINTS native v lokální k-loop větvi přepíná IBZINT z 2 na 6; nelze to zaměňovat s uživatelským běžným TB workflow jen podle čísla. Definice nastavuje mpi=False; paralelní řetězec nebyl certifikován, nejde o prokázanou MPI chybu.
