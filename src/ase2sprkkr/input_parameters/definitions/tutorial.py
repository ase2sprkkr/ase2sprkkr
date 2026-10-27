"""Complete native SPR-KKR 9.7 TUTORIAL dispatch and input parameters."""
import numpy as np
from ._sections import CONTROL, TASK, TAU, ENERGY, SITES, MODE, CPA, STRCONST
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError
from ..task_definition import TaskDefinition
from ..input_parameters_definitions import (
    InputParametersDefinition, InputSectionDefinition as Section,
    InputValueDefinition as V,
)

TUTORIAL_TASKS = {
    "YLM": "Plot a real spherical harmonic, with positive and negative parts separated",
    "JLX": "Plot spherical Bessel j_l(x) and Neumann n_l(x) functions",
    "JLX-ASYMPT": "Compare spherical Bessel/Neumann functions with small- or large-x asymptotes",
    "CHR-POT": "Plot stored potential, charge/spin density and FULLPOT shape functions",
    "POT-DECOMP": "Decompose an ASA potential; input XC must be VBH, MJW or VWN",
    "POT-XC": "Compare VBH, MJW, VWN and PBE XC potentials on the same stored ASA density",
    "RWF-L": "Compare regular scalar-relativistic radial wave functions for different l at EMIN",
    "RWF-E": "Compare regular scalar-relativistic radial wave functions across energies for one l",
    "CORE-ALL": "Plot nonmagnetic core states of potential type 1; ASA with occupied core shells required",
    "DOS-FEG": "Single-site DOS with constant V=0.0001 Ry and B=0 (free-electron demonstration)",
    "DOS-SS": "Single-site DOS using the supplied potential",
    "DOS-CLU": "DOS using an explicitly configured real-space cluster",
    "DOS-BULK": "DOS using an explicitly configured periodic BZ integration",
    "GEN": "General DOS/observable calculation with the supplied potential and TAU settings",
}


def _subtask(option):
    return option._container._container.TASK["TUT-TASK"]()


def _tutorial_ne(option):
    task = _subtask(option)
    return [10 if task.startswith("RWF-") else 300 if task.startswith("DOS-") or task == "GEN" else 1]


def _tutorial_grid(option):
    task = _subtask(option)
    return ["3" if task.startswith(("RWF-", "DOS-")) or task == "GEN" else "0"]


def _tutorial_emin(option):
    return 1e-7 if _subtask(option).startswith("RWF-") else -.2


def _tutorial_ime(option):
    task = _subtask(option)
    return .01 if task.startswith("DOS-") or task == "GEN" else 0.0


def _tutorial_mode(option):
    return "SREL" if _subtask(option).startswith("RWF-") else "FREL"


def _tutorial_bzint(option):
    tau = option._container
    return "CLUSTER" if tau.CLUSTER() or tau.MOL() else "POINTS"


def _tutorial_kkrmode(option):
    return "STANDARD"


def _validate_tutorial(p, _container, why):
    if why == "set":
        return None
    t, e = p.TASK, p.ENERGY
    task = t["TUT-TASK"]()
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if task == "YLM" and abs(t.M()) > t.L():
        error("YLM requires abs(M)<=L; native silently clamps M.")
    if task in ("JLX", "JLX-ASYMPT"):
        if not np.isfinite(t.XMIN()) or not np.isfinite(t.XMAX()) or t.XMIN() < .001 or t.XMAX() <= t.XMIN():
            error("JLX requires finite XMAX>XMIN>=0.001.")
        if task == "JLX-ASYMPT" and t["X-RANGE"]() == "SMALL" and t.XMIN() != .001:
            error("JLX-ASYMPT SMALL forces XMIN=0.001.")
    if task in ("POT-DECOMP", "POT-XC", "CORE-ALL", "RWF-L", "RWF-E") and p.SCF.FULLPOT():
        error(f"{task} requires ASA in this definition; native FULLPOT is unsupported or inconsistent.")
    if task.startswith("RWF-"):
        if p.MODE.MODE() not in ("NREL", "SREL", "SP-SREL"):
            error("RWF uses WAVFUN_READ_SRA; choose NREL, SREL or SP-SREL.")
        if len(e.NE()) != 1 or e.NE()[0] < 2 or e.SPLITSS():
            error("RWF requires one mesh with NE>=2 and no SPLITSS (native EPATH(3)).")
        if e.EMINEV.is_set() or e.EMAXEV.is_set():
            error("RWF rereads absolute EMIN/EMAX in Ry and ignores EMINEV/EMAXEV.")
        if e.EMIN() is None or e.EMAX() is None or not np.isfinite(e.EMIN()) or not np.isfinite(e.EMAX()) or e.EMAX() <= e.EMIN():
            error("RWF requires finite absolute EMAX>EMIN.")
        nl = p.SITES.NL()
        if task == "RWF-E" and nl is not None and t.L() >= max(nl):
            error("RWF-E requires L<NL; native clamps L to NL-1.")
    if task.startswith("DOS-") or task == "GEN":
        if len(e.NE()) != 1 or len(e.GRID()) != 1 or e.SPLITSS():
            error("Use a single TUTORIAL DOS mesh without SPLITSS.")
        if e.GRID()[0] in ("2", "3") and e.NE()[0] < 2:
            error("DOS GRID=2/3 requires NE>=2.")
        if not np.isfinite(e.ImE()) or e.ImE() < 0:
            error("DOS ImE must be finite and nonnegative.")
        for name in ("EMIN", "EMAX", "EMINEV", "EMAXEV"):
            value = e[name]()
            if value is not None and not np.isfinite(value):
                error(f"ENERGY.{name} must be finite.")
    cluster = p.TAU.CLUSTER() or p.TAU.MOL() or p.TAU.BZINT() == "CLUSTER"
    if task in ("DOS-FEG", "DOS-SS"):
        if cluster or p.TAU.KKRMODE() != "STANDARD":
            error("DOS-FEG/DOS-SS require STANDARD KKR without a cluster to use native tau=t.")
    if task == "DOS-CLU":
        if not cluster:
            error("DOS-CLU needs TAU.CLUSTER, MOL or BZINT=CLUSTER; its name alone does not configure a cluster.")
        if p.TAU.CLURAD() is None and p.TAU.NSHLCLU() is None:
            error("DOS-CLU requires TAU.CLURAD or NSHLCLU.")
    if task == "DOS-BULK" and (cluster or p.TAU.BZINT() not in ("POINTS", "WEYL")):
        error("DOS-BULK requires periodic POINTS or WEYL integration.")
    return errors or None


def input_parameters_definition():
    task = TASK("TUTORIAL", add=[
        V("TUT-TASK", Keyword(TUTORIAL_TASKS), "YLM", is_always_added=True,
          info="Select the native educational or diagnostic calculation."),
        V("L", Integer(min=0), 2, info="YLM angular momentum; RWF-E selected l (must be below NL)."),
        V("M", int, 0, info="YLM real harmonic index, abs(M)<=L."),
        V("LMAX", Integer(min=0, max=4), 2, info="Largest l in JLX plots."),
        V("NX", Integer(min=2), 300, info="JLX mesh point count."),
        V("XMIN", .001, info="JLX lower bound; SMALL asymptote forces .001."),
        V("XMAX", 5.0, info="JLX upper bound."),
        V("X-RANGE", Keyword("SMALL", "LARGE"), "SMALL", info="JLX-ASYMPT limit; ignored by plain JLX."),
        V("ITSEL", Integer(min=1), 1, info="Potential type used by RWF-L/RWF-E; must exist in POTFIL."),
        V("NORMALISATION", Keyword("1", "MST"), "1",
          info="RWF: 1 normalizes the plotted real radial part; MST retains multiple-scattering normalization."),
    ])
    scf = Section("SCF", [V("FULLPOT", False,
        info="Native potential format takes precedence; POT-DECOMP/POT-XC/CORE-ALL and audited RWF require ASA.")],
        is_optional=True, is_expert=True)
    mode = MODE.copy(
        add=[
            MODE["MODE"].copy(type=Keyword("NREL", "SREL", "SP-SREL", "FREL"), default_value=_tutorial_mode),
        ],
        is_expert=False, is_optional=False,
    )
    tau = TAU.copy(
        add=[
            TAU["BZINT"].copy(default_value=_tutorial_bzint),
            TAU["KKRMODE"].copy(default_value=_tutorial_kkrmode, is_always_added=True),
        ],
    )
    energy = ENERGY(emin=(_tutorial_emin, "Absolute lower energy in Ry; RWF ignores relative eV bounds", None),
        emax=(1.0, "Absolute upper energy in Ry", None),
        defaults={"GRID": _tutorial_grid, "NE": _tutorial_ne, "ImE": _tutorial_ime})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParametersDefinition(
        "tutorial",
        [CONTROL("TUTORIAL"), tau, energy, SITES, task, scf, mode, CPA, STRCONST],
        validators=_validate_tutorial,
        info="Educational and diagnostic KKR calculations selected by TASK.TUT-TASK: real spherical harmonics "
             "(YLM), spherical Bessel/Neumann functions and asymptotes (JLX, JLX-ASYMPT), stored potential and "
             "densities (CHR-POT), potential decomposition (POT-DECOMP), XC comparison on one density (POT-XC), "
             "radial valence wave functions across l or energy (RWF-L, RWF-E), nonmagnetic core states of "
             "type 1 (CORE-ALL), and free-electron, single-site, cluster or bulk DOS (DOS-FEG, DOS-SS, "
             "DOS-CLU, DOS-BULK). GEN runs the general observable calculation. All branches use native "
             "kkrgen startup and require a potential file. Several branches have native limitations; see their descriptions.",
        description="POT-DECOMP requires ASA with XC=VBH/MJW/VWN. POT-XC compares VBH/MJW/VWN/PBE on the same density, "
                    "without separate SCF runs. CORE-ALL uses only type 1 and forces B=0. RWF uses scalar-relativistic "
                    "wave functions and absolute Ry bounds. DOS-CLU requires explicit TAU cluster geometry; DOS-BULK "
                    "and GEN retain the supplied TAU settings. Native diagnostic branches are not certified for MPI."
    )
    schema._task_definition = TaskDefinition("TUTORIAL", experimental=True, executable="kkrgen", mpi=False)
    return schema
