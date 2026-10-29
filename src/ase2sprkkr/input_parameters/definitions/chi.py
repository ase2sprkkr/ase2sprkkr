"""Serial ASA magnetic susceptibility; native iterations are disabled."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ._sections import CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V, InputSectionDefinition as Section
from ..task_definition import TaskDefinition


def _validate_chi(p, _values, why):
    if why == "set":
        return None
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if p.MODE.MODE() in ("NREL", "SREL", "SP-SREL"):
        error("This CHI definition covers FREL only.")
    if p.SCF.FULLPOT():
        error("Public SPR-KKR does not support CHI with FULLPOT.")
    if p.TAU.BZINT() != "POINTS":
        error("This periodic CHI definition requires BZINT=POINTS.")
    if p.ENERGY.GRID().tolist() != ["5"] or len(p.ENERGY.NE()) != 1 or p.ENERGY.SPLITSS():
        error("This CHI definition uses one GRID=5 contour without SPLITSS.")
    if p.MODE.OP() not in (None, "NONE", "BROOKS"):
        error("This CHI definition covers OP=NONE or BROOKS only.")
    return errors or None


def input_parameters_definition():
    energy = ENERGY(defaults={"GRID": [5], "NE": [30], "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    mode = MODE.copy([V("PRINT", Integer(min=0), is_optional=True,
                         info="SUSC print control; overrides the incoming CHIPRINT when supplied."),
        V("NCHIITER", Integer(min=1, max=1), 1, is_always_added=True,
               info="One pass only: SPR-KKR 9.7 resets NCHIITER to 1 before the convergence loop.")], is_expert=False, is_optional=False)
    schema = InputParameters(
        "chi",
        [CONTROL("CHI"), TAU, mode, CPA, STRCONST, energy,
            Section("SCF", [V("FULLPOT", False)], is_optional=True, is_expert=True),
            TASK("CHI", add=[
                V("CHIPRINT", Integer(min=0), is_optional=True, info="Response preparation print control."),
                V("NOSYMCHI", False, info="Disable response symmetry; also disables global symmetry in native INITVAR."),
            ]), SITES],
        validators=_validate_chi,
        info="Magnetic susceptibility in the serial ASA workflow",
        description="Periodic ASA FREL CHI, one pass only. Native MPI terminates as under construction; "
                    "native magnetization iteration is forcibly disabled. FULLPOT is unavailable in the public version. "
                    "Converge energy and k meshes; this definition does not expose restart, AF signs or manual NEF arrays."
    )
    schema._task_definition = TaskDefinition("CHI", experimental=True, executable="kkrchi", mpi=False, output_files=(
            "{DATASET_ADSI}*_dens.agr",
            "{DATASET_ADSI}*_chifmag.agr",
            "{DATASET_ADSI}_*_GAMMA-OP.dat",
            "{DATASET_ADSI}_*_CHI-r.dat",
            "{DATASET_ADSI}_*_*_I-spin.dat",
        ), cache_files=("{DATASET_ADSI}.tau", "TKTK_tables"))
    return schema
