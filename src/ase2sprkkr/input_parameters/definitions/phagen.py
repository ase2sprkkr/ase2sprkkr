"""PHAGEN input-parameter definition."""

from ._file_manifests import resolve_file
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from .sections import CONTROL, TASK, TAU, ENERGY, SCF, MODE
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition


class PhagenInputParametersDefinition(InputParametersDefinition):
    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        potfil = parameters.CONTROL.POTFIL()
        if potfil:
            files.extend(
                (
                    resolve_file(directory, f"{potfil}_ex_sprkkr"),
                    resolve_file(directory, f"{potfil}.only_coulomb"),
                )
            )
        return tuple(files)


def _validate_phagen(parameters, _values, why):
    if why == "set":
        return None
    if parameters.SCF.FULLPOT():
        return DataValidityError("PHAGEN supports ASA only.")
    if parameters.SCF.VXC() != "VWN":
        return DataValidityError("PHAGEN requires a VWN starting potential.")
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        return DataValidityError("PHAGEN does not support DMFT or LDA+U.")
    if parameters.ENERGY.SPLITSS() or len(parameters.ENERGY.GRID()) != 1 or len(parameters.ENERGY.NE()) != 1:
        return DataValidityError("PHAGEN Hartree SCF step requires a single energy path.")
    return None



def input_parameters_definition():
    energy = ENERGY()
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = PhagenInputParametersDefinition(
        "phagen",
        [CONTROL("PHAGEN"), TASK("PHAGEN"), TAU, energy, SCF, MODE],
        info="Export radial potential, charge density and Coulomb potential for PHAGEN/MsSpec.",
        description="Requires an existing non-spin-polarized ASA VWN potential with NT=1, no CPA/NLCPA "
                    "and no thermal fluctuations. Native export runs one Hartree SCF step with NITER=1 "
                    "and MIX=1 regardless of the input values, writing POTFIL_ex_sprkkr and Coulomb files. "
                    "PHAGEN_CALCTYPE is hardcoded EXP; no phase-shift solver is invoked.",
        validators=_validate_phagen,
    )
    schema._task_definition = TaskDefinition(
        "PHAGEN",
        name="PHAGEN",
        executable="kkrscf",
        mpi=False,
    )
    return schema
