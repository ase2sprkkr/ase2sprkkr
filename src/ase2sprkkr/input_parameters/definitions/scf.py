"""SCF input-parameter definition."""

import numpy as np
from ...common.grammar_types import Integer, SetOf, Keyword
from ...common.warnings import DataValidityError
from ._file_manifests import resolve_file
from .sections import CONTROL, TAU, ENERGY, SCF, SITES, STRCONST, CPA, MODE, TASK
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition, DeclaredFile


class SCFInputParametersDefinition(InputParametersDefinition):
    """SCF schema with files depending on the concrete parameters."""

    def cache_files(self, parameters, directory):
        files = list(super().cache_files(parameters, directory))
        if parameters.SCF.FULLPOT():
            sfnfil = parameters.CONTROL.SFNFIL()
            if sfnfil:
                files.append(resolve_file(directory, sfnfil))
        return tuple(files)

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        potfil = parameters.CONTROL.POTFIL()
        if potfil:
            files.append(resolve_file(directory, f"{potfil}_new"))
        files.append(DeclaredFile(resolve_file(directory, "*_SCFSTART.dos"), optional=True))
        return tuple(files)


def _validate_scf(p, _values, why):
    if why == "set":
        return None
    errors = []
    for name in ("MIX", "MIXOP", "TOL", "EFGUESS"):
        value = p.SCF[name]()
        if value is not None and (not np.isfinite(value) or (name != "EFGUESS" and value <= 0)):
            errors.append(DataValidityError(f"SCF.{name} must be finite" + (" and positive." if name != "EFGUESS" else ".")))
    if p.SCF.FULLPOT() and p.SCF.ALG() == "TCHEBY":
        errors.append(DataValidityError("Native TCHEBY is not available with FULLPOT."))
    grid, ne = p.ENERGY.GRID().tolist(), p.ENERGY.NE().tolist()
    if grid[0] in ("2", "3") and ne[0] < 2:
        errors.append(DataValidityError("GRID=2/3 requires NE>=2 to avoid native one-point integration defects."))
    return errors or None


def input_parameters_definition():
    scf = SCF.copy()
    scf.copy_member("ALG", type=Keyword({"BROYDEN2": "Broyden's second method", "TCHEBY": "Tchebychev mixing", "ANDERSON": "Anderson mixing"}), default_value="BROYDEN2")
    scf.copy_member("NITER", type=Integer(min=1))
    scf.copy_member("ISTBRY", type=Integer(min=1))
    scf.copy_member("ITDEPT", type=Integer(min=3, max=40), info="Broyden/Anderson history depth; native resets values below 3 to 40 and caps values above 40.")
    scf.copy_member("MIXOP", info="Orbital-polarization mixing; omitted native value is initialized before reading MIX and is not necessarily the user-supplied MIX.")
    scf.copy_member("EFGUESS", info="Initial Fermi-energy guess in Ry; skips the initial search when building a starting potential, not all later EF updates.")
    energy = ENERGY()
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = SCFInputParametersDefinition(
        "scf",
        [CONTROL("SCF"), TASK("SCF"), TAU, energy, scf, SITES, STRCONST, CPA, MODE],
        info="SCF - Calculate a self-consistent potential",
        validators=_validate_scf,
    )
    schema._task_definition = TaskDefinition(
        "SCF",
        name="SCF",
        executable="kkrscf",
        mpi=True,
        output_reader="SCF",
    )
    return schema
