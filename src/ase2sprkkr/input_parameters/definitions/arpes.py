"""ARPES task input parameters definition"""

import numpy as np
from ...common.grammar_types import Integer, SetOf, DefKeyword, Keyword, Range, Array, flag
from .sections import CONTROL, TAU, ENERGY, SITES
from ..input_parameters_definitions import (
    InputParametersDefinition as InputParameters,
    InputValueDefinition as V,
    InputSectionDefinition as Section,
)
from ...sprkkr.sprkkr_grammar_types import Site
from ...common.warnings import DataValidityError


class ArpesInputParametersDefinition(InputParameters):
    """Shared SPEC_INPUT schema for ARPES-family tasks."""

    task_selectors = ("ARPES",)

    def task_definition_name(self, parameters):
        return str(parameters.TASK.TASK() or "ARPES").upper()


def _task_default(option):
    """Default the shared TASK selector from the requested task context."""
    parameters = option._get_root_container()
    return parameters.requested_task_name or "ARPES"


def _electron_angle_default(option):
    """SPEC_INPUT defaults differ for photoemission and SPLEED."""
    angle = (45.0 if option.name == 'THETA' else 270.0) if option._container._container.TASK.TASK() == 'SPLEED' else 0.0
    return [angle, angle]


def _electron_spol_default(option):
    return 1 if option._container._container.TASK.TASK() == 'SPLEED' else 2


def angular_scan_errors(parameters):
    """Return inconsistencies rejected by SPEC_INPUT, using its unset defaults."""
    task, electron = parameters.TASK, parameters.SPEC_EL
    if task.INPVER() == 0 or task.TASK() != 'ARPES' or electron.K1() is not None:
        return []
    errors = []
    for angle, count in (('THETA', 'NT'), ('PHI', 'NP')):
        bounds = electron[angle]()
        span = 0. if bounds is None else abs(bounds[-1] - bounds[0])
        points = electron[count]()
        points = 1 if points is None else points
        if (span > 1e-5) != (points > 1):
            if angle == 'PHI' and span <= 1e-5 and 'BETA1' in electron and electron.BETA1() is not None:
                continue
            errors.append(f'SPEC_EL.{angle} and SPEC_EL.{count} must specify a consistent angular scan: '
                          'a nonzero range requires more than one point, and vice versa.')
    return errors


def _validate_angular_scan(parameters, _container, why):
    if why != 'set':
        errors = [DataValidityError(message) for message in angular_scan_errors(parameters)]
        if parameters.TASK.INPVER() == 0:
            return errors
        for section, names in ((parameters.ENERGY, ("EWORK_EV", "IMV_INI_EV", "IMV_FIN_EV")), (parameters.SPEC_PH, ("EPHOT",))):
            for name in names:
                value = section[name]()
                if not np.isfinite(value) or value < 0 or (name == "EPHOT" and value == 0):
                    errors.append(DataValidityError(f"{section.name}.{name} must be finite and nonnegative (EPHOT strictly positive)."))
        electron = parameters.SPEC_EL
        if parameters.TASK.TASK() == "ARPES" and any(electron[name]() is not None for name in ("K3", "K4", "NK3", "NK4")):
            errors.append(DataValidityError("SPR-KKR 9.7 SPEC_INPUT does not read K3/K4/NK3/NK4; use K1/K2."))
        if electron.K1() is not None:
            for vector, count in (("K1", "NK1"), ("K2", "NK2")):
                if electron[vector]() is not None and electron[count]() is None:
                    errors.append(DataValidityError(f"SPEC_EL.{vector} requires {count}."))
        for name in ("THETA", "PHI", "KA", "K1", "K2"):
            value = electron[name]()
            if value is not None and not np.all(np.isfinite(value)):
                errors.append(DataValidityError(f"SPEC_EL.{name} must be finite."))
        return errors


def input_parameters_definition():
    return ArpesInputParametersDefinition(
        "arpes",
        [
            CONTROL("ARPES"),
            TAU,
            ENERGY(
                emin=(None, "Minimum of the energy window in Ry (absolute energy)", -8.0),
                emax=(None, "Maximum of the energy window in Ry (absolute energy)", 5.0),
                add=[
                    V("EWORK_EV", 4.2, info="Inner potential of the bulk crystal in eV"),
                    V(
                        "IMV_INI_EV", 0.05, info="Imaginary part of the potential in eV (initial state)"
                    ),
                    V(
                        "IMV_FIN_EV", 2.0, info="Imaginary part of the potential in eV (final state)"
                    ),
                ],
                defaults={"GRID": 1, "NE": 300},
            ),
            SITES.copy(defaults={"NL": 4}),
            Section(
                "TASK",
                [
                    V(
                        "TASK",
                        DefKeyword(
                            {
                                "ARPES": "Angle resolved photoemission spectroscopy",
                                "AIPES": "Angle integrated photoemission spectroscopy",
                                "SPLEED": "Spin polarized LEED (experimental feature)",
                                "BAND": "band structure calculations (experimental feature)",
                            }
                        ),
                        _task_default,
                        is_required=True,
                        is_optional=False,
                        name_in_grammar=False,
                        info="Type of the calculation",
                    ),
                    V("IQ_AT_SURF", Site.I, 1),
                    V("MILLER_HKL", SetOf(int, length=3), [0, 0, 1]),
                    V("CRYS_VEC", True, info="Miller indices with respect to crystalographic primitive vectors"),
                    V(
                        "STRVER",
                        Integer(min=0, max=1),
                        1,
                        is_expert=True,
                        is_always_added=True,
                        info="Set to 0 to supply the ARPES input file 'struc.inp' manually (and do not generate it).",
                    ),
                    V(
                        "INPVER",
                        Integer(min=0, max=1),
                        1,
                        is_expert=True,
                        is_always_added=True,
                        info="Set to 0 to use an old input.inp from old rslab",
                    ),
                ],
                is_optional=False,
            ),
            Section(
                "SPEC_PH",
                [
                    V("THETA", 45.0, info="Direction of the photon (the polar coordinate)"),
                    V("PHI", 0.0, info="Direction of the photon (the azimuth coordinate)"),
                    V("POL_P", DefKeyword("P", "S", "C+", "C-"), info="Polarization of the light"),
                    V("EPHOT", 25.0, info="Photon energy in eV"),
                    V("ALQ", expert=45.0, info="Alignment of polarization vector or pol.ellipsis"),
                    V(
                        "DELQ",
                        expert=0.0,
                        info="Phase shift between real and imaginary part of e-vector, delq=90 for circular polarized light",
                    ),
                    V(
                        "NPOL",
                        Keyword(
                            {
                                0: "unpolarized and p-s dichroism for the calculation",
                                1: "p-pol or rcp or elliptical (depends on icirc, etc.)",
                                2: "s-pol or lcp or elliptical (depends on icirc, etc.)",
                                3: "dichroism (ddad, ldad)",
                            }
                        ),
                        expert=1,
                        info="Controls the polarization and dichroism",
                    ),
                    V(
                        "ICIRC",
                        Keyword(
                            {
                                0: "elliptically pol. light: alq, delq arbitrary",
                                1: "linear pol. light: alq arbitrary, delq=0",
                                2: "circular pol. light: alq=45, delq = 90",
                            }
                        ),
                        expert=1,
                        info="controls the polarization and dichroism",
                    ),
                    V(
                        "IDREH",
                        Keyword(
                            {
                                0: "linearly polarized (equals icirc=1)",
                                1: "right circular polarization",
                                -1: "left circular polarization",
                            }
                        ),
                        expert=0,
                        info="Helicity of the photons",
                    ),
                    V(
                        "IFSP",
                        Keyword({0: "fixed", 1: "variable"}),
                        is_expert=True,
                        is_required=False,
                        info="Photon azimuth angle type",
                    ),
                    V("THETA_FIX", float, is_expert=True, info="Light and electrons are at fixed polarization angle"),
                ],
                info="",
            ),
            Section(
                "SPEC_EL",
                [
                    V("THETA", Range(float), _electron_angle_default, info="Scattering angle", is_required=False),
                    V("PHI", Range(float), _electron_angle_default, info="Scattering angle", is_required=False),
                    V(
                        "NT",
                        Integer(min=1),
                        1,
                        info="Number of angular values for a rotation in polar coordinate.",
                        is_required=False,
                    ),
                    V(
                        "NP",
                        Integer(min=1),
                        1,
                        info="Number of angular values for a rotation in azimuth coordinate.",
                        is_required=False,
                    ),
                    V("KA", Range(float), info="Scatering in momentum space ", is_required=False),
                    V(
                        "K1",
                        Range(float),
                        info="Translating vector of the scatering in momentum space ",
                        is_required=False,
                    ),
                    V("NK1", Integer(min=1), info="Number of momentum steps for the integration", is_required=False),
                    V(
                        "K2",
                        Range(float),
                        info="Translating vector 2 of the scatering in momentum space ",
                        is_required=False,
                    ),
                    V("NK2", Integer(min=1), info="Number of momentum steps 2 for the integration", is_required=False),
                    V(
                        "K3",
                        Range(float),
                        info="Translating vector 3 of the scatering in momentum space ",
                        is_required=False,
                    ),
                    V("NK3", Integer(min=1), info="Number of momentum steps 3 for the integration", is_required=False),
                    V(
                        "K4",
                        Range(float),
                        info="Translating vector 4 of the scatering in momentum space ",
                        is_required=False,
                    ),
                    V("NK4", Integer(min=1), info="Number of momentum steps 4 for the integration", is_required=False),
                    V("POL_E", DefKeyword("PZ")),
                    V("SPOL", int, _electron_spol_default, is_required=False),
                    V("PSPIN", SetOf(float, length=3), is_required=False),
                    V("BETA1", float, is_required=False, info="Begin of the rotation"),
                    V("BETA2", float, is_required=False, info="End of the rotation"),
                    V("ROTAXIS", SetOf(int, length=3), is_required=False, info="Axis of the rotation"),
                    V(
                        "TYP",
                        Keyword(
                            {
                                0: "i(e) diagram",
                                1: "rotation diagram -> phi scan",
                                2: "scattering-angle diagram -> theta scan",
                                3: "orthonormal projection",
                                4: "stereographic projection",
                            },
                            description="3,4 only for angular resolved\npe (ups, xps) note: nt=np-> nx,ny",
                        ),
                        expert=1,
                        info="Crystal coordinats in splout, xpsrun, or upsrun",
                    ),
                    V("ISTR", SetOf(int, length=2), expert=[0, 0], info="beam number (h,k)"),
                    V("POL0", SetOf(int, length=3), expert=[0, 0, 0], info="initial pol."),
                    V("POL0L", SetOf(int, length=3), expert=[0, 0, 0], info="initial pol. in the laboratory system"),
                    V(
                        "Q1",
                        complex,
                        expert=1.0 + 0.0j,
                        info="Amplitude 1 of the photoelectron used in spin polarized calculations",
                    ),
                    V(
                        "Q2",
                        complex,
                        expert=0.0 + 0.0j,
                        info="Amplitude 2 of the photoelectron used in spin polarized calculations",
                    ),
                    V(
                        "Q3",
                        complex,
                        expert=0.0 + 0.0j,
                        info="Amplitude 3 of the photoelectron used in spin polarized calculations",
                    ),
                    V(
                        "Q4",
                        complex,
                        expert=1.0 + 0.0j,
                        info="Amplitude 4 of the photoelectron used in spin polarized calculations",
                    ),
                ],
            ),
            Section(
                "SPEC_STR",
                [
                    V("N_LAYDBL", SetOf(Integer(min=0), length=2), [10, 10]),
                    V("NLAT_G_VEC", Integer(min=1), 57),
                    V("N_LAYER", Integer(min=1), 50),
                    V("SURF_BAR", SetOf(float), [0.25, 0.25]),
                    V("TRANSP_BAR", flag, False),
                ],
            ),
        ],
        validators=_validate_angular_scan,
        info="Angle resolved photoemission spectroscopy",
    )


""" ARPES task input parameters definition"""
