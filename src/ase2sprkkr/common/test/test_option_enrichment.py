from ase2sprkkr.common.generated_configuration_definitions import NumpyViewDefinition
from ase2sprkkr.common.options import Option


def _plot(_option, **_kwargs):
    return None


def test_generated_option_enrichment_does_not_require_value_type():
    option = Option(NumpyViewDefinition("VIEW", "DATA", plot=_plot))

    assert option.plot() is None
    assert option.number_of_plots() == 1
