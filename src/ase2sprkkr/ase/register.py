#First, register the io format by the old way
try:
    from ase.utils.plugins import ExternalIOFormat
    SPRKKRFormat = ExternalIOFormat(
        desc="SPRKKR potential file",
        module="ase2sprkkr.ase.io",
        code="+F",  # read & write, file-based
        ext=["pot", "pot_new"],
    )
except ImportError:
    pass

#Then the new way
try:
    from ase._4.plugins.calculator import CalculatorPlugin
    from ase._4.plugins.io import IOFormatPlugin

    def calculator():
        from ..sprkkr.calculator import SPRKKR  # NOQA: E402
        return SPRKKR

    def sprkkr_reader():
        from ..ase.io import read_sprkkr
        return read_sprkkr

    def sprkkr_writer():
        from ..ase.io import write_sprkkr
        return write_sprkkr
    from ase2sprkkr.ase import io
    import ase2sprkkr

    __ase_plugins__ = {
        CalculatorPlugin(
            "SPRKKR",
            module=ase2sprkkr,
            implementation = calculator,
            bibtex="""@misc{sprkkr_software,
                            author       = {Ebert, H.},
                            title        = {The Munich {SPRKKR} package, version 8.6},
                            howpublished = {University of Munich},
                            year         = {2024},
                            url          = {https://uni-muenchen.de}
                          }""",
        ),
        IOFormatPlugin(
            name = "SPRKKR",
            description="SPRKKR potential file",
            bibtex="""@misc{sprkkr_software,
                            author       = {Ebert, H.},
                            title        = {The Munich {SPRKKR} package, version 8.6},
                            howpublished = {University of Munich},
                            year         = {2024},
                            url          = {https://uni-muenchen.de}
                          }""",
            module=io,
            code="+F",  # read & write, file-based
            extensions=["pot", "pot_new"],
            reader=sprkkr_reader,
            writer=sprkkr_writer
        )
    }

    def register():
        pass

    registered = True

except ImportError:

    registered = False

    def register():
        global registered

        if registered:
            return

        # fallback to the old way
        else:
            from ase.calculators.calculator import register_calculator_class
            from ase2sprkkr import SPRKKR

            registered = True
            register_calculator_class("sprkkr", SPRKKR)
