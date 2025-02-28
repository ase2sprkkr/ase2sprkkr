registered = False

def register():
    global registered

    if registered:
        return

    # fallback to the old way
    else:
        from ase.calculators.calculator import register_calculator_class
        from ..sprkkr.calculator import SPRKKR  # NOQA: E402
        registered=True
        register_calculator_class('sprkkr', SPRKKR)


try:

    def ase_register(plugin=None):
        global registered
        if registered:
            return
        registered = True
        if plugin:
            rc=plugin.register_calculator
            rio=plugin.register_io_format
        else:
            from ase.plugins.register import \
                register_io_format as fio,   \
                register_calculator as rc

        rc('ase2sprkkr.SPRKKR')
        rio('ase2sprkkr.ase.io', 'SPRKKR potential file',
                            '1F', name='sprkkr', ext='pot')


except ImportError:

    pass
