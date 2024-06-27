try:
    from ase.plugins import register_calculator, register_io_format

    def ase_register():
        register_calculator('ase2sprkkr.SPRKKR')
        register_io_format('ase2sprkkr.ase.io', 'SPRKKR potential file',
                            '1F', name='sprkkr', ext='pot')


except ImportError:

    pass
