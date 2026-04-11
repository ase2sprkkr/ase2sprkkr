#!/usr/bin/env python3
"""Convert SPRKKR output to UppASD input files and Jij plots."""
import sys
from pathlib import Path
from ...common.tools import main  # NOQA

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'
sys.path.append(str(Path(__file__).resolve().parents[3]))

help = 'Convert SPRKKR output to UppASD input files and plots the computed Jij (experimental)'
description = 'Given .pot, .jxc and .dmi file, plot the Jij values.'

def parser(parser):
    parser.add_argument('-p', '--pot-file', type=str, default=None,
                        help='Input potential file (.pot or .pot_new)')
    parser.add_argument('-j', '--jxc-file', type=str, default=None,
                        help='Input exchange interaction file (_XCPLTEN_Jij.dat)')
    parser.add_argument('-d', '--dmi-file', type=str, default=None,
                        help='Input DMI file (_DMIVEC_Dij.dat)')
    parser.add_argument('-o', '--output-dir', type=str, default='.',
                        help='Output directory for files')
    parser.add_argument('-n', '--no-write', action='store_true',
                        help='Skip writing output files')
    parser.add_argument('--plot', action='store_true',
                        help='Generate exchange interaction plots')
    parser.add_argument('--no-plot', action='store_true',
                        help='Skip generating plots')
    parser.add_argument('--separate-plots', action='store_true',
                        help='Generate one plot file per site type')
    parser.add_argument('--axis', type=str, choices=('all', 'x', 'y', 'z'), default='all',
                        help='DMI component to plot: all, x, y, or z')
    parser.add_argument('-r', '--exchange-radius', type=float, default=4.0,
                        help='Maximum distance for exchange interaction plots')
    parser.add_argument('-c', '--cartesian', action='store_const', dest='coordinates', default="lattice",
                        const="cartesian", help='Use cartesian coordinates for interaction instead of lattice ones')
    parser.add_argument('-e', '--exclude', type=str, nargs='*', default=None,
                        help='Comma-separated site selectors to exclude')
    parser.add_argument('-i', '--include', type=str, nargs='*', default=None,
                        help='Comma-separated site selectors to include')
    parser.add_argument('--include-vacuum', action='store_true',
                        help='Include vacuum sites in selector matching and exported outputs')
    parser.add_argument('-f', '--font-size', type=int, default=14,
                        help='Font size for plots')

def run(args):
    import glob

    from ...bindings.uppasd import Coordinates, write_dmfile, write_jfile, write_mom_file, write_pos_file  # NOQA
    from ...output_files.output_files import OutputFile  # NOQA
    from ...output_files.definitions.jxc import JXCOutputFile  # NOQA
    from ...potentials.potentials import Potential  # NOQA

    def _load_jxc_output(filename, potential):
        output = OutputFile.from_file(filename, try_only='JXC', unknown=False)
        output.potential = potential
        return output

    def _plot_exchange_interactions(output, output_dir, selector, exchange_radius, font_size,
                                    separate_plots=False, axis='all'):
        if output.is_Jij():
            plot_name = 'Jij'
        elif axis == 'all':
            plot_name = 'Dij'
        else:
            plot_name = f'Dij_{axis}'

        filename = output_dir / ('{name}.pdf' if separate_plots else f'{plot_name}.pdf')
        out = output.plot(
            layout=2,
            filename=str(filename),
            show=False,
            selector=selector,
            exchange_radius=exchange_radius,
            font_size=font_size,
            axis=axis,
            separate_plots=separate_plots,
        )
        if out:
            print(f'Created plot files in: {output_dir}')
        else:
            print('No interaction data to plot')
        return True

    def find_files_by_pattern(args):
        if args.pot_file is None:
            pot_files = glob.glob('*.pot_new')
            if pot_files:
                args.pot_file = pot_files[0]
            else:
                pot_files = glob.glob('*.pot')
                if pot_files:
                    args.pot_file = pot_files[0]

        if args.jxc_file is None:
            jxc_files = glob.glob('*_XCPLTEN_Jij.dat')
            if jxc_files:
                args.jxc_file = jxc_files[0]

        if args.dmi_file is None:
            dmi_files = glob.glob('*_DMIVEC_Dij.dat')
            if dmi_files:
                args.dmi_file = dmi_files[0]

        return args

    args = find_files_by_pattern(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.pot_file is None:
        print('Error: No potential file found!')
        sys.exit(1)

    try:
        potential = Potential.from_file(args.pot_file)
        atoms = potential.atoms

        if args.jxc_file:
            jxc_output = _load_jxc_output(args.jxc_file, potential)
        else:
            jxc_output = None
            print('No JXC file found - skipping exchange processing')

        if args.dmi_file:
            dmi_output = _load_jxc_output(args.dmi_file, potential)
        elif args.jxc_file:
            print('No DMI file found - skipping DMI export')
            dmi_output = None

        selector_obj =  jxc_output or dmi_output or JXCOutputFile.from_atoms(atoms)
        selector = selector_obj.create_selector(iq=args.include, it=args.include, exclude_it=args.exclude, exclude_vc=not args.include_vacuum)
        coordinates = getattr(Coordinates, args.coordinates)
        print(f'Potential: {args.pot_file}')
        print(f'Number of atoms: {len(atoms.sites)}')
        if selector.it is not ...:
            print(f'Selected IQs: {", ".join(map(lambda i: selector_obj.it_labels[i], selector.it))}')

        if not args.no_write and dmi_output is not None:
            dmi_output.write_uppasd_file(
                output_dir / 'dmfile.dat',
                selector=selector,
                exchange_radius=args.exchange_radius,
                coordinates=coordinates,
            )

        if not args.no_write and jxc_output is not None:
            jxc_output.write_uppasd_file(
                output_dir / 'jfile.dat',
                selector=selector,
                exchange_radius=args.exchange_radius,
                coordinates=coordinates,
            )

        if not args.no_write:
            write_pos_file(
                atoms,
                output_dir / 'posfile.dat',
                selector=selector,
            )
            if not write_mom_file(
                    atoms,
                    output_dir / 'momfile.dat',
                    selector=selector,
                ):
                print('Warning: no selected site type has moments. Skipping momfile.dat.')

        if args.plot and not args.no_plot and jxc_output is not None:
            _plot_exchange_interactions(
                jxc_output,
                output_dir,
                selector=selector,
                exchange_radius=args.exchange_radius,
                font_size=args.font_size,
                separate_plots=args.separate_plots,
                axis=args.axis,
            )
        if args.plot and not args.no_plot and args.dmi_file:
            _plot_exchange_interactions(
                dmi_output,
                output_dir,
                selector=selector,
                exchange_radius=args.exchange_radius,
                font_size=args.font_size,
                separate_plots=args.separate_plots,
                axis=args.axis,
            )

    except FileNotFoundError as exc:
        print(f'Error: {exc}')
        sys.exit(1)
    except Exception as exc:
        print(f'Unexpected error: {exc}')
        sys.exit(1)


if __name__ == '__main__':
    main(globals())
