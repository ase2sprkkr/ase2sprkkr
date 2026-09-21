#!/usr/bin/env python
"""
Plotting the values in SPRKKR output files
"""

from pathlib import Path
import sys
import argparse

if not __package__:
    __package__ = "ase2sprkkr.tools.commands"
sys.path.append(str(Path(__file__).resolve().parents[3]))

from ...common.tools import parse_tuple_function, parse_named_option, append_id_to_filename, parse_inches, main  # NOQA
from ...common.lazy_string import LazyString  # NOQA


@LazyString
def description():
    from ...output_files.output_files import OutputFile

    out = "The type of the file is guessed from the content of the file and from the extension. The currently supported files are: \n"
    defs = OutputFile.definitions.items()
    out += "\n".join(map(lambda x: f"    {x[0].upper()}: {x[1].info()}", defs))
    return out


help = "Plot SPR-KKR output files."

SFN_PLOT_MODES = ("mesh", "cell", "periodic", "shape", "slice", "radial")
SFN_PLOT_MODE_HELP = {
    "mesh": "Show the distinct Voronoi meshes (the SFN default).",
    "cell": "Place each distinct mesh in the crystallographic cell.",
    "periodic": "Show the periodic tessellation clipped to one cell.",
    "shape": "Reconstruct the S=0.5 shape-function isosurface.",
    "slice": "Show a planar section of the reconstructed shape function.",
    "radial": "Show the radial spherical-harmonic coefficients.",
}


def parser(parser):
    def parse_layer(value):
        out = parse_tuple_function(int, 1, 2)(value)
        if len(out) == 1:
            return out[0] - 1
        else:
            return slice(out[0] - 1, out[1])

    parser.add_argument("output", help="SPR-KKR output file name (see the supported files above).", nargs="+")
    parser.add_argument(
        "-o",
        "--output_filename",
        dest="filename",
        type=str,
        help="The plot will be saved to a file with given name, instead of showing it on the screen. If more values are given (see -v), their name will be added to the filename.",
        default=None,
        required=False,
    )
    parser.add_argument(
        "-v",
        "--value",
        dest="value",
        type=str,
        help="Only the value of the given name will be plotted (option can be repeated)",
        action="append",
        default=[],
        required=False,
    )
    parser.add_argument(
        "-V",
        "--show-values",
        dest="show_values",
        help="Show, which values can be plotted.",
        action="store_true",
        default=False
    )

    parser.add_argument(
        "-s",
        "--plot_size",
        dest="figsize",
        default=argparse.SUPPRESS,
        type=parse_tuple_function(parse_inches, 2),
        help='The plot size. Example: "5cm,5cm", "6,4". The default units are inches.',
        required=False,
    )
    parser.add_argument("-c", "--colormap", dest="colormap", type=str, help="Matplotlib colormap", required=False)
    parser.add_argument(
        "-d", "--dpi", dest="dpi", type=float, help="DPI of the resulting image (default 600)", required=False
    )
    parser.add_argument(
        "-n",
        "--norm",
        dest="norm",
        choices=["lin", "log"],
        help="Matplotlib colormap will use linear or logarithmic scale (the default behavior depends on the plotted data)",
        required=False,
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-l",
        "--use_latex",
        dest="latex",
        action="store_true",
        help="Force use LaTex for generating captions. Default is to use the default matplotlib settings (can be configured per user).",
        required=False,
    )
    group.add_argument(
        "-L",
        "--do_not_use_latex",
        dest="latex",
        action="store_false",
        help="Do not use LaTex for generating captions",
        required=False,
    )
    parser.add_argument(
        "-S",
        "--set",
        dest="args",
        type=lambda x: parse_named_option(x, True),
        help="Given a value of the format name=value, pass the value to the plotting function. You can so override various options that matplotlib plotting functions accept (e.g. vmin or vmax for pcolormesh), or the values that can be set using set_<something> functions (e.g. title or (x|y)label). This option can be repeated.",
        action="append",
        default=[],
        required=False,
    )
    parser.set_defaults(latex=None)  # or True/False if you want a default
    group.add_argument(
        "--separate_plots",
        help="Plot each value from file in a separate window or/and file.",
        action="store_true",
        required=False,
    )

    group = parser.add_argument_group("BSF specific options")
    group.add_argument(
        "--sites",
        help="Select a site to plot (for BSF). Either number or two comma delimited numbers from,to. Numbering starts from 1. ",
        type=parse_layer,
        default=argparse.SUPPRESS,
        required=False,
    )
    group.add_argument(
        "--fermi",
        help="Draw a line at Fermi energy. Optional float specifies line width.",
        nargs="?",
        const=True,
        type=float,
        default=argparse.SUPPRESS,
        required=False,
    )

    group = parser.add_argument_group("SFN specific options")
    modes = group.add_mutually_exclusive_group()
    for mode in SFN_PLOT_MODES:
        modes.add_argument(
            f"--{mode}",
            dest="what",
            action="store_const",
            const=mode,
            help=SFN_PLOT_MODE_HELP[mode],
        )
    parser.set_defaults(what=None)


def run(args, global_args):
    from ...output_files.output_files import OutputFile
    from ...output_files.definitions.sfn import SFNOutputFile
    from pathlib import Path
    from ...gui.plot import Multiplot

    kwargs = vars(args).copy()

    outputs = kwargs.pop("output")
    value = kwargs.pop("value")
    filename = kwargs.pop("filename")
    generic_options = dict(kwargs.pop("args"))
    if "what" in generic_options:
        raise ValueError(
            "SFN plot mode cannot be set using '-S what=...'; use one of: "
            + ", ".join(f"--{mode}" for mode in SFN_PLOT_MODES)
        )
    kwargs.update(generic_options)
    show_values = kwargs.pop("show_values")
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    for output in outputs:

        def plot(fn, of, value=None):
            try:
                fn()
            except Exception as e:
                value_label = f" value '{value}'" if value else ""
                raise ValueError(
                    f"File '{of}'{value_label} can not be plotted: {e}"
                ) from e

        of = OutputFile.from_file(output, unknown=False)

        if "what" in kwargs:
            what = kwargs["what"]
            if what not in SFN_PLOT_MODES:
                raise ValueError(
                    f"Unknown SFN plot mode '{what}'; choose one of: "
                    f"{', '.join(SFN_PLOT_MODES)}"
                )
            if not isinstance(of, SFNOutputFile):
                raise ValueError(
                    "SFN plot mode options are only valid for SFN files; "
                    f"'{output}' was read as {type(of).__name__}."
                )
            if value:
                raise ValueError(
                    "An SFN plot mode selects a whole-file plot and "
                    "cannot be combined with '--value'."
                )

        if show_values:
            print(f"In file {output}, there are values to plot:")
            for i in of:
                if hasattr(i, "plot") and callable(i.plot):
                   print(f"   {i.name:<15} {i._definition.info(False)}")
            continue

        for x in ["layer"]:
            if x in kwargs and x not in of.plot_parameters:
                raise ValueError(f"Argument '--{x}' is not valid for {of}")

        if filename:
            if len(outputs) > 1:
                fn = Path(filename)
                fn = str(fn.with_suffix("")) + f"_{Path(output).stem}{fn.suffix}"
            else:
                fn = filename
        else:
            fn = None

        if value:

            def vals():
                for name in value:
                    try:
                        val = of[name.upper()]
                    except KeyError as exc:
                        raise ValueError(
                            f"There is no value named '{name.upper()}' "
                            f"in the output file '{output}'."
                        ) from exc
                    if not hasattr(val, "plot"):
                        raise ValueError(
                            f"Value '{name.upper()}' does not know how it should be plotted."
                        )
                    yield val

            vals = [ i for i in vals() ]
            ln = sum( i.number_of_plots() for i in vals )
            with Multiplot(**kwargs, number_of_plots=len(value), filename=fn, layout=min(ln, 2)) as mp:
                for val in vals:
                    plot(lambda: mp.plot(val), output, value)
        else:
            if not hasattr(of, "plot"):
                raise ValueError(f"File '{of}' does not know, how it should be plotted.")
            plot(lambda: of.plot(filename=fn, **kwargs), output)

if __name__ == "__main__":
    main(globals())
