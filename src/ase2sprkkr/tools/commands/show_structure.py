#!/usr/bin/env python
"""
Display or export atomic structures from an SPR-KKR potential and an optional
in_struct.inp file.
"""

from pathlib import Path

if not __package__:
    __package__ = "ase2sprkkr.tools.commands"
import sys

sys.path.append(str(Path(__file__).resolve().parents[3]))

help = "Visualise potential (and optionally in_struct.inp) files."
description = (
    "Display the structure in ASE, export it using an ASE format, or both. "
    "Without an output filename or an explicit viewing option, display the "
    "structure (the potential structure if --input is omitted)."
)
aliases = ["view"]


def parser(parser):
    parser.add_argument("potential_file", help="SPR-KKR potential file (.pot)")
    parser.add_argument(
        "-i", "--input",
        help="Optional in_struct.inp file; otherwise use the potential structure.",
    )
    parser.add_argument(
        "-o", "--output",
        help="Export the input structure, or the potential structure if --input is omitted.",
    )
    parser.add_argument(
        "-O", "--potential-output", nargs="?", const=True,
        help="Export the potential structure; without a filename, derive <output>_pot.<ext> (requires -o).",
    )
    parser.add_argument(
        "-f", "--format", default="cif",
        help="ASE export format (default: cif).",
    )
    parser.add_argument(
        "-a", "--visualise", action="store_true",
        help="Show the input structure (or the potential structure without --input).",
    )
    parser.add_argument(
        "-A", "--visualise-potential", action="store_true",
        help="Show the potential structure.",
    )
    parser.add_argument(
        "--close_after", "--close-after",
        type=float,
        default=None,
        help="Close each ASE viewer after this many seconds (blocking mode).",
    )
    parser.add_argument(
        "--no-view", action="store_true",
        help="Do not launch ASE viewers, including automatic viewing (useful for headless runs).",
    )
    parser.add_argument(
        "-s", "--scale-radii", type=float, default=0.5,
        help="ASE viewer atomic radius scale (default: 0.5).",
    )
    parser.add_argument(
        "-v", "--vac", dest="vacuum_height", type=float, default=10.0,
        help="Vacuum height in angstroms (default: 10.0).",
    )
    parser.add_argument(
        "-b", "--nbulk", type=int, default=2,
        help="Number of bulk-unit repetitions (default: 2).",
    )


def run(args, global_args=None):
    """Build the requested structures, export selected files and/or open ASE."""
    from ...potentials.potentials import Potential
    from ...sprkkr.structure import structure_file_to_atoms
    from ...ase.visualise import view

    potential_output = args.potential_output
    if potential_output is True:
        if not args.output:
            raise ValueError("--potential-output without a filename requires --output")
        output = Path(args.output)
        potential_output = str(output.with_name(f"{output.stem}_pot{output.suffix}"))

    potential = Potential.from_file(args.potential_file)
    potential_atoms = potential.atoms

    # The primary structure is the input structure, when present, or the
    # structure carried by the potential itself.
    if args.input:
        atoms = structure_file_to_atoms(
            args.input, potential, n_bulk=args.nbulk,
            vacuum_height=args.vacuum_height,
        )
    else:
        atoms = potential_atoms

    if args.output:
        atoms.write(args.output, format=args.format)
    if potential_output:
        potential_atoms.write(potential_output, format=args.format)

    if args.no_view:
        return

    # Explicit viewing options select exactly what to show. Otherwise show
    # the primary structure when it has not been exported to a file.
    explicit_view = args.visualise or args.visualise_potential
    primary_exported = bool(args.output or (not args.input and potential_output))
    show_primary = args.visualise or (not explicit_view and not primary_exported)

    viewer_options = {
        "scale_radii": args.scale_radii,
        "blocking": args.close_after is not None,
        "close_after": args.close_after,
    }
    if args.visualise_potential:
        view(potential_atoms, **viewer_options)
    if show_primary and (args.input or not args.visualise_potential):
        view(atoms, **viewer_options)


if __name__ == "__main__":
    from ...common.tools import main

    main(globals())
