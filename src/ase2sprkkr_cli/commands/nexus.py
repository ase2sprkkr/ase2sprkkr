"""NeXus conversion of native SPR-KKR ARPES spectra."""

from pathlib import Path

from ase2sprkkr_cli.errors import CLIUsageError


help = "Convert SPR-KKR ARPES SPC spectra to NeXus."
description = """Export a native ARPES .spc file to a NeXus (.nxs) file.
The SPC coordinates, intensity values, and axes are preserved by the
existing ARPES NeXus writer. Optionally include SPR-KKR input provenance.

Install the optional exporter dependencies with: pip install 'ase2sprkkr[nexus]'"""


def parser(parser):
    parser.add_argument("spc", type=Path, help="Native SPR-KKR ARPES .spc file.")
    parser.add_argument(
        "-o", "--output", type=Path, metavar="FILE",
        help="Output NeXus filename (default: source filename with .nxs suffix).",
    )
    parser.add_argument(
        "-i", "--input-file", type=Path, metavar="FILE",
        help="Optional SPR-KKR input file (.inp/.in) for experiment settings; "
             "an unambiguous companion input is detected automatically.",
    )


def run(args, global_args=None):
    source = args.spc
    target = args.output if args.output is not None else source.with_suffix(".nxs")

    if source.resolve() == target.resolve():
        raise CLIUsageError("NeXus output must not overwrite the source SPC file.")
    if not source.is_file():
        raise CLIUsageError(f"SPC file not found: {source}")
    if args.input_file is not None and not args.input_file.is_file():
        raise CLIUsageError(f"SPR-KKR input file not found: {args.input_file}")

    # Import only for actual exports: help and unrelated CLI commands must
    # work even when the optional NeXus dependencies are not installed.
    try:
        import h5py  # noqa: F401
    except ImportError as exc:
        raise CLIUsageError(
            "NeXus export requires optional h5py. "
            "Install it with: pip install 'ase2sprkkr[nexus]'"
        ) from exc

    from ase2sprkkr.result_files import ResultFile
    import pyparsing as pp

    try:
        spectrum = ResultFile.from_file(source, try_only="spc")
    except (ValueError, pp.ParseBaseException) as exc:
        raise CLIUsageError(f"Cannot parse ARPES SPC file {source}: {exc}") from exc

    if not callable(getattr(spectrum, "to_nexus", None)):
        raise CLIUsageError(f"The SPC file {source} does not support NeXus export.")

    # Use the canonical writer: no second reshape/transposition of (NP, NT).
    result = spectrum.to_nexus(target, input_file=args.input_file)
    print(f"NeXus exported to {result}")
