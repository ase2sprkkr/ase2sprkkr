"""Analyse crystal symmetry and optionally export transformed unit cells."""

import json
from pathlib import Path


help = "Analyse crystal symmetry and transform unit cells."
description = """Determine space groups, equivalent sites, Wyckoff positions, and
optionally construct a primitive or standardized cell using spglib.
Input/output structures use ASE formats (including registered SPR-KKR POT)."""


def parser(parser):
    parser.add_argument("input", help="Structure file readable by ASE")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--primitive", action="store_true", help="Find a primitive cell")
    mode.add_argument("--standardize", action="store_true", help="Find the standard conventional cell")
    parser.add_argument(
        "--equivalence", choices=("site-types", "occupation"),
        default="site-types", help="Site-equivalence rule (default: site-types)"
    )
    parser.add_argument("--symprec", type=float, default=1e-5, help="Position tolerance in Å")
    parser.add_argument("--angle-tolerance", type=float, default=-1.0, help="Angle tolerance in degrees")
    parser.add_argument("--idealize", action="store_true", help="Idealize in the original basis (also works without cell conversion)")
    parser.add_argument("--ase-arrays", choices=("require-equal", "discard", "choose-first"),
                        default="require-equal", help="Transfer policy for additional ASE arrays")
    parser.add_argument("--sites", action="store_true", help="Show equivalent sites and Wyckoff positions")
    parser.add_argument("--operations", action="store_true", help="Show rotations and translations")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show detailed symmetry information")
    parser.add_argument("--view", action="store_true", help="Visualize the resulting structure")
    parser.add_argument("-o", "--output", help="Write transformed structure to this file")
    parser.add_argument("-f", "--format", help="Override ASE output format")
    parser.add_argument("--report", help="Write JSON analysis report")


def _report(analysis, *, sites=False, operations=False):
    data = analysis.dataset
    result = {
        "spacegroup_number": int(data.number),
        "international": str(data.international),
        "hall": str(data.hall),
        "pointgroup": str(data.pointgroup),
        "number_of_operations": len(data.rotations),
        "number_of_sites": len(analysis.site_tags),
        "number_of_primitive_sites": len(set(map(int, data.mapping_to_primitive))),
        "equivalence": analysis.equivalence,
        "symprec": analysis.symprec,
        "angle_tolerance": analysis.angle_tolerance,
        "transformation_matrix": data.transformation_matrix.tolist(),
        "origin_shift": data.origin_shift.tolist(),
    }
    if sites:
        result["sites"] = [
            {
                "index": i + 1,
                "tag": int(analysis.site_tags[i]),
                "equivalent_atom": int(data.equivalent_atoms[i]) + 1,
                "primitive_atom": int(data.mapping_to_primitive[i]) + 1,
                "wyckoff": str(data.wyckoffs[i]),
                "site_symmetry": str(data.site_symmetry_symbols[i]),
            }
            for i in range(len(analysis.site_tags))
        ]
    if operations:
        result["operations"] = [
            {"rotation": rot.tolist(), "translation": translation.tolist()}
            for rot, translation in zip(data.rotations, data.translations)
        ]
    return result


def run(args, global_args=None):
    from ase.io import read, write
    from ase2sprkkr.sprkkr.symmetry import analyze_symmetry

    if args.output and not (args.primitive or args.standardize or args.idealize):
        raise ValueError("--output requires --primitive, --standardize, or --idealize")
    if args.format and not args.output:
        raise ValueError("--format requires --output")
    if args.output and Path(args.output).resolve() == Path(args.input).resolve():
        raise ValueError("Input file cannot be overwritten")

    atoms = read(args.input)
    analysis = analyze_symmetry(
        atoms,
        equivalence=args.equivalence,
        symprec=args.symprec,
        angle_tolerance=args.angle_tolerance,
    )
    data = analysis.dataset
    result = None
    if args.primitive:
        result = analysis.find_primitive(atoms, idealize=args.idealize, ase_arrays=args.ase_arrays)
    elif args.standardize:
        result = analysis.standardize(atoms, idealize=args.idealize, ase_arrays=args.ase_arrays)
    elif args.idealize:
        result = analysis.idealize(atoms)

    output = _report(
        analysis,
        sites=args.sites or args.verbose,
        operations=args.operations or args.verbose,
    )
    if result is not None:
        output["transformed_sites"] = len(result)
        output["source_sites"] = len(atoms)

    print(f"Space group: {data.international} (No. {data.number})")
    print(f"Hall: {data.hall}; Point group: {data.pointgroup}")
    print(f"Sites: {len(atoms)}; primitive sites: {output['number_of_primitive_sites']}")
    print(f"Symmetry operations: {len(data.rotations)}; equivalence: {args.equivalence}")
    if result is not None:
        print(f"Transformed cell: {len(result)} sites")
    if args.sites or args.verbose:
        print("Index  Species  Equivalent  Primitive  Wyckoff  Site symmetry")
        for site in output["sites"]:
            print(
                f"{site['index']:>5}  {site['tag']:>7}  "
                f"{site['equivalent_atom']:>10}  {site['primitive_atom']:>9}  "
                f"{site['wyckoff']:>7}  {site['site_symmetry']}"
            )
    if args.operations or args.verbose:
        for index, op in enumerate(output["operations"], 1):
            print(f"Operation {index}: rotation={op['rotation']}, translation={op['translation']}")

    if args.report:
        Path(args.report).write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    if args.output:
        kwargs = {"format": args.format} if args.format else {}
        write(args.output, result, **kwargs)
    if args.view:
        from ase.visualize import view
        view(result if result is not None else atoms)
