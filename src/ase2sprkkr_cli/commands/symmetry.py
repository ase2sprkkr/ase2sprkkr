"""Analyse crystal symmetry and optionally export transformed unit cells."""

import json
from pathlib import Path


help = "Analyse crystal symmetry and transform unit cells."
description = """Determine space groups, equivalent sites, Wyckoff positions, and
optionally construct a primitive or standardized cell using spglib.
Input/output structures use ASE formats (including registered SPR-KKR POT)."""


# The tags are taken from the selected CIF data block, not inferred from ASE
# coordinates or normalized to the conventional spglib setting.
_CIF_SYMMETRY_TAGS = {
    "number": (
        "_space_group.it_number", "_space_group_it_number",
        "_symmetry_int_tables_number",
    ),
    "hermann_mauguin": (
        "_space_group.name_h-m_alt", "_space_group_name_h-m_alt",
        "_symmetry_space_group_name_h-m", "_space_group.name_h-m_ref",
        "_space_group.name_h-m_full",
    ),
    "hall": (
        "_space_group.name_hall", "_space_group_name_hall",
        "_symmetry_space_group_name_hall",
    ),
    "crystal_system": (
        "_space_group.crystal_system", "_space_group_crystal_system",
        "_symmetry_cell_setting",
    ),
    "operations": (
        "_space_group_symop.operation_xyz", "_space_group_symop_operation_xyz",
        "_symmetry_equiv_pos_as_xyz",
    ),
}


def _declared_cif_symmetry(info):
    """Extract CIF-declared symmetry, preserving its values and not inferring P1."""
    tags = {key.lower(): value for key, value in info.items() if isinstance(key, str)}
    declared = {}
    for field, candidates in _CIF_SYMMETRY_TAGS.items():
        for name in candidates:
            value = tags.get(name)
            if value is None or (isinstance(value, str) and value.strip() in ("?", ".")):
                continue
            if field == "operations":
                value = [value] if isinstance(value, str) else list(value)
                if not value:
                    continue
            declared[field] = value
            break
    return declared or None


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
    from ase.io.formats import filetype
    from ase2sprkkr.sprkkr.symmetry import analyze_symmetry

    if args.output and not (args.primitive or args.standardize or args.idealize):
        raise ValueError("--output requires --primitive, --standardize, or --idealize")
    if args.format and not args.output:
        raise ValueError("--format requires --output")
    if args.output and Path(args.output).resolve() == Path(args.input).resolve():
        raise ValueError("Input file cannot be overwritten")

    is_cif = filetype(args.input) == "cif"
    atoms = read(args.input, **({"store_tags": True} if is_cif else {}))
    declared = _declared_cif_symmetry(atoms.info) if is_cif else None
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
    if declared:
        output["cif_declared_symmetry"] = declared

        symbol = declared.get("hermann_mauguin")
        number = declared.get("number")
        description = str(symbol) if symbol is not None else "not named"
        if number is not None:
            description += f" (No. {number})"
        print(f"CIF-declared space group: {description}")
        if "hall" in declared:
            print(f"CIF-declared Hall: {declared['hall']}")
        if "crystal_system" in declared:
            print(f"CIF-declared crystal system: {declared['crystal_system']}")
        if "operations" in declared:
            print(f"CIF-declared symmetry operations: {len(declared['operations'])}")

    print(
        f"Space group{' (calculated)' if declared else ''}: "
        f"{data.international} (No. {data.number})"
    )
    print(f"Hall: {data.hall}; Point group: {data.pointgroup}")
    if declared and "number" in declared:
        try:
            declared_number = int(declared["number"])
        except (TypeError, ValueError):
            pass  # Retain and print the original CIF value as-is.
        else:
            if declared_number != int(data.number):
                print("Note: CIF-declared and calculated space-group numbers differ.")
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
        if declared and "operations" in declared:
            for index, operation in enumerate(declared["operations"], 1):
                print(f"CIF operation {index}: {operation}")
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
