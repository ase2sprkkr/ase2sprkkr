"""Downsample native SPR-KKR BSF/SPC files to reproducible parser/plot fixtures.

All channel blocks are sliced in their actual native ordering. Unsupported
headers are rejected rather than producing truncated but plausible files.
Do not use these reduced spectra for quantitative scientific comparisons.
"""
from __future__ import annotations

import argparse
import base64
from pathlib import Path
import re

import numpy as np


def field(header, name):
    values = re.findall(rf"(?m)^\s*{re.escape(name)}\s+(\d+)\s*$", header)
    if not values or len(set(values)) != 1:
        raise ValueError(f"Missing/inconsistent {name}: {values}")
    return int(values[0])


def set_field(header, name, n):
    pattern = rf"(?m)^(\s*{re.escape(name)}\s+)\d+(\s*)$"
    header, hits = re.subn(pattern, lambda m: m[1] + str(n) + m[2], header)
    if hits == 0:
        raise ValueError(f"Cannot rewrite {name}")
    return header


def regular_subset(n, max_points, *, endpoints=False):
    if max_points < 1:
        raise ValueError('Sampling limit must be positive')
    if endpoints and n > 1 and max_points < 2:
        raise ValueError('Endpoint-preserving grids require at least two points')
    if n <= max_points:
        return np.arange(n)
    if endpoints:
        # Only an exact uniform sub-grid preserves the original CONST-E
        # VECK1/VECK2 endpoint vectors without modifying physical geometry.
        count = max((k for k in range(2, max_points + 1)
                     if (n - 1) % (k - 1) == 0), default=2)
        return np.arange(0, n, (n - 1) // (count - 1))
    if max_points == 1:
        return np.array([0])
    step = max(1, (n - 1) // (max_points - 1))
    return np.arange(0, n, step)[:max_points]


def compact_bsf(text, max_e=8, max_k=12):
    lines = text.splitlines(keepends=True)
    marker = max((i for i, v in enumerate(lines) if v.strip() == "#" * 80), default=-1)
    if marker < 0:
        raise ValueError("No BSF raw-data separator")
    header = "".join(lines[:marker + 1])
    numbers = np.fromstring("".join(lines[marker + 1:]), sep=" ")
    kind = re.search(r"(?m)^KEYWORD\s+(\S+)", header)
    mode = re.search(r"(?m)^MODE\s+(\S+)", header)
    if not kind or not mode or kind[1] not in ("BSF", "BSF-SPOL", "BSF-SPN"):
        raise ValueError("Unknown BSF header")
    keyword = kind[1]
    nq = field(header, "NQ_eff")
    channels = 3 if keyword == "BSF" else 4
    if mode[1] == "EK-REL":
        ne, nk = field(header, "NE"), field(header, "NK")
        if numbers.size != channels * ne * nk * nq:
            raise ValueError(f"Wrong raw size {numbers.size} vs {channels * ne * nk * nq}")
        energy = regular_subset(ne, max_e)
        # Preserve every K-direction junction, including the last point of
        # each segment, when reducing a multi-segment band path.
        boundaries = [int(v) for v in re.findall(
            r"(?m)^\s*INDKDIR\s+(\d+)\s*$", header)]
        ndir = field(header, "NKDIR")
        if len(boundaries) != ndir or not boundaries or boundaries[-1] != nk:
            raise ValueError("Inconsistent K-direction boundaries")
        if any(a >= b for a, b in zip((0, *boundaries[:-1]), boundaries)):
            raise ValueError("Nonincreasing K-direction boundaries")
        junctions = {0, nk - 1}
        for begin, end in zip((0, *boundaries[:-1]), boundaries):
            junctions.update((begin, end - 1))
        k = np.array(sorted(set(regular_subset(nk, max_k)) | junctions))
        if keyword == "BSF":
            first = numbers[:2 * ne * nq * nk].reshape(ne, 2, nq, nk)
            second = numbers[2 * ne * nq * nk:].reshape(ne, nq, nk)
            raw = np.concatenate((first[np.ix_(energy, range(2), range(nq), k)].ravel(),
                                  second[np.ix_(energy, range(nq), k)].ravel()))
        else:
            raw = numbers.reshape(ne, 4, nq, nk)[
                np.ix_(energy, range(4), range(nq), k)].ravel()
        header = set_field(set_field(header, "NE", len(energy)), "NK", len(k))
        entries = header.splitlines(keepends=True)
        after_nk = max(i for i, v in enumerate(entries)
                       if re.match(r"^\s*NK\s+\d+\s*$", v)) + 1
        k_lines = entries[after_nk:-1]
        if len(k_lines) != nk:
            raise ValueError(f"K-grid has {len(k_lines)} rows, expected {nk}")
        entries[after_nk:-1] = [k_lines[i] for i in k]
        header = "".join(entries)
        positions = iter(int(np.searchsorted(k, edge - 1, side="right"))
                         for edge in boundaries)
        header = re.sub(
            r"(?m)^(\s*INDKDIR\s+)\d+(\s*)$",
            lambda match: match[1] + str(next(positions)) + match[2],
            header,
        )
        m = re.search(r"(?m)^(EMIN EMAX\s+)(\S+)(\s+)(\S+)(\s*)$", header)
        if not m:
            raise ValueError("No EMIN/EMAX header")
        emin, emax = float(m[2]), float(m[4])
        upper = emin + (emax - emin) * energy[-1] / (ne - 1) if ne > 1 else emax
        header = header[:m.start()] + m[1] + m[2] + m[3] + f"{upper:.10g}" + m[5] + header[m.end():]
    elif mode[1] == "CONST-E":
        nk1, nk2 = field(header, "NK1"), field(header, "NK2")
        if numbers.size != channels * nk1 * nq * nk2:
            raise ValueError(f"Wrong raw size {numbers.size} vs {channels * nk1 * nq * nk2}")
        first = regular_subset(nk1, max_k, endpoints=True)
        second = regular_subset(nk2, max_k, endpoints=True)
        if keyword == "BSF":
            spin = numbers[:2 * nk1 * nq * nk2].reshape(2, nk1, nq, nk2)
            total = numbers[2 * nk1 * nq * nk2:].reshape(nk1, nq, nk2)
            raw = np.concatenate((spin[np.ix_(range(2), first, range(nq), second)].ravel(),
                                  total[np.ix_(first, range(nq), second)].ravel()))
        else:
            raw = numbers.reshape(4, nk1, nq, nk2)[
                np.ix_(range(4), first, range(nq), second)].ravel()
        header = set_field(set_field(header, "NK1", len(first)), "NK2", len(second))
    else:
        raise ValueError(f"Unknown BSF MODE {mode[1]}")
    return header + "\n".join(f"{v:.9g}" for v in raw) + "\n"


def compact_spc(text, max_e=12, max_k=16):
    lines = text.splitlines(keepends=True)
    marker = next((i for i, line in enumerate(lines)
                   if line.startswith("#######################")), None)
    if marker is None:
        raise ValueError("Missing SPC data marker")
    header = "".join(lines[:marker + 1])
    first = lines[marker + 1].strip().split()
    width = len(first)
    if width not in (7, 8):
        raise ValueError(f"SPC width {width} is not 7 or 8")
    values = np.fromstring("".join(lines[marker + 1:]), sep=" ")
    ne, np_, nt = (field(header, n) for n in ("NE", "NP", "NT"))
    if width == 8:
        if values.size != ne * np_ * nt * 8:
            raise ValueError("Malformed energy SPC matrix")
        e, p, t = (regular_subset(n, mx) for n, mx in
                   ((ne, max_e), (np_, max_k), (nt, max_k)))
        data = values.reshape(ne, np_, nt, 8)[np.ix_(e, p, t, range(8))]
        header = set_field(header, "NE", len(e))
    else:
        if values.size != np_ * nt * 7:
            raise ValueError("Malformed kx_ky SPC matrix")
        p, t = regular_subset(np_, max_k), regular_subset(nt, max_k)
        data = values.reshape(np_, nt, 7)[np.ix_(p, t, range(7))]
    header = set_field(set_field(header, "NP", len(p)), "NT", len(t))
    rows = data.reshape(-1, width)
    # The native SPC reader expects columnar lines with a leading blank
    # and scientific notation (as written by SPR-KKR).
    return header + "".join("".join(f" {v:16.5E}" for v in row) + "\n" for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--energies", type=int, default=8)
    parser.add_argument("--kpoints", type=int, default=12)
    parser.add_argument("--emit-base64", action="store_true",
                        help="Emit generated fixture as base64 for review/transfer")
    args = parser.parse_args()
    source = args.input.read_text()
    if args.input.suffix == ".bsf":
        result = compact_bsf(source, args.energies, args.kpoints)
    elif args.input.suffix == ".spc":
        result = compact_spc(source, args.energies, args.kpoints)
    else:
        raise ValueError("Only .bsf/.spc inputs supported")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(result)
    print(f"COMPACT {args.input.name}: {len(source)} -> {len(result)} bytes", flush=True)
    if args.emit_base64:
        print(f"COMPACT_BEGIN {args.output.name}")
        print(base64.b64encode(result.encode()).decode())
        print(f"COMPACT_END {args.output.name}")


if __name__ == "__main__":
    main()
