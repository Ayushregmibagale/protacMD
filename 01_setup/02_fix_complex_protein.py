#!/usr/bin/env python3
"""
Fix complex_protein.pdb for POIs where poi_min=104 (TNIK, MAP4K5, SLK).

pRosettaC combined output contains a duplicate VHL fragment in a second chain-A
segment. split_chains.py keeps it as chain A since res < 104. This script
rebuilds complex_protein.pdb with the duplicate removed.

Usage:
    python fix_complex_protein.py <POI> [rank] [--poi-min 104]
    python fix_complex_protein.py TNIK
    python fix_complex_protein.py MAP4K5 1 --poi-min 104
"""

import sys
import argparse
import shutil
from pathlib import Path


def rename_chain(lines, ch):
    return [l[:21] + ch + l[22:] for l in lines]


def fix(poi, rank, poi_min):
    root     = Path("~/protac-md").expanduser()
    strucdir = root / f"simulations/{poi}/structure/rank{rank}"
    src      = strucdir / f"start_rank{rank}.pdb"
    dest     = strucdir / "complex_protein.pdb"

    if not src.exists():
        sys.exit(f"ERROR: source PDB not found: {src}")
    if not dest.exists():
        sys.exit(f"ERROR: complex_protein.pdb not found (run split_chains.py first): {dest}")

    shutil.copy(dest, dest.with_suffix('.pdb.bak'))
    print(f"[{poi}] Backed up → {dest.with_suffix('.pdb.bak')}")

    # Parse chain-A segments (split by TER) and other chains
    chain_a_segments = []
    other_chains = {'B': [], 'C': []}
    current_segment = []
    in_chain_a = False

    with open(src) as f:
        raw = f.readlines()

    for line in raw:
        rec = line[:6].strip()
        if rec in ('ATOM', 'HETATM'):
            chain  = line[21]
            if chain == 'X':
                continue
            elif chain == 'A':
                current_segment.append(line)
                in_chain_a = True
            elif chain in ('B', 'C'):
                other_chains[chain].append(line)
                in_chain_a = False
        elif rec == 'TER':
            if in_chain_a and current_segment:
                chain_a_segments.append(current_segment)
                current_segment = []
            in_chain_a = False

    if current_segment:
        chain_a_segments.append(current_segment)

    print(f"[{poi}] Found {len(chain_a_segments)} chain-A segments:")
    for i, seg in enumerate(chain_a_segments):
        res_nums = [int(l[22:26]) for l in seg]
        print(f"        Segment {i+1}: {len(seg)} atoms, res {min(res_nums)}-{max(res_nums)}")

    # VHL = the segment that starts at residue 1 and stays within res 1-103
    vhl_seg = next(
        (s for s in chain_a_segments
         if int(s[0][22:26]) == 1 and all(int(l[22:26]) <= 103 for l in s)),
        None,
    )
    if vhl_seg is None:
        sys.exit(f"[{poi}] ERROR: could not identify VHL segment (expected chain-A segment starting at res 1)")

    # Kinase = all atoms with res >= poi_min from non-VHL segments
    kin_lines = []
    dup_count = 0
    for seg in chain_a_segments:
        if seg is vhl_seg:
            continue
        for l in seg:
            rn = int(l[22:26])
            if rn >= poi_min:
                kin_lines.append(l)
            else:
                dup_count += 1

    print(f"[{poi}] VHL (chain A):    {len(vhl_seg)} atoms (res 1-103)")
    print(f"[{poi}] Kinase (chain D): {len(kin_lines)} atoms (res {poi_min}+)")
    print(f"[{poi}] Duplicate (skip): {dup_count} atoms")

    out_lines = []
    out_lines += rename_chain(vhl_seg, 'A');   out_lines.append('TER\n')
    out_lines += rename_chain(kin_lines, 'D'); out_lines.append('TER\n')
    out_lines += other_chains['B'];            out_lines.append('TER\n')
    out_lines += other_chains['C'];            out_lines.append('TER\n')
    out_lines.append('END\n')

    with open(dest, 'w') as f:
        f.writelines(out_lines)

    chain_res = {}
    for l in out_lines:
        if l[:4] == 'ATOM' or l[:6] == 'HETATM':
            ch = l[21]; rn = int(l[22:26])
            chain_res.setdefault(ch, []).append(rn)

    print(f"[{poi}] Fixed complex_protein.pdb:")
    for ch in sorted(chain_res):
        rns = chain_res[ch]
        print(f"        Chain {ch}: {len(rns)} atoms, res {min(rns)}-{max(rns)}")
    print(f"[{poi}] Total: {sum(len(v) for v in chain_res.values())} atoms → {dest}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('poi',  help='POI name, e.g. TNIK')
    parser.add_argument('rank', nargs='?', type=int, default=1, help='Rank (default 1)')
    parser.add_argument('--poi-min', type=int, default=104,
                        help='First kinase residue number (default 104)')
    args = parser.parse_args()
    fix(args.poi, args.rank, args.poi_min)


if __name__ == '__main__':
    main()
