#!/usr/bin/env python3
"""
Fix complex_protein.pdb for POIs where poi_min=104 (TNIK, MAP4K5, SLK).

pRosettaC combined output contains a duplicate VHL fragment in a second chain-A
segment. split_chains.py keeps it as chain A since res < 104. This script
rebuilds complex_protein.pdb with the duplicate removed.

Usage:
    python fix_complex_protein.py <POI> [rank]
    python fix_complex_protein.py TNIK
    python fix_complex_protein.py MAP4K5 1
"""

import sys
import shutil
from pathlib import Path

POI  = sys.argv[1]
RANK = int(sys.argv[2]) if len(sys.argv) > 2 else 1
POI_MIN = 104  # first kinase residue (poi_min used in split_chains.py)

ROOT     = Path("~/protac-md").expanduser()
STRUCDIR = ROOT / f"simulations/{POI}/structure/rank{RANK}"
SRC      = STRUCDIR / f"start_rank{RANK}.pdb"
DEST     = STRUCDIR / "complex_protein.pdb"

shutil.copy(DEST, DEST.with_suffix('.pdb.bak'))
print(f"[{POI}] Backed up → {DEST.with_suffix('.pdb.bak')}")

# Parse chain-A segments (split by TER) and other chains
chain_a_segments = []
other_chains = {'B': [], 'C': []}
current_segment = []

with open(SRC) as f:
    raw = f.readlines()

in_chain_a = False
for line in raw:
    rec = line[:6].strip()
    if rec in ('ATOM', 'HETATM'):
        chain  = line[21]
        resnum = int(line[22:26])
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

print(f"[{POI}] Found {len(chain_a_segments)} chain-A segments:")
for i, seg in enumerate(chain_a_segments):
    res_nums = [int(l[22:26]) for l in seg]
    print(f"        Segment {i+1}: {len(seg)} atoms, res {min(res_nums)}–{max(res_nums)}")

# VHL = the segment that starts at residue 1
vhl_seg = next(s for s in chain_a_segments if int(s[0][22:26]) == 1
               and all(int(l[22:26]) <= 103 for l in s))

# Kinase = all atoms with res >= POI_MIN from non-VHL segments
kin_lines = []
dup_count = 0
for seg in chain_a_segments:
    if seg is vhl_seg:
        continue
    for l in seg:
        rn = int(l[22:26])
        if rn >= POI_MIN:
            kin_lines.append(l)
        else:
            dup_count += 1

print(f"[{POI}] VHL (chain A):   {len(vhl_seg)} atoms (res 1-103)")
print(f"[{POI}] Kinase (chain D): {len(kin_lines)} atoms (res {POI_MIN}+)")
print(f"[{POI}] Duplicate (skip): {dup_count} atoms")

def rename_chain(lines, ch):
    return [l[:21] + ch + l[22:] for l in lines]

out_lines = []
out_lines += rename_chain(vhl_seg, 'A');   out_lines.append('TER\n')
out_lines += rename_chain(kin_lines, 'D'); out_lines.append('TER\n')
out_lines += other_chains['B'];            out_lines.append('TER\n')
out_lines += other_chains['C'];            out_lines.append('TER\n')
out_lines.append('END\n')

with open(DEST, 'w') as f:
    f.writelines(out_lines)

# Summary
chain_res = {}
for l in out_lines:
    if l[:4] == 'ATOM' or l[:6] == 'HETATM':
        ch = l[21]; rn = int(l[22:26])
        chain_res.setdefault(ch, []).append(rn)

print(f"[{POI}] Fixed complex_protein.pdb:")
for ch in sorted(chain_res):
    rns = chain_res[ch]
    print(f"        Chain {ch}: {len(rns)} atoms, res {min(rns)}–{max(rns)}")
print(f"[{POI}] Total: {sum(len(v) for v in chain_res.values())} atoms → {DEST}")
