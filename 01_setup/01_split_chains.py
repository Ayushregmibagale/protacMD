#!/usr/bin/env python3
"""
Split pRosettaC combined PDB into GROMACS-ready files.

pRosettaC concatenates VHL (chain A, res 1-103) and kinase (chain A, res poi_min+)
on the same chain. GROMACS needs them as separate chains.

Outputs (written to the same directory as the input PDB):
  complex_protein.pdb  — chains A (VHL), B (ElonginB), C (ElonginC), D (kinase)
  SJF8240_protac.pdb   — chain X (PROTAC, 79 atoms) for GAFF2 parameterisation

Usage:
    python split_chains.py <combined_pdb> <vhl_max_res> <poi_min_res>

Example for DDR2:
    python split_chains.py start_rank1.pdb 103 553
"""

import sys
import os
from pathlib import Path


def split(pdb_path, vhl_max, poi_min):
    pdb_path = Path(pdb_path)
    out_dir = pdb_path.parent

    protein_lines = []
    protac_lines = []
    prev_chain = None

    with open(pdb_path) as f:
        raw = f.readlines()

    for line in raw:
        rec = line[:6].strip()

        if rec in ('ATOM', 'HETATM'):
            chain = line[21]
            resnum = int(line[22:26])

            if chain == 'X':
                # PROTAC — extract separately
                protac_lines.append(line)
                continue

            if chain == 'A':
                if resnum <= vhl_max:
                    new_chain = 'A'   # VHL stays as chain A
                elif resnum >= poi_min:
                    new_chain = 'D'   # kinase becomes chain D
                else:
                    continue          # gap residues — skip

                # Insert TER when chain changes
                if prev_chain is not None and new_chain != prev_chain:
                    protein_lines.append(f'TER\n')

                line = line[:21] + new_chain + line[22:]

            else:
                new_chain = chain     # B, C unchanged

            if prev_chain is not None and new_chain != prev_chain and rec == 'ATOM':
                if not protein_lines or not protein_lines[-1].startswith('TER'):
                    protein_lines.append('TER\n')

            protein_lines.append(line)
            prev_chain = new_chain

        elif rec == 'TER':
            if protein_lines and not protein_lines[-1].startswith('TER'):
                protein_lines.append('TER\n')
        elif rec == 'END':
            pass   # will add at end

    protein_lines.append('TER\n')
    protein_lines.append('END\n')
    protac_lines.append('END\n')

    protein_out = out_dir / 'complex_protein.pdb'
    protac_out = out_dir / 'SJF8240_protac.pdb'

    with open(protein_out, 'w') as f:
        f.writelines(protein_lines)

    with open(protac_out, 'w') as f:
        f.writelines(protac_lines)

    # Report
    chain_residues = {}
    for line in protein_lines:
        if line.startswith('ATOM'):
            c = line[21]
            r = int(line[22:26])
            chain_residues.setdefault(c, set()).add(r)

    print(f'Wrote: {protein_out}')
    for c in sorted(chain_residues):
        rr = sorted(chain_residues[c])
        label = {'A': 'VHL', 'B': 'ElonginB', 'C': 'ElonginC', 'D': 'DDR2/kinase'}.get(c, c)
        print(f'  Chain {c} ({label}): res {rr[0]}-{rr[-1]}  ({len(rr)} residues)')

    print(f'Wrote: {protac_out}')
    print(f'  PROTAC atoms: {sum(1 for l in protac_lines if l.startswith("HETATM"))}')

    return protein_out, protac_out


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    split(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
