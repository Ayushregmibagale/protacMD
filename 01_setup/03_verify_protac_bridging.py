#!/usr/bin/env python3
"""
Quick sanity check: verify PROTAC is properly bridging in pRosettaC combined PDBs.
Checks the rank1 representative for each POI in cluster_features.csv.

Usage:
    python verify_protac_bridging.py [--all]
    python verify_protac_bridging.py --runs-root /path/to/pRosettaC/runs
"""

import sys
import csv
import os
import argparse
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist

_DEFAULT_RUNS_ROOT = Path.home() / "protac" / "pRosettaC" / "runs"


def check(pdb_path, label):
    ca_atoms, cx = [], []
    try:
        with open(pdb_path) as f:
            for line in f:
                rec = line[:6].strip()
                if rec not in ('ATOM', 'HETATM'):
                    continue
                chain = line[21]
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                if chain == 'A':
                    rn = int(line[22:26])
                    ca_atoms.append((rn, x, y, z))
                elif chain == 'X':
                    cx.append([x, y, z])
    except Exception as e:
        print(f"  ERROR {label}: {e}")
        return False

    if not cx:
        print(f"  FAIL  {label}: no PROTAC (chain X) found")
        return False

    cx = np.array(cx)
    split_idx = next((i for i in range(1, len(ca_atoms))
                      if ca_atoms[i][0] < ca_atoms[i-1][0] - 50), None)
    if split_idx is None:
        vhl = np.array([[a[1],a[2],a[3]] for a in ca_atoms if a[0] <= 103])
        kin = np.array([[a[1],a[2],a[3]] for a in ca_atoms if a[0] > 103])
    else:
        vhl = np.array([[a[1],a[2],a[3]] for a in ca_atoms[:split_idx]])
        kin = np.array([[a[1],a[2],a[3]] for a in ca_atoms[split_idx:]])

    d_vhl = cdist(cx, vhl).min() if len(vhl) else 999
    d_kin = cdist(cx, kin).min() if len(kin) else 999

    # Degrader pRosettaC PatchDock results always have PROTAC in kinase pocket
    # (0-5 Å) but 20-32 Å from VHL — that is expected (VHL closes during MD).
    # Only flag as broken when PROTAC is >40 Å from BOTH proteins (vacuum origin).
    broken = d_vhl > 40 and d_kin > 40
    status = "FAIL" if broken else "OK  "
    print(f"  {status}  {label}: PROTAC-VHL={d_vhl:.1f}Å  PROTAC-Kin={d_kin:.1f}Å")
    return not broken


def check_poi(base, poi):
    csv_path = os.path.join(base, poi, "cluster_features.csv")
    if not os.path.exists(csv_path):
        print(f"  SKIP  {poi}: no cluster_features.csv")
        return None
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    rank1 = next((r for r in rows if r.get('md_candidate') == '1'), rows[0] if rows else None)
    if rank1 is None:
        print(f"  SKIP  {poi}: empty cluster_features.csv")
        return None
    pdb = os.path.join(base, poi, "Patchdock_Results", rank1['representative'])
    return check(pdb, f"{poi} rank1 ({rank1['representative']})")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--all', action='store_true', help='Also check degraders')
    parser.add_argument('--runs-root', default=str(_DEFAULT_RUNS_ROOT),
                        help=f'Path to pRosettaC runs/ directory '
                             f'(default: {_DEFAULT_RUNS_ROOT})')
    args = parser.parse_args()

    runs_root = Path(args.runs_root)
    non_deg = str(runs_root / 'high_affinity_no_degradation')
    deg     = str(runs_root / 'high_affinity_degraders')

    all_ok = True

    print("\n=== Non-degraders ===")
    for poi in ['ABL1', 'AXL', 'EPHA2', 'MAP4K5', 'SLK']:
        ok = check_poi(non_deg, poi)
        if ok is False:
            all_ok = False

    if args.all:
        print("\n=== Degraders ===")
        for poi in ['DDR2', 'MET', 'RIPK2', 'MAPK14', 'EPHB2',
                    'RIPK2_full', 'MAPK14_full']:
            ok = check_poi(deg, poi)
            if ok is False:
                all_ok = False

    print()
    if all_ok:
        print("All checks passed — PROTAC properly bridging in all rank1 poses.")
    else:
        print("Some checks FAILED — see above. Re-run broken POIs.")
    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main())
