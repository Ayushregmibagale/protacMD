# protacMD — GROMACS MD Setup for PROTAC Ternary Complexes

Scripts for preparing and validating GROMACS-ready input files from pRosettaC ternary complex poses. This is the structure-preparation stage that bridges the pRosettaC docking output and the GROMACS MD simulation pipeline.

---

## Background

After pRosettaC samples ternary complex poses (kinase–PROTAC–VHL), the raw combined PDB files need to be split and validated before GROMACS can use them:

- pRosettaC writes both VHL and the kinase on the **same chain A**, which GROMACS cannot handle.
- For N-lobe–truncated kinases (TNIK, MAP4K5, SLK), a duplicate VHL fragment can appear in the combined PDB and must be removed.
- The PROTAC (chain X) needs to be extracted separately for GAFF2 parameterisation.

These three scripts handle that preparation.

---

## Repository Structure

```
protacMD/
└── 01_setup/
    ├── 01_split_chains.py          # Split combined PDB → complex_protein.pdb + SJF8240_protac.pdb
    ├── 02_fix_complex_protein.py   # Remove duplicate VHL fragment (poi_min=104 kinases)
    └── 03_verify_protac_bridging.py # Sanity-check PROTAC bridging in all rank1 poses
```

---

## Prerequisites

| Dependency | Purpose |
|---|---|
| Python ≥ 3.9 | All scripts |
| NumPy | Distance calculations |
| SciPy | `cdist` in the bridging verifier |

```bash
pip install numpy scipy
```

The scripts expect the pRosettaC run directories to live under:
```
~/protac/pRosettaC/runs/
```
This path can be overridden where relevant (see usage below).

---

## Usage

### Step 1 — Split chains

Separates the pRosettaC combined PDB into a protein file and a PROTAC file for GROMACS.

```bash
python 01_setup/01_split_chains.py <combined_pdb> <vhl_max_res> <poi_min_res>

# Examples
python 01_split_chains.py start_rank1.pdb 103 553    # DDR2
python 01_split_chains.py start_rank1.pdb 103 473    # AXL
python 01_split_chains.py start_rank1.pdb 103 104    # TNIK / MAP4K5 / SLK
```

**Output (written next to the input PDB):**

| File | Contents |
|---|---|
| `complex_protein.pdb` | Chain A = VHL (res 1–103), D = kinase, B = ElonginB, C = ElonginC |
| `SJF8240_protac.pdb` | Chain X — PROTAC atoms only (79 heavy atoms for SJF8240) |

**Chain layout after splitting:**

```
Chain A  — VHL (res 1–103)
Chain B  — ElonginB
Chain C  — ElonginC
Chain D  — kinase (original or renumbered residues)
Chain X  — PROTAC (extracted to SJF8240_protac.pdb)
```

---

### Step 2 — Fix complex_protein.pdb (N-lobe truncated kinases only)

For POIs where `poi_min = 104` (TNIK, MAP4K5, SLK), pRosettaC can write a duplicate VHL fragment on chain A. This script removes it and reassigns the correct chains.

```bash
python 01_setup/02_fix_complex_protein.py <POI> [rank] [--poi-min 104]

# Examples
python 02_fix_complex_protein.py TNIK
python 02_fix_complex_protein.py MAP4K5 1
python 02_fix_complex_protein.py SLK 2 --poi-min 104
```

The original `complex_protein.pdb` is backed up as `complex_protein.pdb.bak` before overwriting.

> **When to run this:** only needed for kinases whose domain starts at residue 104, i.e., those where the N-lobe numbering (1–103) collides with VHL's residue range.

---

### Step 3 — Verify PROTAC bridging

Quick sanity check that the PROTAC is geometrically reasonable in each rank1 representative pose.

```bash
python 01_setup/03_verify_protac_bridging.py                  # non-degraders only
python 01_setup/03_verify_protac_bridging.py --all            # include degraders
python 01_setup/03_verify_protac_bridging.py --runs-root /path/to/pRosettaC/runs
```

Output per POI:

```
=== Non-degraders ===
  OK    ABL1  rank1 (combined_12_3_0001.pdb): PROTAC-VHL=8.3Å  PROTAC-Kin=4.1Å
  OK    AXL   rank1 (combined_45_6_0001.pdb): PROTAC-VHL=9.7Å  PROTAC-Kin=3.8Å
  FAIL  SLK   rank1 (combined_7_1_0001.pdb):  PROTAC-VHL=52.1Å PROTAC-Kin=48.3Å
```

**Verdict logic:**

| Status | Criterion |
|---|---|
| `OK` | PROTAC is within 40 Å of at least one protein |
| `FAIL` | PROTAC > 40 Å from **both** VHL and kinase (likely a vacuum-origin artefact) |

> Note: distances of 20–32 Å between PROTAC and VHL are normal in degrader poses — the VHL arm closes around the PROTAC during MD equilibration. Only the extreme `> 40 Å` case indicates a broken pose.

---

## Typical Workflow

```
pRosettaC output
       │
       ▼
01_split_chains.py          ← run for each rank1/2/3 pose
       │
       ▼  (TNIK/MAP4K5/SLK only)
02_fix_complex_protein.py   ← remove duplicate VHL fragment
       │
       ▼
03_verify_protac_bridging.py ← confirm geometry before investing in MD
       │
       ▼
GROMACS pdb2gmx / GAFF2 parameterisation → MD simulation
```

---

## Notes

- All scripts expect the pRosettaC combined PDB convention: VHL on chain A (res 1–103), kinase continuing on chain A, ElonginB on chain B, ElonginC on chain C, PROTAC on chain X.
- `01_split_chains.py` skips residues that fall in the gap between `vhl_max` and `poi_min` (e.g., the VHL–kinase linker residues in the combined file).
- For renumbered full-kinase runs (RIPK2_full, MAPK14_full, residues offset to 200+), use the actual start residue as `poi_min`, e.g. `200`.
