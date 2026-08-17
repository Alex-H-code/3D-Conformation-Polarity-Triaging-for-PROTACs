"""
filter_top_n.py
================
Post-hoc filter: keeps only the N highest-weighted conformers per molecule
from an already-produced SDF, without recomputing anything. Useful for
trimming a weighted/analysis SDF that was written before --top-n filtering
was applied (or applied with the wrong N) instead of rerunning the whole
descriptor calculation.

Two-pass streamed, matching the pattern in analyse_conformers.py: a cheap
first pass (sanitize=False, removeHs=True) reads just Molecule_Index and the
weight property per record to decide what to keep, then a second pass writes
out only the kept records at full fidelity -- no need to hold every
conformer in memory at once.

Usage:
    python filter_top_n.py --input-sdf chunk_0_analysis.sdf --output-sdf chunk_0_analysis_top50.sdf --top-n 50
"""
import argparse
import os

from rdkit import Chem


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-sdf", required=True)
    ap.add_argument("--output-sdf", required=True)
    ap.add_argument("--top-n", type=int, required=True)
    ap.add_argument("--weight-prop", default="Boltzmann_Weight")
    ap.add_argument("--molecule-prop", default="Molecule_Index")
    return ap.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.input_sdf):
        raise FileNotFoundError(f"Input SDF not found: {args.input_sdf}")

    groups = {}
    n_records = 0
    for i, mol in enumerate(Chem.SDMolSupplier(args.input_sdf, removeHs=True, sanitize=False)):
        n_records = i + 1
        if mol is None or not mol.HasProp(args.molecule_prop) or not mol.HasProp(args.weight_prop):
            continue
        groups.setdefault(mol.GetIntProp(args.molecule_prop), []).append(
            (i, mol.GetDoubleProp(args.weight_prop)))

    keep_idxs = set()
    for mol_idx, pairs in groups.items():
        ranked = sorted(pairs, key=lambda p: p[1], reverse=True)
        keep_idxs.update(i for i, _ in ranked[:args.top_n])

    print(f"Loaded {n_records} conformers across {len(groups)} molecules from {args.input_sdf}")
    print(f"Keeping top {args.top_n}/molecule -> {len(keep_idxs)} conformers")

    os.makedirs(os.path.dirname(args.output_sdf) or ".", exist_ok=True)
    writer = Chem.SDWriter(args.output_sdf)
    kept = 0
    for i, mol in enumerate(Chem.SDMolSupplier(args.input_sdf, removeHs=False, sanitize=True)):
        if mol is None or i not in keep_idxs:
            continue
        writer.write(mol)
        kept += 1
    writer.close()

    print(f"Wrote {kept} conformers to {args.output_sdf}")


if __name__ == "__main__":
    main()
