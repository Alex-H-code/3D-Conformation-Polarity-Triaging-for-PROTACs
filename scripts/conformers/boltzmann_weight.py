"""
boltzmann_weight.py
====================
Computes per-conformer Boltzmann weights from MMFF energies and writes them
back as a 'Boltzmann_Weight' SDF property. This is the step analyse_conformers.py
was silently assuming already happened -- nothing in the original pipeline
produced it.

w_i = exp(-E_i / RT) / sum_j exp(-E_j / RT), grouped by Molecule_Index,
at a fixed temperature (default 298 K).

With --top-n, only the N highest-weighted conformers per molecule are
written out (weights are computed from the full conformer set first, so
top-N selection doesn't change the weight values themselves).

Usage:
    python boltzmann_weight.py --input-sdf conformers.sdf --output-sdf conformers_weighted.sdf
    python boltzmann_weight.py --input-sdf conformers.sdf --output-sdf conformers_weighted.sdf --top-n 50
"""
import argparse
import os

import numpy as np
from rdkit import Chem
from rdkit.Chem import SDWriter

R_KCAL = 0.0019872041  # kcal / (mol*K)


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-sdf", required=True)
    ap.add_argument("--output-sdf", required=True)
    ap.add_argument("--temperature", type=float, default=298.0, help="Kelvin")
    ap.add_argument("--energy-prop", default="MMFF_Energy",
                     help="SDF property holding the per-conformer energy (kcal/mol)")
    ap.add_argument("--top-n", type=int, default=None,
                     help="Keep only the N highest-weighted conformers per molecule "
                          "(default: keep all)")
    return ap.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.input_sdf):
        raise FileNotFoundError(f"Input SDF not found: {args.input_sdf}")

    supplier = Chem.SDMolSupplier(args.input_sdf, removeHs=False, sanitize=True)
    mols = [m for m in supplier if m is not None]
    print(f"Loaded {len(mols)} conformers from {args.input_sdf}")

    missing = [m for m in mols if not m.HasProp("Molecule_Index") or not m.HasProp(args.energy_prop)]
    if missing:
        raise ValueError(
            f"{len(missing)} conformers are missing 'Molecule_Index' or '{args.energy_prop}' -- "
            "was this SDF produced by generate_conformers.py?"
        )

    # Group conformer indices by Molecule_Index
    groups = {}
    for i, mol in enumerate(mols):
        groups.setdefault(mol.GetIntProp("Molecule_Index"), []).append(i)

    kT = R_KCAL * args.temperature
    keep_idxs = set()
    for mol_idx, idxs in groups.items():
        energies = np.array([mols[i].GetDoubleProp(args.energy_prop) for i in idxs])
        rel_energies = energies - energies.min()  # numerically stable softmin
        boltz = np.exp(-rel_energies / kT)
        weights = boltz / boltz.sum()
        for i, w in zip(idxs, weights):
            mols[i].SetDoubleProp("Boltzmann_Weight", float(w))

        if args.top_n is not None:
            ranked = sorted(zip(idxs, weights), key=lambda pair: pair[1], reverse=True)
            keep_idxs.update(i for i, _ in ranked[:args.top_n])
        else:
            keep_idxs.update(idxs)

    kept_mols = [mol for i, mol in enumerate(mols) if i in keep_idxs]

    os.makedirs(os.path.dirname(args.output_sdf) or ".", exist_ok=True)
    writer = SDWriter(args.output_sdf)
    for mol in kept_mols:
        writer.write(mol)
    writer.close()

    if args.top_n is not None:
        print(f"Weighted {len(mols)} conformers across {len(groups)} molecules "
              f"(T={args.temperature} K), kept top {args.top_n}/molecule "
              f"-> {len(kept_mols)} conformers -> {args.output_sdf}")
    else:
        print(f"Weighted {len(mols)} conformers across {len(groups)} molecules "
              f"(T={args.temperature} K) -> {args.output_sdf}")


if __name__ == "__main__":
    main()
