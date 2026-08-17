"""
psa_crosscheck.py
==================
Independent validation of the pipeline's PyMOL-computed PSA_3D_Angstrom2:
recomputes the same polar-atom (O, N, and H bonded to O/N) solvent-accessible
surface area using RDKit's rdFreeSASA wrapper (a separate SASA implementation,
Bondi vdW radii, same 1.4 A probe radius) for every conformer in a chunk, and
reports the difference against the value PyMOL already stored on each record.

Does not touch or recompute anything else -- purely a read-only cross-check
of an already-produced chunk_N_analysis_top50.sdf file.

Usage:
    python psa_crosscheck.py --input-sdf chunk_0_analysis_top50.sdf --output-csv chunk_0_psa_crosscheck.csv
"""
import argparse
import csv
import os

from rdkit import Chem
from rdkit.Chem import rdFreeSASA

# Standard Bondi-ish vdW radii (Angstrom) by element -- FreeSASA's built-in
# classifyAtoms() only recognizes standard protein residue/atom names and
# silently returns all-zero radii for an arbitrary small molecule, so radii
# must be supplied explicitly here.
VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
       "S": 1.80, "Cl": 1.75, "Br": 1.85, "P": 1.80}


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-sdf", required=True)
    ap.add_argument("--output-csv", required=True)
    ap.add_argument("--probe-radius", type=float, default=1.4)
    return ap.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.input_sdf):
        raise FileNotFoundError(f"Input SDF not found: {args.input_sdf}")

    polar_query = Chem.MolFromSmarts("[#7,#8,$([#1]~[#7,#8])]").GetAtomWithIdx(0)
    opts = rdFreeSASA.SASAOpts()
    opts.probeRadius = args.probe_radius

    rows = []
    n_records = 0
    for mol in Chem.SDMolSupplier(args.input_sdf, removeHs=False, sanitize=True):
        if mol is None or not mol.HasProp("PSA_3D_Angstrom2") or not mol.HasProp("Source_ID"):
            continue
        n_records += 1
        radii = [VDW.get(atom.GetSymbol(), 1.70) for atom in mol.GetAtoms()]
        free_psa = rdFreeSASA.CalcSASA(mol, radii, confIdx=-1, query=polar_query, opts=opts)
        pymol_psa = mol.GetDoubleProp("PSA_3D_Angstrom2")
        pct_diff = 100 * (pymol_psa - free_psa) / free_psa if free_psa else None

        rows.append({
            "Source_ID": mol.GetProp("Source_ID"),
            "Molecule_Index": mol.GetIntProp("Molecule_Index") if mol.HasProp("Molecule_Index") else "",
            "Conformer_Rank": mol.GetIntProp("Conformer_Rank") if mol.HasProp("Conformer_Rank") else "",
            "Boltzmann_Weight": mol.GetDoubleProp("Boltzmann_Weight") if mol.HasProp("Boltzmann_Weight") else "",
            "pymol_psa": round(pymol_psa, 3),
            "freesasa_psa": round(free_psa, 3),
            "pct_diff": round(pct_diff, 3) if pct_diff is not None else "",
        })

        if n_records % 1000 == 0:
            print(f"  {n_records} conformers processed...")

    if rows:
        os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
        with open(args.output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    diffs = [r["pct_diff"] for r in rows if r["pct_diff"] != ""]
    print(f"\nProcessed {n_records} conformers -> {args.output_csv}")
    if diffs:
        print(f"Mean %% diff: {sum(diffs)/len(diffs):.2f}")
        print(f"Min/Max %% diff: {min(diffs):.2f} / {max(diffs):.2f}")


if __name__ == "__main__":
    main()
