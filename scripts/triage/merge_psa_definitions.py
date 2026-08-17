"""
merge_psa_definitions.py
=========================
Combines the PyMOL-computed 3D-PSA (already Boltzmann-weighted per molecule
in chunk_N_features.csv) with an independent RDKit/FreeSASA recomputation
(per-conformer, from psa_crosscheck.py) into three per-molecule PSA/Shielding
Index values: pymol-only, freesasa-only, and their average.

Boltzmann weights are renormalized within whatever conformer set is present
in the crosscheck CSV (e.g. the top-50 conformers/molecule), since that set's
weights don't sum to 1 on their own -- see boltzmann_average() in
analyse_conformers.py for the same renormalization applied there.

Usage:
    python merge_psa_definitions.py \\
        --crosscheck-csv psa_crosscheck_all.csv \\
        --features-csv chunk_0_features.csv chunk_1_features.csv ... \\
        --output-csv molecules_with_psa_variants.csv
"""
import argparse
import csv
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crosscheck-csv", required=True,
                     help="Output of psa_crosscheck.py (per-conformer PyMOL + FreeSASA PSA)")
    ap.add_argument("--features-csv", required=True, nargs="+",
                     help="One or more chunk_N_features.csv files (for TPSA_2D and other per-molecule context)")
    ap.add_argument("--output-csv", required=True)
    return ap.parse_args()


def main():
    args = parse_args()

    groups = defaultdict(list)
    with open(args.crosscheck_csv, newline="") as f:
        for row in csv.DictReader(f):
            groups[row["Source_ID"]].append((
                float(row["Boltzmann_Weight"]), float(row["pymol_psa"]), float(row["freesasa_psa"])
            ))
    print(f"Loaded conformer-level PSA data for {len(groups)} compounds from {args.crosscheck_csv}")

    boltz = {}
    for sid, entries in groups.items():
        wsum = sum(w for w, _, _ in entries)
        pymol_avg = sum(w * p for w, p, _ in entries) / wsum
        free_avg = sum(w * f for w, _, f in entries) / wsum
        boltz[sid] = (pymol_avg, free_avg)

    rows = []
    for features_path in args.features_csv:
        with open(features_path, newline="") as f:
            for row in csv.DictReader(f):
                sid = row["Source_ID"]
                if sid not in boltz:
                    continue
                pymol_psa, freesasa_psa = boltz[sid]
                avg_psa = (pymol_psa + freesasa_psa) / 2
                tpsa_2d = float(row["TPSA_2D_Angstrom2"])
                rows.append({
                    "Source_ID": sid,
                    "Rgyr_Angstrom_boltzmann": row["Rgyr_Angstrom_boltzmann"],
                    "Intramolecular_HBonds_boltzmann": row["Intramolecular_HBonds_boltzmann"],
                    "TPSA_2D_Angstrom2": round(tpsa_2d, 3),
                    "PSA_3D_pymol": round(pymol_psa, 3),
                    "PSA_3D_freesasa": round(freesasa_psa, 3),
                    "PSA_3D_averaged": round(avg_psa, 3),
                    "Shielding_Index_pymol": round(tpsa_2d - pymol_psa, 3),
                    "Shielding_Index_freesasa": round(tpsa_2d - freesasa_psa, 3),
                    "Shielding_Index_averaged": round(tpsa_2d - avg_psa, 3),
                    "Effective_Sample_Size": row["Effective_Sample_Size"],
                })

    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} compounds to {args.output_csv}")


if __name__ == "__main__":
    main()
