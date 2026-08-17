"""
compute_weighted_score.py
==========================
Computes a continuous triage score per compound from three z-scored 3D
descriptors: score = [w_s*z(Shielding_Index) - w_r*z(Rgyr) + w_h*z(IMHB)]
/ (w_s + w_r + w_h), z-scored against the full input population. Weights
default to 0.5/0.3/0.2 (shielding-dominant), matching the composite used
in the interactive triage tool -- adjust via CLI flags for sensitivity
analysis.

Input is expected to be the output of merge_psa_definitions.py (for a
specific PSA definition) or any CSV with the required columns.

Usage:
    python compute_weighted_score.py \\
        --input-csv molecules_with_psa_variants.csv \\
        --psa-definition averaged \\
        --output-csv scored_compounds.csv \\
        --top-n 200 --top-n-ids-out top200_ids.txt
"""
import argparse
import csv
import statistics


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--psa-definition", choices=["pymol", "freesasa", "averaged"], default="averaged",
                     help="Which Shielding_Index_<def> column to use as the shielding term")
    ap.add_argument("--weight-shielding", type=float, default=0.5)
    ap.add_argument("--weight-rgyr", type=float, default=0.3)
    ap.add_argument("--weight-hbonds", type=float, default=0.2)
    ap.add_argument("--output-csv", required=True)
    ap.add_argument("--top-n", type=int, default=None, help="Also report a top-N cutoff")
    ap.add_argument("--top-n-ids-out", default=None, help="Optional: write top-N Source_IDs, one per line")
    return ap.parse_args()


def zscore(values):
    m = statistics.mean(values)
    sd = statistics.pstdev(values)
    return [(v - m) / sd if sd else 0.0 for v in values], m, sd


def main():
    args = parse_args()

    with open(args.input_csv, newline="") as f:
        rows = list(csv.DictReader(f))
    print(f"Loaded {len(rows)} compounds from {args.input_csv}")

    shield_col = f"Shielding_Index_{args.psa_definition}"
    shield = [float(r[shield_col]) for r in rows]
    rgyr = [float(r["Rgyr_Angstrom_boltzmann"]) for r in rows]
    hbonds = [float(r["Intramolecular_HBonds_boltzmann"]) for r in rows]

    z_shield, *_ = zscore(shield)
    z_rgyr, *_ = zscore(rgyr)
    z_hbonds, *_ = zscore(hbonds)

    total_w = args.weight_shielding + args.weight_rgyr + args.weight_hbonds
    for i, r in enumerate(rows):
        raw = (args.weight_shielding * z_shield[i]
               - args.weight_rgyr * z_rgyr[i]
               + args.weight_hbonds * z_hbonds[i])
        r["Score"] = round(raw / total_w, 4) if total_w else 0.0

    ranked = sorted(rows, key=lambda r: -r["Score"])
    for i, r in enumerate(ranked, start=1):
        r["Rank"] = i

    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(ranked[0].keys()))
        writer.writeheader()
        writer.writerows(ranked)
    print(f"Wrote {len(ranked)} scored/ranked compounds to {args.output_csv}")

    if args.top_n:
        threshold = ranked[min(args.top_n, len(ranked)) - 1]["Score"]
        print(f"Score at rank {args.top_n}: {threshold}")
        if args.top_n_ids_out:
            with open(args.top_n_ids_out, "w") as f:
                for r in ranked[:args.top_n]:
                    f.write(r["Source_ID"] + "\n")
            print(f"Wrote top-{args.top_n} Source_IDs to {args.top_n_ids_out}")


if __name__ == "__main__":
    main()
