"""
filter_by_thresholds.py
========================
Hard-threshold triage filter on Boltzmann-weighted Rgyr, 3D-PSA, and
(optionally) intramolecular H-bonds. Default bounds (Rgyr 4.5-6.0 A,
3D-PSA 135-230 A2) follow the permeability-associated window reported by
Kim, Sheridan, Zhang, Barros, Johnston, and Xiao (J. Chem. Inf. Model.
2025) for heterobifunctional degraders.

Input is expected to be the output of merge_psa_definitions.py (or
compute_weighted_score.py, which passes the same columns through).

Usage:
    python filter_by_thresholds.py \\
        --input-csv molecules_with_psa_variants.csv \\
        --psa-definition averaged \\
        --rgyr-min 4.5 --rgyr-max 6.0 \\
        --psa-min 135 --psa-max 230 \\
        --hbonds-min 1 \\
        --output-csv filtered_compounds.csv
"""
import argparse
import csv


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--psa-definition", choices=["pymol", "freesasa", "averaged"], default="averaged")
    ap.add_argument("--rgyr-min", type=float, default=4.5)
    ap.add_argument("--rgyr-max", type=float, default=6.0)
    ap.add_argument("--psa-min", type=float, default=135.0)
    ap.add_argument("--psa-max", type=float, default=230.0)
    ap.add_argument("--hbonds-min", type=float, default=None,
                     help="Optional strictly-greater-than cutoff on Boltzmann-weighted IMHB (unset = no IMHB filter)")
    ap.add_argument("--output-csv", required=True)
    return ap.parse_args()


def main():
    args = parse_args()
    psa_col = f"PSA_3D_{args.psa_definition}"

    with open(args.input_csv, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    print(f"Loaded {len(rows)} compounds from {args.input_csv}")

    def passes(r):
        if not (args.rgyr_min <= float(r["Rgyr_Angstrom_boltzmann"]) <= args.rgyr_max):
            return False
        if not (args.psa_min <= float(r[psa_col]) <= args.psa_max):
            return False
        if args.hbonds_min is not None and not (float(r["Intramolecular_HBonds_boltzmann"]) > args.hbonds_min):
            return False
        return True

    passed = [r for r in rows if passes(r)]

    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(passed)

    print(f"Passed: {len(passed)}/{len(rows)} ({100*len(passed)/len(rows):.1f}%) -> {args.output_csv}")
    print(f"  Rgyr in [{args.rgyr_min}, {args.rgyr_max}] A")
    print(f"  3D-PSA ({args.psa_definition}) in [{args.psa_min}, {args.psa_max}] A2")
    if args.hbonds_min is not None:
        print(f"  IMHB > {args.hbonds_min}")


if __name__ == "__main__":
    main()
