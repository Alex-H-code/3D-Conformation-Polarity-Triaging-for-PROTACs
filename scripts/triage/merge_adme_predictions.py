"""
merge_adme_predictions.py
==========================
Joins a compound shortlist onto its ADME model predictions by compound ID,
pulling every 'pred(...)' column from the predictions file. Optionally adds
linear-scale companion columns (ER = 10^LogER, CLint = 10^LogCLint, etc.)
alongside any log-scale prediction column requested, since several ADME
endpoints (MDCK-MDR1 efflux ratio, microsomal clearance) are reported in
log10 space by convention.

Usage:
    python merge_adme_predictions.py \\
        --compounds-csv filtered_compounds.csv \\
        --predictions-csv all_predictions_merged.csv \\
        --id-column Source_ID --predictions-id-column compound_id \\
        --output-csv filtered_compounds_with_adme.csv \\
        --linearize "pred(MDCK-MDR1_LogER)" "pred(rLM LogCLint)" "pred(hLM LogCLint)"
"""
import argparse
import csv


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--compounds-csv", required=True)
    ap.add_argument("--id-column", default="Source_ID", help="ID column name in --compounds-csv")
    ap.add_argument("--predictions-csv", required=True)
    ap.add_argument("--predictions-id-column", default="compound_id", help="ID column name in --predictions-csv")
    ap.add_argument("--output-csv", required=True)
    ap.add_argument("--linearize", nargs="*", default=[],
                     help="Log-scale pred(...) column names to also emit as 10^x linear-scale companions")
    return ap.parse_args()


def linear_column_name(log_col):
    # e.g. "pred(MDCK-MDR1_LogER)" -> "pred(MDCK-MDR1_ER)"; "pred(rLM LogCLint)" -> "pred(rLM CLint)"
    return log_col.replace("Log", "", 1)


def main():
    args = parse_args()

    with open(args.compounds_csv, newline="") as f:
        reader = csv.DictReader(f)
        base_fields = reader.fieldnames
        rows = list(reader)
    print(f"Loaded {len(rows)} compounds from {args.compounds_csv}")

    with open(args.predictions_csv, newline="") as f:
        reader = csv.DictReader(f)
        pred_cols = [c for c in reader.fieldnames if c.startswith("pred(")]
        by_id = {row[args.predictions_id_column]: row for row in reader}
    print(f"Loaded {len(by_id)} compounds x {len(pred_cols)} ADME predictions from {args.predictions_csv}")

    linear_cols = [linear_column_name(c) for c in args.linearize]
    out_fields = base_fields + pred_cols + linear_cols

    missing = []
    for row in rows:
        pred_row = by_id.get(row[args.id_column])
        if pred_row is None:
            missing.append(row[args.id_column])
            for c in pred_cols + linear_cols:
                row[c] = ""
            continue
        for c in pred_cols:
            row[c] = pred_row[c]
        for log_col, lin_col in zip(args.linearize, linear_cols):
            row[lin_col] = 10 ** float(pred_row[log_col])

    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} compounds to {args.output_csv}")
    if missing:
        print(f"WARNING: {len(missing)} compounds had no ADME match: {missing[:10]}"
              + (" ..." if len(missing) > 10 else ""))


if __name__ == "__main__":
    main()
