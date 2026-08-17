"""
cluster_compounds.py
=====================
Clusters a compound shortlist by structural similarity using Morgan
(ECFP4-equivalent) fingerprints and Butina clustering on Tanimoto distance.

Can cluster on either the full degrader SMILES or an isolated fragment
(e.g. linker_smiles) pulled from the same predictions file used by
merge_adme_predictions.py. For combinatorially-assembled libraries sharing
large invariant scaffolds (e.g. a fixed E3/POI warhead pair), whole-molecule
clustering needs a much higher similarity threshold than usual (~0.85) to
avoid collapsing nearly everything into one cluster -- see README.

Usage:
    python cluster_compounds.py \\
        --compounds-csv filtered_compounds.csv --id-column Source_ID \\
        --predictions-csv all_predictions_merged.csv --predictions-id-column compound_id \\
        --smiles-column degrader_smiles --similarity-threshold 0.85 \\
        --output-csv filtered_compounds_clustered.csv
"""
import argparse
import csv

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, DataStructs
from rdkit.ML.Cluster import Butina


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--compounds-csv", required=True)
    ap.add_argument("--id-column", default="Source_ID")
    ap.add_argument("--predictions-csv", required=True, help="File providing the SMILES column")
    ap.add_argument("--predictions-id-column", default="compound_id")
    ap.add_argument("--smiles-column", default="degrader_smiles",
                     help="e.g. degrader_smiles (whole molecule) or linker_smiles (linker only)")
    ap.add_argument("--similarity-threshold", type=float, default=0.85,
                     help="Tanimoto similarity above which two compounds join a cluster")
    ap.add_argument("--radius", type=int, default=2)
    ap.add_argument("--n-bits", type=int, default=2048)
    ap.add_argument("--output-csv", required=True)
    return ap.parse_args()


def main():
    args = parse_args()

    with open(args.compounds_csv, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    print(f"Loaded {len(rows)} compounds from {args.compounds_csv}")

    with open(args.predictions_csv, newline="") as f:
        reader = csv.DictReader(f)
        smiles_by_id = {row[args.predictions_id_column]: row[args.smiles_column] for row in reader}

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=args.radius, fpSize=args.n_bits)
    mols_ok = []
    missing = []
    for row in rows:
        smi = smiles_by_id.get(row[args.id_column])
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            missing.append(row[args.id_column])
            continue
        mols_ok.append((row, gen.GetFingerprint(mol)))
    print(f"Fingerprinted {len(mols_ok)}/{len(rows)} compounds (missing/unparseable: {len(missing)})")

    fps = [fp for _, fp in mols_ok]
    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])

    clusters = Butina.ClusterData(dists, n, 1 - args.similarity_threshold, isDistData=True)
    sizes = sorted((len(c) for c in clusters), reverse=True)
    print(f"{len(clusters)} clusters at similarity > {args.similarity_threshold} "
          f"(largest={sizes[0] if sizes else 0}, singletons={sum(1 for s in sizes if s == 1)})")

    for cid, cluster in enumerate(clusters, start=1):
        for idx in cluster:
            row, _ = mols_ok[idx]
            row["cluster_id"] = cid
            row["cluster_size"] = len(cluster)
    for sid in missing:
        for row in rows:
            if row[args.id_column] == sid:
                row["cluster_id"] = ""
                row["cluster_size"] = ""

    out_fields = fieldnames + ["cluster_id", "cluster_size"]
    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (-(r["cluster_size"] or 0), r.get("cluster_id", ""))))

    print(f"Wrote {args.output_csv}")


if __name__ == "__main__":
    main()
