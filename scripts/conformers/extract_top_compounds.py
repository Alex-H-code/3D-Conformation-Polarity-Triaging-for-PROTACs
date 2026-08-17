"""
extract_top_compounds.py
=========================
Pulls every conformer belonging to a given list of compound IDs out of the
per-chunk analysis SDFs, and writes them out as one combined SDF (full 3D
structures + all descriptor properties) and one combined CSV (flattened,
one row per conformer) -- no recomputation, just filtering already-produced
output down to a compound shortlist (e.g. a top-N triage list).

Matches on the 'Source_ID' SDF property written by generate_conformers.py /
carried through boltzmann_weight.py and analyse_conformers.py unchanged.

With --top-n, only the N highest-Boltzmann_Weight conformers per compound are
kept (a cheap first pass ranks candidates by weight before the real pass
writes them out -- same two-pass pattern as filter_top_n.py).

Usage:
    python extract_top_compounds.py \\
        --ids-file top200_ids.txt \\
        --chunk-dir /path/to/3D-Descriptors/new_compounds_chunks \\
        --chunk-pattern "chunk_{i}_analysis_top50.sdf" \\
        --num-chunks 23 \\
        --output-sdf top200_conformers.sdf \\
        --output-csv top200_conformers.csv \\
        --top-n 5
"""
import argparse
import csv
import os

from rdkit import Chem


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids-file", required=True, help="Newline-delimited Source_ID list")
    ap.add_argument("--chunk-dir", required=True)
    ap.add_argument("--chunk-pattern", default="chunk_{i}_analysis_top50.sdf",
                     help="Filename pattern with {i} as the chunk index placeholder")
    ap.add_argument("--num-chunks", type=int, required=True)
    ap.add_argument("--output-sdf", required=True)
    ap.add_argument("--output-csv", required=True)
    ap.add_argument("--top-n", type=int, default=None,
                     help="Keep only the N highest Boltzmann_Weight conformers per compound (default: keep all matches)")
    return ap.parse_args()


PROPS = [
    "Source_ID", "Molecule_Index", "Conformer_Rank", "Conformer",
    "MMFF_Energy_kcal_mol", "Boltzmann_Weight",
    "Rgyr_Angstrom", "PSA_3D_Angstrom2", "TPSA_2D_Angstrom2", "Intramolecular_HBonds",
]


def main():
    args = parse_args()

    with open(args.ids_file) as f:
        wanted = set(line.strip() for line in f if line.strip())
    print(f"Loaded {len(wanted)} target compound IDs from {args.ids_file}")

    chunk_paths = {}
    for i in range(args.num_chunks):
        chunk_path = os.path.join(args.chunk_dir, args.chunk_pattern.format(i=i))
        if os.path.exists(chunk_path):
            chunk_paths[i] = chunk_path
        else:
            print(f"  [chunk {i}] missing, skipping: {chunk_path}")

    # -- optional top-N ranking pass: cheap (properties only, no sanitize),
    # groups candidate conformers by Source_ID and ranks by Boltzmann_Weight
    # before the real pass decides what to write.
    keep_keys = None
    if args.top_n is not None:
        groups = {}
        for i, chunk_path in chunk_paths.items():
            for j, mol in enumerate(Chem.SDMolSupplier(chunk_path, removeHs=True, sanitize=False)):
                if mol is None or not mol.HasProp("Source_ID") or not mol.HasProp("Boltzmann_Weight"):
                    continue
                sid = mol.GetProp("Source_ID")
                if sid not in wanted:
                    continue
                groups.setdefault(sid, []).append((i, j, mol.GetDoubleProp("Boltzmann_Weight")))
        keep_keys = set()
        for sid, entries in groups.items():
            ranked = sorted(entries, key=lambda e: e[2], reverse=True)
            keep_keys.update((i, j) for i, j, _ in ranked[:args.top_n])
        print(f"Top-N filter: keeping up to {args.top_n} conformers/compound "
              f"-> {len(keep_keys)} candidates across {len(groups)} compounds found")

    os.makedirs(os.path.dirname(args.output_sdf) or ".", exist_ok=True)
    sdf_writer = Chem.SDWriter(args.output_sdf)
    csv_rows = []
    found_ids = set()
    n_records = 0

    for i, chunk_path in chunk_paths.items():
        n_chunk_hits = 0
        for j, mol in enumerate(Chem.SDMolSupplier(chunk_path, removeHs=False, sanitize=True)):
            if mol is None or not mol.HasProp("Source_ID"):
                continue
            sid = mol.GetProp("Source_ID")
            if sid not in wanted:
                continue
            if keep_keys is not None and (i, j) not in keep_keys:
                continue
            sdf_writer.write(mol)
            n_chunk_hits += 1
            n_records += 1
            found_ids.add(sid)
            row = {}
            for prop in PROPS:
                if mol.HasProp(prop):
                    val = mol.GetProp(prop)
                    row[prop] = val
                else:
                    row[prop] = ""
            csv_rows.append(row)
        print(f"  [chunk {i}] {n_chunk_hits} matching conformers")

    sdf_writer.close()

    if csv_rows:
        with open(args.output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=PROPS)
            writer.writeheader()
            writer.writerows(csv_rows)

    missing = wanted - found_ids
    print(f"\nWrote {n_records} conformers ({len(found_ids)}/{len(wanted)} compounds) to:")
    print(f"  {args.output_sdf}")
    print(f"  {args.output_csv}")
    if missing:
        print(f"\nWARNING: {len(missing)} requested compound IDs were not found in any chunk:")
        for sid in sorted(missing):
            print(f"  {sid}")


if __name__ == "__main__":
    main()
