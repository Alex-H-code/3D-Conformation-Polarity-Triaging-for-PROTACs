"""
generate_conformers.py
=======================
GPU-accelerated conformer generation (ETKDGv3 embedding + MMFF optimization)
via NVMolKit.

Accepts either an SDF of pre-built structures (--input-sdf) or a CSV of
SMILES + an ID column (--input-csv, e.g. with 'smiles' and 'Title' columns).

Usage:
    python generate_conformers.py \\
        --input-sdf new_compounds.sdf \\
        --output-sdf /vast/scratch/users/$USER/3D-Descriptors/new_compounds-conformers.sdf \\
        --confs-per-molecule 1000

    python generate_conformers.py \\
        --input-csv new_compounds.csv --smiles-column smiles --id-column Title \\
        --confs-per-molecule 1000 \\
        --output-sdf /vast/scratch/users/$USER/3D-Descriptors/new_compounds-conformers.sdf
"""
import argparse
import csv
import os
import time

from rdkit import Chem
from rdkit.Chem.rdDistGeom import ETKDGv3
from rdkit.Chem import SDWriter
from nvmolkit.embedMolecules import EmbedMolecules as nvMolKitEmbed
from nvmolkit.types import HardwareOptions
from nvmolkit.mmffOptimization import MMFFOptimizeMoleculesConfs as nvMolKitMMFFOptimize


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    input_group = ap.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input-sdf", help="Input SDF of pre-built 2D/3D structures")
    input_group.add_argument("--input-csv", help="Input CSV of SMILES + an ID column")
    ap.add_argument("--smiles-column", default="smiles", help="Used with --input-csv")
    ap.add_argument("--id-column", default="Title", help="Used with --input-csv")
    ap.add_argument("--output-sdf", required=True, help="Output SDF path for generated conformers")
    ap.add_argument("--start-index", type=int, default=0,
                     help="Skip this many molecules before reading --max-molecules -- for chunked/array runs "
                          "over a large input, so each chunk's Molecule_Index stays globally unique and chunk "
                          "outputs can be safely concatenated later")
    ap.add_argument("--max-molecules", type=int, default=None,
                     help="Cap on number of molecules to process, starting from --start-index (default: all)")
    ap.add_argument("--confs-per-molecule", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--mmff-max-iters", type=int, default=200)
    ap.add_argument("--embed-preprocessing-threads", type=int, default=2)
    ap.add_argument("--embed-batch-size", type=int, default=25)
    ap.add_argument("--embed-batches-per-gpu", type=int, default=2)
    ap.add_argument("--mmff-preprocessing-threads", type=int, default=4)
    ap.add_argument("--mmff-batch-size", type=int, default=0,
                     help="0 = unbatched (NVMolKit processes all conformers in one pass). At large conformer "
                          "counts this has been observed to crash with 'CUDA error 700: illegal memory access' "
                          "-- try a bounded value (e.g. 200) if that happens")
    return ap.parse_args()


def load_from_sdf(path, start_index, max_molecules):
    supplier = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
    end_index = None if max_molecules is None else start_index + max_molecules
    molecules = []
    for i, mol in enumerate(supplier):
        if i < start_index:
            continue
        if end_index is not None and i >= end_index:
            break
        if mol is None:
            continue
        # Carry the SDF title-line ID through the pipeline (Molecule_Index below
        # is only a positional array index -- merge_features.py needs a stable
        # ID to join conformer-derived features back onto a training/prospective CSV).
        mol.SetProp("Source_ID", mol.GetProp("_Name") if mol.HasProp("_Name") else str(i))
        mol.RemoveAllConformers()
        molecules.append(mol)
    return molecules


def load_from_csv(path, smiles_column, id_column, start_index, max_molecules):
    end_index = None if max_molecules is None else start_index + max_molecules
    molecules = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i < start_index:
                continue
            if end_index is not None and i >= end_index:
                break
            smi = row[smiles_column]
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                print(f"WARNING: could not parse SMILES for {row[id_column]!r}: {smi!r} -- skipping")
                continue
            mol = Chem.AddHs(mol)
            mol.SetProp("Source_ID", str(row[id_column]))
            molecules.append(mol)
    return molecules


def main():
    args = parse_args()

    params = ETKDGv3()
    params.randomSeed = args.seed
    params.useRandomCoords = True

    # -- Step 1: Load molecules -------------------------------------------------
    if args.input_sdf:
        if not os.path.exists(args.input_sdf):
            raise FileNotFoundError(f"SDF file not found: {args.input_sdf}")
        molecules = load_from_sdf(args.input_sdf, args.start_index, args.max_molecules)
        source_desc = args.input_sdf
    else:
        if not os.path.exists(args.input_csv):
            raise FileNotFoundError(f"CSV file not found: {args.input_csv}")
        molecules = load_from_csv(args.input_csv, args.smiles_column, args.id_column,
                                   args.start_index, args.max_molecules)
        source_desc = args.input_csv

    print(f"Successfully loaded {len(molecules)} molecules from {source_desc} "
          f"(rows {args.start_index}-{args.start_index + len(molecules) - 1})")

    # -- Step 2: GPU-accelerated conformer generation ----------------------------
    hardware_opts = HardwareOptions(
        preprocessingThreads=args.embed_preprocessing_threads,
        batchSize=args.embed_batch_size,
        batchesPerGpu=args.embed_batches_per_gpu,
    )

    start_time = time.time()
    nvMolKitEmbed(
        molecules=molecules,
        params=params,
        confsPerMolecule=args.confs_per_molecule,
        maxIterations=-1,
        hardwareOptions=hardware_opts,
    )
    embedding_time = time.time() - start_time
    total_conformers = sum(mol.GetNumConformers() for mol in molecules)

    print(f"Conformer generation completed in {embedding_time:.2f} seconds")
    print(f"Generated {total_conformers} total conformers")
    print(f"Rate: {total_conformers / embedding_time:.1f} conformers/second")

    # -- Step 3: MMFF optimization ------------------------------------------------
    mmff_hardware_opts = HardwareOptions(
        preprocessingThreads=args.mmff_preprocessing_threads,
        batchSize=args.mmff_batch_size,
    )

    start_time = time.time()
    energies = nvMolKitMMFFOptimize(
        molecules=molecules,
        maxIters=args.mmff_max_iters,
        nonBondedThreshold=100.0,
        hardwareOptions=mmff_hardware_opts,
    )
    optimization_time = time.time() - start_time

    print(f"MMFF optimization completed in {optimization_time:.2f} seconds")
    print(f"Rate: {total_conformers / optimization_time:.1f} conformers/second")

    # -- Step 4: Save output with energies + stable per-molecule/conformer IDs ---
    os.makedirs(os.path.dirname(args.output_sdf) or ".", exist_ok=True)
    writer = SDWriter(args.output_sdf)

    for mol_idx, mol in enumerate(molecules):
        for conf_id in range(mol.GetNumConformers()):
            if energies is not None and mol_idx < len(energies):
                energy = energies[mol_idx][conf_id]
                mol.SetDoubleProp("MMFF_Energy", float(energy))
            # Global index (not just this run's local mol_idx) so chunked/array
            # runs over a large input produce non-colliding Molecule_Index values
            # when their output SDFs are later concatenated.
            mol.SetIntProp("Molecule_Index", args.start_index + mol_idx)
            mol.SetIntProp("Conformer_Rank", conf_id)
            writer.write(mol, confId=conf_id)

    writer.close()
    print(f"Saved {total_conformers} conformers to {args.output_sdf}")


if __name__ == "__main__":
    main()
