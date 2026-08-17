"""
analyse_conformers.py
======================
Per-conformer 3D descriptor calculation (Rgyr, 3D-PSA, intramolecular H-bonds)
plus Boltzmann-weighted per-molecule summary. Consumes the output of
scripts/conformers/boltzmann_weight.py (needs 'Boltzmann_Weight',
'Molecule_Index', 'Conformer_Rank', 'MMFF_Energy' SDF properties).

3D-PSA follows Kihlberg 2024 (O, N, and H bonded to O/N; PyMOL solvent-accessible
surface with a 1.4 A probe). Intramolecular H-bonds use a standard geometric
criterion: donor-H...acceptor distance < 2.5 A and D-H...A angle > 120 deg,
excluding donor/acceptor pairs within 2 bonds of each other (1-2/1-3).

Also reports a shielding index (2D TPSA - Boltzmann-weighted 3D-PSA): how much
polar surface area gets hidden by intramolecular folding, relative to the fully
extended/topological value. 2D TPSA is purely topological (RDKit CalcTPSA
ignores 3D coordinates), so it's identical across a molecule's conformers.

With --top-n, only the N highest Boltzmann_Weight conformers per molecule get
the (expensive, PyMOL-based) descriptor calculation -- the rest are skipped
entirely, not just excluded from the summary. Weights themselves must already
be present (from scripts/conformers/boltzmann_weight.py); this only decides
which conformers to spend compute on here.

Usage:
    python analyse_conformers.py \\
        --input-sdf conformers_weighted.sdf \\
        --output-csv conformer_analysis.csv \\
        --output-sdf conformer_analysis.sdf \\
        --summary-csv conformer_features.csv
    python analyse_conformers.py \\
        --input-sdf conformers_weighted.sdf \\
        --output-csv conformer_analysis.csv \\
        --output-sdf conformer_analysis.sdf \\
        --summary-csv conformer_features.csv \\
        --top-n 50
"""
import argparse
import csv
import os
import tempfile
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors3D, rdMolDescriptors
import pymol2


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-sdf", required=True, help="Boltzmann-weighted conformer SDF")
    ap.add_argument("--output-csv", required=True, help="Per-conformer descriptor table")
    ap.add_argument("--output-sdf", required=True, help="SDF with descriptors written back as properties")
    ap.add_argument("--summary-csv", required=True,
                     help="Per-molecule Boltzmann-weighted summary (feature file for merge_features.py)")
    ap.add_argument("--probe-radius", type=float, default=1.4, help="Angstrom, PyMOL SASA probe")
    ap.add_argument("--dot-density", type=int, default=3, help="PyMOL surface sampling density")
    ap.add_argument("--hbond-dist-cutoff", type=float, default=2.5,
                     help="Angstrom, max H...acceptor distance")
    ap.add_argument("--hbond-angle-cutoff", type=float, default=120.0,
                     help="Degrees, min D-H...A angle")
    ap.add_argument("--hbond-min-topological-sep", type=int, default=3,
                     help="Min donor-acceptor bond-path separation (excludes 1-2/1-3 pairs)")
    ap.add_argument("--top-n", type=int, default=None,
                     help="Only run descriptor calculation on the N highest "
                          "Boltzmann_Weight conformers per molecule (default: all)")
    return ap.parse_args()


# -- Boltzmann weighted average ------------------------------------------------
def boltzmann_average(values, weights):
    values = np.array(values)
    weights = np.array(weights)
    return np.sum(weights * values) / np.sum(weights)


# -- Intramolecular H-bond counting --------------------------------------------
def count_intramolecular_hbonds(mol, dist_cutoff, angle_cutoff, min_topological_sep):
    """Geometric H-bond count for a single-conformer mol (conf id 0)."""
    donor_h_pairs = []  # (H idx, donor heavy-atom idx)
    acceptor_idxs = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ("N", "O"):
            acceptor_idxs.append(atom.GetIdx())
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == "H":
                    donor_h_pairs.append((nbr.GetIdx(), atom.GetIdx()))

    if not donor_h_pairs or not acceptor_idxs:
        return 0

    dist_matrix = Chem.GetDistanceMatrix(mol)  # topological (bond-count) distances
    conf = mol.GetConformer()

    count = 0
    for h_idx, d_idx in donor_h_pairs:
        h_pos = np.array(conf.GetAtomPosition(h_idx))
        d_pos = np.array(conf.GetAtomPosition(d_idx))
        for a_idx in acceptor_idxs:
            if a_idx == d_idx or dist_matrix[d_idx][a_idx] < min_topological_sep:
                continue
            a_pos = np.array(conf.GetAtomPosition(a_idx))
            dist_ha = np.linalg.norm(h_pos - a_pos)
            if dist_ha > dist_cutoff:
                continue
            v1, v2 = d_pos - h_pos, a_pos - h_pos
            cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
            if angle >= angle_cutoff:
                count += 1
    return count


def main():
    args = parse_args()

    if not os.path.exists(args.input_sdf):
        raise FileNotFoundError(f"Input SDF not found: {args.input_sdf}")

    # Streamed, not `mols = [m for m in supplier ...]` -- holding every conformer
    # in memory at once (then a second pass over zip(mols, results) to write the
    # output SDF) is what OOM'd on production-scale chunks (18GB for just 40,000
    # conformers). len(supplier) counts records without parsing/holding them.
    supplier = Chem.SDMolSupplier(args.input_sdf, removeHs=False, sanitize=True)
    n_records = len(supplier)
    print(f"Loaded {n_records} conformer records from {args.input_sdf}")

    # -- Optional top-N filter: rank conformers by their already-computed
    # Boltzmann_Weight and skip the expensive PyMOL/RDKit descriptor calc for
    # everything outside the top N per molecule. Only reads Molecule_Index /
    # Boltzmann_Weight off each record (cheap relative to the SASA calc below),
    # not held as full parsed mols, so this doesn't reintroduce the OOM the
    # streamed design above was written to avoid.
    keep_idxs = None
    if args.top_n is not None:
        groups = {}
        for i, mol in enumerate(Chem.SDMolSupplier(args.input_sdf, removeHs=True, sanitize=False)):
            if mol is None or not mol.HasProp("Molecule_Index") or not mol.HasProp("Boltzmann_Weight"):
                continue
            groups.setdefault(mol.GetIntProp("Molecule_Index"), []).append(
                (i, mol.GetDoubleProp("Boltzmann_Weight")))
        keep_idxs = set()
        for mol_idx, pairs in groups.items():
            ranked = sorted(pairs, key=lambda p: p[1], reverse=True)
            keep_idxs.update(i for i, _ in ranked[:args.top_n])
        print(f"Top-N filter: keeping {len(keep_idxs)}/{n_records} conformers "
              f"(top {args.top_n}/molecule across {len(groups)} molecules)")

    # -- Calculate descriptors, writing each conformer's output immediately ---------
    results = []
    start_time = time.time()

    os.makedirs(os.path.dirname(args.output_sdf) or ".", exist_ok=True)
    sdf_writer = Chem.SDWriter(args.output_sdf)

    with pymol2.PyMOL() as pymol:
        pymol.cmd.set("dot_solvent", 1)
        pymol.cmd.set("dot_density", args.dot_density)
        pymol.cmd.set("solvent_radius", args.probe_radius)

        for i, mol in enumerate(supplier):
            if mol is None:
                continue
            if keep_idxs is not None and i not in keep_idxs:
                continue
            try:
                # Rgyr directly from RDKit
                rgyr = Descriptors3D.RadiusOfGyration(mol)

                # 2D TPSA -- topological, identical for every conformer of this
                # molecule; used downstream to compute the shielding index.
                tpsa_2d = rdMolDescriptors.CalcTPSA(mol)

                # Intramolecular H-bonds (geometric criterion, RDKit only)
                n_hbonds = count_intramolecular_hbonds(
                    mol, args.hbond_dist_cutoff, args.hbond_angle_cutoff,
                    args.hbond_min_topological_sep,
                )

                # Write single conformer to temp SDF for PyMOL PSA. mkstemp
                # (not f"/tmp/conf_{i}.sdf") because `i` is only unique within
                # this task's own chunk -- concurrent array tasks on the same
                # node would otherwise race on the same path and load each
                # other's molecule mid-write.
                fd, tmp_sdf = tempfile.mkstemp(suffix=".sdf", prefix=f"conf_{i}_")
                os.close(fd)
                obj = f"conf_{i}"
                try:
                    writer = Chem.SDWriter(tmp_sdf)
                    writer.write(mol)
                    writer.close()

                    pymol.cmd.load(tmp_sdf, obj)

                    # 3D PSA -- O, N and H bonded to O/N, matching Kihlberg
                    pymol.cmd.select("polar_atoms",
                        f"{obj} and (elem O or elem N or "
                        f"(elem H and (neighbor elem O or neighbor elem N)))")
                    psa_3d = pymol.cmd.get_area("polar_atoms", load_b=0)
                finally:
                    # Always clean up, even on failure -- otherwise a PyMOL
                    # error here leaks the loaded object and temp file for
                    # every such conformer, growing session memory and /tmp
                    # usage over a 300k-conformer chunk.
                    try:
                        pymol.cmd.delete(obj)
                    except Exception:
                        pass
                    if os.path.exists(tmp_sdf):
                        try:
                            os.remove(tmp_sdf)
                        except OSError:
                            pass

                energy = mol.GetDoubleProp("MMFF_Energy") if mol.HasProp("MMFF_Energy") else None
                boltzmann = mol.GetDoubleProp("Boltzmann_Weight") if mol.HasProp("Boltzmann_Weight") else None
                conf_rank = mol.GetIntProp("Conformer_Rank") if mol.HasProp("Conformer_Rank") else None
                mol_idx = mol.GetIntProp("Molecule_Index") if mol.HasProp("Molecule_Index") else None
                source_id = mol.GetProp("Source_ID") if mol.HasProp("Source_ID") else str(mol_idx)

                result = {
                    "Conformer": i + 1,
                    "Source_ID": source_id,
                    "Molecule_Index": mol_idx,
                    "Conformer_Rank": conf_rank,
                    "MMFF_Energy_kcal_mol": round(energy, 3) if energy is not None else None,
                    "Boltzmann_Weight": round(boltzmann, 6) if boltzmann is not None else None,
                    "Rgyr_Angstrom": round(rgyr, 3),
                    "PSA_3D_Angstrom2": round(psa_3d, 3),
                    "TPSA_2D_Angstrom2": round(tpsa_2d, 3),
                    "Intramolecular_HBonds": n_hbonds,
                }
                results.append(result)

                # Write this conformer's descriptors back and save it immediately --
                # streamed, so `mol` can be freed rather than held for a second pass.
                mol.SetIntProp("Conformer", result["Conformer"])
                mol.SetIntProp("Molecule_Index", result["Molecule_Index"] or 0)
                mol.SetIntProp("Conformer_Rank", result["Conformer_Rank"] or 0)
                mol.SetDoubleProp("MMFF_Energy_kcal_mol", result["MMFF_Energy_kcal_mol"] or 0.0)
                mol.SetDoubleProp("Boltzmann_Weight", result["Boltzmann_Weight"] or 0.0)
                mol.SetDoubleProp("Rgyr_Angstrom", result["Rgyr_Angstrom"])
                mol.SetDoubleProp("PSA_3D_Angstrom2", result["PSA_3D_Angstrom2"])
                mol.SetDoubleProp("TPSA_2D_Angstrom2", result["TPSA_2D_Angstrom2"])
                mol.SetIntProp("Intramolecular_HBonds", result["Intramolecular_HBonds"])
                sdf_writer.write(mol)

                if (i + 1) % 1000 == 0 or (i + 1) == n_records:
                    elapsed = time.time() - start_time
                    print(f"  Progress: {i+1}/{n_records} conformers "
                          f"({elapsed:.1f}s elapsed, {(i+1)/elapsed:.1f} conformers/sec) -- "
                          f"last: Mol {mol_idx}, Rank {conf_rank}, Rgyr={rgyr:.3f} A, "
                          f"PSA={psa_3d:.1f} A2, HBonds={n_hbonds}")

            except Exception as e:
                print(f"  Conformer {i+1}: failed - {e}")
                continue

    sdf_writer.close()
    total_time = time.time() - start_time
    print(f"\nDescriptor calculation completed in {total_time:.2f} seconds "
          f"({len(results)}/{n_records} conformers, {len(results) / total_time:.1f} conformers/second)")
    print(f"Saved SDF to {args.output_sdf}")

    # -- Write per-conformer CSV ----------------------------------------------------
    if results:
        os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
        with open(args.output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"Saved per-conformer analysis to {args.output_csv}")

    # -- Boltzmann-weighted summary per molecule, written to CSV --------------------
    # Grouped by Source_ID (the original SDF title/ID) rather than Molecule_Index
    # (a positional array index) so merge_features.py can join this back onto the
    # training/prospective CSV by a stable identifier.
    summary_rows = []
    source_ids = sorted(set(r["Source_ID"] for r in results if r["Source_ID"] is not None))

    print("\n-- Boltzmann-weighted summary per molecule --")
    for sid in source_ids:
        mol_results = [r for r in results if r["Source_ID"] == sid]
        weights = np.array([r["Boltzmann_Weight"] for r in mol_results])
        rgyrs = np.array([r["Rgyr_Angstrom"] for r in mol_results])
        psas = np.array([r["PSA_3D_Angstrom2"] for r in mol_results])
        hbonds = np.array([r["Intramolecular_HBonds"] for r in mol_results])

        weighted_rgyr = boltzmann_average(rgyrs, weights)
        weighted_psa = boltzmann_average(psas, weights)
        weighted_hbonds = boltzmann_average(hbonds, weights)
        # 2D TPSA is identical across conformers of the same molecule -- no
        # weighting needed, just read it off any conformer.
        tpsa_2d = mol_results[0]["TPSA_2D_Angstrom2"]
        shielding_index = tpsa_2d - weighted_psa

        # Effective sample size (Kish's ESS): how many conformers are actually
        # contributing to the weighted average, vs. it being dominated by one.
        # weights sums to 1, so ESS = 1 / sum(w_i^2): ESS=1 means one conformer
        # has ~all the weight; ESS=n_conformers means the weight is uniform.
        ess = 1.0 / np.sum(weights ** 2)
        max_weight = float(weights.max())

        summary_rows.append({
            "Source_ID": sid,
            "Molecule_Index": mol_results[0]["Molecule_Index"],
            "n_conformers": len(mol_results),
            "Rgyr_Angstrom_boltzmann": round(weighted_rgyr, 3),
            "PSA_3D_Angstrom2_boltzmann": round(weighted_psa, 3),
            "Intramolecular_HBonds_boltzmann": round(weighted_hbonds, 3),
            "TPSA_2D_Angstrom2": round(tpsa_2d, 3),
            "Shielding_Index_Angstrom2": round(shielding_index, 3),
            "Effective_Sample_Size": round(ess, 2),
            "Max_Conformer_Weight": round(max_weight, 4),
        })

        print(f"\nMolecule {sid}:")
        print(f"  Conformers:                    {len(mol_results)}")
        print(f"  Boltzmann-weighted Rgyr:       {weighted_rgyr:.3f} A")
        print(f"  Boltzmann-weighted PSA:        {weighted_psa:.3f} A2")
        print(f"  Boltzmann-weighted H-bonds:    {weighted_hbonds:.3f}")
        print(f"  2D TPSA:                       {tpsa_2d:.3f} A2")
        print(f"  Shielding index (2D-3D PSA):   {shielding_index:.3f} A2")
        print(f"  Effective sample size:         {ess:.2f} / {len(mol_results)} conformers")
        print(f"  Max single-conformer weight:   {max_weight:.4f}"
              + ("  *** DOMINATED BY ONE CONFORMER ***" if max_weight > 0.9 else ""))
        print(f"  Rgyr range:                    {rgyrs.min():.2f} - {rgyrs.max():.2f} A")
        print(f"  PSA range:                     {psas.min():.1f} - {psas.max():.1f} A2")

    os.makedirs(os.path.dirname(args.summary_csv) or ".", exist_ok=True)
    with open(args.summary_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"\nSaved per-molecule feature summary to {args.summary_csv}")


if __name__ == "__main__":
    main()
