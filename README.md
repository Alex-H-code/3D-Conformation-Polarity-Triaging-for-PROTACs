# Conformational Triage Pipeline for Heterobifunctional Degraders

Reduces a large virtual PROTAC library to a purchasable shortlist using
Boltzmann-weighted 3D conformer descriptors (radius of gyration, solvent-
accessible 3D polar surface area), cross-validated against an independent
SASA implementation, and filtered against literature-derived permeability
windows.

Full methodology and rationale are written up in `docs/methods.md` (the
Methods section text this README summarizes into runnable commands).

## Requirements

This pipeline was built and run on a Slurm-managed HPC cluster with GPU
nodes. **Paths, node names, and partition names throughout the `*_array.sh`
scripts are specific to the cluster this was developed on and need editing
before use elsewhere** (search for `/vast/scratch/`, `gpu-a30-n`, and
`--partition=` to find them).

Two separate conda environments are required — their RDKit versions
genuinely differ, which matters for exact reproduction:

| Environment | Used by | Key packages |
|---|---|---|
| `nvmolkit` | conformer generation, Boltzmann weighting, cross-validation, extraction, clustering | RDKit 2026.03.1, nvmolkit 0.5.1, Python 3.12.13 |
| `pymol_env` | 3D-PSA / intramolecular H-bond calculation | RDKit 2026.03.5, pymol-open-source 3.1.0, Python 3.12.13 |

Environment files: `environment_nvmolkit.yml`, `environment_pymol.yml`
(generate with `conda env export --no-builds > <name>.yml` from each
activated environment). CUDA 12.8 is loaded as a separate cluster module
(`module load CUDA/12.8`), outside conda.

Each `*.sh` script activates its conda environment by full path, defaulting
to a placeholder (`/path/to/condaenvs/<env>`) — override at submission time
with `--export=NVMOLKIT_CONDA_ENV=...` or `PYMOL_CONDA_ENV=...` (or edit the
default in the script) to point at your own environment locations.

**Reproducibility caveat**: conformer embedding uses a fixed random seed
(42), but GPU floating-point operations are not guaranteed bit-identical
across different hardware/driver versions — expect statistical, not exact,
reproduction of individual conformer coordinates.

## Pipeline stages

Run in this order. GPU steps run on the `nvmolkit` environment on GPU
nodes; CPU steps can run on general compute nodes.

### 1. Conformer generation (GPU)

```bash
sbatch --array=0-N%3 --export=INPUT_CSV=compounds.csv,CONFS_PER_MOLECULE=3000,CHUNK_SIZE=100,OUTPUT_DIR=$OUT conformer_gen_array.sh
```

Generates 3,000 MMFF94-minimized conformers per compound (vacuum, no
solvent model) via GPU-batched distance-geometry embedding. **Batch size
must be 200** for the MMFF minimization step — 0 (unbatched) causes an
illegal-memory-access CUDA error, and 1,000 exceeds GPU memory.

### 2. Boltzmann weighting (CPU)

```bash
sbatch --array=0-N --export=CHUNK_DIR=$OUT boltzmann_weight_array.sh
```

Weights each conformer by relative MMFF energy at 298 K (see `docs/methods.md`
for the exact formula).

### 3. Descriptor calculation (CPU, `pymol_env`)

```bash
sbatch --array=0-N --export=CHUNK_DIR=$OUT analyse_conformers_array.sh
```

Computes Rgyr (RDKit), 3D-PSA (PyMOL SASA), and per-compound Boltzmann-
weighted averages, writing `chunk_N_features.csv` (one row per compound) and
`chunk_N_analysis.sdf` (one record per conformer).

### 4. Trim to top-N conformers per compound (CPU, storage/compute efficiency)

```bash
python scripts/conformers/filter_top_n.py --input-sdf chunk_N_analysis.sdf --output-sdf chunk_N_analysis_top50.sdf --top-n 50
```

Retains only the highest-Boltzmann-weight conformers per compound in the
output SDF. Negligible information loss given how sharply peaked most
compounds' weight distributions are (see step 5).

### 5. Cross-validate 3D-PSA against an independent SASA implementation (CPU)

```bash
sbatch --array=0-N --export=CHUNK_DIR=$OUT psa_crosscheck_array.sh
```

Recomputes PyMOL's 3D-PSA with RDKit's FreeSASA wrapper (explicit Bondi
radii — `rdFreeSASA.classifyAtoms()` silently returns zero for non-protein
molecules, so this step assigns radii by element itself). Combine per-chunk
outputs into `psa_crosscheck_all.csv` before the next step.

### 6. Merge PSA definitions per compound

```bash
python scripts/triage/merge_psa_definitions.py \
    --crosscheck-csv psa_crosscheck_all.csv \
    --features-csv chunk_0_features.csv chunk_1_features.csv ... \
    --output-csv molecules_with_psa_variants.csv
```

Produces `PSA_3D_pymol`, `PSA_3D_freesasa`, and `PSA_3D_averaged` (and
matching `Shielding_Index_*`) per compound. The averaged value is what the
filtering step below uses by default, hedging against the ~7% systematic
offset between the two SASA implementations (see `docs/methods.md`).

### 7. Filter to a shortlist

```bash
python scripts/triage/filter_by_thresholds.py \
    --input-csv molecules_with_psa_variants.csv --psa-definition averaged \
    --rgyr-min 4.5 --rgyr-max 6.0 --psa-min 135 --psa-max 230 \
    --output-csv shortlist.csv
```

Default bounds follow Kim, Sheridan, Zhang, Barros, Johnston, and Xiao
(*J. Chem. Inf. Model.* 2025). Add `--hbonds-min 1` for their stricter
H-bonding-inclusive window.

### 8. (Optional) Continuous ranking instead of / alongside hard thresholds

```bash
python scripts/triage/compute_weighted_score.py \
    --input-csv molecules_with_psa_variants.csv --psa-definition averaged \
    --output-csv scored.csv --top-n 200 --top-n-ids-out top200_ids.txt
```

### 9. (Optional) Attach ADME model predictions

```bash
python scripts/triage/merge_adme_predictions.py \
    --compounds-csv shortlist.csv --id-column Source_ID \
    --predictions-csv all_predictions_merged.csv --predictions-id-column compound_id \
    --output-csv shortlist_with_adme.csv \
    --linearize "pred(MDCK-MDR1_LogER)" "pred(rLM LogCLint)" "pred(hLM LogCLint)"
```

Reported alongside the 3D-descriptor filtering as an independent metric —
deliberately not merged into the same composite score, since 3D shielding
and the model's own permeability endpoints capture overlapping signal.

### 10. (Optional) Structural diversity check

```bash
python scripts/triage/cluster_compounds.py \
    --compounds-csv shortlist.csv --id-column Source_ID \
    --predictions-csv all_predictions_merged.csv --predictions-id-column compound_id \
    --smiles-column degrader_smiles --similarity-threshold 0.85 \
    --output-csv shortlist_clustered.csv
```

**Important for combinatorially-assembled libraries**: whole-molecule
clustering needs a much higher similarity threshold (~0.85) than the
conventional default (~0.6) — a shared invariant warhead scaffold across the
library collapses nearly everything into one cluster at low thresholds. Run
again with `--smiles-column linker_smiles --similarity-threshold 0.7` to
cluster on the linker alone instead.

### 11. (Optional) Summary figure

```bash
python scripts/triage/plot_molecular_properties.py \
    --input-csv shortlist_with_adme.csv \
    --columns Rgyr_Angstrom_boltzmann TPSA_2D_Angstrom2 PSA_3D_averaged \
    --labels Rgyr "TPSA (2D)" "3D-PSA (avg)" \
    --output-svg shortlist_properties.svg
```

## Directory structure

```
new_compounds_pipeline/
├── conformer_gen_array.sh, boltzmann_weight_array.sh,
│   analyse_conformers_array.sh, psa_crosscheck_array.sh   # Slurm array wrappers
├── analyse_conformers.py                                  # descriptor calculation
├── scripts/
│   ├── conformers/       # generation, weighting, filtering, cross-validation, extraction
│   └── triage/           # downstream shortlist scoring/filtering/clustering/plotting
├── docs/methods.md        # full Methods write-up
├── environment_nvmolkit.yml, environment_pymol.yml
└── README.md              # this file
```

## Known issues already fixed here (documented so they don't recur)

- **`PIPELINE_DIR` must derive from `$SLURM_SUBMIT_DIR`, not `$BASH_SOURCE`** — Slurm spools batch scripts to `/var/spool/slurmd/...` before running them, so path-of-self tricks resolve to the spool copy, not the submitted script's real location.
- **`rdFreeSASA.classifyAtoms()` silently fails on small molecules** — it's parameterized for protein residues and returns all-zero radii rather than erroring; supply explicit per-element van der Waals radii instead.
- **MMFF batch size must be 200** for GPU conformer minimization — see step 1.

## License / Citation

<!-- Add a LICENSE file (MIT/BSD are typical for tooling like this) and a
     CITATION.cff if you want this repo citable via GitHub's "Cite this
     repository" feature. -->
