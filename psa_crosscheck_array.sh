#!/bin/bash
#SBATCH --job-name=psa_crosscheck_array
#SBATCH --partition=regular
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# Array version of psa_crosscheck.py: task N cross-checks
# $CHUNK_DIR/chunk_N_analysis_top50.sdf's PyMOL-computed PSA_3D_Angstrom2
# against an independently-computed RDKit/FreeSASA value, for every
# conformer in that chunk (read-only -- doesn't touch or regenerate anything).
#
# Usage (array range must match your chunk count):
#   sbatch --array=0-22 --export=CHUNK_DIR=/vast/scratch/users/$USER/3D-Descriptors/new_compounds_chunks psa_crosscheck_array.sh

set -euo pipefail

# Slurm spools the batch script to /var/spool/slurmd/... before running it,
# so $BASH_SOURCE/$0 point at that spool copy, not this file's real location
# -- can't derive the pipeline dir from them. Use SLURM_SUBMIT_DIR (the dir
# sbatch was invoked from) instead; override with PIPELINE_DIR= at
# submission time if you're not submitting from the repo root.
PIPELINE_DIR="${PIPELINE_DIR:-$SLURM_SUBMIT_DIR}"

module load miniconda3/latest
CONDA_SH_PATH="${CONDA_SH_PATH:-/path/to/miniconda3/etc/profile.d/conda.sh}"
source "$CONDA_SH_PATH"
NVMOLKIT_CONDA_ENV="${NVMOLKIT_CONDA_ENV:-/path/to/condaenvs/nvmolkit}"
conda activate "$NVMOLKIT_CONDA_ENV"

CHUNK_DIR="${CHUNK_DIR:?Set CHUNK_DIR at submission time -- must match OUTPUT_DIR from conformer_gen_array.sh}"
INPUT_SDF="$CHUNK_DIR/chunk_${SLURM_ARRAY_TASK_ID}_analysis_top50.sdf"
OUTPUT_CSV="$CHUNK_DIR/chunk_${SLURM_ARRAY_TASK_ID}_psa_crosscheck.csv"

echo "Task $SLURM_ARRAY_TASK_ID: $INPUT_SDF -> $OUTPUT_CSV"

python "$PIPELINE_DIR/scripts/conformers/psa_crosscheck.py" \
    --input-sdf "$INPUT_SDF" \
    --output-csv "$OUTPUT_CSV"
