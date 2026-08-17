#!/bin/bash
#SBATCH --job-name=boltzmann_weight_array
#SBATCH --partition=regular
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# Array version of boltzmann_weight.sh: task N processes
# $CHUNK_DIR/conformers_chunk_N.sdf (from conformer_gen_array.sh).
#
# Usage (array range must match however many chunks conformer_gen_array.sh
# produced for this compound set, and CHUNK_DIR must match its OUTPUT_DIR).
# TOP_N keeps only the N highest-weighted conformers per molecule in the
# output (unset/empty keeps all -- weights are still computed from every
# conformer either way, TOP_N only affects what gets written out):
#   sbatch --array=0-N --export=CHUNK_DIR=/vast/scratch/users/$USER/3D-Descriptors/new_compounds_chunks,TOP_N=50 boltzmann_weight_array.sh

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
TOP_N="${TOP_N:-}"
INPUT_SDF="$CHUNK_DIR/conformers_chunk_${SLURM_ARRAY_TASK_ID}.sdf"
OUTPUT_SDF="$CHUNK_DIR/conformers_chunk_${SLURM_ARRAY_TASK_ID}_weighted.sdf"

echo "Task $SLURM_ARRAY_TASK_ID: $INPUT_SDF -> $OUTPUT_SDF (top_n=${TOP_N:-all})"

python "$PIPELINE_DIR/scripts/conformers/boltzmann_weight.py" \
    --input-sdf "$INPUT_SDF" \
    --output-sdf "$OUTPUT_SDF" \
    ${TOP_N:+--top-n "$TOP_N"}
