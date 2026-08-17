#!/bin/bash
#SBATCH --job-name=analyse_conformers_array
#SBATCH --partition=regular
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# Array version of analyse_conformers.sh: task N processes
# $CHUNK_DIR/conformers_chunk_N_weighted.sdf (from boltzmann_weight_array.sh).
# --time=6:00:00 comfortably covers ~300,000 conformers/chunk at the ~40
# conformers/sec PyMOL throughput measured previously on this cluster --
# rerun a small chunk first to check if your compounds/chunk size differ
# a lot from that.
#
# Usage (run AFTER boltzmann_weight_array.sh finishes for ALL chunks --
# --dependency=afterok waits for the whole array job, not just one task).
# TOP_N runs the descriptor calculation on only the N highest Boltzmann_Weight
# conformers per molecule (unset/empty processes all conformers, matching the
# original behavior):
#   sbatch --array=0-N --export=CHUNK_DIR=/vast/scratch/users/$USER/3D-Descriptors/new_compounds_chunks,TOP_N=50 \
#     --dependency=afterok:<boltzmann_array_job_id> analyse_conformers_array.sh

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
PYMOL_CONDA_ENV="${PYMOL_CONDA_ENV:-/path/to/condaenvs/pymol_env}"
conda activate "$PYMOL_CONDA_ENV"

CHUNK_DIR="${CHUNK_DIR:?Set CHUNK_DIR at submission time -- must match OUTPUT_DIR from conformer_gen_array.sh}"
TOP_N="${TOP_N:-}"
INPUT_SDF="$CHUNK_DIR/conformers_chunk_${SLURM_ARRAY_TASK_ID}_weighted.sdf"
OUTPUT_CSV="$CHUNK_DIR/chunk_${SLURM_ARRAY_TASK_ID}_analysis.csv"
OUTPUT_SDF="$CHUNK_DIR/chunk_${SLURM_ARRAY_TASK_ID}_analysis.sdf"
SUMMARY_CSV="$CHUNK_DIR/chunk_${SLURM_ARRAY_TASK_ID}_features.csv"

echo "Task $SLURM_ARRAY_TASK_ID: $INPUT_SDF -> $OUTPUT_CSV / $SUMMARY_CSV (top_n=${TOP_N:-all})"

python "$PIPELINE_DIR/analyse_conformers.py" \
    --input-sdf "$INPUT_SDF" \
    --output-csv "$OUTPUT_CSV" \
    --output-sdf "$OUTPUT_SDF" \
    --summary-csv "$SUMMARY_CSV" \
    ${TOP_N:+--top-n "$TOP_N"}
