#!/bin/bash
#SBATCH --job-name=nvmolkit_conformers_array
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:A30:1
#SBATCH --nodelist=gpu-a30-n01,gpu-a30-n02,gpu-a30-n03
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# Chunked conformer generation for a large compound set, split across Slurm
# array tasks: each task processes CHUNK_SIZE compounds starting at
# SLURM_ARRAY_TASK_ID * CHUNK_SIZE. Chunking means a single failed/timed-out
# task doesn't lose the whole run (just resubmit that one array index), and
# tasks can run concurrently across the 3 named A100 nodes.
#
# Sizing notes carried over from prior calibration on this cluster (worth
# re-checking with a small run of your own if compound sizes differ a lot):
#   - CHUNK_SIZE=300 with 1000 conformers/molecule used ~24GB peak RAM
#     (32G ceiling; this account's QOS rejects 64G+ for GPU jobs outright --
#     "QOSMaxGRESPerUser").
#   - MMFF_BATCH_SIZE=200 is required, not just a memory optimization: 0
#     (unbatched) crashed with "CUDA error 700: illegal memory access" at
#     ~300 compounds x 1000 conformers, and 1000 crashed with "CUDA error 2:
#     out of memory" (real GPU VRAM). 200 worked cleanly.
#   - CHUNK_SIZE=100 with 3000 conformers/molecule targets the same
#     ~300k-conformer-per-task envelope as the calibrated 300x1000 case
#     (extrapolated, not yet measured directly -- worth confirming with one
#     task before trusting a big --array submission).
#
# Usage:
#   # Array range = ceil(compound_count / CHUNK_SIZE) - 1. E.g. 500 compounds
#   # at CHUNK_SIZE=300 -> ceil(500/300)=2 chunks -> --array=0-1.
#   # %3 caps concurrency at the 3 available A100 nodes.
#   sbatch --array=0-N%3 --export=INPUT_CSV=new_compounds.csv,CHUNK_SIZE=300,CONFS_PER_MOLECULE=1000,MMFF_BATCH_SIZE=200,OUTPUT_DIR=/vast/scratch/users/$USER/3D-Descriptors/new_compounds_chunks conformer_gen_array.sh
#
# After ALL tasks finish, concatenate the per-chunk SDFs -- plain text
# concatenation is valid for SDF, each molecule block is self-delimited by
# '$$$$':
#   cat $OUTPUT_DIR/conformers_chunk_*.sdf > $OUTPUT_DIR/conformers_full.sdf

set -euo pipefail

# Slurm spools the batch script to /var/spool/slurmd/... before running it,
# so $BASH_SOURCE/$0 point at that spool copy, not this file's real location
# -- can't derive the pipeline dir from them. Use SLURM_SUBMIT_DIR (the dir
# sbatch was invoked from) instead; override with PIPELINE_DIR= at
# submission time if you're not submitting from the repo root.
PIPELINE_DIR="${PIPELINE_DIR:-$SLURM_SUBMIT_DIR}"

module load miniconda3/latest
module load CUDA/12.8
CONDA_SH_PATH="${CONDA_SH_PATH:-/path/to/miniconda3/etc/profile.d/conda.sh}"
source "$CONDA_SH_PATH"
conda activate nvmolkit

INPUT_CSV="${INPUT_CSV:?Set INPUT_CSV (or adapt this script for --input-sdf) at submission time}"
SMILES_COLUMN="${SMILES_COLUMN:-smiles}"
ID_COLUMN="${ID_COLUMN:-Title}"
CONFS_PER_MOLECULE="${CONFS_PER_MOLECULE:-3000}"
CHUNK_SIZE="${CHUNK_SIZE:-100}"
MMFF_BATCH_SIZE="${MMFF_BATCH_SIZE:-200}"
SEED="${SEED:-42}"
OUTPUT_DIR="${OUTPUT_DIR:?Set OUTPUT_DIR at submission time -- pick a fresh directory, do not reuse another run chunk output}"

mkdir -p "$OUTPUT_DIR" logs

START_INDEX=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
OUTPUT_SDF="$OUTPUT_DIR/conformers_chunk_${SLURM_ARRAY_TASK_ID}.sdf"

echo "Task $SLURM_ARRAY_TASK_ID: rows $START_INDEX..$((START_INDEX + CHUNK_SIZE - 1)) -> $OUTPUT_SDF (seed=$SEED)"

python "$PIPELINE_DIR/scripts/conformers/generate_conformers.py" \
    --input-csv "$INPUT_CSV" \
    --smiles-column "$SMILES_COLUMN" \
    --id-column "$ID_COLUMN" \
    --start-index "$START_INDEX" \
    --max-molecules "$CHUNK_SIZE" \
    --confs-per-molecule "$CONFS_PER_MOLECULE" \
    --mmff-batch-size "$MMFF_BATCH_SIZE" \
    --seed "$SEED" \
    --output-sdf "$OUTPUT_SDF"
