#!/bin/bash
#SBATCH --job-name=boltzmann_weight
#SBATCH --partition=regular
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Defaults sized for production-scale chunks (~300,000 conformers, multi-GB
# SDF): a 40,000-conformer chunk OOM'd at the original 4G default -- 32G
# matches what generate_conformers.py needed for the same conformer count.
#
# Usage:
#   sbatch --export=INPUT_SDF=/path/conformers.sdf,OUTPUT_SDF=/path/conformers-weighted.sdf boltzmann_weight.sh

set -euo pipefail

module load miniconda3/latest
CONDA_SH_PATH="${CONDA_SH_PATH:-/path/to/miniconda3/etc/profile.d/conda.sh}"
source "$CONDA_SH_PATH"
NVMOLKIT_CONDA_ENV="${NVMOLKIT_CONDA_ENV:-/path/to/condaenvs/nvmolkit}"
conda activate "$NVMOLKIT_CONDA_ENV"

INPUT_SDF="${INPUT_SDF:?Set INPUT_SDF at submission time, e.g. sbatch --export=INPUT_SDF=...,OUTPUT_SDF=... boltzmann_weight.sh}"
OUTPUT_SDF="${OUTPUT_SDF:?Set OUTPUT_SDF at submission time}"

python scripts/conformers/boltzmann_weight.py \
    --input-sdf "$INPUT_SDF" \
    --output-sdf "$OUTPUT_SDF"
