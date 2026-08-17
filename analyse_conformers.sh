#!/bin/bash
#SBATCH --job-name=analyse_conformers
#SBATCH --partition=regular
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Defaults sized for production-scale chunks (~300,000 conformers): the
# original 4G/15min smoke-test defaults were never validated past 200
# conformers. --time is a wide-margin placeholder until the chunk-10
# calibration run (40,000 conformers) gives a real per-conformer rate for
# the PyMOL step, which is untested at this scale.
#
# Usage:
#   sbatch --export=INPUT_SDF=/path/conformers-weighted.sdf,\
#     OUTPUT_CSV=/path/conformer_analysis.csv,OUTPUT_SDF=/path/conformer_analysis.sdf,\
#     SUMMARY_CSV=/path/conformer_features.csv \
#     analyse_conformers.sh

set -euo pipefail

module load miniconda3/latest
CONDA_SH_PATH="${CONDA_SH_PATH:-/path/to/miniconda3/etc/profile.d/conda.sh}"
source "$CONDA_SH_PATH"
PYMOL_CONDA_ENV="${PYMOL_CONDA_ENV:-/path/to/condaenvs/pymol_env}"
conda activate "$PYMOL_CONDA_ENV"

INPUT_SDF="${INPUT_SDF:?Set INPUT_SDF at submission time}"
OUTPUT_CSV="${OUTPUT_CSV:?Set OUTPUT_CSV at submission time}"
OUTPUT_SDF="${OUTPUT_SDF:?Set OUTPUT_SDF at submission time}"
SUMMARY_CSV="${SUMMARY_CSV:?Set SUMMARY_CSV at submission time}"

python analyse_conformers.py \
    --input-sdf "$INPUT_SDF" \
    --output-csv "$OUTPUT_CSV" \
    --output-sdf "$OUTPUT_SDF" \
    --summary-csv "$SUMMARY_CSV"
