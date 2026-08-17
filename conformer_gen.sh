#!/bin/bash
#SBATCH --job-name=nvmolkit_conformers
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:A100:1
#SBATCH --nodelist=gpu-a100-n01,gpu-a100-n02,gpu-a100-n03
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#
# Single-shot conformer generation for a small/moderate compound set. For a
# large set, use conformer_gen_array.sh instead (chunked, resilient to
# individual task failures/timeouts, runs concurrently across GPU nodes).
#
# Set exactly one of INPUT_SDF / INPUT_CSV. Examples:
#
#   # SDF input:
#   sbatch --export=INPUT_SDF=new_compounds.sdf,\
#     OUTPUT_SDF=/vast/scratch/users/$USER/3D-Descriptors/new_compounds-conformers.sdf \
#     conformer_gen.sh
#
#   # CSV input, small smoke test first:
#   sbatch --export=INPUT_CSV=new_compounds.csv,MAX_MOLECULES=10,CONFS_PER_MOLECULE=20,\
#     OUTPUT_SDF=/vast/scratch/users/$USER/3D-Descriptors/smoketest-conformers.sdf \
#     conformer_gen.sh

set -euo pipefail

module load miniconda3/latest
module load CUDA/12.8

CONDA_SH_PATH="${CONDA_SH_PATH:-/path/to/miniconda3/etc/profile.d/conda.sh}"
source "$CONDA_SH_PATH"
conda activate nvmolkit

INPUT_SDF="${INPUT_SDF:-}"
INPUT_CSV="${INPUT_CSV:-}"
SMILES_COLUMN="${SMILES_COLUMN:-smiles}"
ID_COLUMN="${ID_COLUMN:-Title}"
OUTPUT_SDF="${OUTPUT_SDF:?Set OUTPUT_SDF at submission time}"
CONFS_PER_MOLECULE="${CONFS_PER_MOLECULE:-1000}"
MAX_MOLECULES="${MAX_MOLECULES:-}"

if [ -z "$INPUT_SDF" ] && [ -z "$INPUT_CSV" ]; then
    echo "Set exactly one of INPUT_SDF or INPUT_CSV at submission time (see header comment)." >&2
    exit 1
fi

ARGS=(
    --output-sdf "$OUTPUT_SDF"
    --confs-per-molecule "$CONFS_PER_MOLECULE"
)
if [ -n "$INPUT_CSV" ]; then
    ARGS+=(--input-csv "$INPUT_CSV" --smiles-column "$SMILES_COLUMN" --id-column "$ID_COLUMN")
else
    ARGS+=(--input-sdf "$INPUT_SDF")
fi
[ -n "$MAX_MOLECULES" ] && ARGS+=(--max-molecules "$MAX_MOLECULES")

python scripts/conformers/generate_conformers.py "${ARGS[@]}"
