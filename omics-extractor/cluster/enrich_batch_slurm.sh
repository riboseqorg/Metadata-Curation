#!/bin/bash
#SBATCH --job-name=riboseq-enrich
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/enrich_%A_%a.out
#SBATCH --error=logs/enrich_%A_%a.err
#SBATCH --array=0-60  # 61 low-completeness projects

# This script enriches one project using a local LLM on GPU
# Launched as a job array to process multiple projects in parallel

set -e

# Setup directories
# CODE_DIR: Where the omics-extractor code is cloned
# WORK_DIR: Where data, outputs, and config live (and venv)
CODE_DIR=/hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor
WORK_DIR=/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation

# Activate environment from work directory (where venv lives)
source $WORK_DIR/.venv-vllm/bin/activate

# Set config path (in work directory)
export OMICS_EXTRACTOR_CONFIG=$WORK_DIR/config/model_paths.yaml

# Change to work directory for data access
cd $WORK_DIR

# Get list of low-completeness projects from batch summary
PROJECTS_FILE="riboseq_batch_output/low_completeness_projects.txt"

# Get the project for this array task
PROJECT_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $PROJECTS_FILE)

echo "================================================"
echo "SLURM Array Task: $SLURM_ARRAY_TASK_ID"
echo "Processing project: $PROJECT_ID"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "Node: $SLURM_NODELIST"
echo "================================================"

# Check if project already enriched (for resume capability)
INPUT_FILE="riboseq_batch_output/${PROJECT_ID}_traceable.json"
OUTPUT_FILE="riboseq_batch_output/${PROJECT_ID}_enriched.json"

if [ -f "$OUTPUT_FILE" ]; then
    echo "Project $PROJECT_ID already enriched, skipping..."
    exit 0
fi

# Run enrichment with local model
python -m omics_extractor.cli enrich \
    "$INPUT_FILE" \
    --output "$OUTPUT_FILE" \
    --only-if-missing \
    --verbose

echo "✓ Project $PROJECT_ID enriched successfully"
echo "Output: $OUTPUT_FILE"
