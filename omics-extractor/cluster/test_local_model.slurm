#!/bin/bash
#SBATCH --job-name=test-local-model
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/logs/test_model_%j.out
#SBATCH --error=/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/logs/test_model_%j.err

# Test job for local model metadata extraction
# Submits a GPU job to test the local model setup

set -e

CODE_DIR="/hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor"

echo "Starting local model test job..."
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo ""

# Run test script
bash "$CODE_DIR/cluster/test_local_model.sh" "${1:-PRJNA1176138}"

echo ""
echo "Test job complete!"
