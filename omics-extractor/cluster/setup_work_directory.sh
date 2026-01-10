#!/bin/bash
# Setup work directory for cluster enrichment
# Run this script from your work directory on the cluster

set -e

echo "Setting up metadata curation work directory..."

# Create directory structure
echo "Creating directories..."
mkdir -p riboseq_batch_output
mkdir -p logs
mkdir -p config
mkdir -p models  # Optional: for local model downloads

echo "✓ Directory structure created"

# Copy project list from code directory
CODE_DIR=/hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor

if [ -f "$CODE_DIR/riboseq_batch_output/low_completeness_projects.txt" ]; then
    echo "Copying project list..."
    cp "$CODE_DIR/riboseq_batch_output/low_completeness_projects.txt" riboseq_batch_output/
    echo "✓ Project list copied ($(wc -l < riboseq_batch_output/low_completeness_projects.txt) projects)"
else
    echo "⚠ Project list not found. You'll need to copy batch results manually:"
    echo "  rsync -avz local:Metadata-Curation/omics-extractor/riboseq_batch_output/*.json riboseq_batch_output/"
fi

# Copy example config
if [ ! -f "config/model_paths.yaml" ]; then
    echo "Copying example model config..."
    cp "$CODE_DIR/config/model_paths.example.yaml" config/model_paths.yaml
    echo "✓ Config copied to config/model_paths.yaml"
    echo ""
    echo "⚠ IMPORTANT: Edit config/model_paths.yaml with your model path!"
    echo "  nano config/model_paths.yaml"
else
    echo "✓ Config already exists: config/model_paths.yaml"
fi

echo ""
echo "================================================"
echo "Work directory setup complete!"
echo "================================================"
echo ""
echo "Note: Virtual environment (.venv-vllm) will be created here by setup_vllm.sh"
echo ""
echo "Next steps:"
echo "1. Run VLLM setup from code directory (creates venv here):"
echo "   cd /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor"
echo "   bash cluster/setup_vllm.sh"
echo ""
echo "2. Edit config/model_paths.yaml with your model path"
echo ""
echo "3. Copy batch results if not already done:"
echo "   rsync -avz local:Metadata-Curation/omics-extractor/riboseq_batch_output/*.json riboseq_batch_output/"
echo ""
echo "4. Test in interactive session:"
echo "   srun --gres=gpu:1 --cpus-per-task=4 --mem=32G --time=02:00:00 --pty bash"
echo ""
