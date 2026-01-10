#!/bin/bash
# Setup script for VLLM on cluster with A100 GPUs
#
# This script downloads and sets up local LLM models for testing
# Run this once on your cluster to prepare the environment

set -e  # Exit on error

echo "=================================="
echo "VLLM Setup for Cluster"
echo "=================================="
echo

# Check for GPU (optional - can setup on login node)
if command -v nvidia-smi &> /dev/null; then
    echo "✓ Found NVIDIA GPU(s):"
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
    echo
else
    echo "ℹ No GPU detected (running on login node)"
    echo "  VLLM will be installed but requires GPU to run"
    echo
fi

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | cut -d. -f1,2)
echo "Python version: $PYTHON_VERSION"

if (( $(echo "$PYTHON_VERSION < 3.9" | bc -l) )); then
    echo "❌ Python 3.9+ required"
    exit 1
fi

# Create virtual environment if it doesn't exist
if [ ! -d ".venv-vllm" ]; then
    echo "Creating virtual environment..."
    python3 -m venv .venv-vllm
fi

# Activate virtual environment
source .venv-vllm/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install VLLM (this will take a few minutes)
echo "Installing VLLM..."
echo "  (This may take 5-10 minutes on first install)"
pip install vllm

# Install other dependencies
echo "Installing dependencies..."
pip install torch transformers accelerate

# Install omics-extractor package in development mode
echo "Installing omics-extractor package..."
pip install -e .

echo
echo "=================================="
echo "✅ VLLM Setup Complete!"
echo "=================================="
echo
echo "Virtual environment created: .venv-vllm"
echo
echo "To activate in future sessions:"
echo "  source .venv-vllm/bin/activate"
echo
echo "Next steps:"
echo "1. Setup work directory:"
echo "   cd /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation"
echo "   bash /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor/cluster/setup_work_directory.sh"
echo
echo "2. Copy batch results from local machine"
echo
echo "3. Configure model path in work directory:"
echo "   nano /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/config/model_paths.yaml"
echo
echo "4. Test in interactive GPU session (see cluster/SETUP_INSTRUCTIONS.md)"
echo

