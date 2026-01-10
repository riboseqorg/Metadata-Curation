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
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d. -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d. -f2)
echo "Python version: $PYTHON_VERSION"

if [ "$PYTHON_MAJOR" -lt 3 ] || ([ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -lt 9 ]); then
    echo "❌ Python 3.9+ required (found $PYTHON_VERSION)"
    exit 1
fi
echo "✓ Python version OK"
echo

# Virtual environment location (in work directory for space)
VENV_DIR="/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/.venv-vllm"

# Create virtual environment if it doesn't exist
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment in work directory..."
    echo "Location: $VENV_DIR"
    python3 -m venv "$VENV_DIR"
else
    echo "✓ Virtual environment already exists: $VENV_DIR"
fi

# Activate virtual environment
echo "Activating virtual environment..."
source "$VENV_DIR/bin/activate"

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
echo "Virtual environment location:"
echo "  $VENV_DIR"
echo
echo "To activate in future sessions:"
echo "  source $VENV_DIR/bin/activate"
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

