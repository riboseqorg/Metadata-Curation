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

# Install huggingface-hub for model downloads
echo "Installing huggingface-hub..."
pip install huggingface-hub[cli]

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
echo "=================================="
echo "Model Download Options"
echo "=================================="
echo
echo "Would you like to download a model now?"
echo "This is required for testing but can be done later."
echo
read -p "Download model now? (Y/n): " download_choice

if [[ ! "$download_choice" =~ ^[Nn]$ ]]; then
    # Get script directory
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    echo
    bash "$SCRIPT_DIR/download_model.sh"
    echo
    echo "Next steps:"
    echo "1. Test in interactive GPU session (see cluster/SETUP_INSTRUCTIONS.md)"
else
    echo
    echo "To download a model later, run:"
    echo "  bash omics-extractor/cluster/download_model.sh"
    echo
    echo "Or non-interactively:"
    echo "  bash omics-extractor/cluster/download_model.sh --model llama-70b"
    echo "  bash omics-extractor/cluster/download_model.sh --model llama-8b --skip-existing"
    echo
    echo "Next steps:"
    echo "1. Download a model (see above)"
    echo "2. Test in interactive GPU session (see cluster/SETUP_INSTRUCTIONS.md)"
fi
echo

