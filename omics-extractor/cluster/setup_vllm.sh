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

# Check for GPU
if ! command -v nvidia-smi &> /dev/null; then
    echo "❌ nvidia-smi not found. This script requires NVIDIA GPUs."
    exit 1
fi

echo "✓ Found NVIDIA GPU(s):"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo

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

echo
echo "=================================="
echo "✅ VLLM Setup Complete!"
echo "=================================="
echo
echo "Next steps:"
echo "1. Download model weights (or use existing paths)"
echo "2. Update config/model_paths.yaml with your model paths"
echo "3. Run test_local_models.py to verify setup"
echo
echo "To download models (example for Llama 3.3 70B):"
echo "  huggingface-cli download meta-llama/Llama-3.3-70B-Instruct --local-dir /path/to/models/Llama-3.3-70B-Instruct"
echo
echo "Or use existing model paths if already downloaded"
echo

