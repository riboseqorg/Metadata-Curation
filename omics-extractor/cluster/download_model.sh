#!/bin/bash
# Download a local LLM model for metadata extraction
# This downloads to the work directory where we have space

set -e

WORK_DIR="/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation"
MODEL_DIR="$WORK_DIR/models"

echo "=================================="
echo "Model Download for Local LLM"
echo "=================================="
echo ""

# Create models directory
mkdir -p "$MODEL_DIR"

echo "Models will be downloaded to: $MODEL_DIR"
echo ""

# Recommended model options (in order of size/performance)
cat <<EOF
Available models for metadata extraction:

1. Llama-3.3-70B-Instruct (RECOMMENDED)
   - Size: ~140GB
   - Best quality for metadata extraction
   - Requires: 1x A100-80GB or 2x A100-40GB
   - Download time: ~1-2 hours

2. Llama-3.1-8B-Instruct (SMALLER, FASTER)
   - Size: ~16GB
   - Good quality, faster inference
   - Requires: 1x GPU with 24GB+ VRAM
   - Download time: ~10-15 minutes

Which model would you like to download?
EOF

read -p "Enter choice (1 or 2): " choice

case $choice in
    1)
        MODEL_NAME="meta-llama/Llama-3.3-70B-Instruct"
        MODEL_DIR_NAME="Llama-3.3-70B-Instruct"
        ;;
    2)
        MODEL_NAME="meta-llama/Llama-3.1-8B-Instruct"
        MODEL_DIR_NAME="Llama-3.1-8B-Instruct"
        ;;
    *)
        echo "Invalid choice"
        exit 1
        ;;
esac

MODEL_PATH="$MODEL_DIR/$MODEL_DIR_NAME"

# Check if model already exists
if [ -d "$MODEL_PATH" ]; then
    echo ""
    echo "⚠ Model already exists at: $MODEL_PATH"
    read -p "Re-download? (y/N): " confirm
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
        echo "Using existing model"
        exit 0
    fi
    rm -rf "$MODEL_PATH"
fi

echo ""
echo "Downloading: $MODEL_NAME"
echo "To: $MODEL_PATH"
echo ""
echo "This may take a while..."
echo ""

# Check for huggingface-cli
if ! command -v huggingface-cli &> /dev/null; then
    echo "Installing huggingface-hub..."
    pip install huggingface-hub[cli]
fi

# Download model
echo "Starting download..."
huggingface-cli download "$MODEL_NAME" \
    --local-dir "$MODEL_PATH" \
    --local-dir-use-symlinks False

echo ""
echo "=================================="
echo "✅ Model Download Complete!"
echo "=================================="
echo ""
echo "Model location: $MODEL_PATH"
echo ""
echo "Next step: Configure model path"
echo "  nano $WORK_DIR/config/model_paths.yaml"
echo ""
echo "Set the path to:"
echo "  path: \"$MODEL_PATH\""
echo ""
