#!/bin/bash
# Test local model on a single project
# This script:
# 1. Extracts metadata from NCBI for a project
# 2. Enriches it using local LLM
# 3. Shows before/after comparison

set -e

WORK_DIR="/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation"
CODE_DIR="/hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor"

# Default test project (RiboSeq study)
PROJECT_ID="${1:-PRJNA1176138}"
OUTPUT_DIR="$WORK_DIR/test_output"

echo "================================================"
echo "Testing Local Model on Project: $PROJECT_ID"
echo "================================================"
echo ""

# Check for GPU
if ! command -v nvidia-smi &> /dev/null; then
    echo "❌ No GPU detected. This script must run on a GPU node."
    echo ""
    echo "Request GPU with:"
    echo "  srun --gres=gpu:1 --cpus-per-task=4 --mem=32G --time=02:00:00 --pty bash"
    exit 1
fi

echo "✓ GPU detected:"
nvidia-smi --query-gpu=name,memory.free --format=csv,noheader
echo ""

# Activate environment
echo "Activating environment..."
source "$WORK_DIR/.venv-vllm/bin/activate"

# Set config path
export OMICS_EXTRACTOR_CONFIG="$WORK_DIR/config/model_paths.yaml"

# Change to work directory
cd "$WORK_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "================================================"
echo "Step 1: Extract baseline metadata (no LLM)"
echo "================================================"
echo ""

python -m omics_extractor.cli extract "$PROJECT_ID" \
    --output "$OUTPUT_DIR/${PROJECT_ID}_baseline.json" \
    --verbose

echo ""
echo "✓ Baseline extraction complete"

# Show baseline stats
echo ""
echo "Baseline statistics:"
python -c "
import json
data = json.load(open('$OUTPUT_DIR/${PROJECT_ID}_baseline.json'))
stats = data['extraction_statistics']
print(f\"  Samples: {stats['sample_count']}\")
print(f\"  Completeness: {stats['completeness']['percentage']:.1f}%\")
print(f\"  High-confidence: {stats['field_statistics']['high_confidence']}/{stats['field_statistics']['total']}\")
"

echo ""
echo "================================================"
echo "Step 2: Enrich with local LLM"
echo "================================================"
echo ""

python -m omics_extractor.cli enrich \
    "$OUTPUT_DIR/${PROJECT_ID}_baseline.json" \
    --output "$OUTPUT_DIR/${PROJECT_ID}_enriched.json" \
    --only-if-missing \
    --verbose

echo ""
echo "✓ Enrichment complete"

# Show enriched stats
echo ""
echo "Enriched statistics:"
python -c "
import json
data = json.load(open('$OUTPUT_DIR/${PROJECT_ID}_enriched.json'))
stats = data['extraction_statistics']
print(f\"  Samples: {stats['sample_count']}\")
print(f\"  Completeness: {stats['completeness']['percentage']:.1f}%\")
print(f\"  High-confidence: {stats['field_statistics']['high_confidence']}/{stats['field_statistics']['total']}\")
"

echo ""
echo "================================================"
echo "Comparison: Baseline vs Enriched"
echo "================================================"
echo ""

python -c "
import json
baseline = json.load(open('$OUTPUT_DIR/${PROJECT_ID}_baseline.json'))
enriched = json.load(open('$OUTPUT_DIR/${PROJECT_ID}_enriched.json'))

b_stats = baseline['extraction_statistics']
e_stats = enriched['extraction_statistics']

print(f\"Completeness:      {b_stats['completeness']['percentage']:.1f}% → {e_stats['completeness']['percentage']:.1f}%\")
print(f\"High-confidence:   {b_stats['field_statistics']['high_confidence']} → {e_stats['field_statistics']['high_confidence']}\")
print(f\"Improvement:       +{e_stats['completeness']['percentage'] - b_stats['completeness']['percentage']:.1f}% completeness\")
"

echo ""
echo "================================================"
echo "✅ Test Complete!"
echo "================================================"
echo ""
echo "Output files:"
echo "  Baseline:  $OUTPUT_DIR/${PROJECT_ID}_baseline.json"
echo "  Enriched:  $OUTPUT_DIR/${PROJECT_ID}_enriched.json"
echo ""
echo "To view full output:"
echo "  cat $OUTPUT_DIR/${PROJECT_ID}_enriched.json | jq ."
echo ""
