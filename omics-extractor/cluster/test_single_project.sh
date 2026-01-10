#!/bin/bash
# Quick test script to extract metadata for a single project
# Uses Claude API (no local model needed)

set -e

# Check if project ID provided
if [ -z "$1" ]; then
    echo "Usage: $0 <PROJECT_ID>"
    echo ""
    echo "Example:"
    echo "  $0 PRJNA1176138"
    echo ""
    echo "This will extract metadata using Claude API (no local model needed)"
    exit 1
fi

PROJECT_ID=$1
OUTPUT_DIR="test_output"

echo "================================================"
echo "Testing metadata extraction for: $PROJECT_ID"
echo "================================================"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check for API key
if [ -z "$ANTHROPIC_API_KEY" ]; then
    echo "⚠ Warning: ANTHROPIC_API_KEY not set"
    echo "  Set it with: export ANTHROPIC_API_KEY=your-key-here"
    echo "  Or extraction will use baseline (no LLM enrichment)"
    echo ""
fi

# Run extraction with traceable output
echo "Extracting metadata..."
python -m omics_extractor.cli extract "$PROJECT_ID" \
    --output "$OUTPUT_DIR/${PROJECT_ID}_traceable.json" \
    --verbose

echo ""
echo "================================================"
echo "✅ Extraction complete!"
echo "================================================"
echo ""
echo "Output file: $OUTPUT_DIR/${PROJECT_ID}_traceable.json"
echo ""

# Show quick summary
echo "Quick summary:"
python -c "
import json
data = json.load(open('$OUTPUT_DIR/${PROJECT_ID}_traceable.json'))
stats = data['extraction_statistics']
print(f\"  Samples: {stats['sample_count']}\")
print(f\"  Runs: {stats['run_count']}\")
print(f\"  Completeness: {stats['completeness']['percentage']:.1f}%\")
print(f\"  High-confidence fields: {stats['field_statistics']['high_confidence']}\")
print(f\"  Total fields: {stats['field_statistics']['total']}\")
"

echo ""
echo "To view full details:"
echo "  cat $OUTPUT_DIR/${PROJECT_ID}_traceable.json | jq ."
echo ""
