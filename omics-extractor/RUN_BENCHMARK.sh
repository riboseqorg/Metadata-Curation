#!/bin/bash
# Quick benchmark runner - compares baseline vs LLM extraction

echo "======================================================================"
echo "OMICS METADATA EXTRACTION BENCHMARK"
echo "======================================================================"
echo ""

# Check for API key
if [ -z "$ANTHROPIC_API_KEY" ]; then
    echo "❌ ERROR: ANTHROPIC_API_KEY not set"
    echo ""
    echo "To run this benchmark:"
    echo "1. Get your API key from: https://console.anthropic.com/settings/keys"
    echo "2. Run:"
    echo "   export ANTHROPIC_API_KEY=sk-ant-your-key-here"
    echo "   ./RUN_BENCHMARK.sh"
    echo ""
    exit 1
fi

echo "✓ API key found"
echo "✓ Running benchmark on real RiboSeq data..."
echo ""

# Run the benchmark
~/.local/bin/uv run python scripts/quick_benchmark.py > benchmark_results.txt

echo ""
echo "======================================================================"
echo "BENCHMARK COMPLETE"
echo "======================================================================"
