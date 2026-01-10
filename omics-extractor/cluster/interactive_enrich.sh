#!/bin/bash
# Alternative: Interactive SLURM session for testing/debugging
# Use this for:
# - Testing the enrichment pipeline
# - Debugging model loading
# - Processing 1-3 projects manually

set -e

echo "Starting interactive GPU session..."
echo "This will request:"
echo "  - 1 GPU"
echo "  - 4 CPUs"
echo "  - 32GB RAM"
echo "  - 2 hour time limit"
echo ""

# Request interactive session
srun --partition=gpu \
     --gres=gpu:1 \
     --cpus-per-task=4 \
     --mem=32G \
     --time=02:00:00 \
     --pty bash

# After session starts, you'll be on a GPU node
# Then run:
#
# cd /path/to/Metadata-Curation/omics-extractor
# source .venv-vllm/bin/activate
#
# # Test model loading
# python cluster/test_local_models.py --model llama-3.3-70b
#
# # Enrich a single project
# python -m omics_extractor.cli enrich \
#     riboseq_batch_output/PRJDB10544_traceable.json \
#     --output riboseq_batch_output/PRJDB10544_enriched.json \
#     --only-if-missing \
#     --verbose
