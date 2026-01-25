# Cluster Setup

Quick setup for local models on cluster.

## Setup

```bash
# 1. Setup VLLM + download model
bash omics-extractor/cluster/setup_vllm.sh

# 2. Setup config
cp omics-extractor/config/model_paths.example.yaml config/model_paths.yaml
```

## Usage

```bash
# Extract metadata (CPU)
omics-extract extract PRJNA1170270

# Enrich with local model (GPU)
omics-extract enrich PRJNA1170270_metadata.json --model qwen-2.5-7b

# Batch process many projects (GPU)
omics-extract batch-enrich *_metadata.json --output-dir enriched/ --model mistral-7b
```

## Models

**Free (no auth):**
- `qwen-2.5-7b` - 16GB, single GPU, best quality
- `mistral-7b` - 14GB, single GPU, good alternative

Download: `bash omics-extractor/cluster/download_model.sh --model <name>`

## Notes

- GPU count auto-detected
- `tensor_parallel_size` auto-adjusted based on available GPUs
- VLLM batch processing is automatic
