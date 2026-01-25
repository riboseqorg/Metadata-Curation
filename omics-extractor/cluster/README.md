# Cluster Setup

Quick setup for local models on cluster.

## Setup

```bash
# 1. Run setup (installs VLLM + downloads model)
bash omics-extractor/cluster/setup_vllm.sh

# 2. Test it
python omics-extractor/cluster/test_local_models.py --model qwen-2.5-7b
```

## Models

**Free (no auth):**
- `qwen-7b` - Download: `bash omics-extractor/cluster/download_model.sh --model qwen-7b`
- `mistral` - Download: `bash omics-extractor/cluster/download_model.sh --model mistral`

**Gated (requires HF login):**
- `llama-8b` - Requires: `huggingface-cli login`

## Config

Copy: `cp omics-extractor/config/model_paths.example.yaml config/model_paths.yaml`

GPU count auto-detected. If config says `tensor_parallel_size: 2` but you have 1 GPU, it auto-adjusts.

## Troubleshooting

**"World size (2) > available GPUs (1)"**: Fixed automatically by GPU auto-detection.

**"Cannot access gated repo"**: Use free model or run `huggingface-cli login`

**"CUDA out of memory"**: Use smaller model (qwen-7b instead of qwen-72b)
