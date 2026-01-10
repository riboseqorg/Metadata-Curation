### Cluster Setup for Local Model Testing

This directory contains scripts for testing local LLM models on your A100 cluster.

## Quick Start

### 1. Setup VLLM (One-time)

On your cluster, run:

```bash
# Navigate to omics-extractor
cd /path/to/omics-extractor

# Run setup script (installs VLLM and dependencies)
./cluster/setup_vllm.sh
```

This will:
- Create a virtual environment (`.venv-vllm`)
- Install VLLM with CUDA support
- Install required dependencies

### 2. Configure Model Paths

Copy the example config and update with your model paths:

```bash
cp config/model_paths.example.yaml config/model_paths.yaml
```

Edit `config/model_paths.yaml`:

```yaml
local_models:
  llama-3.3-70b:
    path: "/data/models/Llama-3.3-70B-Instruct"  # Your actual path
    provider: "vllm"
    tensor_parallel_size: 2
    gpu_memory_utilization: 0.9
```

### 3. Test Models

Activate the environment and run tests:

```bash
# Activate VLLM environment
source .venv-vllm/bin/activate

# Test all models
python cluster/test_local_models.py

# Or test specific model
python cluster/test_local_models.py --model llama-3.3-70b

# Or test on specific project
python cluster/test_local_models.py --model qwen-2.5-72b --project PRJNA1234567
```

## Getting Model Weights

### Option A: Download from HuggingFace

If you don't have models downloaded yet:

```bash
# Install huggingface_hub
pip install huggingface_hub

# Download Llama 3.3 70B (requires ~140GB disk space)
huggingface-cli download meta-llama/Llama-3.3-70B-Instruct \
  --local-dir /data/models/Llama-3.3-70B-Instruct

# Download Qwen 2.5 72B
huggingface-cli download Qwen/Qwen2.5-72B-Instruct \
  --local-dir /data/models/Qwen2.5-72B-Instruct

# Download Mixtral 8x22B
huggingface-cli download mistralai/Mixtral-8x22B-Instruct-v0.1 \
  --local-dir /data/models/Mixtral-8x22B-Instruct-v0.1
```

**Note**: You may need HuggingFace tokens for some models (especially Llama).
Get one from: https://huggingface.co/settings/tokens

```bash
huggingface-cli login
```

### Option B: Use Existing Models

If models are already on your cluster, just point to them in `config/model_paths.yaml`.

## Recommended Models for Testing

### Llama 3.3 70B Instruct (Top Choice)
- **Best balance** of quality and speed
- **Memory**: ~140GB disk, ~40GB GPU memory
- **Speed**: ~400-600ms per sample
- **Quality**: Very close to Claude

### Qwen 2.5 72B Instruct
- **Excellent** at structured extraction
- **Memory**: ~144GB disk, ~40GB GPU memory
- **Speed**: ~500-700ms per sample
- **Quality**: Sometimes better than Llama for metadata

### Mixtral 8x22B Instruct
- **Very fast** with good quality
- **Memory**: ~150GB disk, ~40GB GPU memory
- **Speed**: ~300-500ms per sample
- **Quality**: Good but slightly lower than Llama/Qwen

## GPU Requirements

### For 70B Models:
- **Minimum**: 1x A100 80GB
- **Recommended**: 2x A100 40GB (tensor parallelism)
- Set `tensor_parallel_size: 2` in config

### For 8B Models (testing only):
- **Minimum**: 1x A100 40GB
- Can also run on single GPU

## What the Test Script Does

1. **Extracts baseline metadata** using structured fields + ontologies
2. **Identifies samples** that need LLM enrichment (missing critical fields)
3. **Tests each local model** on those samples
4. **Measures speed** (latency per sample)
5. **Shows results** (what the LLM added/filled)

### Example Output:

```
==================================================
Testing: llama-3.3-70b
==================================================

Configuration:
  path: /data/models/Llama-3.3-70B-Instruct
  provider: vllm
  tensor_parallel_size: 2

1. Extracting baseline metadata from PRJNA1170270...
   ✓ Extracted 1 samples

2. Gap analysis:
   Total samples: 1
   Need LLM: 0

   ✓ All samples complete - no LLM needed!
```

Or if samples need enrichment:

```
5. Running LLM extraction...
   ✓ Extraction complete!
   Time: 0.58s
   Tokens: 156

6. Results:
   Before LLM:
     tissue: testis
     cell_type: [MISSING]
     strain: C57BL/6

   After LLM:
     tissue: testis
     cell_type: spermatocyte
     strain: C57BL/6
```

## Comparing with Claude

To compare local models with Claude API baseline:

1. Set `ANTHROPIC_API_KEY` environment variable
2. Run the same test with Claude:

```bash
export ANTHROPIC_API_KEY=sk-ant-your-key
python scripts/demo_extraction_modes.py
```

3. Compare results:
   - Which model matches Claude's output?
   - Which is faster?
   - Which fills more gaps correctly?

## Troubleshooting

### "CUDA out of memory"
- Reduce `gpu_memory_utilization` in config (try 0.8)
- Increase `tensor_parallel_size` to use more GPUs
- Try a smaller model first (8B)

### "Model not found"
- Check path in `config/model_paths.yaml` is correct
- Ensure model weights are actually downloaded
- Check permissions on model directory

### "Import error: vllm not found"
- Make sure you activated the vllm environment:
  ```bash
  source .venv-vllm/bin/activate
  ```
- Or re-run setup: `./cluster/setup_vllm.sh`

### Slow inference
- Check you're using `tensor_parallel_size: 2` for 70B models
- Ensure CUDA is working: `nvidia-smi`
- Try VLLM optimizations in config

## Production Deployment

After testing, to use local models in production:

1. Choose the best-performing model from tests
2. Update main extraction pipeline to use it:
   ```python
   from omics_extractor.extraction.llm_providers import create_provider

   provider = create_provider(
       "vllm",
       model_path="/data/models/Llama-3.3-70B-Instruct",
       tensor_parallel_size=2
   )
   ```
3. Run at scale on your datasets

## Cost Savings

Using local 70B models instead of Claude API:

- **Claude**: $0.01-0.03 per sample
- **Local 70B**: $0 per sample (just electricity + GPU amortization)

For 10,000 samples:
- Claude: $100-300
- Local: ~$10-20 in compute costs

**Breakeven**: After ~1,000-5,000 samples, local models are cheaper!
