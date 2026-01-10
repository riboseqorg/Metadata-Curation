# How to Test Local Models on Your Cluster

## TL;DR - Copy/Paste This

```bash
# On your cluster:
cd /path/to/Metadata-Curation/omics-extractor

# 1. Setup VLLM (one-time, ~10 minutes)
./cluster/setup_vllm.sh

# 2. Configure your model paths
cp config/model_paths.example.yaml config/model_paths.yaml
# Edit model_paths.yaml with your actual paths

# 3. Activate environment and test
source .venv-vllm/bin/activate
python cluster/test_local_models.py --model llama-3.3-70b
```

## Prerequisites

- ✅ A100 GPU(s) available
- ✅ CUDA installed
- ✅ Python 3.9+
- ✅ Model weights downloaded (or use script below)

## Step-by-Step

### Step 1: Get Model Weights (if you don't have them)

```bash
# Install huggingface CLI
pip install huggingface_hub

# Login (needed for Llama)
huggingface-cli login
# Paste your token from https://huggingface.co/settings/tokens

# Download Llama 3.3 70B (~140GB, takes 1-2 hours)
huggingface-cli download meta-llama/Llama-3.3-70B-Instruct \
  --local-dir /scratch/models/Llama-3.3-70B-Instruct

# Or use existing models if already on your cluster
# Ask your sysadmin where shared models are stored
```

### Step 2: Setup VLLM

```bash
cd /path/to/Metadata-Curation/omics-extractor

# Run setup script
chmod +x cluster/setup_vllm.sh
./cluster/setup_vllm.sh

# This creates .venv-vllm with VLLM installed
```

### Step 3: Configure Model Paths

```bash
# Copy example config
cp config/model_paths.example.yaml config/model_paths.yaml

# Edit with your paths
nano config/model_paths.yaml
```

Update the paths to match your setup:

```yaml
local_models:
  llama-3.3-70b:
    path: "/scratch/models/Llama-3.3-70B-Instruct"  # ← Your path here
    provider: "vllm"
    tensor_parallel_size: 2  # Use 2 GPUs if you have them
    gpu_memory_utilization: 0.9
```

### Step 4: Test

```bash
# Activate environment
source .venv-vllm/bin/activate

# Test single model
python cluster/test_local_models.py --model llama-3.3-70b

# Or test all models
python cluster/test_local_models.py
```

## Expected Output

```
==================================================
LOCAL MODEL TESTING ON CLUSTER
==================================================

Config: config/model_paths.yaml
Project: PRJNA1170270
Models to test: 1

==================================================
Testing: llama-3.3-70b
==================================================

Configuration:
  path: /scratch/models/Llama-3.3-70B-Instruct
  provider: vllm
  tensor_parallel_size: 2

1. Extracting baseline metadata from PRJNA1170270...
Fetching runs for PRJNA1170270...
Found 4 runs
Building study metadata...
Building run and sample metadata...
  Processing run 1/4: SRR30910340
  Processing run 2/4: SRR30910341
  Processing run 3/4: SRR30910342
  Processing run 4/4: SRR30910343
Complete! Study: 1, Samples: 1, Runs: 4
   ✓ Extracted 1 samples

2. Gap analysis:
   Total samples: 1
   Need LLM: 0

   ✓ All samples complete - no LLM needed!

==================================================
✅ Testing Complete!
==================================================
```

## Testing on a Sample That Needs LLM

```bash
# Test on a project with missing metadata
python cluster/test_local_models.py --model llama-3.3-70b --project PRJNA1234567
```

This will show:
- What baseline extraction found
- What the LLM added
- How fast it was (latency)
- Token usage

## Troubleshooting

### "CUDA out of memory"
```yaml
# Reduce GPU memory usage in config
gpu_memory_utilization: 0.7  # Try lower value
```

### "Model path not found"
```bash
# Check your path is correct
ls /scratch/models/Llama-3.3-70B-Instruct
# Should show: config.json, pytorch_model.bin, etc.
```

### "vllm not found"
```bash
# Make sure you activated the environment
source .venv-vllm/bin/activate

# Check vllm is installed
pip list | grep vllm
```

## After Testing

Once you know which model works best:

1. Note the speed (latency per sample)
2. Note the quality (what it correctly fills)
3. Use that model in production

See `CLUSTER_SETUP_COMPLETE.md` for production usage.

## Getting Help

- Check `cluster/README.md` for detailed documentation
- See `CLUSTER_SETUP_COMPLETE.md` for complete overview
- See `docs/PRACTICAL_LLM_STRATEGY.md` for system architecture

That's it! You're ready to test local models on your cluster.
