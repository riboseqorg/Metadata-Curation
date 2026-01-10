# Quick Start: Local Model Testing on Cluster

Complete workflow to test local LLM-based metadata extraction on the cluster.

## Overview

This will:
1. Install Python dependencies
2. Download a local LLM (Llama 3.1 8B or Llama 3.3 70B)
3. Test extraction + enrichment on a single project
4. Show before/after comparison

## Step 1: Setup Environment (5-10 minutes)

```bash
# Pull latest code
cd /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor
git pull

# Setup work directory
cd /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation
bash /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor/cluster/setup_work_directory.sh

# Install Python packages (creates venv in work directory)
cd /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor
bash cluster/setup_vllm.sh
```

**Wait for installation to complete** (~5-10 minutes)

## Step 2: Download Model (10 minutes - 2 hours)

```bash
cd /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation
bash /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor/cluster/download_model.sh
```

**Choose a model:**
- Option 1: Llama-3.3-70B (~140GB, best quality, ~1-2 hour download)
- Option 2: Llama-3.1-8B (~16GB, faster, ~10-15 min download)

**Recommendation for testing:** Start with Option 2 (8B model) for faster testing.

## Step 3: Configure Model Path

After download completes:

```bash
cd /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation
nano config/model_paths.yaml
```

**For 8B model:**
```yaml
local_models:
  llama-3.1-8b:
    path: "/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/models/Llama-3.1-8B-Instruct"
    provider: "vllm"
    tensor_parallel_size: 1
    gpu_memory_utilization: 0.9

default: "llama-3.1-8b"
```

**For 70B model:**
```yaml
local_models:
  llama-3.3-70b:
    path: "/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/models/Llama-3.3-70B-Instruct"
    provider: "vllm"
    tensor_parallel_size: 2  # Use 2 GPUs
    gpu_memory_utilization: 0.9

default: "llama-3.3-70b"
```

## Step 4: Test with Interactive Session

```bash
# Request GPU
srun --gres=gpu:1 --cpus-per-task=4 --mem=32G --time=01:00:00 --pty bash

# Once on GPU node, run test
cd /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor
bash cluster/test_local_model.sh PRJNA1176138
```

This will:
1. Extract metadata from NCBI (baseline)
2. Enrich with local LLM
3. Show before/after comparison

## Step 5: Or Submit as SLURM Job

```bash
# Create logs directory
mkdir -p /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/logs

# Submit test job
cd /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor
sbatch cluster/test_local_model.slurm PRJNA1176138

# Monitor job
squeue -u $USER
tail -f /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/logs/test_model_*.out
```

## Expected Output

```
================================================
Comparison: Baseline vs Enriched
================================================

Completeness:      35.0% → 87.5%
High-confidence:   14 → 35
Improvement:       +52.5% completeness

✅ Test Complete!
```

## Troubleshooting

### Model doesn't load
```bash
# Check GPU
nvidia-smi

# Verify model path exists
ls /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/models/

# Check config
cat /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/config/model_paths.yaml
```

### Out of memory
- Use smaller model (8B instead of 70B)
- Or request more memory: `--mem=64G`
- Or use 2 GPUs for 70B: `--gres=gpu:2`

### Download interrupted
Re-run download script - it will resume from where it left off.

## Next Steps

After successful test:

1. **Test more projects**: Try different project IDs
2. **Run batch enrichment**: Use the batch scripts to enrich all 61 low-completeness projects
3. **Compare models**: Test both 8B and 70B to see quality difference

## Common Projects to Test

RiboSeq projects known to exist:
- PRJNA1176138 (good test case)
- PRJNA449378
- PRJNA1170270
- PRJDB10544

Usage:
```bash
bash cluster/test_local_model.sh PRJNA449378
```
