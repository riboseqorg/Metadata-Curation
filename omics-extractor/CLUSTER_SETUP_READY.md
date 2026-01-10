# Cluster Setup - Ready to Run

## Summary

The cluster enrichment pipeline is ready to run on 61 low-completeness RiboSeq projects (1,273 samples).

## Recommendation: Batch Job Array ✅

**Use batch submission** rather than interactive for these reasons:

### Why Batch is Better

1. **Scale:** 61 projects is too many for interactive
2. **Parallelization:** Can use multiple GPUs simultaneously
3. **Fault Tolerance:** Each project is independent - failures don't cascade
4. **Resume Capability:** Automatically skips completed projects
5. **Better Scheduling:** Cluster optimizes GPU allocation
6. **Logging:** Each job has separate output/error logs

### Quick Start

```bash
# 1. Prepare batch (already done!)
cd /path/to/Metadata-Curation/omics-extractor
python cluster/prepare_enrichment_batch.py

# 2. Edit SLURM script with your paths
nano cluster/enrich_batch_slurm.sh
# Change: WORK_DIR=/path/to/Metadata-Curation/omics-extractor

# 3. Configure model paths
cp config/model_paths.example.yaml config/model_paths.yaml
nano config/model_paths.yaml
# Set path to your Llama-3.3-70B model

# 4. Create log directory
mkdir -p logs

# 5. Submit jobs (limit to 10 concurrent GPUs)
sbatch --array=0-60%10 cluster/enrich_batch_slurm.sh

# 6. Monitor
squeue -u $USER
tail -f logs/enrich_*.out
```

---

## When to Use Interactive Instead

Use `cluster/interactive_enrich.sh` only for:

- **Testing:** Verify model loads correctly
- **Debugging:** Troubleshoot errors with direct interaction
- **Small Scale:** Processing 1-3 projects manually

Not recommended for the full 61-project batch.

---

## What's Been Prepared

### Files Created

1. **`cluster/enrich_batch_slurm.sh`** - SLURM job array script
2. **`cluster/prepare_enrichment_batch.py`** - Batch preparation tool
3. **`cluster/interactive_enrich.sh`** - Interactive session script (for testing)
4. **`cluster/CLUSTER_ENRICHMENT_GUIDE.md`** - Complete documentation

### Generated Files (from preparation)

1. **`riboseq_batch_output/low_completeness_projects.txt`** - 61 project IDs
2. **`riboseq_batch_output/enrichment_config.json`** - Resource estimates

---

## The 61 Projects

**Total scope:**
- 61 projects with <50% completeness
- 1,273 samples needing enrichment
- 1,521 sequencing runs

**Estimated resources:**
- ~21 samples per project (average)
- ~0.7 minutes per project (estimated with LLM)
- 2 hours allocated per job (safety margin)

**Expected improvement:**
- Current: 0-50% completeness
- Target: 70-90% completeness after LLM enrichment
- Average gain: ~40-50% increase in completeness

---

## Resource Configuration

### Recommended (Single GPU per Job)

```bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1           # 1 GPU
#SBATCH --cpus-per-task=4      # 4 CPUs
#SBATCH --mem=32G              # 32GB RAM
#SBATCH --time=02:00:00        # 2 hours
#SBATCH --array=0-60%10        # 61 jobs, max 10 concurrent
```

**Requirements:**
- A100 GPU (40GB or 80GB VRAM)
- 70B model with quantization or tensor parallelism
- VLLM for fast inference

### Alternative (Dual GPU for Full Precision)

```bash
#SBATCH --gres=gpu:2           # 2 GPUs
#SBATCH --cpus-per-task=8      # 8 CPUs
#SBATCH --mem=64G              # 64GB RAM
```

**Advantages:**
- Full precision 70B model
- Faster inference
- Better quality (no quantization loss)

---

## Expected Timeline

### Parallel (10 GPUs)
- **Wall time:** 1-2 hours
- **GPU hours:** 10-20 hours total
- **Cost:** ~$20-40 (at $2/GPU-hour)

### Sequential (1 GPU)
- **Wall time:** 10-20 hours
- **GPU hours:** 10-20 hours total
- **Cost:** ~$20-40 (at $2/GPU-hour)

### Comparison to Claude API
- **Claude cost:** 1,273 samples × $0.02 = ~$25
- **Similar cost, but local gives:**
  - Privacy (no external API calls)
  - Reusability (can re-run anytime)
  - Customization (tune models/prompts)
  - Benchmarking (compare different models)

---

## Step-by-Step on Cluster

### 1. SSH to Cluster

```bash
ssh your-cluster.edu
```

### 2. Setup Project

```bash
cd /path/to/Metadata-Curation/omics-extractor

# Setup VLLM environment (if not done)
./cluster/setup_vllm.sh

# Activate
source .venv-vllm/bin/activate
```

### 3. Configure Model

```bash
# Copy and edit model config
cp config/model_paths.example.yaml config/model_paths.yaml
nano config/model_paths.yaml
```

Set your model path:
```yaml
local_models:
  llama-3.3-70b:
    path: "/scratch/models/Llama-3.3-70B-Instruct"  # YOUR PATH HERE
    provider: "vllm"
    tensor_parallel_size: 1  # Use 2 for dual-GPU
    gpu_memory_utilization: 0.9
```

### 4. Edit SLURM Script

```bash
nano cluster/enrich_batch_slurm.sh
```

Change line 16:
```bash
WORK_DIR=/scratch/username/Metadata-Curation/omics-extractor  # YOUR PATH
```

### 5. Create Log Directory

```bash
mkdir -p logs
```

### 6. Submit Jobs

```bash
# Submit with concurrency limit of 10
sbatch --array=0-60%10 cluster/enrich_batch_slurm.sh
```

### 7. Monitor

```bash
# Check queue
squeue -u $USER

# Watch first job
tail -f logs/enrich_*_0.out

# Check progress
ls riboseq_batch_output/*_enriched.json | wc -l
# Should count up to 61
```

### 8. Resume if Needed

If some jobs fail, just resubmit:
```bash
sbatch --array=0-60%10 cluster/enrich_batch_slurm.sh
```

Script automatically skips completed projects.

---

## Verification

### Before Submission

```bash
# Check batch preparation worked
ls riboseq_batch_output/low_completeness_projects.txt
wc -l riboseq_batch_output/low_completeness_projects.txt
# Should show: 61

# Verify model path exists
ls /path/to/your/model/
```

### After First Job Completes

```bash
# Check first enriched file exists
ls riboseq_batch_output/*_enriched.json | head -1

# Compare before/after
PROJECT=$(ls riboseq_batch_output/*_enriched.json | head -1 | xargs basename | sed 's/_enriched.json//')
echo "Project: $PROJECT"

# Show improvement
python -c "
import json
before = json.load(open('riboseq_batch_output/${PROJECT}_traceable.json'))
after = json.load(open('riboseq_batch_output/${PROJECT}_enriched.json'))
print(f'Completeness: {before[\"extraction_statistics\"][\"completeness\"][\"percentage\"]}% → {after[\"extraction_statistics\"][\"completeness\"][\"percentage\"]}%')
"
```

---

## Troubleshooting

See `cluster/CLUSTER_ENRICHMENT_GUIDE.md` for:
- Model loading issues
- Out of memory errors
- Job timeout solutions
- Performance optimization

Quick fixes:
```bash
# Test model loads
srun --partition=gpu --gres=gpu:1 --pty bash
source .venv-vllm/bin/activate
python cluster/test_local_models.py --model llama-3.3-70b

# Check logs for errors
grep -i error logs/enrich_*.err

# Find which jobs failed
sacct -u $USER --format=JobID,State,ExitCode | grep FAILED
```

---

## Summary

✅ **Batch preparation complete** - 61 projects ready
✅ **Scripts ready** - SLURM job array configured
✅ **Documentation complete** - Full guide available
✅ **Resource estimates** - ~1-2 hours on 10 GPUs

**Next action:** Configure model paths and submit jobs on cluster

**Files to customize:**
1. `cluster/enrich_batch_slurm.sh` - Set WORK_DIR
2. `config/model_paths.yaml` - Set model path

**Command to run:**
```bash
sbatch --array=0-60%10 cluster/enrich_batch_slurm.sh
```
